"""Synthetic tests of recon/recon_wave_tse_from_twix_nifti.py helpers."""

import argparse
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import nibabel as nib
import numpy as np
import pypulseq as pp
import torch

import recon_wave_tse_from_twix_nifti as cli
from utils.l1_noise_model import colored_noise
from utils.mb2_cg_sense_wave import cg_sense_mb2_wave_2d, mb2_wave_forward
from utils.wave_cg_sense_precondition import fftc_dim


def unit_rss_sensitivities(ncoil, nx, nx_os, ny, seed):
    """Return random smooth unit-RSS sensitivities on the native readout rows.

    Returns:
        torch.Tensor: complex64 (coil, Nx_os, Ny), zero outside the centered
        native rows.
    """
    rng = np.random.default_rng(seed)
    x, y = np.mgrid[:nx, :ny]
    maps = np.stack([np.exp(-((x - rng.uniform(0, nx))**2 + (y - rng.uniform(0, ny))**2) / 50)
                     * np.exp(1j * rng.uniform(0, 2 * np.pi)) for _ in range(ncoil)])
    maps /= np.sqrt(np.sum(abs(maps)**2, axis=0))
    sens = np.zeros((ncoil, nx_os, ny), np.complex64)
    sens[:, cli.native_rows(nx, nx_os)] = maps
    return torch.as_tensor(sens)


def synthetic_mb2_group(ncoil=8, nx=32, nx_os=96, ny=128, seed=10, psf=None):
    """Return ``group_terms``-like MB2 Ry2 terms, the receiver profile and the CAIPI-encoded truth.

    Broad smooth unit-RSS coil maps cover the native rows; band 1 carries a
    FOV/4 CAIPI shift. The data are the encoded truth plus colored noise with
    correlated coils and a readout profile that is flat (1.3) inside and
    rolls off to 0.3 at both readout edges.

    Returns:
        tuple: ``(terms, profile)``.
    """
    rng = np.random.default_rng(seed)
    x, y = np.mgrid[:nx, :ny] / np.array([nx, ny])[:, None, None]
    sens = np.zeros((2, ncoil, nx_os, ny), np.complex64)
    for band in range(2):
        maps = np.stack([np.exp(-((x - rng.uniform()) ** 2 + (y - rng.uniform()) ** 2) / .08)
                         * np.exp(2j * np.pi * (rng.uniform() + .3 * rng.uniform() * (x + y)))
                         for _ in range(ncoil)])
        sens[band, :, cli.native_rows(nx, nx_os)] = maps / np.sqrt(np.sum(abs(maps) ** 2, axis=0))
    phase = np.ones((2, ny), np.complex64)
    phase[1] = np.exp(2j * np.pi * np.arange(ny) / 4)
    mask = np.zeros(ny, np.float32)
    mask[::2] = 1
    truth = np.zeros((2, nx_os, ny), np.complex64)
    truth[:, cli.native_rows(nx, nx_os)] = 20 * (np.cos(4 * x) + 1j * np.sin(3 * y))[None]
    a = rng.standard_normal((ncoil, ncoil)) + 1j * rng.standard_normal((ncoil, ncoil))
    cov = a @ a.conj().T / ncoil + np.eye(ncoil)
    u = abs(np.arange(nx_os) - nx_os / 2) / (nx_os / 2)
    profile = np.where(u < .6, 1.3, 1.3 - (u - .6) / .4)
    terms = {"sens": torch.as_tensor(sens), "phase": torch.as_tensor(phase),
             "psf": torch.ones((2, nx_os, ny), dtype=torch.complex64) if psf is None else psf,
             "mask": torch.as_tensor(mask[None, None]), "raw_sli": 0, "physical": np.array([0, 1]),
             "lines": np.flatnonzero(mask)}
    clean = mb2_wave_forward(torch.as_tensor(truth), terms["sens"], terms["phase"], terms["psf"], terms["mask"])
    terms["data"] = (clean + colored_noise(cov, profile, (ncoil, nx_os, ny), seed + 1)) * terms["mask"]
    return terms, profile


def lambda0_fit(terms):
    """Return the lambda-0 solution and its predicted k-space for ``terms``."""
    solution = cg_sense_mb2_wave_2d(terms["data"], terms["sens"], terms["phase"], terms["psf"], terms["mask"],
                                    n_iter=500, tol=1e-8, verbose=False)
    return solution, mb2_wave_forward(solution, terms["sens"], terms["phase"], terms["psf"], terms["mask"])


def fake_slice_worker(physical, twix, labels, folders, basis):
    """Stand in for the ESPIRiT worker: save the slice index and the worker's OMP_NUM_THREADS.

    Returns:
        float: Zero seconds.
    """
    np.savez(folders["arrays"] / f"csm_diagnostics_{physical:02d}.npz", threads=os.environ.get("OMP_NUM_THREADS"))
    np.save(folders["cache"] / f"csm_slice_{physical:02d}.npy", np.full(2, physical))
    return 0.


def translated_excitation_sequence(positions, sign):
    """Return a Pulseq sequence whose excitations encode ``positions`` in the RF waveform.

    Every excitation plays the same envelope times ``exp(sign*2i*pi*z*m(t))``
    on a constant slice gradient, followed by one ADC block; ``sign=-1`` is the
    convention accepted by the reconstruction.

    Returns:
        pypulseq.Sequence: In-memory sequence with one ADC per position.
    """
    system = pp.Opts(rf_dead_time=100e-6, rf_ringdown_time=30e-6, adc_dead_time=10e-6)
    dwell, n = 10e-6, 200
    t = (np.arange(n) + .5) * dwell
    envelope = np.sinc(np.linspace(-3, 3, n)) * np.hanning(n)
    gz = pp.make_trapezoid("z", system=system, amplitude=2e5, flat_time=n * dwell + 2 * system.rf_dead_time,
                           rise_time=200e-6)
    adc = pp.make_adc(num_samples=16, dwell=10e-6, delay=system.adc_dead_time, system=system)
    seq = pp.Sequence(system)
    for z in positions:
        rf = pp.make_arbitrary_rf(envelope.astype(complex), np.pi / 2, dwell=dwell, delay=gz.rise_time + 50e-6,
                                  system=system, use="excitation", no_signal_scaling=True)
        moment = gz.amplitude * (t - rf.center)
        rf.signal = envelope * np.exp(sign * 2j * np.pi * z * moment)
        seq.add_block(rf, gz)
        seq.add_block(adc)
    return seq


class GeometryHelperTests(unittest.TestCase):
    """Check readout rows, PE coordinates, PSF tags and the Wave PSF convention."""

    def test_native_rows_and_phase_coordinates(self):
        """The nominal FOV is centered in the oversampled readout; PE coordinates are centered."""
        self.assertEqual(cli.native_rows(440, 1320), slice(440, 880))
        with self.assertRaises(ValueError):
            cli.native_rows(440, 1000)
        y = cli.phase_coordinates_m(220, .22)
        self.assertAlmostEqual(y[110], 0)
        self.assertAlmostEqual(y[0], -.11)
        self.assertAlmostEqual(y[1] - y[0], .001)

    def test_psf_tags(self):
        """Wave PSF settings map to distinct file tags."""
        self.assertEqual(cli.psf_tag(None, False), "auto")
        self.assertEqual(cli.psf_tag([1120, 1520], False), "range1120-1520")
        self.assertEqual(cli.psf_tag(None, True), "auto_cfallback")

    def test_wave_psf_uses_negative_sidecar_z(self):
        """The PSF phase combines theory and fitted terms at z = -sidecar position."""
        nx_os, ny = 6, 4
        rng = np.random.default_rng(2)
        meta = {"fov": [.22, .22], "psfFovZ": .16, "slicePositions": np.array([-.03, .045])}
        cal = {"delta_y": rng.normal(size=nx_os), "delta_z": rng.normal(size=nx_os)}
        fit = {key: rng.normal(size=nx_os) for key in ("a_processed", "b_processed", "c_processed")}
        psf = cli.wave_psf(np.array([0, 1]), meta, cal, fit, ny)
        self.assertEqual(psf.shape, (2, nx_os, ny))
        y = cli.phase_coordinates_m(ny, .22)
        band, k, line = 1, 3, 0
        z = -.045
        angle = (-2 * np.pi * (cal["delta_y"][k] * y[line] + cal["delta_z"][k] * z)
                 + fit["a_processed"][k] * y[line] / .22 + fit["b_processed"][k] * z / .16
                 + fit["c_processed"][k])
        self.assertAlmostEqual(np.angle(psf[band, k, line] * np.exp(-1j * angle)), 0, places=5)
        np.testing.assert_allclose(abs(psf), 1, atol=1e-6)


class SliceMappingTests(unittest.TestCase):
    """Check the scalar-carrier and RF-waveform slice-group mappings."""

    def test_scalar_carrier_mapping(self):
        """Carrier over slice gradient selects the nearest group center within 3 um."""
        centers = np.array([-.003, 0., .003])
        gradient = np.full(4, 2e5)
        carrier = centers[[2, 0, 1, 2]] * gradient
        actual, report = cli.map_excitations_to_groups(None, None, carrier, gradient, {}, centers)
        np.testing.assert_array_equal(actual, [2, 0, 1, 2])
        self.assertTrue(report["method"].startswith("scalar"))
        with self.assertRaises(ValueError):
            cli.map_excitations_to_groups(None, None, carrier + 2e5 * 1e-5, gradient, {}, centers)

    def test_waveform_translation_sign(self):
        """Envelopes translated with exp(-2i*pi*z*m) pass; the opposite sign is rejected."""
        positions = np.array([-.036, .012, 0.03, -.006])
        report = cli.verify_waveform_slice_translation(translated_excitation_sequence(positions, -1), positions)
        self.assertLess(report["max_envelope_error"], 1e-6)
        self.assertGreater(report["opposite_sign_median_envelope_error"], .1)
        with self.assertRaises(ValueError):
            cli.verify_waveform_slice_translation(translated_excitation_sequence(positions, 1), positions)

    def test_zero_carrier_uses_group_by_slc(self):
        """Without a carrier, the sidecar groupBySlc mapping is used after the waveform check."""
        centers = np.array([-.03, 0., .03])
        group_by_slc = np.array([2, 3, 1])  # one-based, as written by MATLAB
        slc = np.array([0, 1, 2, 1])
        seq = translated_excitation_sequence(centers[group_by_slc[slc] - 1], -1)
        actual, report = cli.map_excitations_to_groups(
            seq, {"SLC": slc}, np.zeros(4), np.full(4, 2e5), {"groupBySlc": group_by_slc}, centers)
        np.testing.assert_array_equal(actual, group_by_slc[slc] - 1)
        self.assertTrue(report["method"].startswith("RF waveform"))


class CalibrationMapTests(unittest.TestCase):
    """Check the parallel ESPIRiT scheduling (the ESPIRiT call itself is replaced)."""

    def test_parallel_workers_fill_only_missing_slices(self):
        """Spawned single-threaded workers fill the missing slices; diagnostics without a map are refused."""
        meta = {"csmAdcIndices": np.zeros((4, 1), int)}
        threads = os.environ.get("OMP_NUM_THREADS")
        with tempfile.TemporaryDirectory() as tmp:
            folders = cli.output_folders(tmp)
            np.save(folders["cache"] / "csm_slice_01.npy", np.zeros(1))
            with mock.patch.object(cli, "estimate_slice_map_worker", fake_slice_worker), \
                    mock.patch.object(cli.LOGGER, "disabled", True):
                cli.estimate_maps({"twix": None, "labels": None}, None, None, meta, None, folders, -1, 2)
            self.assertEqual(os.environ.get("OMP_NUM_THREADS"), threads)
            for p in (0, 2, 3):
                np.testing.assert_array_equal(np.load(folders["cache"] / f"csm_slice_{p:02d}.npy"), [p, p])
                with np.load(folders["arrays"] / f"csm_diagnostics_{p:02d}.npz") as saved:
                    self.assertEqual(str(saved["threads"]), "1")
            np.testing.assert_array_equal(np.load(folders["cache"] / "csm_slice_01.npy"), [0])
            (folders["cache"] / "csm_slice_03.npy").unlink()
            with self.assertRaises(ValueError):
                cli.estimate_maps({"twix": None, "labels": None}, None, None, meta, None, folders, -1, 2)


class SolverTests(unittest.TestCase):
    """Check the exact MB1 solution and the MB2 operator gate."""

    def test_mb1_least_squares_is_exact_on_a_complete_grid(self):
        """A complete Cartesian grid gives the image back exactly on the native rows."""
        ncoil, nx, nx_os, ny = 4, 8, 24, 12
        sens = unit_rss_sensitivities(ncoil, nx, nx_os, ny, seed=0)[None]
        rng = np.random.default_rng(1)
        image = np.zeros((nx_os, ny), np.complex64)
        image[cli.native_rows(nx, nx_os)] = rng.normal(size=(nx, ny)) + 1j * rng.normal(size=(nx, ny))
        phase = torch.as_tensor(np.exp(1j * rng.uniform(0, 2 * np.pi, (1, ny))).astype(np.complex64))
        data = fftc_dim(fftc_dim(sens[0] * torch.as_tensor(image), 1), 2) * phase[0][None, None]
        result, prediction = cli.mb1_least_squares(data, sens, phase)
        np.testing.assert_allclose(result[0].numpy(), image, atol=1e-4)
        self.assertLess(cli.relative_norm(prediction - data, data), 1e-5)

    def test_operator_tests_pass_and_record_once(self):
        """A consistent MB2 encoding passes the adjoint and recovery gates and writes one report."""
        ncoil, nx, nx_os, ny = 4, 8, 24, 16
        sens = torch.stack([unit_rss_sensitivities(ncoil, nx, nx_os, ny, seed=s) for s in (3, 4)])
        rng = np.random.default_rng(5)
        phase = torch.as_tensor(np.exp(1j * rng.uniform(0, 2 * np.pi, (2, ny))).astype(np.complex64))
        psf = torch.as_tensor(np.exp(1j * rng.uniform(-1, 1, (2, nx_os, ny))).astype(np.complex64))
        mask = torch.ones((1, 1, ny))
        with tempfile.TemporaryDirectory() as tmp:
            folders = cli.output_folders(tmp)
            cli.operator_tests(sens, phase, psf, mask, folders, 400, 1e-6)
            report = json.loads((folders["reports"] / "operator_tests.json").read_text())
            self.assertLess(report["adjoint_error"], 1e-5)
            self.assertTrue(report["solver"]["converged"])
            cli.operator_tests(sens, phase, psf, mask, folders, 400, 1e-6)


class L1Tests(unittest.TestCase):
    """Check the L1 settings, the per-scan noise profile and the per-group L1 solve."""

    def test_tags_and_support(self):
        """Tags encode C and non-default FISTA settings; MB1 and indivisible matrices are refused."""
        self.assertEqual(cli.l1_tag(.5, 3, 300, 1e-4), "l1c0p5")
        self.assertEqual(cli.l1_tag(1, 2, 600, 1e-5), "l1c1_lev2_it600_tol1e-05")
        meta = {"Nx": 440, "Ny": 220}
        cli.check_l1_support({"MB": 2}, meta, 2)
        for audit, levels in (({"MB": 1}, 2), ({"MB": 2}, 3)):
            with self.assertRaises(ValueError):
                cli.check_l1_support(audit, meta, levels)

    def test_noise_profile_corrects_the_fit_leverage(self):
        """The corrected profile recovers the FOV noise level that the raw lambda-0 residual underestimates."""
        terms, profile = synthetic_mb2_group(ny=256)
        rows = cli.native_rows(32, 96)
        _, prediction = lambda0_fit(terms)
        cov, whitener, whitening = cli.l1_group_whitening(terms, prediction)
        self.assertLess(whitening["margin_leakage"], 1e-6)
        result = cli.noise_profile_from_fit(terms, prediction, cov, whitener, rows, 500, 1e-8)
        margin = cli.far_margin_rows(96, cli.L1_MARGIN_FRACTION)
        expected = profile[rows].mean() / profile[margin].mean()
        level = whitening["margin_variance_mean"]
        # Two band unknowns per four coil samples remain at Ry2 with 8 coils: the
        # unwhitened fit of correlated-coil noise leaves at least half of it.
        remaining = result["synthetic_power"][rows].mean() / level
        self.assertTrue(.5 <= remaining < .8, remaining)
        self.assertAlmostEqual(result["relative_variance"][rows].mean() / expected, 1, delta=.05)
        self.assertAlmostEqual(result["measured_power"][rows].mean() / level / expected, remaining, delta=.05)
        self.assertAlmostEqual(result["relative_variance"][margin].mean(), 1, delta=.05)

    def test_noise_profile_gates(self):
        """Wave spread into the margins, too few data per unknown and too few outer lines are refused."""
        rng = np.random.default_rng(11)
        spread = torch.as_tensor(np.exp(1j * rng.uniform(-3, 3, (2, 96, 128))).astype(np.complex64))
        terms, _ = synthetic_mb2_group(psf=spread)
        with self.assertRaises(ValueError):
            cli.l1_group_whitening(terms, lambda0_fit(terms)[1])
        rows = cli.native_rows(32, 96)
        terms, _ = synthetic_mb2_group(ncoil=4)
        _, prediction = lambda0_fit(terms)
        cov, whitener, _ = cli.l1_group_whitening(terms, prediction)
        with self.assertRaisesRegex(ValueError, "too few data per unknown|did not converge"):
            cli.noise_profile_from_fit(terms, prediction, cov, whitener, rows, 500, 1e-8)
        terms, _ = synthetic_mb2_group()
        _, prediction = lambda0_fit(terms)
        cov, whitener, _ = cli.l1_group_whitening(terms, prediction)
        terms["lines"] = terms["lines"][abs(terms["lines"] - 64) < 50]
        with self.assertRaisesRegex(ValueError, "outer"):
            cli.noise_profile_from_fit(terms, prediction, cov, whitener, rows, 500, 1e-8)

    def test_l1_group_uses_the_noise_null(self):
        """lambda_b is C * sigma_b * diagonal_b; FISTA stays finite, native-row only and fits less tightly."""
        terms, profile = synthetic_mb2_group()
        rows = cli.native_rows(32, 96)
        start, prediction = lambda0_fit(terms)
        cov, whitener, _ = cli.l1_group_whitening(terms, prediction)
        regions = np.ones((2, 32, 128), bool)
        regions[:, :4] = False
        relative = profile / profile[cli.far_margin_rows(96, cli.L1_MARGIN_FRACTION)].mean()
        null = cli.l1_noise_null(3, terms, cov, whitener, relative, regions, rows, 500, 1e-8)
        self.assertEqual(null["seed"], cli.L1_NULL_SEED + 3)
        self.assertTrue(all(s > 0 for s in null["sigma_lambda0_noise_rms"]))
        self.assertEqual(null, cli.l1_noise_null(3, terms, cov, whitener, relative, regions, rows, 500, 1e-8))
        args = argparse.Namespace(l1_threshold_multiple=.5, l1_wavelet_levels=2, l1_max_iterations=60,
                                  l1_tolerance=1e-4)
        x, report = cli.solve_l1_group(terms, start, whitener, null, rows, args)
        np.testing.assert_allclose(report["lambda"], [.5 * s * d for s, d in zip(
            null["sigma_lambda0_noise_rms"], null["normal_diagonal"])])
        self.assertTrue(torch.isfinite(x).all())
        self.assertTrue(torch.all(x[:, :rows.start] == 0) and torch.all(x[:, rows.stop:] == 0))
        measured = terms["data"] * terms["mask"]
        self.assertGreater(report["relative_data_residual"], cli.relative_norm(prediction - measured, measured))

    def test_reconstruct_l1_writes_once_and_resumes(self):
        """The L1 stage writes every product once, resumes from it and shares noise nulls across C."""
        groups = [synthetic_mb2_group(seed=seed)[0] for seed in (40, 41, 42)]
        rows = cli.native_rows(32, 96)
        meta = {"Nx": 32, "Ny": 128, "Nx_os": 96, "Nslices": 6}
        encoding = {"pairs": np.array([[0, 3], [1, 4], [2, 5]]),
                    "mask": np.stack([t["mask"][0, 0].numpy() for t in groups])}
        baseline = np.zeros((32, 128, 6), np.complex64)
        for terms, pair in zip(groups, encoding["pairs"]):
            baseline[:, :, pair] = lambda0_fit(terms)[0][:, rows].numpy().transpose(1, 2, 0)
        yaps = {("sSliceArray", "asSlice", "0", "dReadoutFOV"): 16.0,
                ("sSliceArray", "asSlice", "0", "dPhaseFOV"): 128.0,
                ("sSliceArray", "asSlice", "0", "dThickness"): 3.0,
                ("sSliceArray", "asSlice", "0", "sNormal", "dTra"): 1.0,
                ("sSliceArray", "asSlice", "0", "sPosition", "dTra"): 0.0}
        scan = {"hdr": {"MeasYaps": yaps, "Meas": {}, "Protocol": {}}}
        audit = {"MB": 2, "voxel_mm": [.5, 1., 3.], "center_sct_mm": [0., 0., 0.]}
        args = argparse.Namespace(l1_threshold_multiple=.5, l1_wavelet_levels=2, l1_max_iterations=40,
                                  l1_tolerance=1e-4, max_iterations=500, tolerance=1e-8, device="cpu",
                                  min_free_gib=0)
        rng = np.random.default_rng(43)
        with tempfile.TemporaryDirectory() as tmp:
            calibration_folders = cli.output_folders(Path(tmp) / "calibration")
            folders = cli.output_folders(Path(tmp) / "tse")
            for p in range(6):
                cli.save_npz(calibration_folders["arrays"] / f"csm_diagnostics_{p:02d}.npz",
                             calibration_rss=rng.uniform(size=(32, 128)).astype(np.float32))
            cli.save_npy(folders["arrays"] / "lambda0_complex_physical_slices.npy", baseline)

            def run(threshold_multiple):
                """Run the L1 stage with group data taken from ``groups``."""
                args.l1_threshold_multiple = threshold_multiple
                with mock.patch.object(cli, "group_terms", side_effect=lambda scan, group, *_: groups[group]) \
                        as loader, mock.patch.object(cli, "l1_noise_null", wraps=cli.l1_noise_null) as null, \
                        mock.patch.object(cli.LOGGER, "disabled", True):
                    image, reports = cli.reconstruct_l1(scan, meta, None, calibration_folders, folders, encoding,
                                                        None, None, audit, args, baseline)
                return image, reports, loader.call_count, null.call_count

            image, reports, loads, nulls = run(.5)
            self.assertEqual((loads, nulls), (4, 3))
            tag = "l1c0p5_lev2_it40"
            for name in ("cache/l1_noise_profile.npz", "reports/l1_noise_profile.json",
                         "figures/l1_noise_profile.png", "reports/l1_noise_null_group_02.json",
                         f"arrays/{tag}_group_02.npy", f"arrays/{tag}_complex_physical_slices.npy",
                         f"nifti/{tag}_magnitude.nii.gz", f"reports/{tag}_completed.json"):
                self.assertTrue((Path(tmp) / "tse" / name).exists(), name)
            profile = json.loads((folders["reports"] / "l1_noise_profile.json").read_text())
            self.assertEqual(profile["group_0based"], 1)
            np.testing.assert_array_equal(image[:, :, [1, 4]],
                                          np.load(folders["arrays"] / f"{tag}_group_01.npy").transpose(1, 2, 0))
            again, _, loads, nulls = run(.5)
            np.testing.assert_array_equal(again, image)
            self.assertEqual((loads, nulls), (0, 0))
            stronger, reports_c1, loads, nulls = run(1.)
            self.assertEqual((loads, nulls), (3, 0))
            self.assertTrue(np.allclose([r["lambda"][0] for r in reports_c1],
                                        [2 * r["lambda"][0] for r in reports]))
            self.assertFalse(np.array_equal(stronger, image))
            report_path = folders["reports"] / f"{tag}_group_00.json"
            report = json.loads(report_path.read_text())
            report["noise_profile_sha256"] = "0" * 64
            report_path.write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, "another noise profile"):
                run(.5)


class InputOutputTests(unittest.TestCase):
    """Check sidecar discovery, folder creation and exclusive writes."""

    def test_resolve_inputs_finds_one_sidecar(self):
        """Sidecars are found in metadata/ next to v151/; a second candidate is ambiguous."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("calib.seq", "calib_geometry_labels.mat", "calib.dat", "tse.dat",
                         "v151/tse.seq", "metadata/tse_PEorder.mat"):
                (root / name).parent.mkdir(parents=True, exist_ok=True)
                (root / name).write_text("x")
            args = argparse.Namespace(calibration_twix=root / "calib.dat", calibration_seq=root / "calib.seq",
                                      calibration_labels=None, tse_twix=root / "tse.dat",
                                      tse_seq=root / "v151/tse.seq", tse_peorder=None)
            inputs = cli.resolve_inputs(args)
            self.assertEqual(inputs["tse"]["peorder"], (root / "metadata/tse_PEorder.mat").resolve())
            self.assertEqual(inputs["calibration"]["labels"], (root / "calib_geometry_labels.mat").resolve())
            (root / "tse_PEorder.mat").write_text("x")
            with self.assertRaises(ValueError):
                cli.resolve_inputs(args)

    def test_output_folders_reject_escaping_symlinks(self):
        """Typed subfolders are created; a subfolder symlink that leaves the root is refused."""
        with tempfile.TemporaryDirectory() as tmp:
            folders = cli.output_folders(Path(tmp) / "out")
            self.assertEqual(sorted(folders), sorted(cli.ARTIFACT_FOLDERS))
            other = Path(tmp) / "elsewhere"
            other.mkdir()
            root = Path(tmp) / "bad"
            root.mkdir()
            os.symlink(other, root / "arrays")
            with self.assertRaises(ValueError):
                cli.output_folders(root)

    def test_exclusive_and_provenance_writes(self):
        """Reports are never overwritten; provenance must match; evidence keeps the first value."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.json"
            cli.ensure_json(path, {"value": np.float64(1.5), "shape": (2, 3)})
            cli.ensure_json(path, {"value": 1.5, "shape": [2, 3]})
            with self.assertRaises(ValueError):
                cli.ensure_json(path, {"value": 2.0, "shape": [2, 3]})
            with self.assertRaises(FileExistsError):
                cli.save_json(path, {})
            evidence = Path(tmp) / "evidence.json"
            cli.save_json_if_absent(evidence, {"error": 1e-7})
            cli.save_json_if_absent(evidence, {"error": 2e-7})
            self.assertEqual(json.loads(evidence.read_text()), {"error": 1e-7})


class ExportTests(unittest.TestCase):
    """Check the NIfTI export conventions on a synthetic TWIX header."""

    def test_export_flips_axes_and_reads_back(self):
        """Magnitude/phase NIfTI hold the PE- and slice-flipped volume with the requested voxel size."""
        yaps = {("sSliceArray", "asSlice", "0", "dReadoutFOV"): 4.0,
                ("sSliceArray", "asSlice", "0", "dPhaseFOV"): 6.0,
                ("sSliceArray", "asSlice", "0", "dThickness"): 3.0,
                ("sSliceArray", "asSlice", "0", "sNormal", "dTra"): 1.0,
                ("sSliceArray", "asSlice", "0", "sPosition", "dTra"): 30.0}
        scan = {"hdr": {"MeasYaps": yaps, "Meas": {}, "Protocol": {}}}
        rng = np.random.default_rng(6)
        image = (rng.normal(size=(8, 6, 4)) + 1j * rng.normal(size=(8, 6, 4))).astype(np.complex64)
        audit = {"voxel_mm": [.5, 1., 3.], "center_sct_mm": [0., 0., 30.]}
        with tempfile.TemporaryDirectory() as tmp:
            folders = cli.output_folders(tmp)
            paths = cli.export_volume(image, scan, audit, folders, "test")
            magnitude = nib.load(paths["magnitude"])
            expected = np.flip(abs(image), axis=(1, 2)).astype(np.float32)
            np.testing.assert_array_equal(np.asarray(magnitude.dataobj), expected)
            np.testing.assert_allclose(magnitude.header.get_zooms(), (.5, 1., 3.), atol=1e-6)
            self.assertEqual(paths, cli.export_volume(image, scan, audit, folders, "test"))


if __name__ == "__main__":
    unittest.main()
