"""Synthetic tests of recon/recon_wave_tse_from_twix_nifti.py helpers."""

import argparse
import json
import os
from pathlib import Path
import tempfile
import unittest

import nibabel as nib
import numpy as np
import pypulseq as pp
import torch

import recon_wave_tse_from_twix_nifti as cli
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
