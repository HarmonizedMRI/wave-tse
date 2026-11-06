"""Synthetic tests of the calibration appended to a TSE and of MB1 Wave support."""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest import mock

import numpy as np
import pypulseq as pp
import sigpy as sp
import torch

import recon_wave_tse_from_twix_nifti as cli
from utils import appended_calibration as appended
from utils import wave_psf_coefficient_processing as processing
from utils.mb2_cg_sense_wave import mb2_wave_forward


def center_block(n, width):
    """Return the one-based centered index block of ``seq/utils/designFlashCalibration.m``."""
    center = n // 2 + 1
    return np.arange(center - width // 2, center - width // 2 + width)


def appended_sidecar(ny=24, nx=8, ro_os=3, nacs=(12, 12), grid_nz=20, spacing=.003, nslices=10,
                     ncalib_y=12, ncalib_z=16, psf_fov_z=.048, tse_adcs=100):
    """Return a TSE sidecar with an appended calibration in the MATLAB layout.

    Parts 0/1 (y projections), 2 (no-Wave z projection) and 4 (ACS) are
    present; part 3 is absent. Sizes are small synthetic stand-ins.

    Returns:
        dict: Sidecar as ``loadmat(simplify_cells=True)`` returns it.
    """
    ky_calib, ky_center, ky_acs = center_block(ny, ncalib_y), center_block(ny, 1), center_block(ny, nacs[0])
    kz_center, kz_acs = center_block(ncalib_z, 1), center_block(grid_nz, nacs[1])
    areas_y = (np.arange(ny) - ny // 2) / .22
    areas_psf = (np.arange(ncalib_z) - ncalib_z // 2) / psf_fov_z
    areas_acs = (np.arange(grid_nz) - grid_nz // 2) / (grid_nz * spacing)
    parts = [(0, ky_calib, kz_center, areas_psf), (1, ky_calib, kz_center, areas_psf),
             (2, ky_center, np.arange(1, ncalib_z + 1), areas_psf), (4, ky_acs, kz_acs, areas_acs)]
    table, psf, acs = [], np.zeros((4, max(ncalib_y, ncalib_z)), int), np.zeros(nacs, int)
    for part, ky_list, kz_list, kz_areas in parts:
        for i, ky in enumerate(ky_list):
            for j, kz in enumerate(kz_list):
                table.append({"setID": part, "mode": 1, "isRef": 1, "kyPhysicalIndex": int(ky),
                              "kzPhysicalIndex": int(kz), "kyArea": areas_y[ky - 1], "kzArea": kz_areas[kz - 1],
                              "linLocal": i, "parLocal": j, "rawLin0": len(table)})
                if part < 4:
                    psf[part, max(i, j)] = len(table)
                else:
                    acs[i, j] = len(table)
    positions = (np.arange(nslices) - (nslices - 1) / 2) * spacing
    calibration = {
        "storageLayoutVersion": "appended_ref_lin_v1", "expectedAdcCount": len(table), "acqTable": table,
        "psfAdcIndices": psf, "psfValidMask": (psf > 0).astype(float), "psfValidLengths": (psf > 0).sum(1),
        "acsAdcIndices": acs, "acsKyIndices": ky_acs, "acsKzIndices": kz_acs,
        "acsKzAreas": areas_acs[kz_acs - 1], "acsFovZ": grid_nz * spacing, "acsGridSpacingZ": spacing,
        "acsGridNz": grid_nz, "Nacs": np.array(nacs),
        "psfYCoordinates": (np.arange(ncalib_y) - ncalib_y // 2) * .22 / ncalib_y,
        "psfZCoordinates": (np.arange(ncalib_z) - ncalib_z // 2) * psf_fov_z / ncalib_z, "psfFovZ": psf_fov_z,
        "Nx": nx, "Ny": ny, "Nx_os": nx * ro_os, "fov": np.array([.22, .22]), "slicePositions": positions,
        "sliceThickness": spacing, "Ncycles": 20, "gwave_max": 10., "wavePolarity": 1., "tseAdcCount": tse_adcs}
    return {"isAppendCalibration": 1, "calibration": calibration, "SMSFactor": 1}


def fake_refscan(lin):
    """Return a metadata-only refscan stream with zero storage counters."""
    lin = np.asarray(lin, dtype=float)
    stream = SimpleNamespace(Lin=lin, NAcq=len(lin), NCol=24, NCha=2)
    for name in ("Par", "Sli", "Set", "Eco", "Rep", "Ave", "Seg"):
        setattr(stream, name, np.zeros(len(lin)))
    return stream


def scanner_kz_data(profile, z, kz):
    """Return Fourier data of a z profile with the scanner convention exp(+2i*pi*k*z).

    Args:
        profile: (z sample, coil) values on the positions ``z`` [m].
        z: Sample positions [m].
        kz: Encoding moments [1/m].

    Returns:
        numpy.ndarray: complex (kz line, coil).
    """
    return np.exp(2j * np.pi * np.outer(kz, z)) @ profile


class AppendedMetadataTests(unittest.TestCase):
    """Check the sidecar contract and the refscan index."""

    def test_metadata_is_validated(self):
        """A consistent sidecar loads; broken maps, kz moments or a missing calibration are refused."""
        meta = appended.load_appended_metadata(appended_sidecar())
        self.assertEqual(meta["kind"], "appended")
        self.assertEqual(meta["Nslices"], 10)
        broken = appended_sidecar()
        broken["calibration"]["acsAdcIndices"][0, 0] = broken["calibration"]["acsAdcIndices"][0, 1]
        with self.assertRaises(ValueError):
            appended.load_appended_metadata(broken)
        broken = appended_sidecar()
        broken["calibration"]["acsKzAreas"] = -broken["calibration"]["acsKzAreas"]
        with self.assertRaises(ValueError):
            appended.load_appended_metadata(broken)
        with self.assertRaises(ValueError):
            appended.load_appended_metadata({"isAppendCalibration": 0})

    def test_refscan_index(self):
        """Every compact LIN maps to its refscan record; duplicates are refused."""
        meta = appended.load_appended_metadata(appended_sidecar())
        count = meta["expectedAdcCount"]
        order = np.random.default_rng(0).permutation(count)
        names, local = appended.index_appended_stream({"refscan": fake_refscan(order)}, meta)
        np.testing.assert_array_equal(order[local], np.arange(count))
        self.assertTrue(np.all(names == "refscan"))
        order[1] = order[0]
        with self.assertRaises(ValueError):
            appended.index_appended_stream({"refscan": fake_refscan(order)}, meta)


class ZConventionTests(unittest.TestCase):
    """Check that the ACS planes fall on the physical slice centers."""

    def test_full_kz_grid_places_slices_on_planes(self):
        """With every kz line, point sources at the slice centers return exactly at their planes."""
        meta = appended.load_appended_metadata(appended_sidecar(nacs=(12, 20), grid_nz=20))
        grid = appended.acs_z_grid(meta)
        np.testing.assert_allclose(grid["positions"][grid["slice_planes"]], meta["slicePositions"], atol=1e-12)
        weights = np.random.default_rng(1).normal(size=(10, 3)) + 1j
        data = scanner_kz_data(weights, meta["slicePositions"], meta["acsKzAreas"])  # (kz, coil)
        image = sp.ifft(appended.embed_kz(data.T, meta, grid), axes=(-1,))
        np.testing.assert_allclose(image[:, grid["slice_planes"]].T, np.sqrt(grid["n"]) * weights, atol=1e-5)
        mirrored = scanner_kz_data(weights, -meta["slicePositions"], meta["acsKzAreas"])
        image = sp.ifft(appended.embed_kz(mirrored.T, meta, grid), axes=(-1,))
        self.assertGreater(np.max(abs(image[:, grid["slice_planes"]].T - np.sqrt(grid["n"]) * weights)), 1)

    def test_off_grid_slices_are_refused(self):
        """Slice centers between z planes cannot be sampled exactly."""
        sidecar = appended_sidecar()
        sidecar["calibration"]["slicePositions"] = sidecar["calibration"]["slicePositions"] + .001
        meta = appended.load_appended_metadata(sidecar)
        sidecar["calibration"]["slicePositions"][0] += .0005
        with self.assertRaises(ValueError):
            appended.acs_z_grid(appended.load_appended_metadata(sidecar))
        self.assertEqual(len(appended.acs_z_grid(meta)["slice_planes"]), 10)

    def test_z_projection_check_prefers_the_stated_orientation(self):
        """The ACS ky=0 planes agree with the low-passed z projection only in the stated orientation."""
        sidecar = appended_sidecar(nacs=(12, 16), grid_nz=24, nslices=16, ncalib_z=40, psf_fov_z=.06)
        meta = appended.load_appended_metadata(sidecar)
        # Compact profiles inside both encoded z FOVs (60 and 72 mm), so neither aliases.
        z = np.linspace(-.028, .028, 561)
        taper = np.hanning(len(z))
        profile = np.stack([taper * np.exp(-((z - center) / .012) ** 2) * np.exp(1j * 1000 * center * z)
                            for center in (-.02, -.007, .008, .02)], axis=-1)
        table = meta["acqTable"]

        def read(scan, index, indices):
            """Return Fourier data of the profile at the kz of the requested records (center kx only)."""
            kz = np.array([table[i - 1]["kzArea"] for i in indices])
            data = scanner_kz_data(profile, z, kz) * (z[1] - z[0])
            return np.broadcast_to(data[None], (meta["Nx_os"],) + data.shape).copy()

        with mock.patch.object(appended, "read_compact_records", side_effect=read):
            report = appended.compare_acs_with_z_reference(None, None, meta)
        stated = report["hypotheses"]["1"]["median_supported_correlation"]
        reversed_order = report["hypotheses"]["-1"]["median_supported_correlation"]
        self.assertGreater(stated, .99)
        self.assertGreater(stated - reversed_order, .1)


class AcsEspiritTests(unittest.TestCase):
    """Check slice2d ESPIRiT of the 3D ACS against known smooth coil maps."""

    def test_maps_match_the_coils_at_the_slice_planes(self):
        """Maps from a 12 x 12 ACS follow z-dependent coils at the slice centers, not their mirror."""
        meta = appended.load_appended_metadata(appended_sidecar(ny=24, nacs=(12, 12), grid_nz=20, nslices=10))
        grid = appended.acs_z_grid(meta)
        ny, positions = meta["Ny"], grid["positions"]
        y = (np.arange(ny) - ny // 2) * .22 / ny
        yy, zz = np.meshgrid(y, positions, indexing="ij")
        # Coil-specific z phase slopes make the maps differ from their z mirror.
        coils = np.stack([np.exp(-((yy - cy) ** 2 + (zz - cz) ** 2) / .0015) * np.exp(1j * (3 * yy + slope * zz + p))
                          for cy, cz, slope, p in ((-.05, -.02, -60, 0), (.05, -.02, 20, 1), (-.05, .025, 60, 2),
                                                   (.05, .025, -20, 3))])
        coils /= np.sqrt(np.sum(abs(coils) ** 2, axis=0))
        obj = (1 + .3 * np.cos(30 * yy)) * (abs(yy) < .08) * (abs(zz) < .025)
        image = coils * obj  # (coil, y, z plane)
        # ky: centered FFT as for the TSE; kz: the scanner convention at the ACS moments.
        hybrid_y = sp.fft(image, axes=(1,))[:, meta["acsKyIndices"] - 1]
        kz_data = np.einsum("cyz,kz->cyk", hybrid_y, np.exp(2j * np.pi * np.outer(meta["acsKzAreas"], positions)))
        maps, eigen, rss = appended.espirit_planes(kz_data[None].astype(np.complex64), meta, ny, 12, 0,
                                                   sp.Device(-1))
        self.assertEqual(maps.shape, (1, 4, ny, 10))
        truth = coils[:, :, grid["slice_planes"]]
        support = obj[:, grid["slice_planes"]] > 0
        for candidate, expected in ((truth, True), (truth[:, :, ::-1], False)):
            correlation = abs(np.sum(maps[0].conj() * candidate, axis=0))
            if expected:
                self.assertGreater(np.median(correlation[support]), .97)
            else:
                self.assertLess(np.median(correlation[support]), .9)
        np.testing.assert_allclose(np.sum(abs(maps[0]) ** 2, axis=0), 1, atol=1e-4)
        self.assertGreater(rss[0][support].min(), 0)


class MB1WaveTests(unittest.TestCase):
    """Check the y-only Wave PSF processing and the one-band Wave solver."""

    def test_y_only_coefficient_processing(self):
        """The y-only mode fits a and c = c_y as the two-axis mode does, returns b = 0 and keeps the a/c gates."""
        x = np.arange(512, dtype=float)
        w = 2 * np.pi * 4 / len(x)
        a, b, c = np.sin(w * x), .7 * np.cos(w * x + .2), .3 + .5 * np.sin(w * x + .4)
        result = processing.process_wave_yz_coefficients(a, np.full_like(a, np.nan), c, np.full_like(a, 5.),
                                                          128, 384, expected_cycles=4, axes="y")
        both = processing.process_wave_yz_coefficients(a, b, c / 2, c / 2, 128, 384, expected_cycles=4)
        for key in ("a_processed", "c_processed"):
            np.testing.assert_allclose(result[key].numpy(), both[key].numpy(), atol=1e-6)
        self.assertLess(np.max(abs(result["a_processed"].numpy() - a)[128:384]), .02)
        np.testing.assert_array_equal(result["b_processed"].numpy(), 0)
        report = result["report"]
        self.assertEqual(report["wave_axes"], "y")
        self.assertNotIn("b_frequency_within_10_percent_of_sequence", report["period_validation_gates"])
        with self.assertRaises(ValueError):
            processing.process_wave_yz_coefficients(a, a, c, c, 128, 384, expected_cycles=4, axes="z")

    def test_y_only_wave_psf(self):
        """Without z the PSF holds the y terms only and is the same for every band."""
        rng = np.random.default_rng(3)
        nx_os, ny = 6, 4
        meta = {"fov": [.22, .22], "psfFovZ": .16, "slicePositions": np.array([-.03, .045])}
        cal = {"delta_y": rng.normal(size=nx_os)}
        fit = {key: rng.normal(size=nx_os) for key in ("a_processed", "b_processed", "c_processed")}
        psf = cli.wave_psf(np.array([1]), meta, cal, fit, ny, use_z=False)
        y = cli.phase_coordinates_m(ny, .22)
        angle = -2 * np.pi * cal["delta_y"][:, None] * y[None] + fit["a_processed"][:, None] * y[None] / .22 + \
            fit["c_processed"][:, None]
        np.testing.assert_allclose(np.angle(psf[0] * np.exp(-1j * angle)), 0, atol=1e-5)

    def test_one_band_operator_passes_the_operator_tests(self):
        """The Wave operator and its CG solver accept one band."""
        ncoil, nx, nx_os, ny = 4, 8, 24, 16
        rng = np.random.default_rng(6)
        maps = rng.normal(size=(ncoil, nx, ny)) + 1j * rng.normal(size=(ncoil, nx, ny))
        sens = np.zeros((1, ncoil, nx_os, ny), np.complex64)
        sens[0, :, cli.native_rows(nx, nx_os)] = maps / np.sqrt(np.sum(abs(maps) ** 2, axis=0))
        psf = torch.as_tensor(np.exp(1j * rng.uniform(-1, 1, (1, nx_os, ny))).astype(np.complex64))
        with tempfile.TemporaryDirectory() as tmp:
            folders = cli.output_folders(tmp)
            cli.operator_tests(torch.as_tensor(sens), torch.ones((1, ny), dtype=torch.complex64), psf,
                               torch.ones((1, 1, ny)), folders, 400, 1e-6)
            self.assertTrue(json.loads((folders["reports"] / "operator_tests.json").read_text())["solver"]["converged"])

    def test_mb1_mapping_check_deconvolves_the_wave(self):
        """Wave-encoded complete MB1 slices match their own maps only after the PSF deconvolution."""
        # An even slice count: the reversed partner of every checked slice is another slice.
        ncoil, nx, nx_os, ny, nslices = 4, 8, 24, 16, 6
        rng = np.random.default_rng(7)
        x, yy = np.mgrid[:nx, :ny]
        meta = {"Nx": nx, "Ny": ny, "Nx_os": nx_os, "Nslices": nslices, "fov": [.22, .22], "psfFovZ": .16,
                "slicePositions": np.linspace(-.06, .06, nslices)}
        cal = {"delta_y": 30 * np.sin(np.linspace(0, 4 * np.pi, nx_os))}
        fit = {"a_processed": np.zeros(nx_os), "b_processed": np.zeros(nx_os), "c_processed": np.zeros(nx_os)}
        encoding = {"mapping": np.arange(nslices), "pairs": np.arange(nslices)[:, None],
                    "phase": np.ones((nslices, 1, ny), np.complex64)}
        data = {}
        with tempfile.TemporaryDirectory() as tmp:
            calibration_folders, folders = cli.output_folders(Path(tmp) / "cal"), cli.output_folders(Path(tmp) / "tse")
            np.save(calibration_folders["cache"] / "global_coil_basis.npy", np.eye(ncoil))
            for p in range(nslices):
                centers = rng.uniform(0, [nx, ny], (ncoil, 2))
                maps = np.stack([np.exp(-((x - cx) ** 2 + (yy - cy) ** 2) / 20) for cx, cy in centers]) + 0j
                maps /= np.sqrt(np.sum(abs(maps) ** 2, axis=0))
                np.save(calibration_folders["cache"] / f"csm_slice_{p:02d}.npy", maps.astype(np.complex64))
                sens = np.zeros((1, ncoil, nx_os, ny), np.complex64)
                sens[0, :, cli.native_rows(nx, nx_os)] = maps
                image = np.zeros((1, nx_os, ny), np.complex64)
                image[0, cli.native_rows(nx, nx_os)] = 1 + rng.uniform(size=(nx, ny))
                psf = torch.as_tensor(cli.wave_psf(np.array([p]), meta, cal, fit, ny, use_z=False))
                data[p] = mb2_wave_forward(torch.as_tensor(image), torch.as_tensor(sens), torch.ones((1, ny),
                                           dtype=torch.complex64), psf, torch.ones((1, 1, ny))).numpy()
            with mock.patch.object(cli, "load_group_data", side_effect=lambda scan, group, *a: (data[group], group)):
                cli.verify_mb1_calibration_mapping(None, None, meta, calibration_folders, folders, encoding,
                                                   cal, fit, {"wave": True})
            report = json.loads((folders["reports"] / "mb1_flash_mapping_validation.json").read_text())
            self.assertTrue(report["passed"])
            self.assertLess(max(r["stated_mapping_residual"] for r in report["slices"]), 1e-4)


def appended_test_sequence(meta, path, kz_sign=1):
    """Write a gradient-only Pulseq file with two TSE ADCs and the calibration ADCs of ``meta``.

    Every line plays its ky/kz moment before the ADC and rewinds it after, so
    k-space at the ADC equals the sidecar moments (times ``kz_sign`` along z).
    """
    system = pp.Opts(max_grad=80, grad_unit="mT/m", max_slew=200, slew_unit="T/m/s", adc_dead_time=10e-6)
    seq = pp.Sequence(system)
    adc = pp.make_adc(num_samples=meta["Nx_os"], dwell=10e-6, delay=system.adc_dead_time, system=system)

    def line(ky, kz, labels):
        """Add one encoded ADC line with its rewinders."""
        for sign in (1, None, -1):
            if sign is None:
                seq.add_block(adc, *labels)
                continue
            moves = [pp.make_trapezoid(channel, area=sign * area, system=system)
                     for channel, area in (("y", ky), ("z", kz)) if abs(area) > 0]
            if moves:
                seq.add_block(*moves)

    def labels(values):
        """Return SET labels of a value mapping."""
        return [pp.make_label(label=key, type="SET", value=value) for key, value in values.items()]

    for lin in range(int(meta["tseAdcCount"])):
        line(0, 0, labels({"LIN": lin, "SLC": 0, "REF": False}))
    for record in meta["acqTable"]:
        line(record["kyArea"], kz_sign * record["kzArea"],
             labels({"LIN": record["rawLin0"], "SLC": 0, "PAR": 0, "SET": 0, "ECO": 0, "REF": True}))
    seq.write(str(path))


class SequenceContractTests(unittest.TestCase):
    """Check the sequence-level contract of an appended calibration and the pre-scan tool."""

    def test_contract_accepts_the_stated_grid_and_rejects_a_flipped_kz(self):
        """Correct ACS moments pass; a negated kz encoding is refused."""
        meta = appended.load_appended_metadata(appended_sidecar(tse_adcs=2))
        with tempfile.TemporaryDirectory() as tmp:
            good, flipped = Path(tmp) / "good.seq", Path(tmp) / "flipped.seq"
            appended_test_sequence(meta, good)
            appended_test_sequence(meta, flipped, kz_sign=-1)
            result = cli.appended_sequence_contract(good, meta)
            np.testing.assert_allclose(result["delta_y"], 0, atol=1e-9)
            self.assertNotIn("delta_z", result)
            with self.assertRaisesRegex(ValueError, "kz grid"):
                cli.appended_sequence_contract(flipped, meta)

    def test_pre_scan_tool_reports_the_export(self):
        """The evaluation tool passes a consistent export and reports a broken sidecar."""
        sys_path = str(Path(__file__).resolve().parents[2] / "evaluation")
        import sys
        sys.path.insert(0, sys_path)
        try:
            import check_appended_calibration as tool
        finally:
            sys.path.remove(sys_path)
        sidecar = appended_sidecar(tse_adcs=2)
        meta = appended.load_appended_metadata(appended_sidecar(tse_adcs=2))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tse_cal.seq"
            appended_test_sequence(meta, path)
            with mock.patch.object(tool, "loadmat", return_value=sidecar):
                summary = tool.check_export(path, Path(tmp) / "unused_PEorder.mat")
            self.assertTrue(summary["passed"], summary["error"])
            self.assertEqual(summary["tse_adc_count"], 2)
            sidecar["calibration"]["tseAdcCount"] = 3
            with mock.patch.object(tool, "loadmat", return_value=sidecar):
                summary = tool.check_export(path, Path(tmp) / "unused_PEorder.mat")
            self.assertFalse(summary["passed"])


class AppendedCliTests(unittest.TestCase):
    """Check input resolution, the TSE/calibration ADC split and the 3D-ACS map files."""

    def test_resolve_inputs_without_a_separate_calibration(self):
        """Without --calibration-twix the calibration inputs are the TSE files."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("tse.dat", "v151/tse.seq", "metadata/tse_PEorder.mat"):
                (root / name).parent.mkdir(parents=True, exist_ok=True)
                (root / name).write_text("x")
            args = argparse.Namespace(calibration_twix=None, calibration_seq=None, calibration_labels=None,
                                      tse_twix=root / "tse.dat", tse_seq=root / "v151/tse.seq", tse_peorder=None)
            inputs = cli.resolve_inputs(args)
            self.assertEqual(inputs["calibration_mode"], "appended")
            self.assertEqual(inputs["calibration"]["labels"], inputs["tse"]["peorder"])
            self.assertEqual(inputs["calibration"]["twix"], inputs["tse"]["twix"])

    def test_tse_adc_count(self):
        """TSE ADCs end where the stated REF calibration readouts begin."""
        labels = {"LIN": np.arange(7), "REF": np.array([0, 0, 0, 0, 1, 1, 1])}
        pe = {"isAppendCalibration": 1, "calibration": {"tseAdcCount": 4, "expectedAdcCount": 3}}
        self.assertEqual(cli.tse_adc_count(labels, pe), 4)
        with self.assertRaises(ValueError):
            cli.tse_adc_count(labels, {"isAppendCalibration": 1,
                                       "calibration": {"tseAdcCount": 4, "expectedAdcCount": 2}})
        with self.assertRaises(ValueError):
            cli.tse_adc_count(labels, {})
        self.assertEqual(cli.tse_adc_count({"LIN": np.arange(3), "REF": np.zeros(3)}, {}), 3)
        self.assertEqual(cli.wave_axes({"wave": True}), (True, True))
        self.assertEqual(cli.wave_axes({"wave": True, "wave_y": True, "wave_z": False}), (True, False))

    def test_acs_maps_are_written_per_slice(self):
        """The 3D-ACS path writes unit-RSS maps and diagnostics of every slice and refuses partial folders."""
        meta = appended.load_appended_metadata(appended_sidecar())
        rng = np.random.default_rng(8)
        hybrid = (rng.normal(size=(meta["Nx"], 3, 12, 12)) + 1j * rng.normal(size=(meta["Nx"], 3, 12, 12)))
        with tempfile.TemporaryDirectory() as tmp:
            folders = cli.output_folders(tmp)
            with mock.patch.object(cli, "read_acs", return_value=hybrid.astype(np.complex64)), \
                    mock.patch.object(cli, "ESPIRIT_CALIB_WIDTH", 12), mock.patch.object(cli.LOGGER, "disabled", True):
                cli.estimate_maps({}, None, None, meta, None, folders, -1, 1)
                for p in range(meta["Nslices"]):
                    maps = np.load(folders["cache"] / f"csm_slice_{p:02d}.npy")
                    self.assertEqual(maps.shape, (3, meta["Nx"], meta["Ny"]))
                    np.testing.assert_allclose(np.sum(abs(maps) ** 2, axis=0), 1, atol=1e-4)
                    with np.load(folders["arrays"] / f"csm_diagnostics_{p:02d}.npz") as saved:
                        self.assertEqual(saved["calibration_rss"].shape, (meta["Nx"], meta["Ny"]))
                (folders["cache"] / "csm_slice_03.npy").unlink()
                with self.assertRaises(ValueError):
                    cli.estimate_maps({}, None, None, meta, None, folders, -1, 1)


if __name__ == "__main__":
    unittest.main()
