"""Focused tests for compact LIN routing and bounded record reads."""

from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from utils.compact_calibration_import import index_compact_streams, read_records
from utils.compact_calibration_import import compare_direct_slice_fourier_z


def fake_stream(lin):
    """Return a tiny metadata-only stream with zero scanner counters."""
    lin = np.asarray(lin, dtype=float)
    stream = SimpleNamespace(Lin=lin, NAcq=len(lin), NCol=6, NCha=2)
    for name in ("Par", "Sli", "Set", "Eco", "Rep", "Ave", "Seg"):
        setattr(stream, name, np.zeros(len(lin)))
    return stream


class CompactImportTests(unittest.TestCase):
    """Exercise zero-based Python counters and one-based sidecar boundaries."""

    def test_global_ref_lin_not_trimmed(self):
        """Keep REF global LIN even if the reader's dense shape trims its prefix."""
        scan = {"image": fake_stream([0, 1]), "refscan": fake_stream([4, 2, 3])}
        meta = {"expectedAdcCount": 5, "imageAdcIndices": np.array([1, 2]),
                "refAdcIndices": np.array([3, 4, 5]), "Nx_os": 6}
        names, local = index_compact_streams(scan, meta)
        np.testing.assert_array_equal(local, [0, 1, 1, 2, 0])
        np.testing.assert_array_equal(names, ["image", "image", "refscan", "refscan", "refscan"])

    def test_duplicate_lin_rejected(self):
        """Reject duplicated global records instead of silently averaging them."""
        scan = {"image": fake_stream([0, 1]), "refscan": fake_stream([2, 2, 4])}
        meta = {"expectedAdcCount": 5, "imageAdcIndices": np.array([1, 2]),
                "refAdcIndices": np.array([3, 4, 5]), "Nx_os": 6}
        with self.assertRaises(ValueError):
            index_compact_streams(scan, meta)

    def test_original_record_correction_indices(self):
        """Pass real acquisition indices and preserve arbitrary request order."""
        stream = fake_stream(np.arange(9))
        stream.memPos = np.arange(9) * 100
        calls = []

        def read(mem, **kwargs):
            """Emulate mapVBVD and verify per-record correction alignment."""
            np.testing.assert_array_equal(mem // 100, kwargs["cIxToRaw"])
            calls.append(len(mem))
            return np.broadcast_to(kwargs["cIxToRaw"][None, None],
                                   tuple(kwargs["outSize"])).astype(np.complex64)

        stream.readData = read
        data = read_records(stream, [7, 1, 8, 0, 3], batch_size=2)
        self.assertEqual(data.shape, (6, 5, 2))
        np.testing.assert_array_equal(data[0, :, 0], [7, 1, 8, 0, 3])
        self.assertEqual(calls, [2, 2, 1])

    def test_direct_slice_fourier_direction(self):
        """Recover reversed direct-slice direction from calibration alone."""
        import torch
        from utils.wave_cg_sense_precondition import fftc_dim
        z = (np.arange(80) - 40) * .002
        positions = (np.arange(50) - 24.5) * .003
        profile = np.stack([np.exp(-((z - center) / .025)**2)
                            for center in (-.06, -.02, .02, .06)], axis=-1).astype(np.complex64)
        projection = fftc_dim(torch.from_numpy(profile[None].repeat(6, axis=0)), dim=1).numpy()
        csm_map = np.arange(1, 1601).reshape(50, 32)
        meta = {"Nx_os": 6, "Ny": 220, "slicePositions": positions,
                "psfZCoordinates": z, "csmAdcIndices": csm_map,
                "psfAdcIndices": np.tile(np.arange(1601, 1681), (4, 1)),
                "psfValidLengths": [80, 80, 80, 80],
                "acqTable": [{"kyPhysicalIndex": 95 + index % 32} for index in range(1600)]}

        def read(scan, index, indices):
            """Supply synthetic Fourier or direct-slice coil vectors by record ID."""
            if len(indices) == 80:
                return projection
            physical = (int(indices[0]) - 1) // 32
            vector = np.array([np.interp(-positions[physical], z, profile[:, c].real)
                               for c in range(4)], dtype=np.complex64)
            return np.broadcast_to(vector, (6, 1, 4)).copy()

        with patch("utils.compact_calibration_import.read_compact_records", side_effect=read):
            report = compare_direct_slice_fourier_z(None, None, meta)
        self.assertGreater(report["hypotheses"]["-1"]["median_supported_correlation"], .99999)
        self.assertLess(report["hypotheses"]["1"]["median_supported_correlation"], .7)


if __name__ == "__main__":
    unittest.main()
