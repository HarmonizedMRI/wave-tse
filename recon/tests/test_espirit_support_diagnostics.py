"""Small CPU regression tests for optional ESPIRiT support diagnostics."""

import contextlib
import io
import unittest

import numpy as np
import sigpy as sp
import torch

from utils.coil_sensitivity import estimate_espirit_csm_2d_slicewise


class EspiritSupportDiagnosticsTests(unittest.TestCase):
    """Check backward compatibility and crop-independent calibration evidence."""

    def test_optional_diagnostics_preserve_maps(self):
        """Compare default/diagnostic outputs and the eigenvalue crop relation."""
        rng = np.random.default_rng(17)
        shape = (32, 16, 1, 4)
        data = torch.from_numpy((rng.normal(size=shape) +
                                 1j * rng.normal(size=shape)).astype(np.complex64))
        args = (data, np.eye(4, dtype=np.complex64), (32, 32), 8, 1, sp.Device(-1))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            original = estimate_espirit_csm_2d_slicewise(*args, calib_width=12, crop=.95)
            cropped, diag = estimate_espirit_csm_2d_slicewise(
                *args, calib_width=12, crop=.95, return_diagnostics=True)
            full, full_diag = estimate_espirit_csm_2d_slicewise(
                *args, calib_width=12, crop=0, return_diagnostics=True)
        self.assertEqual(original.shape, (4, 32, 32, 1))
        np.testing.assert_allclose(original, cropped)
        for name in ("eigenvalues", "calibration_rss"):
            self.assertEqual(diag[name].shape, (32, 32, 1))
            self.assertTrue(np.isfinite(diag[name]).all())
            np.testing.assert_allclose(diag[name], full_diag[name])
        np.testing.assert_allclose(cropped, full * (full_diag["eigenvalues"] > .95)[None], atol=1e-5)
        self.assertGreater(float(np.sum(diag["calibration_rss"]**2)), 0.)


if __name__ == "__main__":
    unittest.main()
