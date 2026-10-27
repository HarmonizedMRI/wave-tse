"""Tests for evaluation/compare_reconstructions.py on synthetic volumes."""

import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

import nibabel as nib
import numpy as np


EVALUATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EVALUATION_ROOT))

import compare_reconstructions as cr  # noqa: E402


def run_cli(argv):
    """Run the command line quietly and return its exit status."""
    with contextlib.redirect_stdout(io.StringIO()):
        return cr.main(argv)


class CompareReconstructionsTests(unittest.TestCase):
    """Volume and NIfTI comparisons pass on identical data and fail on changes."""

    def setUp(self):
        """Create a temporary folder and one random complex volume."""
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        rng = np.random.default_rng(0)
        self.volume = (rng.normal(size=(6, 4, 3)) + 1j * rng.normal(size=(6, 4, 3))).astype(np.complex64)

    def tearDown(self):
        """Remove the temporary folder."""
        self.tmp.cleanup()

    def write_nifti(self, name, data, shift=0.):
        """Write ``data`` with a 0.5 x 1 x 3 mm affine, optionally translated."""
        affine = np.diag([.5, 1., 3., 1.])
        affine[2, 3] = shift
        path = self.root / name
        nib.save(nib.Nifti1Image(data.astype(np.float32), affine), path)
        return path

    def test_identical_volumes_pass_and_changes_fail(self):
        """Identical arrays pass; a 1% change or a shape change fails with metrics."""
        same = cr.compare_volumes(self.volume, self.volume.copy(), 1e-5)
        self.assertTrue(same["passed"])
        self.assertEqual(same["relative_l2_difference"], 0)
        changed = self.volume.copy()
        changed[..., 2] *= 1.01
        result = cr.compare_volumes(changed, self.volume, 1e-5)
        self.assertFalse(result["passed"])
        self.assertEqual(result["worst_slice"], 2)
        self.assertFalse(cr.compare_volumes(self.volume[..., :2], self.volume, 1e-5)["passed"])

    def test_nifti_affine_difference_fails(self):
        """Equal data with a translated affine fail the NIfTI comparison."""
        magnitude = abs(self.volume)
        a = self.write_nifti("a.nii.gz", magnitude)
        b = self.write_nifti("b.nii.gz", magnitude)
        c = self.write_nifti("c.nii.gz", magnitude, shift=1.)
        self.assertTrue(cr.compare_nifti(a, b, 1e-5, 1e-4)["passed"])
        result = cr.compare_nifti(a, c, 1e-5, 1e-4)
        self.assertFalse(result["passed"])
        self.assertAlmostEqual(result["affine_max_abs_difference_mm"], 1.)

    def test_command_line_exit_status_and_report(self):
        """The command returns 0 on a match, 1 on a mismatch, and writes a new report."""
        np.save(self.root / "new.npy", self.volume)
        np.save(self.root / "old.npy", self.volume)
        np.save(self.root / "bad.npy", 2 * self.volume)
        a = self.write_nifti("a.nii.gz", abs(self.volume))
        report = self.root / "report.json"
        status = run_cli(["--candidate", str(self.root / "new.npy"), "--reference", str(self.root / "old.npy"),
                          "--nifti", str(a), str(a), "--report", str(report)])
        self.assertEqual(status, 0)
        self.assertTrue(json.loads(report.read_text())["passed"])
        self.assertEqual(run_cli(["--candidate", str(self.root / "bad.npy"),
                                  "--reference", str(self.root / "old.npy")]), 1)


if __name__ == "__main__":
    unittest.main()
