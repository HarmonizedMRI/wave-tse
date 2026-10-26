"""Tests for evaluation/compare_sequence_exports.py on synthetic exports."""

import contextlib
import hashlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import scipy.io as sio


EVALUATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EVALUATION_ROOT))

import compare_sequence_exports as cse  # noqa: E402


# Optional scanned v3 export for the signature test, for example
# 20261007_tse_wave_TE15p28_50sl_skyra/v151/
# tse_m2b7v1_B15_cF2_N20_TRA_440x220x50_F220_R1_os3_E17_ES15p29_TE107_TR15000_skyra.seq;
# the test is skipped unless WAVE_TSE_REFERENCE_SEQ names such a file.
REFERENCE_SEQ = Path(os.environ.get("WAVE_TSE_REFERENCE_SEQ", ""))

DEFINITIONS = [("Name", "tse_test"), ("OutputPath", "/machine/a/out/"), ("TR", "15")]
BLOCKS = ["    1   100   1   0   0   0   0   0", "    2 1e-05   0   0   0   0   1   0"]


def write_signed_seq(path, definitions=DEFINITIONS, blocks=BLOCKS):
    """Write a minimal Pulseq-format file with a valid md5 signature.

    Args:
        path: Output ``.seq`` path; parent folders are created.
        definitions: ``(key, value)`` pairs of the ``[DEFINITIONS]`` section.
        blocks: Text rows of the ``[BLOCKS]`` section.

    Returns:
        Path: ``path``.
    """
    lines = ["# Pulseq sequence file", "# Created by MATLAB mr toolbox", "",
             "[VERSION]", "major 1", "minor 5", "revision 1", "",
             "[DEFINITIONS]", *[f"{key} {value} " for key, value in definitions], "",
             "# Format of blocks:", "[BLOCKS]", *blocks, "",
             "[SHAPES]", "", "shape_id 1", "num_samples 2", "1", "0"]
    text = "\n".join(lines) + "\n"
    digest = hashlib.md5(text.encode("ascii")).hexdigest()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + f"\n[SIGNATURE]\n# synthetic test signature\nType md5\nHash {digest}\n")
    return path


def write_sidecar(path, values=(1.0, 2.0), pulse_mat="/machine/a/rf.mat",
                  output_path="/machine/a/out", dtype=np.float64):
    """Write a small MAT sidecar with nested struct, logical and path fields.

    Args:
        path: Output ``.mat`` path; parent folders are created.
        values: Numeric values stored in ``meta.values``.
        pulse_mat: Value of the path field ``mbPulseMatFile``.
        output_path: Value of the nested path field ``meta.OutputPath``.
        dtype: NumPy dtype of ``meta.values``.

    Returns:
        Path: ``path``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    sio.savemat(path, {
        "PEorder": np.arange(6, dtype=np.float64).reshape(2, 3),
        "mbPulseMatFile": pulse_mat,
        "flags": np.array([True, False]),
        "meta": {"OutputPath": output_path, "values": np.asarray(values, dtype=dtype)},
    })
    return path


def compare_pair(candidate, reference, rtol=0.0, atol=0.0):
    """Run ``compare_export`` with the default ignore sets.

    Args:
        candidate: Candidate ``.seq`` path.
        reference: Reference ``.seq`` path.
        rtol: Relative tolerance.
        atol: Absolute tolerance.

    Returns:
        dict: Comparison result.
    """
    return cse.compare_export(candidate, reference, set(cse.DEFAULT_IGNORED_DEFINITIONS),
                              set(cse.DEFAULT_IGNORED_FIELDS), rtol, atol)


class SequenceComparisonTests(unittest.TestCase):
    """Verify signature checks and section comparison of .seq files."""

    def setUp(self):
        """Create a temporary folder for each test."""
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        """Remove the temporary folder."""
        self.tmp.cleanup()

    def test_identical_exports_match_with_valid_signatures(self):
        """Byte-identical files match as identical text."""
        candidate = write_signed_seq(self.root / "cand" / "a.seq")
        reference = write_signed_seq(self.root / "ref" / "a.seq")
        result = compare_pair(candidate, reference)
        self.assertTrue(result["match"])
        self.assertTrue(result["sequence"]["text_identical"])
        self.assertTrue(result["sequence"]["candidate_signature"]["valid"])

    def test_output_path_is_an_expected_definition_difference(self):
        """OutputPath and the resulting signature change do not fail the match."""
        definitions = [("Name", "tse_test"), ("OutputPath", "/machine/b/out/"), ("TR", "15")]
        candidate = write_signed_seq(self.root / "cand" / "a.seq", definitions=definitions)
        reference = write_signed_seq(self.root / "ref" / "a.seq")
        result = compare_pair(candidate, reference)
        self.assertTrue(result["match"])
        self.assertEqual(result["sequence"]["expected_differences"]["OutputPath"],
                         {"candidate": "/machine/b/out/", "reference": "/machine/a/out/"})

    def test_block_change_is_reported_with_file_line_numbers(self):
        """A changed block duration fails and reports both file line numbers."""
        blocks = ["    1   200   1   0   0   0   0   0", BLOCKS[1]]
        candidate = write_signed_seq(self.root / "cand" / "a.seq", blocks=blocks)
        reference = write_signed_seq(self.root / "ref" / "a.seq")
        result = compare_pair(candidate, reference)
        self.assertFalse(result["match"])
        (section,) = result["sequence"]["sections"]
        self.assertEqual(section["section"], "BLOCKS")
        self.assertEqual(section["numeric_differences"], 1)
        self.assertEqual(section["max_abs_numeric_deviation"], 100.0)
        reference_lines = reference.read_text().splitlines()
        line_number = section["first_differences"][0]["reference_line"]
        self.assertEqual(reference_lines[line_number - 1], BLOCKS[0])

    def test_formatting_only_difference_is_numerically_equal(self):
        """Equal numbers written differently match without identical text."""
        blocks = [BLOCKS[0], "    2 0.00001   0   0   0   0   1   0"]
        candidate = write_signed_seq(self.root / "cand" / "a.seq", blocks=blocks)
        reference = write_signed_seq(self.root / "ref" / "a.seq")
        result = compare_pair(candidate, reference)
        self.assertTrue(result["match"])
        self.assertFalse(result["sequence"]["text_identical"])
        self.assertEqual(result["sequence"]["sections"][0]["status"], "numerically equal")

    def test_edited_content_breaks_the_signature(self):
        """Editing a signed file invalidates its md5 signature and the match."""
        candidate = write_signed_seq(self.root / "cand" / "a.seq")
        reference = write_signed_seq(self.root / "ref" / "a.seq")
        candidate.write_text(candidate.read_text().replace("major 1", "major 1 "))
        result = compare_pair(candidate, reference)
        self.assertFalse(result["sequence"]["candidate_signature"]["valid"])
        self.assertFalse(result["match"])

    @unittest.skipUnless(REFERENCE_SEQ.is_file(), "set WAVE_TSE_REFERENCE_SEQ to a scanned v3 export")
    def test_real_reference_signature_is_valid(self):
        """The md5 rule reproduces the stored hash of a scanned v3 export."""
        signature = cse.read_sequence_file(REFERENCE_SEQ)["signature"]
        self.assertTrue(signature["valid"])


class SidecarComparisonTests(unittest.TestCase):
    """Verify MAT-file and text sidecar comparison."""

    def setUp(self):
        """Create a temporary folder for each test."""
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        """Remove the temporary folder."""
        self.tmp.cleanup()

    def compare_mat(self, candidate, reference, rtol=0.0, atol=0.0):
        """Compare two MAT sidecars with the default ignored fields.

        Args:
            candidate: Candidate MAT path.
            reference: Reference MAT path.
            rtol: Relative tolerance.
            atol: Absolute tolerance.

        Returns:
            dict: Result of ``compare_mat_files``.
        """
        return cse.compare_mat_files(candidate, reference, set(cse.DEFAULT_IGNORED_FIELDS),
                                     rtol, atol)

    def test_path_fields_are_expected_differences(self):
        """Machine paths at the top level and inside structs do not fail."""
        candidate = write_sidecar(self.root / "c.mat", pulse_mat="/machine/b/rf.mat",
                                  output_path="/machine/b/out")
        reference = write_sidecar(self.root / "r.mat")
        result = self.compare_mat(candidate, reference)
        self.assertTrue(result["match"])
        self.assertEqual(set(result["expected_differences"]),
                         {"mbPulseMatFile", "meta.OutputPath"})

    def test_numeric_change_is_reported_by_field_path(self):
        """A changed nested value fails with its dotted path and deviation."""
        candidate = write_sidecar(self.root / "c.mat", values=(1.0, 2.5))
        reference = write_sidecar(self.root / "r.mat")
        result = self.compare_mat(candidate, reference)
        self.assertFalse(result["match"])
        self.assertEqual(result["differences"][0]["path"], "meta.values")
        self.assertEqual(result["differences"][0]["max_abs_deviation"], 0.5)
        self.assertTrue(self.compare_mat(candidate, reference, atol=0.6)["match"])

    def test_dtype_change_is_reported(self):
        """single versus double storage of the same values is a difference."""
        candidate = write_sidecar(self.root / "c.mat", dtype=np.float32)
        reference = write_sidecar(self.root / "r.mat")
        result = self.compare_mat(candidate, reference)
        self.assertFalse(result["match"])
        self.assertIn("dtype", result["differences"][0]["issue"])

    def test_text_report_keys(self):
        """Ignored report keys are expected; other changed lines fail."""
        reference = self.root / "r.txt"
        reference.write_text("Name: a\nOutputPath: /machine/a\nADCCount: 1904\n")
        candidate = self.root / "c.txt"
        candidate.write_text("Name: a\nOutputPath: /machine/b\nADCCount: 1904\n")
        ignored = set(cse.DEFAULT_IGNORED_FIELDS)
        self.assertTrue(cse.compare_text_reports(candidate, reference, ignored)["match"])
        candidate.write_text("Name: a\nOutputPath: /machine/b\nADCCount: 1600\n")
        result = cse.compare_text_reports(candidate, reference, ignored)
        self.assertFalse(result["match"])
        self.assertEqual(result["differences"][0]["line"], 3)


class FolderModeTests(unittest.TestCase):
    """Verify pairing by file name, sidecar discovery and the CLI exit status."""

    def run_cli(self, argv):
        """Run the command-line entry point without printing its report.

        Args:
            argv: Argument list without the program name.

        Returns:
            int: Exit status of ``main``.
        """
        with contextlib.redirect_stdout(io.StringIO()):
            return cse.main(argv)

    def setUp(self):
        """Create candidate and reference folders in the generator layouts."""
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.candidate_dir = root / "regenerated"
        self.reference_dir = root / "reference"
        # MB2 layout: v151/<name>.seq with metadata/<name>_PEorder.mat.
        for folder in (self.candidate_dir, self.reference_dir):
            write_signed_seq(folder / "v151" / "tse_mb2.seq")
            write_sidecar(folder / "metadata" / "tse_mb2_PEorder.mat")
        # Calibration layout: the .seq and its geometry labels side by side.
        for folder, output in ((self.candidate_dir, "/machine/b"), (self.reference_dir, "/m/a")):
            write_signed_seq(folder / "calib.seq")
            write_sidecar(folder / "calib_geometry_labels.mat", output_path=output)
            (folder / "calib_geometry_labels.txt").write_text("Name: calib\n")

    def tearDown(self):
        """Remove the temporary folders."""
        self.tmp.cleanup()

    def test_cli_pairs_exports_and_writes_report(self):
        """Both exports and all sidecars are found; the CLI returns 0."""
        report = self.candidate_dir / "report.json"
        status = self.run_cli(["--candidate-dir", str(self.candidate_dir),
                           "--reference-dir", str(self.reference_dir),
                           "--report", str(report)])
        self.assertEqual(status, 0)
        results = json.loads(report.read_text())
        self.assertEqual([result["name"] for result in results], ["calib.seq", "tse_mb2.seq"])
        self.assertEqual([s["suffix"] for s in results[0]["sidecars"]],
                         ["_geometry_labels.mat", "_geometry_labels.txt"])
        self.assertEqual([s["suffix"] for s in results[1]["sidecars"]], ["_PEorder.mat"])

    def test_missing_candidate_sidecar_fails(self):
        """A sidecar present only in the reference fails the pair."""
        (self.candidate_dir / "metadata" / "tse_mb2_PEorder.mat").unlink()
        status = self.run_cli(["--candidate-dir", str(self.candidate_dir),
                           "--reference-dir", str(self.reference_dir)])
        self.assertEqual(status, 1)

    def test_missing_reference_raises(self):
        """A candidate without a same-named reference is an error."""
        write_signed_seq(self.candidate_dir / "v151" / "extra.seq")
        with self.assertRaises(FileNotFoundError):
            self.run_cli(["--candidate-dir", str(self.candidate_dir),
                      "--reference-dir", str(self.reference_dir)])


if __name__ == "__main__":
    unittest.main()
