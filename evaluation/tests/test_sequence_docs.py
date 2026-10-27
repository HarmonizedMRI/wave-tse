"""Check that the line numbers quoted in docs/sequence.md match the generators."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SEQUENCE_DOC = REPO_ROOT / "docs" / "sequence.md"
GENERATORS = {
    "tse_2d_wave_mb2v_v3.m": REPO_ROOT / "seq" / "tse_2d_wave_mb2v_v3.m",
    "tse_2d_wave_mb1.m": REPO_ROOT / "seq" / "tse_2d_wave_mb1.m",
    "gre_mixed_dim_mb2_calibration.m": REPO_ROOT / "seq" / "gre_mixed_dim_mb2_calibration.m",
}
TABLE_ROW = re.compile(r"^\| (\d+)(?:–(\d+))? \| (.+?) \|")
CODE_NAME = re.compile(r"`([A-Za-z_]\w*)`")
EDIT_REFERENCE = re.compile(r"lines? (\d+)(?:–(\d+))? ((?:`\w+ = [^`]*`\s*)+)")
EDIT_NAME = re.compile(r"`(\w+) = ")


def assignment_lines(path):
    """Map each assigned MATLAB variable name to the lines that assign it.

    Args:
        path: MATLAB source file.

    Returns:
        dict[str, list[int]]: One-based line numbers of ``name = ...``
        statements that start a line, keyed by variable name.
    """
    lines = {}
    for number, text in enumerate(path.read_text().splitlines(), start=1):
        match = re.match(r"\s*(\w+)\s*=(?!=)", text)
        if match:
            lines.setdefault(match.group(1), []).append(number)
    return lines


def documented_parameters(doc_text):
    """Collect the parameter rows of the per-generator tables.

    Args:
        doc_text: Contents of docs/sequence.md.

    Returns:
        list[tuple[str, str, int, int]]: ``(generator, name, first, last)``
        for every variable name in a table row; ``first``/``last`` is the
        quoted line or line range.
    """
    rows = []
    generator = None
    for line in doc_text.splitlines():
        if line.startswith("#"):
            heading = re.match(r"^### `seq/(.+\.m)`$", line)
            generator = heading.group(1) if heading else None
            continue
        match = TABLE_ROW.match(line)
        if generator and match:
            first = int(match.group(1))
            last = int(match.group(2) or first)
            for name in CODE_NAME.findall(match.group(3)):
                rows.append((generator, name, first, last))
    return rows


def documented_edits(doc_text):
    """Collect the ``line N `name = value;``` edits of the reproduction table.

    Args:
        doc_text: Contents of docs/sequence.md.

    Returns:
        list[tuple[str, str, int, int]]: ``(generator, name, first, last)``
        for every edit; the generator follows from the row's export type.
    """
    section = doc_text.split("## Reproducing the reference exports", 1)[1]
    section = section.split("\n## ", 1)[0]
    generator_by_prefix = {
        "| MB2": "tse_2d_wave_mb2v_v3.m",
        "| MB1": "tse_2d_wave_mb1.m",
        "| Calibration": "gre_mixed_dim_mb2_calibration.m",
    }
    edits = []
    for line in section.splitlines():
        generator = next((name for prefix, name in generator_by_prefix.items()
                          if line.startswith(prefix)), None)
        if generator is None:
            continue
        for match in EDIT_REFERENCE.finditer(line):
            first = int(match.group(1))
            last = int(match.group(2) or first)
            for name in EDIT_NAME.findall(match.group(3)):
                edits.append((generator, name, first, last))
    return edits


class SequenceDocLineNumberTest(unittest.TestCase):
    """Every quoted line number must point at the assignment it names."""

    @classmethod
    def setUpClass(cls):
        """Read the document and index the generator assignments once."""
        cls.doc_text = SEQUENCE_DOC.read_text()
        cls.assignments = {name: assignment_lines(path) for name, path in GENERATORS.items()}

    def check_references(self, references):
        """Assert that each reference range contains an assignment of its name.

        Args:
            references: ``(generator, name, first, last)`` tuples.
        """
        self.assertTrue(references)
        for generator, name, first, last in references:
            with self.subTest(generator=generator, name=name, lines=(first, last)):
                found = self.assignments[generator].get(name, [])
                self.assertTrue(any(first <= line <= last for line in found),
                                f"{name} is assigned on lines {found}, not {first}-{last}")

    def test_parameter_tables(self):
        """The parameter tables quote the line of each default."""
        self.check_references(documented_parameters(self.doc_text))

    def test_reproduction_edits(self):
        """The reproduction table quotes the line of each edited parameter."""
        self.check_references(documented_edits(self.doc_text))


if __name__ == "__main__":
    unittest.main()
