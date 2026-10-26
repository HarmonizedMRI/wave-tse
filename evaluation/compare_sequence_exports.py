"""Compare regenerated Pulseq exports and their sidecars with reference exports.

The comparison is exact by default. For every ``.seq`` pair it verifies both
md5 signatures and compares all sections line by line after removing the
definitions that depend on the generating machine (``OutputPath`` by default)
and the ``[SIGNATURE]`` section, whose hash also covers those definitions.
When a section differs in text, a token-wise numeric comparison reports
whether the difference is formatting only or a numeric change. Matching
sidecars (``*_PEorder.mat``, ``*_geometry_labels.mat`` and
``*_geometry_labels.txt``) are compared variable by variable; fields that hold
paths of the generating machine are reported as expected differences.

Run from the repository root, for example::

    python evaluation/compare_sequence_exports.py \
        --candidate-dir /path/to/regenerated/out_path \
        --reference-dir /path/to/reference/protocol_folder \
        --report /path/to/comparison_report.json

The exit status is 0 when every pair matches apart from expected differences
and 1 otherwise.
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio


def main(argv=None):
    """Compare every requested export pair and print the result.

    Args:
        argv: Command-line arguments without the program name; ``None`` reads
            ``sys.argv[1:]``.

    Returns:
        int: 0 when all pairs match apart from expected differences, else 1.

    Raises:
        FileNotFoundError: If an input is missing or a candidate sequence has
            no uniquely named reference.

    Side effects:
        Prints a report and, with ``--report``, writes it as JSON.
    """
    args = parse_args(argv)
    pairs = collect_sequence_pairs(args)
    ignored_definitions = set(DEFAULT_IGNORED_DEFINITIONS) | set(args.ignore_definition)
    ignored_fields = set(DEFAULT_IGNORED_FIELDS) | set(args.ignore_field)

    results = []
    for candidate_seq, reference_seq in pairs:
        result = compare_export(candidate_seq, reference_seq, ignored_definitions,
                                ignored_fields, args.rtol, args.atol)
        print_export_result(result)
        results.append(result)

    n_match = sum(result["match"] for result in results)
    print(f"\n{n_match}/{len(results)} export pairs match apart from expected differences.")
    if args.report is not None:
        args.report.write_text(json.dumps(json_ready(results), indent=2) + "\n")
        print(f"Report written to {args.report}")
    return 0 if n_match == len(results) else 1


# [DEFINITIONS] keys whose values depend on the generating machine.
DEFAULT_IGNORED_DEFINITIONS = ("OutputPath",)

# Sidecar fields (matched at any struct depth) and "Key: value" report keys
# that record paths of the generating machine.
DEFAULT_IGNORED_FIELDS = ("OutputPath", "MBGeometryMatFile", "mbPulseMatFile", "out_path")

# Sidecars written by seq/tse_2d_wave_mb2v_v3.m, seq/tse_2d_wave_mb1.m and
# seq/gre_mixed_dim_mb2_calibration.m.
SIDECAR_SUFFIXES = ("_PEorder.mat", "_geometry_labels.mat", "_geometry_labels.txt")

MAX_REPORTED_DIFFERENCES = 10
SECTION_HEADER = re.compile(r"\[([^\]]+)\]")


def parse_args(argv):
    """Parse the command line.

    Args:
        argv: Argument list without the program name, or ``None``.

    Returns:
        argparse.Namespace: Options with ``pathlib.Path`` values for paths and
            float tolerances ``rtol``/``atol``.
    """
    parser = argparse.ArgumentParser(
        description="Compare regenerated Pulseq exports and sidecars with reference exports.")
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--candidate", type=Path, help="regenerated .seq file")
    inputs.add_argument("--candidate-dir", type=Path,
                        help="folder searched recursively for regenerated .seq files")
    parser.add_argument("--reference", type=Path, help="reference .seq file (with --candidate)")
    parser.add_argument("--reference-dir", type=Path,
                        help="folder searched recursively for reference .seq files "
                             "with the candidates' file names (with --candidate-dir)")
    parser.add_argument("--ignore-definition", action="append", default=[], metavar="NAME",
                        help="additional [DEFINITIONS] key treated as an expected difference")
    parser.add_argument("--ignore-field", action="append", default=[], metavar="NAME",
                        help="additional sidecar field or report key treated as an "
                             "expected difference")
    parser.add_argument("--rtol", type=float, default=0.0,
                        help="relative tolerance for numeric values (default 0: exact)")
    parser.add_argument("--atol", type=float, default=0.0,
                        help="absolute tolerance for numeric values (default 0: exact)")
    parser.add_argument("--report", type=Path, help="optional JSON report path")
    args = parser.parse_args(argv)
    if args.candidate is not None and args.reference is None:
        parser.error("--candidate requires --reference")
    if args.candidate_dir is not None and args.reference_dir is None:
        parser.error("--candidate-dir requires --reference-dir")
    return args


def collect_sequence_pairs(args):
    """Pair regenerated and reference sequence files.

    Args:
        args: Parsed options from ``parse_args``.

    Returns:
        list[tuple[Path, Path]]: ``(candidate, reference)`` pairs, sorted by
            candidate path in folder mode.

    Raises:
        FileNotFoundError: If an input is missing, the candidate folder holds
            no ``.seq`` file, or a candidate name does not occur exactly once
            below the reference folder.
    """
    if args.candidate is not None:
        for path in (args.candidate, args.reference):
            if not path.is_file():
                raise FileNotFoundError(f"Sequence file not found: {path}")
        return [(args.candidate, args.reference)]

    for folder in (args.candidate_dir, args.reference_dir):
        if not folder.is_dir():
            raise FileNotFoundError(f"Folder not found: {folder}")
    candidates = sorted(args.candidate_dir.rglob("*.seq"))
    if not candidates:
        raise FileNotFoundError(f"No .seq file below {args.candidate_dir}")
    references = {}
    for path in args.reference_dir.rglob("*.seq"):
        references.setdefault(path.name, []).append(path)

    pairs = []
    for candidate in candidates:
        matches = references.get(candidate.name, [])
        if len(matches) != 1:
            raise FileNotFoundError(
                f"Expected exactly one reference named {candidate.name} below "
                f"{args.reference_dir}; found {len(matches)}")
        pairs.append((candidate, matches[0]))
    return pairs


def compare_export(candidate_seq, reference_seq, ignored_definitions, ignored_fields, rtol, atol):
    """Compare one sequence pair and all sidecars found for it.

    Args:
        candidate_seq: Regenerated ``.seq`` path.
        reference_seq: Reference ``.seq`` path.
        ignored_definitions: ``[DEFINITIONS]`` keys that may differ.
        ignored_fields: Sidecar fields and report keys that may differ.
        rtol: Relative tolerance for numeric values.
        atol: Absolute tolerance for numeric values.

    Returns:
        dict: ``sequence`` (see ``compare_sequences``), ``sidecars`` (one
            entry per sidecar suffix present on either side) and ``match``.
    """
    sequence = compare_sequences(candidate_seq, reference_seq, ignored_definitions, rtol, atol)
    candidate_sidecars = find_sidecars(candidate_seq)
    reference_sidecars = find_sidecars(reference_seq)

    sidecars = []
    for suffix in SIDECAR_SUFFIXES:
        candidate = candidate_sidecars.get(suffix)
        reference = reference_sidecars.get(suffix)
        if candidate is None and reference is None:
            continue
        if candidate is None or reference is None:
            sidecars.append({
                "suffix": suffix, "match": False,
                "issue": "missing candidate sidecar" if candidate is None
                         else "missing reference sidecar",
                "candidate": str(candidate) if candidate else None,
                "reference": str(reference) if reference else None,
            })
            continue
        if suffix.endswith(".mat"):
            result = compare_mat_files(candidate, reference, ignored_fields, rtol, atol)
        else:
            result = compare_text_reports(candidate, reference, ignored_fields)
        result["suffix"] = suffix
        sidecars.append(result)

    return {
        "name": candidate_seq.name,
        "candidate": str(candidate_seq),
        "reference": str(reference_seq),
        "sequence": sequence,
        "sidecars": sidecars,
        "match": sequence["match"] and all(sidecar["match"] for sidecar in sidecars),
    }


def find_sidecars(seq_path):
    """Locate the sidecars that a generator wrote for a sequence export.

    The MB2 generator writes ``<out_path>/v151/<name>.seq`` and
    ``<out_path>/metadata/<name>_PEorder.mat``, the MB1 generator writes the
    sidecar to ``<out_path>/<name>_PEorder.mat``, and the calibration writes
    ``<name>.seq`` beside ``<name>_geometry_labels.{mat,txt}``. Reference
    folders may also keep sidecars in a sibling ``metadata/`` folder.

    Args:
        seq_path: ``.seq`` file path.

    Returns:
        dict[str, Path]: Existing sidecar paths keyed by suffix.
    """
    folders = [seq_path.parent, seq_path.parent / "metadata",
               seq_path.parent.parent, seq_path.parent.parent / "metadata"]
    found = {}
    for suffix in SIDECAR_SUFFIXES:
        for folder in folders:
            path = folder / f"{seq_path.stem}{suffix}"
            if path.is_file():
                found[suffix] = path
                break
    return found


def read_sequence_file(path):
    """Split a Pulseq file into sections and verify its md5 signature.

    Pulseq hashes the file bytes written before the newline that precedes
    ``[SIGNATURE]``; that newline belongs to the signature block.

    Args:
        path: ``.seq`` file path.

    Returns:
        dict: ``sections`` maps each section name (``""`` for the leading
            comment block) to ``(line_number, text)`` tuples in file order,
            without newline characters. ``signature`` holds the stored type
            and hash, the recomputed md5 and ``valid``.
    """
    data = path.read_bytes()
    signature = {"stored_type": None, "stored_hash": None, "computed_md5": None, "valid": False}
    signature_start = data.rfind(b"\n[SIGNATURE]")
    if signature_start >= 0:
        signature["computed_md5"] = hashlib.md5(data[:signature_start]).hexdigest()
        tail = data[signature_start:].decode("ascii", errors="replace")
        stored_type = re.search(r"^Type\s+(\S+)\s*$", tail, re.M)
        stored_hash = re.search(r"^Hash\s+([0-9a-fA-F]+)\s*$", tail, re.M)
        signature["stored_type"] = stored_type.group(1) if stored_type else None
        signature["stored_hash"] = stored_hash.group(1).lower() if stored_hash else None
        signature["valid"] = (signature["stored_type"] == "md5"
                              and signature["stored_hash"] == signature["computed_md5"])

    sections = {"": []}
    current = ""
    for line_number, line in enumerate(data.decode("utf-8", errors="replace").splitlines(),
                                       start=1):
        header = SECTION_HEADER.fullmatch(line.strip())
        if header:
            current = header.group(1)
            sections.setdefault(current, [])
            continue
        sections[current].append((line_number, line))
    return {"sections": sections, "signature": signature}


def compare_sequences(candidate_path, reference_path, ignored_definitions, rtol, atol):
    """Compare two Pulseq files section by section.

    Args:
        candidate_path: Regenerated ``.seq`` path.
        reference_path: Reference ``.seq`` path.
        ignored_definitions: ``[DEFINITIONS]`` keys that may differ.
        rtol: Relative tolerance of the numeric diagnosis.
        atol: Absolute tolerance of the numeric diagnosis.

    Returns:
        dict: Both signature checks, ``expected_differences`` with the values
            of ignored definitions, ``sections`` listing every section that is
            not text-identical, ``text_identical`` and ``match``. A pair
            matches when both signatures are valid and every section is
            text-identical or numerically equal within tolerance.
    """
    candidate = read_sequence_file(candidate_path)
    reference = read_sequence_file(reference_path)
    expected = {}
    sections = []
    names = list(dict.fromkeys(list(reference["sections"]) + list(candidate["sections"])))
    for name in names:
        if name == "SIGNATURE":
            continue
        candidate_lines = candidate["sections"].get(name)
        reference_lines = reference["sections"].get(name)
        if candidate_lines is None or reference_lines is None:
            sections.append({"section": name, "status": "different",
                             "issue": "missing in candidate" if candidate_lines is None
                                      else "missing in reference"})
            continue
        if name == "DEFINITIONS":
            candidate_lines, candidate_removed = split_ignored_definitions(
                candidate_lines, ignored_definitions)
            reference_lines, reference_removed = split_ignored_definitions(
                reference_lines, ignored_definitions)
            for key in sorted(set(candidate_removed) | set(reference_removed)):
                if candidate_removed.get(key) != reference_removed.get(key):
                    expected[key] = {"candidate": candidate_removed.get(key),
                                     "reference": reference_removed.get(key)}
        result = compare_section_lines(name, candidate_lines, reference_lines, rtol, atol)
        if result["status"] != "identical":
            sections.append(result)

    signatures_valid = candidate["signature"]["valid"] and reference["signature"]["valid"]
    return {
        "candidate_signature": candidate["signature"],
        "reference_signature": reference["signature"],
        "expected_differences": expected,
        "sections": sections,
        "text_identical": not sections,
        "match": signatures_valid and all(s["status"] == "numerically equal" for s in sections),
    }


def split_ignored_definitions(lines, ignored_definitions):
    """Remove ignored keys from ``[DEFINITIONS]`` lines.

    Args:
        lines: ``(line_number, text)`` tuples of the section.
        ignored_definitions: Keys to remove.

    Returns:
        tuple: Remaining ``(line_number, text)`` tuples and a dict mapping each
            removed key to its value text.
    """
    kept, removed = [], {}
    for line_number, text in lines:
        fields = text.split(maxsplit=1)
        if fields and fields[0] in ignored_definitions:
            removed[fields[0]] = fields[1].strip() if len(fields) > 1 else ""
        else:
            kept.append((line_number, text))
    return kept, removed


def compare_section_lines(name, candidate_lines, reference_lines, rtol, atol):
    """Compare the lines of one section by position.

    Args:
        name: Section name, used in the result.
        candidate_lines: ``(line_number, text)`` tuples of the candidate.
        reference_lines: ``(line_number, text)`` tuples of the reference.
        rtol: Relative tolerance of the numeric diagnosis.
        atol: Absolute tolerance of the numeric diagnosis.

    Returns:
        dict: ``status`` is ``"identical"`` (same text), ``"numerically
            equal"`` (same line count and every differing line equal token by
            token within tolerance) or ``"different"``; counts, the largest
            numeric deviation and up to ``MAX_REPORTED_DIFFERENCES`` differing
            line pairs with their file line numbers.
    """
    candidate_text = [text for _, text in candidate_lines]
    reference_text = [text for _, text in reference_lines]
    if candidate_text == reference_text:
        return {"section": name, "status": "identical", "lines": len(reference_text)}

    n_lines = max(len(candidate_text), len(reference_text))
    text_differences = 0
    numeric_differences = 0
    max_deviation = 0.0
    first_differences = []
    for index in range(n_lines):
        candidate = candidate_text[index] if index < len(candidate_text) else None
        reference = reference_text[index] if index < len(reference_text) else None
        if candidate == reference:
            continue
        text_differences += 1
        equal, deviation = False, float("nan")
        if candidate is not None and reference is not None:
            equal, deviation = tokens_numerically_equal(candidate, reference, rtol, atol)
            if np.isfinite(deviation):
                max_deviation = max(max_deviation, deviation)
        if equal:
            continue
        numeric_differences += 1
        if len(first_differences) < MAX_REPORTED_DIFFERENCES:
            first_differences.append({
                "candidate_line": candidate_lines[index][0] if candidate is not None else None,
                "reference_line": reference_lines[index][0] if reference is not None else None,
                "candidate": candidate,
                "reference": reference,
            })

    same_length = len(candidate_text) == len(reference_text)
    return {
        "section": name,
        "status": "numerically equal" if same_length and numeric_differences == 0
                  else "different",
        "candidate_lines": len(candidate_text),
        "reference_lines": len(reference_text),
        "text_differences": text_differences,
        "numeric_differences": numeric_differences,
        "max_abs_numeric_deviation": max_deviation,
        "first_differences": first_differences,
    }


def tokens_numerically_equal(candidate, reference, rtol, atol):
    """Compare two whitespace-separated lines token by token.

    Args:
        candidate: Candidate line text.
        reference: Reference line text.
        rtol: Relative tolerance for tokens that parse as numbers.
        atol: Absolute tolerance for tokens that parse as numbers.

    Returns:
        tuple[bool, float]: Whether all tokens agree (numbers within
            tolerance, other tokens exactly) and the largest absolute
            deviation between numeric tokens (NaN when token counts differ).
    """
    candidate_tokens = candidate.split()
    reference_tokens = reference.split()
    if len(candidate_tokens) != len(reference_tokens):
        return False, float("nan")
    equal = True
    max_deviation = 0.0
    for candidate_token, reference_token in zip(candidate_tokens, reference_tokens):
        try:
            candidate_value = float(candidate_token)
            reference_value = float(reference_token)
        except ValueError:
            equal = equal and candidate_token == reference_token
            continue
        if np.isnan(candidate_value) and np.isnan(reference_value):
            continue
        deviation = abs(candidate_value - reference_value)
        max_deviation = max(max_deviation, deviation)
        equal = equal and bool(np.isclose(candidate_value, reference_value,
                                          rtol=rtol, atol=atol))
    return equal, max_deviation


def compare_mat_files(candidate_path, reference_path, ignored_fields, rtol, atol):
    """Compare two MAT-file sidecars variable by variable.

    Args:
        candidate_path: Regenerated sidecar path.
        reference_path: Reference sidecar path.
        ignored_fields: Field names (any depth) that may differ.
        rtol: Relative tolerance for numeric arrays.
        atol: Absolute tolerance for numeric arrays.

    Returns:
        dict: Paths, the number of reference variables, up to
            ``MAX_REPORTED_DIFFERENCES`` differences with dotted field paths,
            the total difference count, ``expected_differences`` and ``match``.
    """
    candidate = load_mat_variables(candidate_path)
    reference = load_mat_variables(reference_path)
    differences, expected = [], {}
    compare_values(candidate, reference, "", ignored_fields, rtol, atol, differences, expected)
    return {
        "candidate": str(candidate_path),
        "reference": str(reference_path),
        "variables": len(reference),
        "difference_count": len(differences),
        "differences": differences[:MAX_REPORTED_DIFFERENCES],
        "expected_differences": expected,
        "match": not differences,
    }


def load_mat_variables(path):
    """Load the user variables of a MAT file (format 5/7).

    Args:
        path: MAT-file path.

    Returns:
        dict: Variable name to value, with MATLAB structs as dicts, cell
            arrays as lists and char arrays as ``str`` (SciPy
            ``simplify_cells``); MAT header entries are dropped. Numeric
            values use their MATLAB class (``mat_dtype``), not the compressed
            storage type, so an integral double is loaded as float64.
    """
    variables = sio.loadmat(path, simplify_cells=True, mat_dtype=True)
    return {name: value for name, value in variables.items() if not name.startswith("__")}


def compare_values(candidate, reference, path, ignored_fields, rtol, atol, differences, expected):
    """Recursively compare two loaded MAT values and record differences.

    Args:
        candidate: Candidate value (dict, list, str, scalar or array).
        reference: Reference value of the same field.
        path: Dotted field path of the value; ``""`` for the top level.
        ignored_fields: Field names whose differences are expected.
        rtol: Relative tolerance for numeric arrays.
        atol: Absolute tolerance for numeric arrays.
        differences: List extended in place with difference records.
        expected: Dict updated in place with expected differences by path.

    Returns:
        None.
    """
    field = re.sub(r"\[.*$", "", path.rsplit(".", 1)[-1]) if path else ""
    if field in ignored_fields:
        if not values_equal(candidate, reference):
            expected[path] = {"candidate": summarize_value(candidate),
                              "reference": summarize_value(reference)}
        return

    if isinstance(candidate, dict) or isinstance(reference, dict):
        if not (isinstance(candidate, dict) and isinstance(reference, dict)):
            differences.append({"path": path, "issue": "struct/non-struct mismatch"})
            return
        for key in sorted(set(candidate) | set(reference)):
            child = f"{path}.{key}" if path else key
            if key not in candidate:
                differences.append({"path": child, "issue": "missing in candidate"})
            elif key not in reference:
                differences.append({"path": child, "issue": "missing in reference"})
            else:
                compare_values(candidate[key], reference[key], child, ignored_fields,
                               rtol, atol, differences, expected)
        return

    if isinstance(candidate, list) or isinstance(reference, list):
        if not (isinstance(candidate, list) and isinstance(reference, list)):
            differences.append({"path": path, "issue": "cell/non-cell mismatch"})
            return
        if len(candidate) != len(reference):
            differences.append({"path": path,
                                "issue": f"length {len(candidate)} vs {len(reference)}"})
            return
        for index, (candidate_item, reference_item) in enumerate(zip(candidate, reference)):
            compare_values(candidate_item, reference_item, f"{path}[{index}]", ignored_fields,
                           rtol, atol, differences, expected)
        return

    if isinstance(candidate, str) or isinstance(reference, str):
        if candidate != reference:
            differences.append({"path": path, "candidate": summarize_value(candidate),
                                "reference": summarize_value(reference)})
        return

    candidate_array = np.asarray(candidate)
    reference_array = np.asarray(reference)
    if candidate_array.dtype == object or reference_array.dtype == object:
        if candidate_array.shape != reference_array.shape:
            differences.append({"path": path, "issue": f"shape {candidate_array.shape} "
                                                       f"vs {reference_array.shape}"})
            return
        for index in np.ndindex(reference_array.shape):
            compare_values(candidate_array[index], reference_array[index],
                           f"{path}[{','.join(map(str, index))}]", ignored_fields,
                           rtol, atol, differences, expected)
        return

    if candidate_array.shape != reference_array.shape:
        differences.append({"path": path, "issue": f"shape {candidate_array.shape} "
                                                   f"vs {reference_array.shape}"})
        return
    if candidate_array.dtype != reference_array.dtype:
        differences.append({"path": path, "issue": f"dtype {candidate_array.dtype} "
                                                   f"vs {reference_array.dtype}"})
        return
    if candidate_array.dtype.kind not in "iufc":
        # Logical, byte and other fixed-size values are compared exactly.
        if not np.array_equal(candidate_array, reference_array):
            differences.append({"path": path, "issue": f"{candidate_array.dtype} values differ"})
        return
    if not np.allclose(candidate_array, reference_array, rtol=rtol, atol=atol, equal_nan=True):
        deviation = np.abs(candidate_array.astype(np.complex128)
                           - reference_array.astype(np.complex128))
        differences.append({"path": path,
                            "max_abs_deviation": float(np.nanmax(deviation)),
                            "differing_elements": int(np.count_nonzero(
                                ~np.isclose(candidate_array, reference_array, rtol=rtol,
                                            atol=atol, equal_nan=True)))})


def values_equal(candidate, reference):
    """Return whether two loaded MAT values are exactly equal.

    Args:
        candidate: Candidate value.
        reference: Reference value.

    Returns:
        bool: True for identical structure and values.
    """
    differences = []
    compare_values(candidate, reference, "", (), 0.0, 0.0, differences, {})
    return not differences


def compare_text_reports(candidate_path, reference_path, ignored_fields):
    """Compare two ``Key: value`` text reports line by line.

    Args:
        candidate_path: Regenerated report path.
        reference_path: Reference report path.
        ignored_fields: Keys whose value lines may differ.

    Returns:
        dict: Paths, line count, up to ``MAX_REPORTED_DIFFERENCES``
            differing lines, the total count, ``expected_differences`` and
            ``match``.
    """
    candidate_lines = candidate_path.read_text().splitlines()
    reference_lines = reference_path.read_text().splitlines()
    differences, expected = [], {}
    for index in range(max(len(candidate_lines), len(reference_lines))):
        candidate = candidate_lines[index] if index < len(candidate_lines) else None
        reference = reference_lines[index] if index < len(reference_lines) else None
        if candidate == reference:
            continue
        candidate_key = candidate.split(":", 1)[0] if candidate is not None else None
        reference_key = reference.split(":", 1)[0] if reference is not None else None
        if candidate_key == reference_key and reference_key in ignored_fields:
            expected[reference_key] = {"candidate": candidate, "reference": reference}
            continue
        differences.append({"line": index + 1, "candidate": candidate, "reference": reference})
    return {
        "candidate": str(candidate_path),
        "reference": str(reference_path),
        "lines": len(reference_lines),
        "difference_count": len(differences),
        "differences": differences[:MAX_REPORTED_DIFFERENCES],
        "expected_differences": expected,
        "match": not differences,
    }


def summarize_value(value):
    """Return a short JSON-friendly description of a loaded MAT value.

    Args:
        value: Loaded value.

    Returns:
        str: The string itself (truncated to 300 characters) or a
            ``type/shape/dtype`` description.
    """
    if isinstance(value, str):
        return value if len(value) <= 300 else value[:297] + "..."
    array = np.asarray(value) if not isinstance(value, (dict, list)) else None
    if array is None:
        return f"{type(value).__name__} with {len(value)} entries"
    if array.size == 1 and array.dtype.kind in "biufc":
        return repr(array.item())
    return f"array shape={array.shape} dtype={array.dtype}"


def print_export_result(result):
    """Print a readable summary of one export comparison.

    Args:
        result: Dict returned by ``compare_export``.

    Returns:
        None.
    """
    sequence = result["sequence"]
    if not sequence["match"]:
        status = "DIFFERENT"
    elif sequence["text_identical"]:
        status = "MATCH (identical text)"
    else:
        status = "MATCH (numerically equal; formatting differs)"
    print(f"\n== {result['name']}")
    print(f"   candidate: {result['candidate']}")
    print(f"   reference: {result['reference']}")
    print(f"   sequence: {status}")
    for side in ("candidate", "reference"):
        signature = sequence[f"{side}_signature"]
        if not signature["valid"]:
            print(f"   {side} md5 signature INVALID: stored {signature['stored_hash']}, "
                  f"computed {signature['computed_md5']}")
    for key, values in sequence["expected_differences"].items():
        print(f"   expected definition difference {key}: candidate {values['candidate']!r}, "
              f"reference {values['reference']!r}")
    for section in sequence["sections"]:
        if "issue" in section:
            print(f"   [{section['section']}] {section['issue']}")
            continue
        print(f"   [{section['section']}] {section['status']}: "
              f"{section['text_differences']} text / {section['numeric_differences']} numeric "
              f"differing lines, {section['candidate_lines']}/{section['reference_lines']} lines, "
              f"max |numeric deviation| {section['max_abs_numeric_deviation']:.3g}")
        for difference in section["first_differences"]:
            print(f"      line {difference['candidate_line']}/{difference['reference_line']}: "
                  f"{difference['candidate']!r} | {difference['reference']!r}")
    for sidecar in result["sidecars"]:
        if "issue" in sidecar:
            print(f"   sidecar {sidecar['suffix']}: DIFFERENT ({sidecar['issue']})")
            continue
        status = "MATCH" if sidecar["match"] else f"DIFFERENT ({sidecar['difference_count']})"
        print(f"   sidecar {sidecar['suffix']}: {status}")
        for difference in sidecar["differences"]:
            print(f"      {difference}")
        for path in sidecar["expected_differences"]:
            print(f"      expected difference: {path}")


def json_ready(value):
    """Convert NumPy scalars, paths and tuples for JSON serialization.

    Args:
        value: Nested result structure.

    Returns:
        The same structure built from JSON-serializable Python types.
    """
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


if __name__ == "__main__":
    sys.exit(main())
