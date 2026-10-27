"""Compare a reconstructed TSE volume and its NIfTI export with a reference result.

The complex volumes are compared as saved by the reconstruction (``.npy``,
shape (Nx, Ny, physical slice)); NIfTI exports are compared by their affine
and data. A regression passes when the relative L2 difference of the complex
volumes and of every NIfTI pair is at most ``--rtol`` and the affines agree
within ``--affine-atol`` mm.

Run from the repository root, for example::

    python evaluation/compare_reconstructions.py \
        --candidate /path/to/new/arrays/lambda0_complex_physical_slices.npy \
        --reference /path/to/old/arrays/complex_physical_slices.npy \
        --nifti /path/to/new/nifti/lambda0_magnitude.nii.gz \
                /path/to/old/nifti/old_magnitude.nii.gz \
        --report /path/to/regression_report.json

The exit status is 0 when every comparison passes and 1 otherwise.
"""

import argparse
import json
from pathlib import Path

import nibabel as nib
import numpy as np


def main(argv=None):
    """Compare the requested volumes and NIfTI pairs and print the result.

    Args:
        argv: Command-line arguments without the program name; ``None`` reads
            ``sys.argv[1:]``.

    Returns:
        int: 0 when all comparisons pass, else 1.

    Side effects:
        Prints a report and, with ``--report``, writes it as new JSON.
    """
    args = parse_args(argv)
    results = []
    if args.candidate is not None:
        results.append(compare_volumes(np.load(args.candidate), np.load(args.reference), args.rtol))
        results[-1].update(kind="complex volume", candidate=str(args.candidate),
                           reference=str(args.reference))
    for candidate, reference in args.nifti:
        results.append(compare_nifti(Path(candidate), Path(reference), args.rtol, args.affine_atol))
    for result in results:
        print_result(result)
    passed = all(result["passed"] for result in results)
    print(f"\n{sum(r['passed'] for r in results)}/{len(results)} comparisons pass.")
    if args.report is not None:
        with args.report.open("x") as handle:
            json.dump({"passed": passed, "rtol": args.rtol, "affine_atol_mm": args.affine_atol,
                       "comparisons": results}, handle, indent=2)
        print(f"Report written to {args.report}")
    return 0 if passed else 1


def parse_args(argv):
    """Return the parsed command-line arguments.

    Raises:
        SystemExit: Through argparse when nothing is compared or a volume
            lacks its counterpart.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--candidate", type=Path, help="New complex volume (.npy).")
    parser.add_argument("--reference", type=Path, help="Reference complex volume (.npy).")
    parser.add_argument("--nifti", nargs=2, action="append", default=[], metavar=("CANDIDATE", "REFERENCE"),
                        help="A NIfTI pair to compare; may be repeated.")
    parser.add_argument("--rtol", type=float, default=1e-5,
                        help="Largest accepted relative L2 difference (default: 1e-5).")
    parser.add_argument("--affine-atol", type=float, default=1e-4,
                        help="Largest accepted affine entry difference [mm] (default: 1e-4).")
    parser.add_argument("--report", type=Path, help="New JSON report file.")
    args = parser.parse_args(argv)
    if (args.candidate is None) != (args.reference is None):
        parser.error("--candidate and --reference must be given together")
    if args.candidate is None and not args.nifti:
        parser.error("Nothing to compare")
    return args


def compare_volumes(candidate, reference, rtol):
    """Compare two arrays of the same shape.

    Args:
        candidate, reference: Real or complex arrays.
        rtol: Largest accepted relative L2 difference.

    Returns:
        dict: Shapes, relative L2 difference, largest absolute difference
        relative to the reference maximum, the worst slice along the last axis
        and the pass flag.
    """
    result = {"candidate_shape": list(candidate.shape), "reference_shape": list(reference.shape)}
    if candidate.shape != reference.shape:
        return {**result, "passed": False, "reason": "shape differs"}
    difference = candidate.astype(np.complex128) - reference.astype(np.complex128)
    reference_norm = np.linalg.norm(reference.astype(np.complex128))
    relative = float(np.linalg.norm(difference) / max(reference_norm, 1e-30))
    per_slice = np.linalg.norm(difference.reshape(-1, difference.shape[-1]), axis=0) / np.maximum(
        np.linalg.norm(reference.astype(np.complex128).reshape(-1, reference.shape[-1]), axis=0), 1e-30)
    finite = bool(np.isfinite(candidate).all())
    return {**result, "relative_l2_difference": relative,
            "max_abs_difference_over_reference_max": float(np.max(abs(difference)) /
                                                           max(float(np.max(abs(reference))), 1e-30)),
            "worst_slice": int(np.argmax(per_slice)), "worst_slice_relative_l2": float(per_slice.max()),
            "candidate_finite": finite, "passed": finite and relative <= rtol}


def compare_nifti(candidate, reference, rtol, affine_atol):
    """Compare the affines and data of two NIfTI files.

    Args:
        candidate, reference: NIfTI paths.
        rtol: Largest accepted relative L2 data difference.
        affine_atol: Largest accepted affine entry difference [mm].

    Returns:
        dict: Data comparison of :func:`compare_volumes` plus the affine
        difference, both orientation codes and the pass flag.
    """
    new, old = nib.load(candidate), nib.load(reference)
    result = compare_volumes(np.asarray(new.dataobj), np.asarray(old.dataobj), rtol)
    affine_difference = float(np.max(abs(new.affine - old.affine)))
    codes = ["".join(nib.aff2axcodes(image.affine)) for image in (new, old)]
    return {**result, "kind": "NIfTI", "candidate": str(candidate), "reference": str(reference),
            "affine_max_abs_difference_mm": affine_difference, "axis_codes": codes,
            "passed": result["passed"] and affine_difference <= affine_atol}


def print_result(result):
    """Print one comparison in a short human-readable form."""
    status = "PASS" if result["passed"] else "FAIL"
    print(f"[{status}] {result['kind']}: {result['candidate']}")
    print(f"       reference: {result['reference']}")
    if "relative_l2_difference" not in result:
        print(f"       {result['reason']}: {result['candidate_shape']} vs {result['reference_shape']}")
        return
    print(f"       relative L2 difference {result['relative_l2_difference']:.3e}, "
          f"max |difference| / max |reference| {result['max_abs_difference_over_reference_max']:.3e}, "
          f"worst slice {result['worst_slice']} ({result['worst_slice_relative_l2']:.3e})")
    if "affine_max_abs_difference_mm" in result:
        print(f"       affine difference {result['affine_max_abs_difference_mm']:.3e} mm, "
              f"axis codes {result['axis_codes'][0]} / {result['axis_codes'][1]}")


if __name__ == "__main__":
    raise SystemExit(main())
