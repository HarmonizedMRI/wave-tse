"""Check a TSE export with an appended FLASH calibration before scanning.

The reconstruction audits the same contracts on the acquired data; this
script runs the sequence-level part on the exported ``.seq`` file and its
``_PEorder.mat`` sidecar alone, so that an export can be checked right after
generation:

- the sidecar's sequence signature equals the ``.seq`` file;
- the calibration REF readouts follow the TSE ADCs with the compact LIN
  labels and zero SLC/PAR/SET/ECO;
- every Wave projection differs from its reference by one line-independent
  Wave trajectory;
- the ACS lines sit on the stated ky/kz grid with one common kx;
- the ACS z grid has a plane at every slice center.

Run from the repository root, for example::

    python evaluation/check_appended_calibration.py \
        --seq /path/to/out_path/v151/tse_..._cal_skyra.seq \
        --report /path/to/check_report.json

The exit status is 0 when every check passes and 1 otherwise.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "recon"))
import recon_wave_tse_from_twix_nifti as cli  # noqa: E402
from utils.appended_calibration import acs_z_grid, load_appended_metadata  # noqa: E402


def main(argv=None):
    """Run the checks and print the summary.

    Args:
        argv: Command-line arguments without the program name; ``None`` reads
            ``sys.argv[1:]``.

    Returns:
        int: 0 when every check passes, else 1.

    Side effects:
        Prints a JSON summary and, with ``--report``, writes it as a new file.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--seq", type=Path, required=True, help="Exported TSE .seq file.")
    parser.add_argument("--peorder", type=Path,
                        help="Its *_PEorder.mat sidecar (default: found as by the reconstruction).")
    parser.add_argument("--report", type=Path, help="New JSON report file.")
    args = parser.parse_args(argv)
    summary = check_export(args.seq, args.peorder)
    text = json.dumps(cli.json_ready(summary), indent=2)
    print(text)
    if args.report is not None:
        with args.report.open("x") as handle:
            handle.write(text + "\n")
    return 0 if summary["passed"] else 1


def check_export(sequence_path, peorder=None):
    """Return the check summary of one export.

    Args:
        sequence_path: Exported TSE ``.seq`` file.
        peorder: Its sidecar, or None to locate it beside the sequence.

    Returns:
        dict: ``passed``, the first ``error`` (or None) and, when the checks
        ran, the TSE/calibration ADC counts, calibration timing, ACS grid and
        slice planes.
    """
    sequence_path = Path(sequence_path).resolve(strict=True)
    if peorder is None:
        args = argparse.Namespace(calibration_twix=None, tse_twix=sequence_path, tse_seq=sequence_path,
                                  tse_peorder=None)
        peorder = cli.resolve_inputs(args)["tse"]["peorder"]
    summary = {"sequence": str(sequence_path), "peorder": str(peorder), "passed": False, "error": None}
    try:
        pe = loadmat(peorder, simplify_cells=True)
        meta = load_appended_metadata(pe)
        seq, labels, _, _, _ = cli.sequence_details(sequence_path)
        if "sequenceMetadata" in pe and seq.signature_value != pe["sequenceMetadata"]["SequenceSignature"]:
            raise ValueError("TSE sequence signature differs from its sidecar")
        tse_adcs = cli.tse_adc_count(labels, pe)
        evidence = cli.appended_sequence_contract(sequence_path, meta)
        grid = acs_z_grid(meta)
    except ValueError as error:
        summary["error"] = str(error)
        return summary
    validation = meta.get("validation", {})
    summary.update({
        "passed": True, "tse_adc_count": tse_adcs, "calibration_adc_count": int(meta["expectedAdcCount"]),
        "psf_valid_lengths": np.asarray(meta["psfValidLengths"]).tolist(),
        "acs_lines_ky_kz": list(np.asarray(meta["acsAdcIndices"]).shape),
        "wave_axes": sorted(key[-1] for key in evidence if key.startswith("delta_")),
        "calibration_te_tr_ms": [float(meta[key]) * 1e3 if key in meta else None for key in ("TE", "TR")],
        "calibration_duration_s": validation.get("duration"),
        "sequence_duration_s": float(seq.duration()[0]),
        "acs_z_grid": {"planes": grid["n"], "spacing_mm": grid["spacing"] * 1e3, "z0_mm": grid["z0"] * 1e3,
                       "slice_planes": grid["slice_planes"].tolist()}})
    return summary


if __name__ == "__main__":
    sys.exit(main())
