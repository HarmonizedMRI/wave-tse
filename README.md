# Wave TSE

Pulseq-based 2D Wave-CAIPI turbo spin echo (TSE) sequence generation for
single-band (MB1) and two-band simultaneous-multislice (MB2) acquisitions,
a mixed-dimensional FLASH calibration, and Python reconstruction from
Siemens TWIX data.

The repository has two intentionally separate parts:

- MATLAB generates the MB2 TSE, the MB1 TSE and the calibration as three
  independent Pulseq sequences.
- Python reconstructs TSE data from the calibration TWIX, the TSE TWIX, and
  the matching `.seq` files and sidecars.

`pyproject.toml` defines only the Python environment. The repository is not an
installable Python package.

## Repository layout

```text
.
├── README.md
├── LICENSE
├── pyproject.toml
├── docs/
│   ├── sequence.md
│   └── calibration_storage.md
├── seq/
│   ├── tse_2d_wave_mb2v_v3.m
│   ├── tse_2d_wave_mb1.m
│   ├── gre_mixed_dim_mb2_calibration.m
│   ├── rf_pulses/
│   │   └── mb2_tse_ls_tbw7_beta150_v1verse_B1max15p000uT.mat
│   └── utils/
├── evaluation/
│   ├── compare_sequence_exports.py
│   └── tests/
└── external/
    └── safe_pns_prediction/   (submodule)
```

The reconstruction (`recon/`) is being assembled from the validated research
pipeline and will be documented in `docs/reconstruction.md`.

## Clone

```bash
git clone --recurse-submodules git@github.com:HarmonizedMRI/wave-tse.git
cd wave-tse
```

## Generate sequences

| Generator | Sequence |
|---|---|
| `seq/tse_2d_wave_mb2v_v3.m` | MB2 TSE, 50 slices in 25 SMS groups, blipped CAIPI, optional y/z Wave |
| `seq/tse_2d_wave_mb1.m` | MB1 TSE, 50 slices, conventional sinc RF, optional y Wave |
| `seq/gre_mixed_dim_mb2_calibration.m` | FLASH calibration: 3D y/z Wave PSF projections and 2D slice-wise CSM |

Open a generator in MATLAB, edit its in-file parameters and run it, for
example:

```matlab
cd seq
tse_2d_wave_mb2v_v3
```

The first run asks for the Pulseq folder, the optional safe_pns_prediction
folder, the output folder and the optional scanner `.asc` file, and saves them
in `seq/wave_tse_path_settings.json` (ignored by Git).

The MB2 RF is the distributed SigPy-optimized 50-slice TBW7 beta150 bundle in
`seq/rf_pulses/`. To use a different MB2 pulse, generate a bundle with the
same variables; see [Sequence generation](docs/sequence.md#rf-pulses).

`evaluation/compare_sequence_exports.py` compares regenerated exports and
their sidecars with reference exports, for example after updating Pulseq:

```bash
python evaluation/compare_sequence_exports.py \
  --candidate-dir /path/to/regenerated/out_path \
  --reference-dir /path/to/reference/protocol_folder
```

## Documentation

- [Sequence generation](docs/sequence.md)
- [Calibration storage](docs/calibration_storage.md)

## Attribution

Author: Yiyun Dong, Athinoula A. Martinos Center for Biomedical Imaging.

The TSE implementation is based on Berkin's Wave-TSE code in the
HarmonizedMRI `wave-haste` repository and the Pulseq TSE example code by
Juergen Hennig and Maxim Zaitsev. The FLASH calibration organization and
path-setting workflow follow the Wave-MPRAGE project.

## License

MIT License. See `LICENSE`. The `external/safe_pns_prediction` submodule keeps
its own license.
