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

`environment.yml` (conda) and `requirements.txt` (venv) define the Python
environment. The repository is not an installable Python package.

## Repository layout

```text
.
├── README.md
├── LICENSE
├── environment.yml
├── requirements.txt
├── docs/
│   ├── sequence.md
│   ├── calibration_storage.md
│   ├── reconstruction.md
│   └── troubleshooting.md
├── seq/
│   ├── tse_2d_wave_mb2v_v3.m
│   ├── tse_2d_wave_mb1.m
│   ├── gre_mixed_dim_mb2_calibration.m
│   ├── rf_pulses/
│   │   └── mb2_tse_ls_tbw7_beta150_v1verse_B1max15p000uT.mat
│   └── utils/
├── recon/
│   ├── recon_wave_tse_from_twix_nifti.py
│   ├── utils/
│   └── tests/
├── evaluation/
│   ├── compare_sequence_exports.py
│   ├── compare_reconstructions.py
│   └── tests/
└── external/
    └── safe_pns_prediction/   (submodule)
```

## Clone

```bash
git clone --recurse-submodules git@github.com:HarmonizedMRI/wave-tse.git
cd wave-tse
```

## Python environment

The Python tools need the packages in `requirements.txt`. Both routes below
install the same list. The code was tested with Python 3.12, NumPy 2.4,
SciPy 1.17, PyTorch 2.13, SigPy 0.1.27, pyMapVBVD 0.6.1 and PyPulseq 1.5.0.

With conda:

```bash
conda env create -f environment.yml
conda activate wave-tse
```

With a Python virtual environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Everything runs on the CPU by default. For GPU computation, PyTorch needs a
CUDA-enabled build; the Linux wheels from PyPI include one. GPU ESPIRiT
through SigPy also needs the CuPy wheel that matches the CUDA driver, for
example `python -m pip install cupy-cuda12x` or `cupy-cuda13x`.

Run the tests from the repository root:

```bash
python -m unittest discover -s evaluation/tests
python -m unittest discover -s recon/tests -t recon
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

## Reconstruct

`recon/recon_wave_tse_from_twix_nifti.py` reconstructs one MB1 or MB2 TSE scan
with the calibration acquired in the same session:

```bash
python recon/recon_wave_tse_from_twix_nifti.py \
  --calibration-twix /path/to/meas_..._pulseq_calib.dat \
  --calibration-seq /path/to/gre_mixed_dim_mb2_calib_....seq \
  --tse-twix /path/to/meas_..._pulseq_tse_....dat \
  --tse-seq /path/to/v151/tse_....seq \
  --calibration-dir /path/to/output/calibration \
  --output-dir /path/to/output/tse_scan
```

The sidecars are found next to their `.seq` files. The command validates the
headers, sequences, sidecars, ADC chronology, slice mapping and trajectories
before it computes the shared calibration (noise whitening, coil compression,
slice-wise ESPIRiT maps, Wave PSF) and the reconstruction: the exact
least-squares solution for complete no-Wave MB1 data, lambda-0 CG-SENSE for
MB2 data. It writes magnitude and phase NIfTI in scanner geometry. Output
files are never overwritten; reruns resume from the saved checkpoints. See
[Reconstruction](docs/reconstruction.md).

`evaluation/compare_reconstructions.py` compares a reconstruction with a
reference result, for example for a regression after a software update.

## Documentation

- [Sequence generation](docs/sequence.md)
- [Calibration storage](docs/calibration_storage.md)
- [Reconstruction](docs/reconstruction.md)
- [Troubleshooting](docs/troubleshooting.md)

## Attribution

Author: Yiyun Dong, Athinoula A. Martinos Center for Biomedical Imaging.

The TSE implementation is based on Berkin's Wave-TSE code in the
HarmonizedMRI `wave-haste` repository and the Pulseq TSE example code by
Juergen Hennig and Maxim Zaitsev. The FLASH calibration organization and
path-setting workflow follow the Wave-MPRAGE project.

## License

MIT License. See `LICENSE`. The `external/safe_pns_prediction` submodule keeps
its own license.
