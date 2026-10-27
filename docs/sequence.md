# Sequence generation

## Entry points

| Generator | Sequence | RF | Writes to `out_path` |
|---|---|---|---|
| `seq/tse_2d_wave_mb2v_v3.m` | MB2 TSE: 50 physical slices in 25 SMS groups, blipped CAIPI, optional two-axis Wave | TBW7 beta150 vMB-VERSE bundle `seq/rf_pulses/mb2_tse_ls_tbw7_beta150_v1verse_B1max15p000uT.mat` | `v151/tse_m2b7v1_*.seq`, `metadata/<name>_PEorder.mat` |
| `seq/tse_2d_wave_mb1.m` | MB1 TSE: 50 single-band slices, optional y-sine Wave | Pulseq sinc pulses (`rfDesign = 'conventional_sinc'`) | `v151/tse_mb1_*.seq`, `<name>_PEorder.mat` |
| `seq/gre_mixed_dim_mb2_calibration.m` | Mixed-dimension spoiled-GRE calibration: 3D y/z Wave PSF projections plus 2D slice-wise CSM lines for the same 50 slices | Pulseq sinc pulses; the TBW7 bundle supplies slice geometry only | `gre_mixed_dim_mb2_calib_*.seq`, `<name>_geometry_labels.mat`, `<name>_geometry_labels.txt` |

The three sequences are independent measurements. A session acquires the
calibration and one or more TSE scans with the same prescription; the
reconstruction reads both TWIX files. An integrated calibration-plus-TSE
sequence is not part of this repository.

MATLAB is the authoritative sequence implementation. The calibration states
the storage contract of its raw data in [Calibration storage](calibration_storage.md).

## Prerequisites

- MATLAB and the Pulseq MATLAB toolbox. The generators write Pulseq v1.5.1
  files; `write_v141_format = true` additionally writes v1.4.1 copies.
- For the optional PNS/CNS and forbidden-frequency checks: the
  `external/safe_pns_prediction` submodule and the Siemens gradient-system
  `.asc` file of the scanner.

```bash
git submodule update --init external/safe_pns_prediction
```

## Running a generator

Open the generator in MATLAB, edit the parameters in the file, and run it, for
example:

```matlab
cd seq
tse_2d_wave_mb2v_v3
```

All parameters are ordinary variables near the top of each file; the
generators read no environment variables.

### Path settings

On the first run, `seq/utils/configurePathSettings.m` prompts for four paths
and saves them in `seq/wave_tse_path_settings.json`, which is shared by the
three generators and ignored by Git:

| Setting | Meaning |
|---|---|
| `pulseq_path` | Pulseq repository root or its `matlab/` folder (must contain `+mr`) |
| `safe_pns_prediction_path` | safe_pns_prediction folder, e.g. `<repo>/external/safe_pns_prediction`; optional |
| `out_path` | Export folder; pressing Enter uses MATLAB's current folder |
| `system_asc_file` | Scanner gradient-system `.asc` file; optional |

Later runs ask whether to reuse the saved settings and which entries to
change. Enter `-` to clear an optional entry. Paths starting with `~` are
expanded to the home directory, and the TSE generators record the expanded
`out_path` in the `OutputPath` definition.

The calibration refuses to overwrite an existing export; select a new
`out_path` for each regeneration.

## RF pulses

The repository distributes one RF bundle,
`seq/rf_pulses/mb2_tse_ls_tbw7_beta150_v1verse_B1max15p000uT.mat`
(sha256 `4bb7239185eaf17f6298ba99f2f6da8ecbea11778e022e5c27c54b126d8f35b5`;
identical to the research bundle except that the provenance fields
`reference_arrays_path` and `v3_excitation_source_path` hold relative paths). It
holds the SigPy-optimized TBW7 beta150 MB2 excitation, first-165 and
repeated-150 refocusing waveforms with their time-varying vMB-VERSE slice
gradients for the 50-slice, 3 mm, gap-free geometry (25 SMS groups, pairs
`[p, p+25]`), designed for a 15 uT B1 limit. The repeated-echo train was
simulated for ETL 17 at ESP 15.2857 ms. The bundle is simulation-selected; it
is not scanner, profile or production approved, and the sequence definitions
record these flags.

- **MB2 TSE.** The generator loads the bundle at line 91 and accepts only this
  file: it checks the required variables (`requiredMbFields`, line 274), the
  candidate tag and reference-array checksum (lines 330–340) and the status
  flags that follow. To use another MB2 pulse, generate a MAT file with the
  same variables and conventions, point `mbPulseMatFile` to it, and update
  the expected tag and checksum deliberately. Product-sequence RF is not
  distributed.
- **MB1 TSE.** The default `rfDesign = 'conventional_sinc'` (line 76) builds
  2.5 ms excitation and 3.1 ms refocusing sinc pulses with Pulseq and needs no
  file. The alternative `rfDesign = 'sigpy_9ms'` needs a single-band template
  MAT that is not distributed.
- **Calibration.** The FLASH pulses are Pulseq sinc pulses. The calibration
  reads slice positions and SMS grouping from the same TBW7 bundle (line 37)
  and checks its tag and checksum (lines 201–205), so CSM slices and TSE
  slices share one authoritative geometry.

## Main parameters

Line numbers refer to the current revision.

### `seq/tse_2d_wave_mb2v_v3.m`

| Line | Parameter | Default | Notes |
|---:|---|---|---|
| 99 | `write_v141_format` | `false` | also write v1.4.1 |
| 100 | `writeSequence` | `true` | `false` runs all checks without writing |
| 102 | `execute_system_check_code` | `false` | PNS/CNS and forbidden-frequency checks |
| 104 | `sys_type` | `'skyra'` | `'prisma'`, `'skyra'`, `'TerraX'`, ... (hardware limits and the expected `.asc`) |
| 108 | `ETL` | `17` | screened prefix 8–17 |
| 109 | `TE1Requested` | `107e-3/7` | echo spacing [s]; 15.2857 ms is the RF-simulated value |
| 110 | `TEeffRequested` | `107e-3` | snapped to the nearest echo |
| 111 | `allowUnvalidatedEchoSpacing` | `true` | permits an ESP outside the RF simulation (timing checks only) |
| 116 | `minimumEchoSpacing` | `14.89e-3` | shortest accepted echo spacing [s]; see [Echo spacing](#echo-spacing) |
| 135 | `TRRequested` | `15000e-3` | raised to the minimum when infeasible |
| 176 | `Ry` | `1` | in-plane acceleration (1–4) |
| 181 | `isUseCaipi` | `true` | balanced slice-axis CAIPI blips |
| 182 | `caipiShiftFractionFullFOV` | `0.5` | 1/2, 1/4, 1/6, 1/8 for Ry 1–4 |
| 187 | `isUseWave_sin` | `false` | y (PE) sine Wave |
| 188 | `isUseWave_cos` | `false` | z (slice) cosine Wave |
| 189 | `gwave_max` | `10` | Wave amplitude [mT/m]; slew-limited at 140 T/m/s |
| 191 | `Ncycles` | `20` | 10 or 20 Wave cycles per readout |

### `seq/tse_2d_wave_mb1.m`

| Line | Parameter | Default | Notes |
|---:|---|---|---|
| 67 | `write_v141_format` | `false` | |
| 68 | `writeSequence` | `true` | |
| 69 | `execute_system_check_code` | `1` | |
| 70 | `sys_type` | `'skyra'` | |
| 76 | `rfDesign` | `'conventional_sinc'` | see RF pulses |
| 80–82 | `ETL`, `TE1Requested`, `TEeffRequested` | `17`, `107e-3/7`, `107e-3` | |
| 86 | `minimumEchoSpacing` | `14.89e-3` | same floor as the MB2 generator |
| 116 | `TRRequested` | `8000e-3` | raised to the minimum when infeasible |
| 169 | `Ry` | `1` | |
| 176 | `Ncycles` | `20` | |
| 177 | `isUseWave_sin` | `false` | MB1 supports the y-sine Wave only |

### `seq/gre_mixed_dim_mb2_calibration.m`

| Line | Parameter | Default | Notes |
|---:|---|---|---|
| 52 | `write_v141_format` | `false` | |
| 53 | `writeSequence` | `true` | |
| 54 | `isUseWaveCalibration` | `true` | PSF + CSM; `false` acquires the CSM only |
| 56 | `execute_system_check_code` | `true` | PNS/CNS above 0.95 is an error |
| 59 | `sys_type` | `'prisma'` | |
| 62–63 | `alphaPsf`, `alphaCsm` | `7`, `20` | flip angles [deg] |
| 80–85 | `Ncalib_y`, `Ncalib_z`, `Nacs_y` | `72`, `80`, `32` | PSF and CSM line counts |
| 98–100 | `gwave_max`, `Ncycles` | `10`, `20` | must match the Wave TSE |

The calibration readout, ADC, y-sine, z-cosine and bridge timing reproduce
`tse_2d_wave_mb2v_v3.m`. Use a calibration exported with the same `Ncycles`
and Wave amplitude as the Wave TSE it calibrates.

### Echo spacing

Both TSE generators accept an echo spacing `TE1Requested` only at or above
`minimumEchoSpacing` = 14.89 ms. This floor is the shortest spacing that the
MB2 timing allows for every Wave mode on every supported system. Skyra with y
and z Wave sets it: the low-PNS readout spoiler takes 1.79 ms, and a 0.37 ms
cosine bridge sits on each side of the readout. Without Wave the MB2 timing
alone allows 14.15 ms on Skyra, the Prisma limits allow shorter spacings, and
the MB1 timing allows 13.52 ms. The shared floor keeps every export matchable
across Wave modes and between MB1 and MB2.

The default 15.2857 ms (`107e-3/7`, effective TE 107 ms at echo 7) is the only
spacing included in the RF echo-train simulation of the TBW7 bundle. Other
spacings pass the sequence timing checks only; the MB2 generator then prints
a warning and records
`EchoSpacingValidationMode = sequence_timing_only_not_rf_simulated`.

## System checks

With `execute_system_check_code` true, the generators run `seq.calcPNS`
(PNS, and CNS on CimaX) and `forbiddenFreqCheck` on the gradient model given by
`system_asc_file`. The file name must match the generator's `sys_type` table,
for example `MP_GPA_K2309_2250V_951A_AS82.asc` for `prisma` and
`MP_GPA_K2309_2250V_793A_GC99.asc` for `skyra`, and
`safe_pns_prediction_path` must be set; otherwise generation stops with an
error rather than skipping the checks. The TSE generators warn above
PNS/CNS 0.95, the calibration stops. Forbidden-frequency results are
warnings. The TSE definitions record `ExecuteSystemCheckCode`.

Passing these checks does not make an export scanner-safe; RF power, SAR and
local scanner checks remain mandatory.

## Outputs

| Generator | Files |
|---|---|
| MB2 | `<out_path>/v151/<name>.seq`, `<out_path>/metadata/<name>_PEorder.mat`, optional `<out_path>/v141/<name>_v141.seq` |
| MB1 | `<out_path>/v151/<name>.seq`, `<out_path>/<name>_PEorder.mat`, optional `<out_path>/v141/<name>_v141.seq` |
| Calibration | `<out_path>/<name>.seq`, `<out_path>/<name>_geometry_labels.mat`, `<out_path>/<name>_geometry_labels.txt` |

The `_PEorder.mat` sidecars hold the PE tables, SMS grouping and mappings,
CAIPI phase tables and RF provenance used by the reconstruction; keep each
sidecar with its `.seq` and the acquired TWIX data.

File names encode the main parameters:

- MB2: `tse_m2b7v1_B15_cF2_wYZ_N20_TRA_440x220x50_F220_R1_os3_E17_ES15p29_TE107_TR15000_skyra`
  — `B` design B1 limit [uT]; `cF<n>` CAIPI shift of FOV/n of the full PE
  FOV (`noC` without CAIPI); `wYZ`, `wY`, `wZ` active Wave axes; `N` Wave
  cycles; matrix and FOV [mm]; `R` acceleration (`R3p73` is Ry 3 with the
  complete 73-line grid); readout oversampling; `E` ETL; `ES` echo spacing
  and `TE` effective TE [ms]; `TR` realized TR [ms]; `sys_type`.
- MB1: `tse_mb1_convEx2p5R3p1_Ncyc20_TRA_440x220x50_FOV220_Ry1_os3_ETL17_TE15p29_TR13666_skyra`
  — RF design, optional `wY`, Wave cycles, matrix, FOV [mm], `Ry`, readout
  oversampling, ETL, echo spacing (`TE` tag) [ms], realized TR [ms], `sys_type`.
- Calibration: `gre_mixed_dim_mb2_calib_waveYZ3dproj_csm2d_Ncyc20_440x220x80_FOV220_Nacs32_os3_prisma`
  — `waveYZ3dproj_csm2d` (PSF + CSM) or `csm2donly`, Wave cycles, PSF matrix
  `Nx x Ny x Ncalib_z`, FOV [mm], CSM lines, readout oversampling, `sys_type`.

## Reproducing the reference exports

The generators were imported unchanged except for their path handling and
the echo-spacing floor, which leaves exports at the default spacing
unchanged. Regenerating with the parameters of an acquired export is
therefore expected to reproduce it; the comparison below tests this. Starting
from the committed defaults, edit only the listed lines:

| Reference export | Edits | System checks |
|---|---|---|
| MB2 `tse_m2b7v1_B15_cF2_N20_..._R1_..._TR15000_skyra` | none | off |
| MB2 `tse_m2b7v1_B15_cF2_N20_..._R1_..._TR12000_skyra` | line 135 `TRRequested = 12000e-3;` | off |
| MB2 `tse_m2b7v1_B15_cF2_wYZ_N20_..._R1_..._TR15000_skyra` | lines 187–188 `isUseWave_sin = true;` `isUseWave_cos = true;` | off |
| MB2 `tse_m2b7v1_B15_cF2_wYZ_N20_..._R1_..._TR12000_skyra` | lines 187–188 as above; line 135 `TRRequested = 12000e-3;` | off |
| MB2 `tse_m2b7v1_B15_cF4_wYZ_N20_..._R2_..._TR12000_skyra` | line 176 `Ry = 2;` line 182 `caipiShiftFractionFullFOV = 0.25;` lines 187–188 as above; line 135 `TRRequested = 12000e-3;` | off |
| MB1 `tse_mb1_convEx2p5R3p1_Ncyc20_..._Ry1_..._TR13666_skyra` | none (TR 8 s is raised to 13.666 s) | on, skyra `.asc` |
| Calibration `..._Nacs32_os3_prisma` | none | on, prisma `.asc` |
| Calibration `..._Nacs32_os3_skyra` | line 59 `sys_type = 'skyra';` | on, skyra `.asc` |

Then compare each export folder with the reference folder that holds the
same file names:

```bash
python evaluation/compare_sequence_exports.py \
  --candidate-dir /path/to/regenerated/out_path \
  --reference-dir /path/to/reference/protocol_folder \
  --report /path/to/comparison_report.json
```

The comparison is exact by default. Two differences are expected and listed
separately: the TSE `OutputPath` definition (and therefore the `[SIGNATURE]`
hash) and the sidecar fields that store paths of the generating machine
(`mbPulseMatFile`, `out_path`, `calibrationMetadata.OutputPath`,
`calibrationMetadata.MBGeometryMatFile`). Every other line of every `.seq`
section and every sidecar variable must agree. When a section differs only in
number formatting, for example after a Pulseq update, the report marks it as
numerically equal.

## Scanner protocol notes

- Prescribe the calibration and the TSE scans back to back with the same
  transverse plane, FOV center, rotation, readout and phase-encoding
  directions, and keep the receive coil, active channels and subject position
  unchanged. Matched array sizes alone do not prove matched geometry; the
  reconstruction validates positions, trajectories and orientation.
- The 50 slices are contiguous 3 mm slices with centers from -73.5 mm to
  +73.5 mm; MB2 excites physical slices `p` and `p+25` together.
- Keep each exported `.seq` and its sidecar with the TWIX data acquired with
  it; reconstruction needs the exact pair of every measurement.
