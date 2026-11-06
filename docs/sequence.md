# Sequence generation

## Entry points

| Generator | Sequence | RF | Writes to `out_path` |
|---|---|---|---|
| `seq/tse_2d_wave_mb2v_v3.m` | MB2 TSE: 50 physical slices in 25 SMS groups, blipped CAIPI, two-axis Wave (on by default), followed by the [appended FLASH calibration](#appended-flash-calibration) | TBW7 beta150 vMB-VERSE bundle `seq/rf_pulses/mb2_tse_ls_tbw7_beta150_v1verse_B1max15p000uT.mat` | `v151/tse_m2b7v1_*.seq`, `metadata/<name>_PEorder.mat` |
| `seq/tse_2d_wave_mb1.m` | MB1 TSE: 50 single-band slices, y-sine Wave (on by default), followed by the appended FLASH calibration | Pulseq sinc pulses (`rfDesign = 'conventional_sinc'`) | `v151/tse_mb1_*.seq`, `<name>_PEorder.mat` |
| `seq/gre_mixed_dim_mb2_calibration.m` | Separate mixed-dimension spoiled-GRE calibration: 3D y/z Wave PSF projections plus 2D slice-wise CSM lines for the same 50 slices | Pulseq sinc pulses; the TBW7 bundle supplies slice geometry only | `gre_mixed_dim_mb2_calib_*.seq`, `<name>_geometry_labels.mat`, `<name>_geometry_labels.txt` |

By default each TSE generator writes one integrated measurement: the TSE and
its 3D slab FLASH calibration in the same `.seq` file, so a session needs a
single scan and the reconstruction a single TWIX file. With
`isAppendCalibration = false` the generators write the TSE alone; such a
scan is reconstructed with a separately acquired
`gre_mixed_dim_mb2_calibration.m` export of the same session. The mixed
calibration also serves to compare its direct-slice coil maps with those of
the appended calibration.

MATLAB is the authoritative sequence implementation. The raw-data storage of
both calibrations is described in [Calibration storage](calibration_storage.md).

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
generators read no environment variables. The scanner system is the one
exception, see [Scanner system](#scanner-system).

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

### Scanner system

After the path settings, each generator asks for the scanner system
(`sys_type`), as the Wave-MPRAGE and Wave-GRE generators do:

```text
Select scanner/system name:
   1) prisma
   2) skyra
   ...
Select 1-9 [default 2: skyra]:
```

Enter a menu number or a name; pressing Enter keeps `sysTypeDefault`
(`'skyra'` in all three generators). The selection sets the gradient and
slew limits, the gradient model expected by the [system checks](#system-checks)
and the `_<sys_type>` suffix of the file names. Select the same system for the
calibration and the TSE scans of a session.

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

- **MB2 TSE.** The generator loads the bundle at line 105 and accepts only this
  file: it checks the required variables (`requiredMbFields`, line 307), the
  candidate tag and reference-array checksum (lines 363–373) and the status
  flags that follow. To use another MB2 pulse, generate a MAT file with the
  same variables and conventions, point `mbPulseMatFile` to it, and update
  the expected tag and checksum deliberately. Product-sequence RF is not
  distributed.
- **MB1 TSE.** The default `rfDesign = 'conventional_sinc'` (line 97) builds
  2.5 ms excitation and 3.1 ms refocusing sinc pulses with Pulseq and needs no
  file. The alternative `rfDesign = 'sigpy_9ms'` needs a single-band template
  MAT that is not distributed.
- **Calibration.** The FLASH pulses are Pulseq sinc pulses. The calibration
  reads slice positions and SMS grouping from the same TBW7 bundle (line 42)
  and checks its tag and checksum (lines 214–218), so CSM slices and TSE
  slices share one authoritative geometry.

## Main parameters

Line numbers refer to the current revision.

### `seq/tse_2d_wave_mb2v_v3.m`

| Line | Parameter | Default | Notes |
|---:|---|---|---|
| 113 | `write_v141_format` | `false` | also write v1.4.1 |
| 114 | `writeSequence` | `true` | `false` runs all checks without writing |
| 116 | `execute_system_check_code` | `true` | PNS/CNS and forbidden-frequency checks; needs `system_asc_file` and safe_pns_prediction |
| 121–124 | `sysTypeDefault`, `sys_type` | `'skyra'` | asked at run time, see [Scanner system](#scanner-system) |
| 129 | `ETL` | `17` | screened prefix 8–17 |
| 130 | `TE1Requested` | `107e-3/7` | echo spacing [s]; 15.2857 ms is the RF-simulated value |
| 131 | `TEeffRequested` | `107e-3` | snapped to the nearest echo |
| 132 | `allowUnvalidatedEchoSpacing` | `true` | permits an ESP outside the RF simulation (timing checks only) |
| 137 | `minimumEchoSpacing` | `14.89e-3` | shortest accepted echo spacing [s]; see [Echo spacing](#echo-spacing) |
| 156 | `TRRequested` | `15000e-3` | raised to the minimum when infeasible |
| 197 | `Ry` | `1` | in-plane acceleration (1–4) |
| 202 | `isUseCaipi` | `true` | balanced slice-axis CAIPI blips |
| 203 | `caipiShiftFractionFullFOV` | `0.5` | 1/2, 1/4, 1/6, 1/8 for Ry 1–4 |
| 210 | `isUseWave_sin` | `true` | y (PE) sine Wave |
| 211 | `isUseWave_cos` | `true` | z (slice) cosine Wave; set both to `false` for a Cartesian export |
| 212 | `gwave_max` | `10` | Wave amplitude [mT/m]; slew-limited at 140 T/m/s |
| 214 | `Ncycles` | `20` | 10 or 20 Wave cycles per readout |
| 228 | `isAppendCalibration` | `true` | append the [3D FLASH calibration](#appended-flash-calibration); `false` writes the TSE alone |
| 229–230 | `calibrationAcsSize`, `calibrationAcsFovZ` | `[32 32]`, `180e-3` | ky x kz ACS lines and ACS kz FOV [m] |

The reconstruction accepts MB2 data with both Wave axes or with neither.

With CAIPI, kz at the ADCs is not zero by design: it alternates between 0
and the blip moment, -1/(2 x 75 mm) = -6.67 1/m for the default FOV/2 shift
and the 75 mm band separation, so a k-space viewer shows every other echo
slightly off the ky axis. The blip is rewound after every readout, and the
slice-axis moment before each refocusing pulse is the same in every echo
interval. The generator prints these values under `Slice-axis kz diagnostics`
and stops when the kz residual left after removing the designed CAIPI moment
exceeds 1% of the PE k-space step (0.045 1/m). With `isUseCaipi = false` and
`caipiShiftFractionFullFOV = 0`, kz is zero at every ADC.

### `seq/tse_2d_wave_mb1.m`

| Line | Parameter | Default | Notes |
|---:|---|---|---|
| 81 | `write_v141_format` | `false` | |
| 82 | `writeSequence` | `true` | |
| 83 | `execute_system_check_code` | `1` | |
| 87–90 | `sysTypeDefault`, `sys_type` | `'skyra'` | asked at run time, see [Scanner system](#scanner-system) |
| 97 | `rfDesign` | `'conventional_sinc'` | see RF pulses |
| 101–103 | `ETL`, `TE1Requested`, `TEeffRequested` | `17`, `107e-3/7`, `107e-3` | |
| 107 | `minimumEchoSpacing` | `14.89e-3` | same floor as the MB2 generator |
| 137 | `TRRequested` | `8000e-3` | raised to the minimum when infeasible |
| 190 | `Ry` | `1` | the reconstruction needs the complete grid of `Ry = 1` |
| 197 | `Ncycles` | `20` | |
| 201 | `isUseWave_sin` | `true` | y-sine Wave; `false` for a Cartesian export. MB1 has no z Wave |
| 215 | `isAppendCalibration` | `true` | append the 3D FLASH calibration |
| 216–217 | `calibrationAcsSize`, `calibrationAcsFovZ` | `[32 32]`, `180e-3` | as for MB2 |

### `seq/gre_mixed_dim_mb2_calibration.m`

| Line | Parameter | Default | Notes |
|---:|---|---|---|
| 57 | `write_v141_format` | `false` | |
| 58 | `writeSequence` | `true` | |
| 59 | `isUseWaveCalibration` | `true` | PSF + CSM; `false` acquires the CSM only |
| 61 | `execute_system_check_code` | `true` | PNS/CNS above 0.95 is an error |
| 68–71 | `sysTypeDefault`, `sys_type` | `'skyra'` | asked at run time, see [Scanner system](#scanner-system) |
| 75–76 | `alphaPsf`, `alphaCsm` | `7`, `20` | flip angles [deg] |
| 93–98 | `Ncalib_y`, `Ncalib_z`, `Nacs_y` | `72`, `80`, `32` | PSF and CSM line counts |
| 111–113 | `gwave_max`, `Ncycles` | `10`, `20` | must match the Wave TSE |

The calibration readout, ADC, y-sine, z-cosine and bridge timing reproduce
`tse_2d_wave_mb2v_v3.m`. Use a calibration exported with the same `Ncycles`
and Wave amplitude as the Wave TSE it calibrates.

### Appended FLASH calibration

With `isAppendCalibration = true`, `seq/utils/designFlashCalibration.m`
designs a 3D slab FLASH calibration from the TSE's own ADC, readout and Wave
events, validates it on its own sequence object, and
`seq/utils/appendFlashCalibration.m` appends it after the last TSE train.
All parts share one RF-spoiled steady state: a 1 ms, TBW 20 sinc pulse with a
7 degree flip angle that excites a 160 mm slab over the 150 mm slice stack,
100 dummy shots before the first part and 10 settling shots before each
further part.

| Part | Content | MB2 | MB1 |
|---:|---|---|---|
| 0, 1 | no-Wave and y-sine Wave projections, 72 ky at kz = 0 | with Wave | with Wave |
| 2 | no-Wave z projection, 80 kz at ky = 0 over 160 mm | always | no |
| 3 | z-cosine Wave projection on the lines of part 2 | with Wave | no |
| 4 | no-Wave 3D ACS, 32 ky x 32 kz, kz FOV 180 mm | always | always |

Parts 0 to 3 reproduce the PSF projections of
`gre_mixed_dim_mb2_calibration.m`. The ACS replaces its 2D slice-wise CSM:
the reconstruction estimates the coil maps by ESPIRiT in (ky, kz) at every
readout position on a 3 mm z grid whose planes are the slice centers (the
grid is shifted by half a plane for the even slice count). The kz FOV covers
the excited slab, and the 32 kz lines give the maps about 5.6 mm resolution
along z. Part 2 lets the reconstruction verify the z orientation of the ACS
planes on every MB2 scan; for MB1 the TSE slices themselves verify it.

The generator states the calibration TE, TR, ADC count and duration when it
runs, about 19 s for MB2 with Wave and 15 s for MB1 with Wave (a TR of about
12 ms). Every calibration ADC is a REF readout with a compact LIN, see
[Calibration storage](calibration_storage.md#appended-calibration-appended_ref_lin_v1).
The sidecar `calibration` struct holds the tables the reconstruction needs.

Before scanning, check an export with

```bash
python evaluation/check_appended_calibration.py \
  --seq /path/to/out_path/v151/tse_..._cal_skyra.seq
```

It runs the sequence-level checks of the reconstruction audit on the `.seq`
file and its sidecar (signature, TSE/calibration ADC split and labels, Wave
pairs, ACS ky/kz grid, z planes at the slice centers) and prints the ADC
counts and the calibration timing. In a k-space viewer the calibration ADCs
follow the TSE as the ky and kz projection lines and the 32 x 32 ky-kz ACS
grid.

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
`system_asc_file`. The model must be the one the generator's `sys_type`
table names, for example `MP_GPA_K2309_2250V_951A_AS82.asc` for `prisma` and
`MP_GPA_K2309_2250V_793A_GC99.asc` for `skyra`. When the saved file belongs
to another system, the generator uses the model of the selected system from
the same folder and prints a note; keep the `.asc` files of all your systems
in one folder. `safe_pns_prediction_path` must be set. Without the matching
model or the submodule, generation stops with an error rather than skipping
the checks. The TSE generators warn above
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
CAIPI phase tables and RF provenance used by the reconstruction and, with an
appended calibration, its tables (`calibration`) and the sequence signature
(`sequenceMetadata`); keep each sidecar with its `.seq` and the acquired TWIX
data.

File names encode the main parameters:

- MB2: `tse_m2b7v1_B15_cF2_wYZ_N20_TRA_440x220x50_F220_R1_os3_E17_ES15p29_TE107_TR15000_skyra`
  — `B` design B1 limit [uT]; `cF<n>` CAIPI shift of FOV/n of the full PE
  FOV (`noC` without CAIPI); `wYZ`, `wY`, `wZ` active Wave axes; `N` Wave
  cycles; matrix and FOV [mm]; `R` acceleration (`R3p73` is Ry 3 with the
  complete 73-line grid); readout oversampling; `E` ETL; `ES` echo spacing
  and `TE` effective TE [ms]; `TR` realized TR [ms]; `_cal` with the appended
  calibration; `sys_type`. The default export is
  `tse_m2b7v1_B15_cF2_wYZ_N20_TRA_440x220x50_F220_R1_os3_E17_ES15p29_TE107_TR15000_cal_skyra`.
- MB1: `tse_mb1_convEx2p5R3p1_Ncyc20_TRA_440x220x50_FOV220_Ry1_os3_ETL17_TE15p29_TR13666_skyra`
  — RF design, optional `wY`, Wave cycles, matrix, FOV [mm], `Ry`, readout
  oversampling, ETL, echo spacing (`TE` tag) [ms], realized TR [ms], `_cal`
  with the appended calibration, `sys_type`.
- Calibration: `gre_mixed_dim_mb2_calib_waveYZ3dproj_csm2d_Ncyc20_440x220x80_FOV220_Nacs32_os3_prisma`
  — `waveYZ3dproj_csm2d` (PSF + CSM) or `csm2donly`, Wave cycles, PSF matrix
  `Nx x Ny x Ncalib_z`, FOV [mm], CSM lines, readout oversampling, `sys_type`.

## Reproducing the reference exports

The generators were imported unchanged except for their path handling, the
system prompt, the Wave defaults, the appended calibration and the
echo-spacing floor. With `isAppendCalibration = false` none of these changes
the exported TSE or its sidecar for the same parameters, so regenerating with
the parameters of an acquired export is expected to reproduce it; the
comparison below tests this. The scanned reference exports have no appended
calibration. Starting from the committed defaults, edit only the listed
lines, and at the system prompt select the system named at the end of the
reference file name:

| Reference export | Edits | System checks |
|---|---|---|
| MB2 `tse_m2b7v1_B15_cF2_N20_..._R1_..._TR15000_skyra` | line 228 `isAppendCalibration = false;` line 116 `execute_system_check_code = false;` lines 210–211 `isUseWave_sin = false;` `isUseWave_cos = false;` | off |
| MB2 `tse_m2b7v1_B15_cF2_N20_..._R1_..._TR12000_skyra` | lines 228, 116, 210–211 as above; line 156 `TRRequested = 12000e-3;` | off |
| MB2 `tse_m2b7v1_B15_cF2_wYZ_N20_..._R1_..._TR15000_skyra` | line 228 `isAppendCalibration = false;` line 116 `execute_system_check_code = false;` | off |
| MB2 `tse_m2b7v1_B15_cF2_wYZ_N20_..._R1_..._TR12000_skyra` | lines 228, 116 as above; line 156 `TRRequested = 12000e-3;` | off |
| MB2 `tse_m2b7v1_B15_cF4_wYZ_N20_..._R2_..._TR12000_skyra` | lines 228, 116 as above; line 197 `Ry = 2;` line 203 `caipiShiftFractionFullFOV = 0.25;` line 156 `TRRequested = 12000e-3;` | off |
| MB1 `tse_mb1_convEx2p5R3p1_Ncyc20_..._Ry1_..._TR13666_skyra` | line 215 `isAppendCalibration = false;` line 201 `isUseWave_sin = false;` (TR 8 s is raised to 13.666 s) | on, skyra `.asc` |
| Calibration `..._Nacs32_os3_prisma` | none (select `prisma`) | on, prisma `.asc` |
| Calibration `..._Nacs32_os3_skyra` | none | on, skyra `.asc` |

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

- With the appended calibration, one measurement holds the TSE and its
  calibration; no separate calibration scan is needed. The calibration runs
  after the last TSE train and adds about 15 to 19 s.
- A separate mixed calibration (for a TSE without the appended calibration,
  or to compare the two coil-map methods) must be prescribed back to back
  with the TSE with the same transverse plane, FOV center, rotation, readout
  and phase-encoding directions, and the receive coil, active channels and
  subject position must stay unchanged. Matched array sizes alone do not
  prove matched geometry; the reconstruction validates positions,
  trajectories and orientation.
- The 50 slices are contiguous 3 mm slices with centers from -73.5 mm to
  +73.5 mm; MB2 excites physical slices `p` and `p+25` together.
- Keep each exported `.seq` and its sidecar with the TWIX data acquired with
  it; reconstruction needs the exact pair of every measurement.
