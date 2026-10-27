# Troubleshooting

## Sequence generation

**`ESP ... is below the 14.89 ms floor`.** `TE1Requested` is shorter than the
shared echo-spacing floor of the TSE generators; see
[Echo spacing](sequence.md#echo-spacing).

**`... was not included in the RF echo-train simulation` (warning).** Any
echo spacing other than 15.2857 ms passes the sequence timing checks only.

**A system check stops generation.** With `execute_system_check_code` true,
`system_asc_file` must name the gradient model of the selected `sys_type` and
`safe_pns_prediction_path` must point to the
`external/safe_pns_prediction` submodule. Run
`git submodule update --init external/safe_pns_prediction` after cloning.

**The calibration refuses to write.** It never overwrites an export; choose a
new `out_path`.

**`compare_sequence_exports.py` reports differences.** Only `OutputPath`, the
signature hash and the sidecar path fields may differ. Check that the
generator parameters match the reference export (see the reproduction table
in [Sequence generation](sequence.md#reproducing-the-reference-exports)) and
that the same Pulseq version was used; formatting-only differences are marked
as numerically equal.

## Reconstruction

**`Existing provenance differs: .../input_manifest.json`.** The folder holds
results of other inputs or settings. Use a new `--output-dir` (or a new
`--calibration-dir` when the calibration inputs, `--virtual-coils` or the
ESPIRiT settings changed). Wave PSF settings do not require a new calibration
folder; their results are tagged.

**`Recorded sequence ... differs`.** The TWIX header names another `.seq`
file. Pass the exact sequence exported for that acquisition.

**`Only transverse prescriptions are supported`.** The slice-direction and
NIfTI conventions are validated for transverse stacks only.

**`Raw ADC chronology is not the sequence order ...`.** The TWIX file and the
sequence do not belong together, or the file ended early. When a scan was
stopped after its last useful readouts, `--allow-missing-final-adcs N` accepts
up to `N` missing final ADCs, provided the acquired ADCs are an exact prefix of
the sequence order. MB1 reconstruction still needs a complete grid.

**`RF/Gz physical mapping error` or `RF waveform slice translation not
verified`.** The excitation frequencies or RF waveforms do not place the slice
groups at the sidecar positions. Check that the sidecar belongs to the
sequence and that the scan used the same calibration geometry.

**`Scanner MDH slice offsets differ ... for both signs`.** The scanner slice
positions do not follow the mapped groups in order and spacing. The raw SLI
counters, sequence labels and sidecar are inconsistent; do not reconstruct.

**`TSE and calibration trajectories do not match`.** The readout, PE or Wave
encoding of the TSE differs from the calibration, for example a calibration
exported with other Wave cycles or amplitude.

**`Calibration does not establish the negative Fourier-z convention`.** The
calibration's direct slices and Fourier-z projection do not agree; the
calibration acquisition cannot be used.

**`Wave PSF validation failed`.** Inspect
`figures/wave_coefficient_fit_<tag>.png` and
`reports/wave_coefficient_fit_<tag>.json` in the calibration folder. If the
automatic interval rejects a calibration whose central readout is clean,
rerun with a reviewed central interval `--psf-fit-range START STOP` read off
that figure (the research pipeline used `1120 1520` for the 2640-sample
readout of the TerraX 0.25 mm scan); the held-out samples outside the interval
must still pass. A setting that failed once fails again immediately; choose
another setting.

**`MB1/FLASH physical mapping not established`.** The MB1 slices match the
slice-reversed FLASH maps at least as well as the stated ones. The slice
mapping or the calibration geometry is wrong; do not reconstruct.

**`MB2 operator tests failed` or `Group ... failed its finite/convergence
gate`.** The failed report is saved with a time stamp in `reports/`. Increase
`--max-iterations` only if CG converges slowly but steadily; nonfinite values
point to invalid sensitivities or data.

**`Partial checkpoint of group ...`.** A previous run stopped between writing
a group array and its report. Move the incomplete pair aside or use a new
output folder.

**`CSM diagnostics without maps for physical slices ...`.** An ESPIRiT run
stopped between the diagnostics and the map of these slices. Move the listed
`arrays/csm_diagnostics_XX.npz` files aside and rerun.

**ESPIRiT takes hours on the CPU.** SigPy's CPU ESPIRiT uses one thread; use
`--csm-workers N` (N parallel processes) or `--csm-device 0`.

**GPU options.** `--device cuda` needs a CUDA-enabled PyTorch;
`--csm-device 0` needs CuPy for SigPy. Without them, run on the CPU.

### L1-wavelet regularization

**`L1-wavelet regularization is implemented for MB2 scans only`.** MB1 data
are reconstructed exactly from a complete grid; omit
`--l1-threshold-multiple`.

**`The native ... matrix does not allow ... periodized wavelet levels`.**
Both native matrix sizes must be divisible by `2**levels`; the 440 x 220
protocol needs `--l1-wavelet-levels 2`.

**`The lambda-0 prediction holds ...% of the far-margin power`.** The Wave
spread reaches the outer 20% of the oversampled readout, so these margins
cannot measure the noise covariance. L1 needs more readout oversampling or a
smaller Wave spread for this protocol.

**`... too few data per unknown to estimate the noise profile` or `Only ...
sampled PE lines lie in the outer 30% of the grid`.** The lambda-0 residual of
the central group cannot measure the readout noise profile: the fit leaves
less than 20% of the noise somewhere in the FOV (too few sampled lines or
virtual coils per band unknown), or fewer than 8 outer PE lines are sampled.

**`Changed noise profile`, `Partial noise profile`, `Noise-only result of
another profile or group` or `L1 checkpoint of group ... belongs to another
noise profile`.** The folder holds L1 results of another run. Move
`cache/l1_noise_profile.npz`, `reports/l1_noise_profile.json`,
`reports/l1_noise_null_group_XX.json` and the `l1c*` files of `arrays/`,
`reports/`, `nifti/` and `figures/` aside; the lambda-0 results can stay.

**`The synthetic noise solve of the noise profile did not converge` or `The
noise-only solve of group ... did not converge`.** CG on noise data reached
`--max-iterations`. A larger budget changes the input manifest and needs a new
output folder.

**`L1 group ... used the full budget of ... FISTA iterations` (warning).**
The result is kept. If `fista_relative_change_last10` in the group report is
still well above `--l1-tolerance`, rerun with a larger `--l1-max-iterations`,
which writes a new tag.
