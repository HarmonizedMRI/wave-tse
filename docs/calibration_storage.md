# Mixed calibration storage: `compact_adc_lin_v1`

This contract applies to `seq/gre_mixed_dim_mb2_calibration.m` (introduced
2026-10-14). Calibration exports written before it, without
`storageLayoutVersion`, keep the older SET/PAR/SLC layout described at the
end.

## Why the scanner counters are compact

The Wave calibration encodes y projection samples, z projection samples,
physical CSM slices and scientific parts. Storing them in LIN, PAR, SLC and
SET would make a rectangular reader allocate the product of the independent
maxima: the PSF alone has 304 acquired readouts but would occupy
72 x 80 x 4 = 23,040 slots, 98.68% empty, and a CSM-only export would leave
four leading empty SET positions. Reducing the receive-channel count can hide
that allocation without fixing it.

## Scanner storage contract

Every acquired ADC uses one unique zero-based LIN, in acquisition order.
PAR, SLC, SET and ECO are always zero. Dummy shots have no ADC. CSM readouts
keep REF=true and PSF readouts REF=false, so the CSM arrives in the TWIX
`refscan` stream and the PSF in the `image` stream. RF spoiling, RF and
gradient waveforms, physical slice locations, Wave amplitude, timing and shot
order are not affected by the storage layout.

| Mode | Acquired ADCs | LIN range | Other numeric counters |
|---|---:|---|---|
| Wave PSF + CSM (`isUseWaveCalibration = true`, default) | 1,904 | 0 through 1,903 | PAR=SLC=SET=ECO=0 |
| CSM only (`isUseWaveCalibration = false`) | 1,600 | 0 through 1,599 | PAR=SLC=SET=ECO=0 |

These are the default 50-slice / 72-y / 80-z / 32-CSM-line counts. Counts and
mappings are derived from the actual configuration, not from filename
numbers. The scientific `acqTable.setID` and `acqTable.slcLabel` fields remain
in the sidecar for compatibility; they no longer equal scanner SET or SLC.
`acqTable.rawLin0` is the scanner LIN value.

The scanner SLC counter does not report physical slice identity. Physical CSM
slices must be recovered from the sidecar. Never FFT or interpolate across the
CSM slice dimension.

## Sidecar maps and index bases

The matching `<name>_geometry_labels.mat` is required for reconstruction. It
contains the scientific variables of the calibration and:

- `storageLayoutVersion`: `compact_adc_lin_v1`.
- `acqTable`: one record per acquired ADC, including physical ky/kz indices,
  RF mode, scientific part ID, physical slice index/position, REF routing and
  zero-based `rawLin0`.
- `psfAdcIndices`: 4 x max(Ncalib_y, Ncalib_z), normally 4 x 80. Entries are
  **one-based global ADC-record indices**; zero marks an unused entry.
- `psfValidMask`, `psfValidLengths`: normally [72 72 80 80]. In CSM-only mode
  the PSF map is all zero and the lengths are [0 0 0 0].
- `csmAdcIndices`: Nslices x Nacs_y, normally 50 x 32. Rows are ascending
  physical slice index; columns are ascending local CSM ky. Entries are
  one-based global ADC-record indices, with no missing entries.
- `imageAdcIndices`, `refAdcIndices`: global record indices of the PSF/image
  and CSM/reference streams.
- `partStart`, `partStop`: one-based global ADC bounds of each enabled
  scientific part, matching `parts` and `acqTable`.
- `storageLabelExtents`: extents in [SET SLC LIN PAR] order. Their product
  equals the acquired ADC count, which rules out Cartesian-product padding.
- `calibrationMetadata`: scientific definitions kept out of the `.seq`,
  including geometry, coordinates, mappings, RF parameters, timing, Wave
  design and emitted amplitude, system-check settings, output location and
  target-TSE provenance. The sequence md5 signature and signature type bind
  the sidecar to its export.

PSF row meanings:

| Row (one-based) | Scientific part ID | Encoding | Valid columns |
|---|---:|---|---:|
| 1 | 0 | y reference, no Wave | Ncalib_y |
| 2 | 1 | y sine Wave | Ncalib_y |
| 3 | 2 | z reference, no Wave | Ncalib_z |
| 4 | 3 | z cosine Wave | Ncalib_z |

Unused PSF positions exist only in this small reconstruction table; the
sequence does not acquire extra ADCs for them. All physical ky/kz coordinates
come from `acqTable`, `parts` and the saved phase-area arrays. The CSM RF is a
normal slice-selective FLASH sinc pulse; the TSE RF MAT supplies geometry and
provenance only. `calibrationMetadata.TargetTseRF*` fields describe the
target TSE, not the emitted FLASH RF.

## Reading the data

Do not slice arrays by SET/PAR/SLC and do not infer physical geometry from
the scanner storage counters. With MATLAB mapVBVD, load image and reference
readouts in acquisition order with `.unsorted()` and place them by their
stored LIN counters (MATLAB mapVBVD `.Lin` already includes +1):

```matlab
meta = load(matchingGeometryLabelsMat);
assert(strcmp(meta.storageLayoutVersion, 'compact_adc_lin_v1'));
raw = complex(zeros(nSamples, nCoils, meta.expectedAdcCount, 'single'));
% For each present stream (image, refscan):
% stream = twix.image or twix.refscan;
% raw(:,:,double(stream.Lin)) = stream.unsorted();
% Verify stream counts, unique LIN values and complete global coverage.

psf = complex(zeros(nSamples, nCoils, size(meta.psfAdcIndices,1), ...
    size(meta.psfAdcIndices,2), 'single'));
for row = 1:size(meta.psfAdcIndices,1)
    for col = find(meta.psfValidMask(row,:))
        psf(:,:,row,col) = raw(:,:,meta.psfAdcIndices(row,col));
    end
end
csm = reshape(raw(:,:,meta.csmAdcIndices(:)), ...
    nSamples, nCoils, size(meta.csmAdcIndices,1), size(meta.csmAdcIndices,2));
% psf: [ADC sample, coil, projection row, projection column]
% csm: [ADC sample, coil, physical slice, local ky]
```

Python pymapvbvd reports zero-based LIN values. Use the MDH LIN values rather
than positions in a sorted or padded array: REF readouts start at global
LIN=304 in the default combined acquisition, and some readers remove leading
empty LIN positions from the returned reference array. The stored counters
still provide the global mapping. CSM-only data have no PSF/image readouts.

Historical calibration files without `storageLayoutVersion` use the old
SET/PAR/SLC contract and need a separate loading branch.

## `.seq` definition allowlist

The calibration source sets only `Name` and `FOV` explicitly. FOV is the
physical envelope [220, 220, 160] mm at the default prescription, not a dense
storage matrix. Pulseq adds its four raster definitions and `TotalDuration`.
No matrix, slice-count or position arrays, scientific part IDs, coordinate
grids, paths, provenance, Wave diagnostics or descriptive text are sent as
scanner definitions; all of them are in the sidecar.

## Generation

`isUseWaveCalibration = true` (the default) acquires PSF + CSM; `false`
acquires the CSM only. `writeSequence = false` runs the diagnostics and builds
the metadata without writing files. Generation refuses to overwrite an
existing `.seq`, metadata MAT, report TXT or requested legacy v1.4.1 export;
choose a new `out_path` in the path settings instead.

## Validation record

MATLAB R2024b validation on 2026-10-14 used isolated temporary exports:

- Wave-enabled mode: 1,904 acquired ADCs; all 18,272 physical blocks exactly
  match the previous layout after removing label fields.
- CSM-only mode: 1,600 acquired ADCs; all 14,800 physical blocks exactly
  match.
- Pulseq timing, ADC label evolution, RF-center TE, shot spacing and paired
  y/z Wave trajectory checks pass. Every sidecar map is a bijection over the
  acquired records, and the scientific acquisition table matches the earlier
  export field for field, with only `rawLin0` added.
- An independent Python/SciPy readback verified the serialized ADC labels,
  sidecar mappings, [72 72 80 80] / [0 0 0 0] valid PSF lengths and the
  seven-key scanner definition allowlist.
- With system checks enabled, Skyra PNS remained within the source's 0.95
  threshold; CNS is not evaluated on Skyra.
- The forbidden-frequency check warned in the 585 Hz (100 Hz width) and
  1,120 Hz (220 Hz width) bands: band RMS relative to the configured maxGrad
  was 7.46% / 3.47% in combined mode and 7.65% / 3.96% in CSM-only mode,
  above the utility's 1% heuristic. The ASC file supplies no explicit absolute
  band-energy limit. Pulseq also warns about shape restoration of the
  arbitrary cosine gradient.

RF power and SAR do not depend on the storage layout. Compact exports were
subsequently acquired and reconstructed (TerraX 2026-10-14, Prisma 2026-10-15
and later TerraX sessions), each with 304 PSF `image` and 1,600 CSM `refscan`
records.
