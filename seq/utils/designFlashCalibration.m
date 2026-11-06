function cal = designFlashCalibration(cfg)
%DESIGNFLASHCALIBRATION Design and validate the 3D slab FLASH calibration of a TSE.
%
%   cal = designFlashCalibration(cfg)
%
% The calibration is appended to a TSE in the same sequence file by
% appendFlashCalibration. It acquires in one RF-spoiled slab steady state:
%   part 0  no-Wave y projection reference, Ncalib_y ky x 1 kz   (y Wave)
%   part 1  y-sine Wave projection on the lines of part 0         (y Wave)
%   part 2  no-Wave z projection reference, 1 ky x Ncalib_z kz   (z Wave or
%                                                                 z reference)
%   part 3  z-cosine Wave projection on the lines of part 2       (z Wave)
%   part 4  no-Wave 3D ACS, Nacs(1) ky x Nacs(2) kz, for the coil maps
% Every ADC is a REF readout with a compact LIN 0...N-1 in acquisition
% order; PAR, SLC, SET and ECO are zero (docs/calibration_storage.md). The
% readout, ADC, y-sine and z-cosine Wave and the cosine bridges reproduce
% the TSE, and the PSF projections reproduce gre_mixed_dim_mb2_calibration.m.
%
% Input:
%   cfg  Struct with the required fields
%          sys_sequence, sys_wave, sys_lowPNS  Pulseq system structs of the TSE
%          physical_slew_max   physical slew limit [T/m/s]
%          ax                  axis names d1 (readout), d2 (PE), d3 (slice)
%          fov                 in-plane FOV [m], 1x2
%          N                   in-plane matrix [Nx Ny]
%          ro_os               readout oversampling factor
%          adc                 ADC event of the TSE readout
%          readoutTime         TSE readout-block duration [s]
%          Ncycles, gwave_max, swave_max, wavePolarity  Wave design of the TSE
%          includeYWave        acquire parts 0 and 1
%          includeZWave        acquire parts 2 and 3
%          includeZReference   acquire part 2 also without z Wave
%          slicePositions      physical slice centers [m], ascending
%          sliceThickness      [m]
%        and the optional fields (defaults in parentheses)
%          tseGyWave, tseGzWave  TSE Wave events; when given, the calibration
%                                Wave waveforms must equal them
%          cosineBridgeDurationBudget (370e-6 s), flipAngle (7 deg),
%          rfDuration (1e-3 s), rfTimeBandwidth (20), rfApodization (0.42),
%          rfSpoilingInc (50 deg), NdummyPsf (100), NsettlePerPart (10),
%          Ncalib_y (72), Ncalib_z (80), psfFovZ (160e-3 m),
%          psfSlabThickness (psfFovZ), psfSlabCenter (0 m), Nacs ([32 32]),
%          acsFovZ (180e-3 m), acsGridSpacingZ (slice center spacing),
%          roSpoil (3), readoutRampTime (250e-6 s), verbose (true).
%
% Output:
%   cal  Struct with the events, timing and label objects used by
%        appendFlashCalibration and addFlashCalibrationShot, and cal.meta,
%        the sidecar description of the calibration: acquisition table,
%        PSF and ACS record maps, geometry, timing and the validation report.
%
% Errors:
%   Raises an assertion error when the configuration is inconsistent, when
%   the calibration Wave differs from the TSE Wave, or when the standalone
%   timing, label, trajectory or TE/TR checks fail.

cfg = applyDefaults(cfg);
required = {'sys_sequence', 'sys_wave', 'sys_lowPNS', 'physical_slew_max', ...
    'ax', 'fov', 'N', 'ro_os', 'adc', 'readoutTime', 'Ncycles', ...
    'gwave_max', 'swave_max', 'wavePolarity', 'includeYWave', ...
    'includeZWave', 'includeZReference', 'slicePositions', 'sliceThickness'};
for iField = 1:numel(required)
    assert(isfield(cfg, required{iField}), ...
        'designFlashCalibration: missing cfg.%s.', required{iField});
end

sys = cfg.sys_sequence;
sysLow = cfg.sys_lowPNS;
sysWave = cfg.sys_wave;
ax = cfg.ax;
Nx = cfg.N(1);
Ny = cfg.N(2);
Nx_os = Nx*cfg.ro_os;
deltakX = 1/cfg.fov(1);
deltakY = 1/cfg.fov(2);
Nacs_y = cfg.Nacs(1);
Nacs_z = cfg.Nacs(2);
includeYWave = logical(cfg.includeYWave);
includeZWave = logical(cfg.includeZWave);
includeZReference = logical(cfg.includeZReference) || includeZWave;

%% Geometry
% The slab covers the slice stack; the ACS kz FOV covers the slab so that
% the excited magnetization does not alias, and its grid spacing places the
% reconstructed z planes on the slice centers.
slicePositions = double(cfg.slicePositions(:).');
assert(numel(slicePositions) >= 1 && all(diff(slicePositions) > 0), ...
    'Slice positions must be strictly ascending.');
stackExtent = slicePositions(end)-slicePositions(1)+cfg.sliceThickness;
assert(cfg.psfSlabThickness >= stackExtent-1e-12 && ...
    cfg.psfFovZ >= cfg.psfSlabThickness-1e-12, ...
    'The slab must cover the slice stack and lie inside the PSF z FOV.');
assert(cfg.acsFovZ >= cfg.psfSlabThickness-1e-12, ...
    'The ACS kz FOV must cover the excited slab.');
acsGridNz = round(cfg.acsFovZ/cfg.acsGridSpacingZ);
assert(abs(acsGridNz*cfg.acsGridSpacingZ-cfg.acsFovZ) < 1e-9, ...
    'acsFovZ must be an integer multiple of acsGridSpacingZ.');
sliceSteps = (slicePositions-slicePositions(1))/cfg.acsGridSpacingZ;
assert(max(abs(sliceSteps-round(sliceSteps))) < 1e-6, ...
    'Slice centers must lie on the ACS z grid spacing.');
assert(all(round([Nacs_y Nacs_z]) == [Nacs_y Nacs_z]) && ...
    Nacs_y <= Ny && Nacs_z <= acsGridNz && min(Nacs_y, Nacs_z) >= 24, ...
    'The ACS needs 24 to Ny ky and 24 to acsGridNz kz lines.');
assert(cfg.Ncalib_y > 1 && cfg.Ncalib_y <= Ny && cfg.Ncalib_z > 1, ...
    'Invalid PSF projection widths.');
assert(cfg.NdummyPsf >= 0 && cfg.NsettlePerPart >= 0 && ...
    cfg.NdummyPsf == round(cfg.NdummyPsf) && ...
    cfg.NsettlePerPart == round(cfg.NsettlePerPart), ...
    'Dummy and settling counts must be nonnegative integers.');

%% RF and readout
[rf, gzRf, gzReph] = mr.makeSincPulse(cfg.flipAngle*pi/180, ...
    'Duration', cfg.rfDuration, 'SliceThickness', cfg.psfSlabThickness, ...
    'apodization', cfg.rfApodization, 'timeBwProduct', cfg.rfTimeBandwidth, ...
    'system', sysLow, 'use', 'excitation');
rfBlockDuration = mr.calcDuration(rf, gzRf);

% The calibration uses the ADC event of the TSE itself.
adc = cfg.adc;
adc.phaseOffset = 0;
adc.freqOffset = 0;
assert(adc.numSamples == Nx_os, 'The TSE ADC does not hold Nx*ro_os samples.');
Tread = adc.numSamples*adc.dwell;
readoutTime = cfg.readoutTime;
assert(abs(adc.delay-sys.adcDeadTime) < 1e-9 && ...
    abs(readoutTime-(Tread+2*sys.adcDeadTime)) < 1e-9, ...
    'The TSE readout block must be the ADC plus one dead time on each side.');

readoutRampTime = cfg.readoutRampTime;
readoutAmplitude = Nx*deltakX/readoutTime;
gxRampUp = mr.makeExtendedTrapezoid(ax.d1, 'times', [0 readoutRampTime], ...
    'amplitudes', [0 readoutAmplitude], 'system', sysLow);
gxRead = mr.makeExtendedTrapezoid(ax.d1, 'times', [0 readoutTime], ...
    'amplitudes', [readoutAmplitude readoutAmplitude], 'system', sysLow);
gxRampDown = mr.makeExtendedTrapezoid(ax.d1, 'times', [0 readoutRampTime], ...
    'amplitudes', [readoutAmplitude 0], 'system', sysLow);

%% Wave events
% Designed exactly as in the TSE; the comparisons below make any difference
% an error instead of a silently mismatched PSF.
gyWave = [];
gzWave = [];
gzCosPreBridge = [];
gzCosPostBridge = [];
gxBridgeHold = [];
waveTiming = struct('actualAmplitude_mTm', 0, 'slewPeak_T_m_s', 0);
cosineWaveTiming = struct('actualAmplitude_mTm', 0, 'slewPeak_T_m_s', 0);
cosineBridgeDuration = 0;
if includeYWave
    [gyWave, ~, waveTiming] = defineSineReadoutWave(ax.d2, Tread, ...
        adc.delay, readoutTime, sysWave, cfg.Ncycles, cfg.wavePolarity, ...
        cfg.gwave_max, cfg.swave_max, cfg.physical_slew_max, adc, ...
        cfg.verbose, cfg.verbose);
    assertSameWave(gyWave, cfg.tseGyWave, 'y-sine');
end
if includeZWave
    [gzWave, waveAreaZ, cosineWaveTiming] = defineCosineReadoutWave(ax.d3, ...
        Tread, adc.delay, readoutTime, sysWave, cfg.Ncycles, ...
        cfg.wavePolarity, cfg.gwave_max, cfg.swave_max, ...
        cfg.physical_slew_max, adc, cfg.verbose, cfg.verbose);
    assertSameWave(gzWave, cfg.tseGzWave, 'z-cosine');
    gzCosPreBridge = makeFixedDurationThreePointBridge(ax.d3, 0, ...
        cosineWaveTiming.edgeAmplitude, -cosineWaveTiming.areaToAdcCenter, ...
        cfg.cosineBridgeDurationBudget, sysLow);
    gzCosPostBridge = makeFixedDurationThreePointBridge(ax.d3, ...
        cosineWaveTiming.edgeAmplitude, 0, -(gzCosPreBridge.area+waveAreaZ), ...
        cfg.cosineBridgeDurationBudget, sysLow);
    cosineBridgeDuration = max(mr.calcDuration(gzCosPreBridge), ...
        mr.calcDuration(gzCosPostBridge));
    assert(abs(gzCosPreBridge.area+cosineWaveTiming.areaToAdcCenter) <= 1e-6 && ...
        abs(gzCosPreBridge.area+waveAreaZ+gzCosPostBridge.area) <= 1e-6, ...
        'Cosine Wave moments do not close.');
    gxBridgeHold = mr.makeExtendedTrapezoid(ax.d1, ...
        'times', [0 cosineBridgeDuration], ...
        'amplitudes', [readoutAmplitude readoutAmplitude], 'system', sysLow);
end

%% Scientific parts and compact acquisition table
MODE_NOWAVE = 1;
MODE_WAVE_Y = 2;
MODE_WAVE_Z = 3;
phaseAreasY = ((0:Ny-1)-floor(Ny/2))*deltakY;
phaseAreasZPsf = ((0:cfg.Ncalib_z-1)-floor(cfg.Ncalib_z/2))/cfg.psfFovZ;
phaseAreasZAcs = ((0:acsGridNz-1)-floor(acsGridNz/2))/cfg.acsFovZ;
kyCalib = centerBlockIndices(Ny, cfg.Ncalib_y);
kyCenter = centerBlockIndices(Ny, 1);
kyAcs = centerBlockIndices(Ny, Nacs_y);
kzCalib = 1:cfg.Ncalib_z;
kzCenterPsf = centerBlockIndices(cfg.Ncalib_z, 1);
kzAcs = centerBlockIndices(acsGridNz, Nacs_z);

parts = struct('id', {}, 'name', {}, 'mode', {}, 'kyList', {}, ...
    'kzList', {}, 'kyAreas', {}, 'kzAreas', {}, 'nDummy', {});
if includeYWave
    parts(end+1) = makePart(0, 'nowave_y_projection_3d', MODE_NOWAVE, ...
        kyCalib, kzCenterPsf, phaseAreasY(kyCalib), phaseAreasZPsf(kzCenterPsf));
    parts(end+1) = makePart(1, 'wave_y_projection_3d', MODE_WAVE_Y, ...
        kyCalib, kzCenterPsf, phaseAreasY(kyCalib), phaseAreasZPsf(kzCenterPsf));
end
if includeZReference
    parts(end+1) = makePart(2, 'nowave_z_projection_3d', MODE_NOWAVE, ...
        kyCenter, kzCalib, phaseAreasY(kyCenter), phaseAreasZPsf(kzCalib));
end
if includeZWave
    parts(end+1) = makePart(3, 'wave_z_projection_3d', MODE_WAVE_Z, ...
        kyCenter, kzCalib, phaseAreasY(kyCenter), phaseAreasZPsf(kzCalib));
end
parts(end+1) = makePart(4, 'nowave_acs_3d', MODE_NOWAVE, kyAcs, kzAcs, ...
    phaseAreasY(kyAcs), phaseAreasZAcs(kzAcs));
% One steady-state preparation before the first part, short settling
% trains between parts; all parts share the RF, TE and TR.
for p = 1:numel(parts)
    if p == 1
        parts(p).nDummy = cfg.NdummyPsf;
    else
        parts(p).nDummy = cfg.NsettlePerPart;
    end
end

acqTable = struct('partArrayIdx', {}, 'setID', {}, 'mode', {}, ...
    'isRef', {}, 'kyPhysicalIndex', {}, 'kzPhysicalIndex', {}, ...
    'kyArea', {}, 'kzArea', {}, 'linLocal', {}, 'parLocal', {}, 'rawLin0', {});
partStart = zeros(1, numel(parts));
partStop = zeros(1, numel(parts));
for p = 1:numel(parts)
    partStart(p) = numel(acqTable)+1;
    for iKy = 1:numel(parts(p).kyList)
        for iKz = 1:numel(parts(p).kzList)
            row.partArrayIdx = p;
            row.setID = parts(p).id;
            row.mode = parts(p).mode;
            row.isRef = 1;
            row.kyPhysicalIndex = parts(p).kyList(iKy);
            row.kzPhysicalIndex = parts(p).kzList(iKz);
            row.kyArea = parts(p).kyAreas(iKy);
            row.kzArea = parts(p).kzAreas(iKz);
            row.linLocal = iKy-1;
            row.parLocal = iKz-1;
            row.rawLin0 = numel(acqTable); % zero-based storage index, not physical ky
            acqTable(end+1) = row; %#ok<AGROW>
        end
    end
    partStop(p) = numel(acqTable);
end
expectedAdcCount = numel(acqTable);
psfAdcIndices = zeros(4, max(cfg.Ncalib_y, cfg.Ncalib_z));
acsAdcIndices = zeros(Nacs_y, Nacs_z);
for iAdc = 1:expectedAdcCount
    row = acqTable(iAdc);
    if row.setID < 4
        column = max(row.linLocal, row.parLocal)+1;
        assert(psfAdcIndices(row.setID+1, column) == 0, ...
            'Duplicate PSF mapping at ADC %d.', iAdc);
        psfAdcIndices(row.setID+1, column) = iAdc;
    else
        assert(acsAdcIndices(row.linLocal+1, row.parLocal+1) == 0, ...
            'Duplicate ACS mapping at ADC %d.', iAdc);
        acsAdcIndices(row.linLocal+1, row.parLocal+1) = iAdc;
    end
end
psfValidMask = psfAdcIndices > 0;
psfValidLengths = sum(psfValidMask, 2).';
assert(all(acsAdcIndices(:) > 0) && isequal(sort([psfAdcIndices(psfValidMask); ...
    acsAdcIndices(:)]), (1:expectedAdcCount).'), ...
    'PSF/ACS tables must map every ADC exactly once.');

%% Prephasers, spoilers and the common TE/TR
% Shared by all parts; the durations follow the largest PE moments of the
% acquisition table, so every shot has the same timing.
allKyAreas = [acqTable.kyArea];
allKzAreas = [acqTable.kzArea];
adcCenterWithinReadout = adc.delay+0.5*adc.numSamples*adc.dwell;
groPreArea = -(0.5*readoutAmplitude*readoutRampTime+ ...
    readoutAmplitude*cosineBridgeDuration+readoutAmplitude*adcCenterWithinReadout);
preDurations = [mr.calcDuration(mr.makeTrapezoid(ax.d1, 'Area', groPreArea, 'system', sysLow)), ...
    mr.calcDuration(mr.makeTrapezoid(ax.d2, 'Area', max(abs(allKyAreas)), 'system', sysLow)), ...
    mr.calcDuration(mr.makeTrapezoid(ax.d3, 'Area', gzReph.area+min(allKzAreas), 'system', sysLow)), ...
    mr.calcDuration(mr.makeTrapezoid(ax.d3, 'Area', gzReph.area+max(allKzAreas), 'system', sysLow))];
preGradDuration = ceil(max(preDurations)/sys.gradRasterTime)*sys.gradRasterTime;
groPre = mr.makeTrapezoid(ax.d1, 'Area', groPreArea, ...
    'Duration', preGradDuration, 'system', sysLow);

readoutArea = readoutAmplitude*(readoutTime+readoutRampTime);
groSpoilArea = cfg.roSpoil*readoutArea/2;
spoilAreaZ = 4/cfg.psfSlabThickness;
postDurations = [mr.calcDuration(mr.makeTrapezoid(ax.d1, 'Area', groSpoilArea, 'system', sysLow)), ...
    mr.calcDuration(mr.makeTrapezoid(ax.d2, 'Area', max(abs(allKyAreas)), 'system', sysLow)), ...
    mr.calcDuration(mr.makeTrapezoid(ax.d3, 'Area', spoilAreaZ-min(allKzAreas), 'system', sysLow)), ...
    mr.calcDuration(mr.makeTrapezoid(ax.d3, 'Area', spoilAreaZ-max(allKzAreas), 'system', sysLow))];
minimumPostDuration = ceil(max(postDurations)/sys.gradRasterTime)*sys.gradRasterTime;

rfCenterToEnd = rfBlockDuration-(rf.delay+mr.calcRfCenter(rf));
readoutLeadAfterPrephase = readoutRampTime+cosineBridgeDuration+adcCenterWithinReadout;
TE = ceil((rfCenterToEnd+preGradDuration+readoutLeadAfterPrephase)/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
delayTE = round((TE-rfCenterToEnd-readoutLeadAfterPrephase)/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
assert(delayTE >= preGradDuration-1e-12, 'The derived TE is too short for the prephasers.');
baseDuration = rfBlockDuration+delayTE+2*readoutRampTime+ ...
    2*cosineBridgeDuration+readoutTime;
shotDuration = ceil((baseDuration+minimumPostDuration)/sys.gradRasterTime)* ...
    sys.gradRasterTime;
spoilDuration = shotDuration-baseDuration;
assert(spoilDuration >= minimumPostDuration-1e-12, ...
    'The derived TR is too short for the rewinders and spoilers.');

%% Assemble the calibration struct
cal = struct();
cal.sys_sequence = sys;
cal.sys_lowPNS = sysLow;
cal.ax = ax;
cal.rf = rf;
cal.gzRf = gzRf;
cal.rfBlockEnvelope = mr.makeDelay(rfBlockDuration);
cal.gzRephArea = gzReph.area;
cal.slabCenter = cfg.psfSlabCenter;
cal.delayTE = delayTE;
cal.preGradDuration = preGradDuration;
cal.groPre = groPre;
cal.gxRampUp = gxRampUp;
cal.gxRead = gxRead;
cal.gxRampDown = gxRampDown;
cal.gxBridgeHold = gxBridgeHold;
cal.gyWave = gyWave;
cal.gzWave = gzWave;
cal.gzCosPreBridge = gzCosPreBridge;
cal.gzCosPostBridge = gzCosPostBridge;
cal.spoilDuration = spoilDuration;
cal.groSpoilArea = groSpoilArea;
cal.spoilAreaZ = spoilAreaZ;
cal.adc = adc;
cal.rfSpoilingInc = cfg.rfSpoilingInc;
cal.MODE_WAVE_Y = MODE_WAVE_Y;
cal.MODE_WAVE_Z = MODE_WAVE_Z;
cal.parts = parts;
cal.acqTable = acqTable;
cal.partStart = partStart;
cal.partStop = partStop;
cal.lblLIN = arrayfun(@(v) mr.makeLabel('SET', 'LIN', v), ...
    0:expectedAdcCount-1, 'UniformOutput', false);
cal.lblZero = {mr.makeLabel('SET', 'PAR', 0), mr.makeLabel('SET', 'SLC', 0), ...
    mr.makeLabel('SET', 'SET', 0), mr.makeLabel('SET', 'ECO', 0), ...
    mr.makeLabel('SET', 'REF', true)};

meta = struct();
meta.storageLayoutVersion = 'appended_ref_lin_v1';
meta.expectedAdcCount = expectedAdcCount;
meta.acqTable = acqTable;
meta.parts = rmfield(parts, {'kyAreas', 'kzAreas'});
meta.partStart = partStart;
meta.partStop = partStop;
meta.psfAdcIndices = psfAdcIndices;
meta.psfValidMask = double(psfValidMask);
meta.psfValidLengths = psfValidLengths;
meta.acsAdcIndices = acsAdcIndices;
meta.acsKyIndices = kyAcs;
meta.acsKzIndices = kzAcs;
meta.acsKzAreas = phaseAreasZAcs(kzAcs);
meta.acsFovZ = cfg.acsFovZ;
meta.acsGridSpacingZ = cfg.acsGridSpacingZ;
meta.acsGridNz = acsGridNz;
meta.Nacs = [Nacs_y Nacs_z];
meta.psfYCoordinates = ((0:cfg.Ncalib_y-1)-floor(cfg.Ncalib_y/2))*cfg.fov(2)/cfg.Ncalib_y;
meta.psfZCoordinates = cfg.psfSlabCenter+ ...
    ((0:cfg.Ncalib_z-1)-floor(cfg.Ncalib_z/2))*cfg.psfFovZ/cfg.Ncalib_z;
meta.psfFovZ = cfg.psfFovZ;
meta.psfSlabThickness = cfg.psfSlabThickness;
meta.psfSlabCenter = cfg.psfSlabCenter;
meta.Ncalib_y = cfg.Ncalib_y;
meta.Ncalib_z = cfg.Ncalib_z;
meta.includeYWave = double(includeYWave);
meta.includeZWave = double(includeZWave);
meta.includeZReference = double(includeZReference);
meta.Nx = Nx;
meta.Ny = Ny;
meta.Nx_os = Nx_os;
meta.ro_os = cfg.ro_os;
meta.fov = cfg.fov(1:2);
meta.dwell = adc.dwell;
meta.Tread = Tread;
meta.slicePositions = slicePositions;
meta.sliceThickness = cfg.sliceThickness;
meta.Ncycles = cfg.Ncycles;
meta.gwave_max = cfg.gwave_max;
meta.swave_max = cfg.swave_max;
meta.wavePolarity = cfg.wavePolarity;
meta.waveAmplitude_mTm = max(waveTiming.actualAmplitude_mTm, ...
    cosineWaveTiming.actualAmplitude_mTm);
meta.cosineBridgeDuration = cosineBridgeDuration;
meta.flipAngle = cfg.flipAngle;
meta.rfDuration = cfg.rfDuration;
meta.rfTimeBandwidth = cfg.rfTimeBandwidth;
meta.rfApodization = cfg.rfApodization;
meta.rfSpoilingInc = cfg.rfSpoilingInc;
meta.NdummyPsf = cfg.NdummyPsf;
meta.NsettlePerPart = cfg.NsettlePerPart;
meta.TE = TE;
meta.TR = shotDuration;
meta.shotCount = expectedAdcCount+sum([parts.nDummy]);
cal.meta = meta;

%% Standalone validation
% The identical construction is checked on its own sequence object, which
% keeps the TSE k-space and label checks independent of the calibration.
cal.meta.validation = validateStandalone(cal, Nx_os, deltakY);
if cfg.verbose
    fprintf(['3D FLASH calibration: %d ADCs (%s), TE/TR %.3f/%.3f ms, ', ...
        '%.2f s, ACS %d x %d over a %.0f mm kz FOV.\n'], expectedAdcCount, ...
        mat2str(psfValidLengths), TE*1e3, shotDuration*1e3, ...
        cal.meta.validation.duration, Nacs_y, Nacs_z, cfg.acsFovZ*1e3);
end
end

%% Local helpers

function cfg = applyDefaults(cfg)
%APPLYDEFAULTS Fill the optional configuration fields with their defaults.
% cfg is the user struct; returns it with every missing optional field set.
    defaults = struct('tseGyWave', [], 'tseGzWave', [], ...
        'cosineBridgeDurationBudget', 370e-6, 'flipAngle', 7, ...
        'rfDuration', 1e-3, 'rfTimeBandwidth', 20, 'rfApodization', 0.42, ...
        'rfSpoilingInc', 50, 'NdummyPsf', 100, 'NsettlePerPart', 10, ...
        'Ncalib_y', 72, 'Ncalib_z', 80, 'psfFovZ', 160e-3, ...
        'psfSlabCenter', 0, 'Nacs', [32 32], 'acsFovZ', 180e-3, ...
        'roSpoil', 3, 'readoutRampTime', 250e-6, 'verbose', true);
    names = fieldnames(defaults);
    for iName = 1:numel(names)
        if ~isfield(cfg, names{iName})
            cfg.(names{iName}) = defaults.(names{iName});
        end
    end
    if ~isfield(cfg, 'psfSlabThickness')
        cfg.psfSlabThickness = cfg.psfFovZ;
    end
    if ~isfield(cfg, 'acsGridSpacingZ')
        spacing = diff(double(cfg.slicePositions(:)));
        assert(~isempty(spacing) && max(abs(spacing-spacing(1))) < 1e-12, ...
            'Set acsGridSpacingZ for nonuniform or single-slice stacks.');
        cfg.acsGridSpacingZ = spacing(1);
    end
end

function part = makePart(id, name, mode, kyList, kzList, kyAreas, kzAreas)
%MAKEPART Describe one scientific calibration part.
% Inputs are the part ID and name, the Wave mode, the one-based ky/kz grid
% indices and their PE moments [1/m]. The output is a scalar struct; the
% dummy count is assigned by the caller.
    part.id = id;
    part.name = name;
    part.mode = mode;
    part.kyList = kyList;
    part.kzList = kzList;
    part.kyAreas = kyAreas;
    part.kzAreas = kzAreas;
    part.nDummy = 0;
end

function idx = centerBlockIndices(Ndim, nCenter)
%CENTERBLOCKINDICES Return centered one-based indices within a dimension.
% Inputs are a full dimension and the requested width. The output is a row
% vector; an oversized or nonpositive request raises an assertion.
    assert(nCenter > 0 && nCenter <= Ndim, ...
        'Requested center block does not fit its dimension.');
    centerIdx = floor(Ndim/2)+1;
    if mod(nCenter, 2) == 0
        idx = (centerIdx-nCenter/2):(centerIdx+nCenter/2-1);
    else
        idx = (centerIdx-floor(nCenter/2)):(centerIdx+floor(nCenter/2));
    end
end

function assertSameWave(calibrationWave, tseWave, name)
%ASSERTSAMEWAVE Require the calibration Wave to equal the TSE Wave.
% calibrationWave and tseWave are Pulseq arbitrary gradients (tseWave may be
% empty, which skips the check); name labels the error. No output.
    if isempty(tseWave)
        return;
    end
    a = calibrationWave.waveform(:);
    b = tseWave.waveform(:);
    assert(strcmp(calibrationWave.channel, tseWave.channel) && ...
        numel(a) == numel(b) && max(abs(a-b)) <= 1e-9*max(abs(b)) && ...
        abs(calibrationWave.first-tseWave.first) <= 1e-9*max(abs(b)) && ...
        abs(calibrationWave.last-tseWave.last) <= 1e-9*max(abs(b)), ...
        'The calibration %s Wave differs from the TSE Wave.', name);
end

function report = validateStandalone(cal, Nx_os, deltakY)
%VALIDATESTANDALONE Check the calibration on its own sequence object.
% cal is the assembled calibration, Nx_os the ADC samples per readout and
% deltakY the PE k-space step [1/m]. Returns a struct with the duration,
% block count, measured TE and the largest trajectory deviations. Raises an
% assertion error when timing, labels, Wave pairs, ACS encoding, TE or the
% shot spacing are wrong.
    seqCal = mr.Sequence(cal.sys_sequence);
    appendFlashCalibration(seqCal, cal);
    [ok, errorReport] = seqCal.checkTiming;
    assert(ok, 'Calibration timing check failed:\n%s', strjoin(errorReport, '\n'));

    count = cal.meta.expectedAdcCount;
    lbl = seqCal.evalLabels('evolution', 'adc');
    assert(isfield(lbl, 'LIN') && isfield(lbl, 'REF') && ...
        numel(lbl.LIN) == count && isequal(lbl.LIN(:), (0:count-1).') && ...
        all(lbl.REF(:) ~= 0), ...
        'Calibration ADC labels differ from the compact REF/LIN contract.');
    names = {'PAR', 'SLC', 'SET', 'ECO'};
    for iName = 1:numel(names)
        if isfield(lbl, names{iName})
            assert(all(lbl.(names{iName})(:) == 0), ...
                'Calibration label %s is not zero.', names{iName});
        end
    end

    warning('off', 'mr:restoreShape');
    [ktrajAdc, tAdc, ~, ~, tExcitation] = seqCal.calculateKspacePP();
    warning('on', 'mr:restoreShape');
    k = reshape(ktrajAdc, 3, Nx_os, count);
    assert(all(squeeze(k(1,end,:)) > squeeze(k(1,1,:))), ...
        'Readout polarity is not consistently positive x.');
    rows = cal.acqTable;
    setID = [rows.setID];
    report = struct();
    % Wave/reference pairs may differ only along their Wave axis, by the
    % same amount on every line.
    % Rows: [reference SET, Wave SET, k-space axis of the Wave (2 = y, 3 = z)].
    pairDefinitions = [0 1 2; 2 3 3];
    for iPair = 1:size(pairDefinitions, 1)
        reference = find(setID == pairDefinitions(iPair, 1));
        wave = find(setID == pairDefinitions(iPair, 2));
        if isempty(wave)
            continue;
        end
        waveAxis = pairDefinitions(iPair, 3);
        assert(numel(wave) == numel(reference), ...
            'Wave SET %d and its reference differ in length.', pairDefinitions(iPair, 2));
        delta = k(:,:,wave)-k(:,:,reference);
        others = setdiff(1:3, waveAxis);
        assert(max(abs(delta(others,:,:)), [], 'all') < 1e-5, ...
            'Wave SET %d differs from its reference outside the Wave axis.', ...
            pairDefinitions(iPair, 2));
        assert(max(abs(delta(waveAxis,:,:)-delta(waveAxis,:,1)), [], 'all') < deltakY*1e-4, ...
            'The Wave difference of SET %d varies across lines.', pairDefinitions(iPair, 2));
    end
    % Every shot encodes its designed ky/kz; no-Wave lines are constant
    % during the readout and all lines share the readout kx.
    noWave = [rows.mode] == 1;
    kyError = squeeze(k(2,:,noWave))-[rows(noWave).kyArea];
    kzError = squeeze(k(3,:,noWave))-[rows(noWave).kzArea];
    kxError = squeeze(k(1,:,:))-k(1,:,1).';
    report.maxNoWaveKyError = max(abs(kyError), [], 'all');
    report.maxNoWaveKzError = max(abs(kzError), [], 'all');
    report.maxKxError = max(abs(kxError), [], 'all');
    assert(report.maxNoWaveKyError < 1e-3 && report.maxNoWaveKzError < 1e-3 && ...
        report.maxKxError < 1e-3, ...
        'Calibration ky/kz encoding or readout kx differs from the design.');

    adcCenters = mean(reshape(tAdc, Nx_os, count), 1);
    excitations = tExcitation(isfinite(tExcitation));
    measuredTE = zeros(1, count);
    for iAdc = 1:count
        measuredTE(iAdc) = adcCenters(iAdc)-excitations(find(excitations < adcCenters(iAdc), 1, 'last'));
    end
    report.measuredTE = median(measuredTE);
    assert(max(abs(measuredTE-cal.meta.TE)) <= max(1e-9, cal.adc.dwell/100), ...
        'Measured calibration TE differs from the design.');
    for p = 1:numel(cal.parts)
        centers = adcCenters(cal.partStart(p):cal.partStop(p));
        if numel(centers) > 1
            assert(max(abs(diff(centers)-cal.meta.TR)) <= max(1e-9, cal.adc.dwell/100), ...
                'ADC spacing is not one TR inside calibration part %d.', cal.parts(p).id);
        end
    end
    [report.duration, report.blocks] = seqCal.duration();
end

function [gWave, area, timing] = defineSineReadoutWave(channel, Tread, adcDelay, ...
        readoutBlockDur, sys, Ncycles, wavePolarity, gwave_max, swave_max, ...
        physical_slew_max, adc, waveInfoFlag, debugFlag)
%DEFINESINEREADOUTWAVE Construct an ADC-aligned full-block sine gradient.
% Copied from gre_mixed_dim_mb2_calibration.m. Inputs describe the axis,
% active ADC timing, complete readout-block duration, polarity, scanner
% limits and print flags. Outputs are the zero-delay arbitrary gradient, its
% zeroth moment [1/m] and realized timing metadata. Assertions reject raster,
% ADC, block-duration, slew or amplitude inconsistencies.
    dt = sys.gradRasterTime;
    assert(ismember(wavePolarity, [-1 1]), 'wavePolarity must be +1 or -1.');
    [G0, w, TreadRaster, actualAmplitude_mTm] = designWaveAmplitude( ...
        Tread, sys, Ncycles, gwave_max, swave_max, physical_slew_max, waveInfoFlag);
    nReadSamples = round(TreadRaster/dt);
    assert(abs(nReadSamples*dt-TreadRaster) < dt/10, ...
        'TreadRaster must contain an integer number of gradient samples.');
    nBlockSamples = round(readoutBlockDur/dt);
    assert(abs(nBlockSamples*dt-readoutBlockDur) < dt/10, ...
        'The readout block duration must be on the gradient raster.');
    nPadPre = round(adcDelay/dt);
    assert(abs(nPadPre*dt-adcDelay) < dt/10, ...
        'The ADC delay must be on the gradient raster for full-block Wave padding.');
    nPadPost = nBlockSamples-nPadPre-nReadSamples;
    assert(nPadPost >= 0 && nReadSamples >= 2, ...
        'The active sine does not fit in the readout block.');

    activePhase = 2*pi*Ncycles*(0:(nReadSamples-1))/(nReadSamples-1);
    activeWave = wavePolarity*G0*sin(activePhase);
    activeWave([1 end]) = 0;
    fullBlockWave = [zeros(1, nPadPre), activeWave, zeros(1, nPadPost)];
    gWave = mr.makeArbitraryGrad(channel, fullBlockWave, 'system', sys, ...
        'first', 0, 'last', 0);
    gWave.delay = 0;
    area = gWave.area;
    assertArbitraryWaveformWithinSystem(fullBlockWave, 0, 0, sys, 'active sine wave');
    slewPeak = max(abs(diff(fullBlockWave)))/dt/sys.gamma;
    assert(slewPeak <= min(swave_max, physical_slew_max)*(1+1e-9), ...
        'Active sine slew exceeds the requested Wave slew limit.');
    assert(abs(mr.calcDuration(gWave)-nBlockSamples*dt) < max(dt/10, sys.adcRasterTime), ...
        'Wave event duration does not match the complete readout block.');
    timing.G0 = G0;
    timing.omega = w;
    timing.actualAmplitude_mTm = actualAmplitude_mTm;
    timing.slewPeak_T_m_s = slewPeak;
    timing.area = area;
    if debugFlag
        fprintf('Calibration y-sine Wave: amplitude=%.6f mT/m, slew=%.6f T/m/s, cycles=%d\n', ...
            actualAmplitude_mTm, slewPeak, Ncycles);
    end
end

function [gWave, area, timing] = defineCosineReadoutWave(channel, Tread, adcDelay, ...
        readoutBlockDur, sys, Ncycles, wavePolarity, gwave_max, swave_max, ...
        physical_slew_max, adc, waveInfoFlag, debugFlag)
%DEFINECOSINEREADOUTWAVE Construct an ADC-aligned full-block cosine gradient.
% Copied from gre_mixed_dim_mb2_calibration.m. Inputs describe the axis,
% ADC/readout timing, Wave design, scanner limits and print flags. Outputs
% are the +G0-endpoint arbitrary gradient, its area [1/m] and realized
% timing/moment metadata, including the moment up to the ADC center.
    dt = sys.gradRasterTime;
    assert(ismember(wavePolarity, [-1 1]), 'wavePolarity must be +1 or -1.');
    [G0, ~, TreadRaster, actualAmplitude_mTm] = designWaveAmplitude( ...
        Tread, sys, Ncycles, gwave_max, swave_max, physical_slew_max, waveInfoFlag);
    Gedge = wavePolarity*G0;
    nReadSamples = round(TreadRaster/dt);
    nBlockSamples = round(readoutBlockDur/dt);
    nPadPre = round(adcDelay/dt);
    nPadPost = nBlockSamples-nPadPre-nReadSamples;
    assert(nReadSamples >= 2 && nPadPost >= 0, ...
        'The active cosine does not fit in the readout block.');

    activePhase = 2*pi*Ncycles*(0:(nReadSamples-1))/(nReadSamples-1);
    activeWave = wavePolarity*G0*cos(activePhase);
    activeWave([1 end]) = Gedge;
    fullBlockWave = [Gedge*ones(1, nPadPre), activeWave, Gedge*ones(1, nPadPost)];
    gWave = mr.makeArbitraryGrad(channel, fullBlockWave, 'system', sys, ...
        'first', Gedge, 'last', Gedge);
    gWave.delay = 0;
    area = gWave.area;
    adcCenterTime = adc.delay+0.5*adc.numSamples*adc.dwell;
    nCellsToAdcCenter = round(adcCenterTime/dt);
    assert(abs(nCellsToAdcCenter*dt-adcCenterTime) < dt/10, ...
        'ADC center must lie on the gradient raster.');
    assertArbitraryWaveformWithinSystem(fullBlockWave, Gedge, Gedge, sys, 'active cosine wave');
    timing.edgeAmplitude = Gedge;
    timing.actualAmplitude_mTm = actualAmplitude_mTm;
    timing.slewPeak_T_m_s = max(abs(diff(fullBlockWave)))/dt/sys.gamma;
    assert(timing.slewPeak_T_m_s <= min(swave_max, physical_slew_max)*(1+1e-9), ...
        'Active cosine slew exceeds the requested Wave slew limit.');
    timing.area = area;
    timing.areaToAdcCenter = sum(fullBlockWave(1:nCellsToAdcCenter))*dt;
    if debugFlag
        fprintf('Calibration z-cosine Wave: amplitude=%.6f mT/m, area to ADC center=%.9g 1/m\n', ...
            actualAmplitude_mTm, timing.areaToAdcCenter);
    end
end

function gBridge = makeFixedDurationThreePointBridge(channel, Gstart, Gend, ...
        areaTarget, duration, sys)
%MAKEFIXEDDURATIONTHREEPOINTBRIDGE Design a low-PNS endpoint bridge.
% Copied from gre_mixed_dim_mb2_calibration.m. Inputs are the channel, start
% and end amplitudes [Hz/m], the target area [1/m], the fixed duration [s]
% and the system limits; the midpoint time minimizes the peak slew. Returns
% the extended trapezoid; an infeasible design raises an error.
    dt = sys.gradRasterTime;
    nIntervals = round(duration/dt);
    Ttotal = nIntervals*dt;
    assert(nIntervals >= 2 && abs(Ttotal-duration) <= dt/10, ...
        'Cosine bridge duration must contain at least two raster intervals.');
    best = [];
    bestSlew = inf;
    for iMid = 1:(nIntervals-1)
        Tmid = iMid*dt;
        Gmid = (2*areaTarget-Tmid*Gstart-(Ttotal-Tmid)*Gend)/Ttotal;
        peakSlew = max(abs(Gmid-Gstart)/Tmid, abs(Gend-Gmid)/(Ttotal-Tmid));
        if max(abs([Gstart Gmid Gend])) <= sys.maxGrad*(1+1e-12) && ...
                peakSlew <= sys.maxSlew*(1+1e-12) && peakSlew < bestSlew
            best = struct('Tmid', Tmid, 'Gmid', Gmid);
            bestSlew = peakSlew;
        end
    end
    if isempty(best)
        error(['No feasible separated three-point bridge was found at ', ...
               'the accepted %.3f ms duration.'], duration*1e3);
    end
    gBridge = mr.makeExtendedTrapezoid(channel, ...
        'times', [0 best.Tmid Ttotal], 'amplitudes', [Gstart best.Gmid Gend], ...
        'system', sys);
    assert(abs(gBridge.area-areaTarget) <= 1e-6, ...
        'Three-point bridge area error %.9g 1/m exceeds tolerance.', ...
        gBridge.area-areaTarget);
end

function [G0_pulseq, w, TreadRaster, actualAmplitude_mTm] = designWaveAmplitude( ...
        Tread, sys, Ncycles, gwave_max, swave_max, physical_slew_max, waveInfoFlag)
%DESIGNWAVEAMPLITUDE Calculate the slew-limited Wave amplitude.
% Copied from gre_mixed_dim_mb2_calibration.m. Inputs specify the duration,
% raster/system limits, cycles and requested limits. Outputs are the Pulseq
% amplitude [Hz/m], angular frequency [rad/s], rasterized duration [s] and
% physical amplitude [mT/m]. Assertions reject invalid cycles or duration.
    assert(ismember(Ncycles, [10 20]), ...
        'Ncycles must be one of the supported protocol options: 10 or 20.');
    wavepoints = round(Tread/sys.gradRasterTime);
    TreadRaster = wavepoints*sys.gradRasterTime;
    assert(abs(TreadRaster-Tread) < sys.gradRasterTime/10, ...
        'Tread must be on the gradient raster.');
    w = 2*pi/(TreadRaster/Ncycles);
    swave_max = min(physical_slew_max, swave_max);
    activePhase = 2*pi*Ncycles*(0:(wavepoints-1))/(wavepoints-1);
    rasterSlewPerUnitAmplitude = max([abs(diff(sin(activePhase))), ...
        abs(diff(cos(activePhase)))])/sys.gradRasterTime;
    amplitude_T_m = min(gwave_max*1e-3, swave_max/rasterSlewPerUnitAmplitude);
    G0_pulseq = amplitude_T_m*sys.gamma;
    actualAmplitude_mTm = amplitude_T_m*1e3;
    if waveInfoFlag
        fprintf('Calibration Wave amplitude: requested/actual %.6f/%.6f mT/m.\n', ...
            gwave_max, actualAmplitude_mTm);
    end
end

function assertArbitraryWaveformWithinSystem(waveform, firstValue, lastValue, system, segmentName)
%ASSERTARBITRARYWAVEFORMWITHINSYSTEM Check an arbitrary gradient against limits.
% waveform holds the raster samples [Hz/m] and firstValue/lastValue the
% block-boundary values; system supplies maxGrad, maxSlew and the gradient
% raster, and segmentName labels the error. The slew includes the half-raster
% steps to the boundary values. No output; an assertion error is raised when
% a limit is exceeded by more than a relative 1e-10.
    waveform = waveform(:);
    slew = [2*(firstValue-waveform(1)); waveform(2:end)-waveform(1:end-1); ...
        2*(waveform(end)-lastValue)]/system.gradRasterTime;
    assert(max(abs(waveform)) <= system.maxGrad*(1+1e-10), ...
        '%s exceeds the gradient limit.', segmentName);
    assert(max(abs(slew)) <= system.maxSlew*(1+1e-10), ...
        '%s exceeds the slew limit.', segmentName);
end
