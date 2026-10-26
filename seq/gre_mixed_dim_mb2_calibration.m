%% Mixed-dimensional spoiled-GRE calibration for yz-Wave MB2 TSE
%
% Scientific part IDs (sidecar metadata only; scanner SET is always zero):
%   0: no-wave, 72 ky x 1 kz, 3D slab y-projection reference
%   1: y-sine only, same 3D projection grid as part 0
%   2: no-wave, 1 ky x 80 kz, 3D slab z-projection reference
%   3: z-cosine only, same 3D projection grid as part 2
%   4: no-wave, 32 ky at each of 50 direct 2D slices, REF=true
%
% Every acquired ADC has a unique sequential LIN; PAR/SLC/SET/ECO are zero.
% The sidecar maps ADC records to a 4-by-80 PSF table and a 50-by-32 CSM
% table. Physical slice positions and ky/kz encoding are not storage counters.
% See docs/calibration_storage.md before reconstruction.
%
% SETs 0--3 share one slab-selective RF, one steady-state preparation, and
% short per-SET settling trains. SET 4 switches to the geometry-matched 2D
% single-band RF and therefore receives its own 250-shot (5-cycle) dummy train.
% The PSF z dimension is true Fourier-encoded kz; the CSM SLC dimension is
% direct image space and must never be Fourier transformed or interpolated.
%
% The x readout, ADC, y-sine, z-cosine, and separated z bridges reproduce
% tse_2d_wave_mb2v_v3.m. The calibration and target therefore share ADC time
% samples, physical gradient axes, Wave amplitude, cycles, polarity, and
% bridge timing. The PSF fit should use physical coordinates in metres and
% store slope coefficients in rad/m before evaluating them on the TSE grid.

clearvars; close all; clc
format long

%% Paths

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
addpath(fullfile(script_dir, 'utils'));
mbPulseMatFile = fullfile(script_dir, 'rf_pulses', ...
    'mb2_tse_ls_tbw7_beta150_v1verse_B1max15p000uT.mat');

% Machine-specific paths are entered interactively and saved beside this
% script in wave_tse_path_settings.json (ignored by Git). The .seq file and
% its _geometry_labels.mat/.txt sidecars are written directly to out_path.
pathSettings = configurePathSettings(script_dir);
pulseq_path = pathSettings.pulseq_path;
safe_pns_prediction_path = pathSettings.safe_pns_prediction_path;
out_path = pathSettings.out_path;
system_asc_file = pathSettings.system_asc_file;

%% Parameters

% ------------------------- write options -------------------------
write_v141_format = false;
writeSequence = true;
isUseWaveCalibration = true;
stopAfterTimingCheck = false;
execute_system_check_code = true;

% ------------------------- system -------------------------
sys_type = 'prisma';

% ------------------------- sequence timing / RF -------------------------
alphaPsf = 7;                               % [deg], repeated slab spoiled GRE
alphaCsm = 20;                               % [deg], geometry-matched 2D CSM
rfDurationPsf = 1e-3;                        % [s]
rfDurationCsm = 2e-3;                        % [s]
rfSpoilingInc = 50;                          % [deg]
NdummyPsf = 100;                             % one 3D slab preparation
NsettlePsfPerSet = 10;                       % short settling before SETs 1--3
NdummyCsm = 250;                             % 5 full 50-slice cycles after RF switch
csmShotSpacingFloor = 15e-3;                 % [s], conservative global RF-duty floor
maximumDesiredSequenceDuration = 40;         % [s], fail before writing if exceeded

% ------------------------- geometry / sampling -------------------------
fov = [220 220]*1e-3;                        % [m], x/y
N = [440 220];
ro_os = 3;
ro_dur = 6600e-6;                            % [s]
ro_spoil = 3;

Ncalib_y = 72;
Ncalib_z = 80;
psfFovZ = 160e-3;                            % [m], 2 mm grid covering 150 mm stack
psfSlabThickness = psfFovZ;                  % [m], no excitation beyond encoded FOV
psfSlabCenter = 0;                           % [m], scanner physical z coordinate
Nacs_y = 32;

Nslices = 50;
sliceThickness = 3e-3;
sliceSpacing = 0;
sliceCenterSpacing = sliceThickness+sliceSpacing;
slOrientation = 'TRA';

ax.d1 = 'x';
ax.d2 = 'y';
ax.d3 = 'z';

% ------------------------- wave -------------------------
gwave_max = 10;                              % [mT/m], matches MB2 TSE
swave_max = 140;                             % [T/m/s], matches MB2 TSE
Ncycles = 20;                                % supported protocol options: 10 or 20
wavePolarity = 1;
waveRepresentation = 'fullblock-arbitrary';
cosineBridgeDurationBudget = 370e-6;          % Skyra low-PNS bridge per readout side

%% System limits

slew_safety_margin = 0.7;
grad_safety_margin = 0.9;
lowPNS_slew_safety_margin = 0.41;
lowPNS_grad_safety_margin = grad_safety_margin;

switch sys_type
    case {'prisma', 'C2_simulate_prisma', 'prisma_XA30A'}
        physical_slew_max = 200; physical_grad_max = 80; B0 = 2.89;
    case 'premier'
        physical_slew_max = 200; physical_grad_max = 70; B0 = 3;
    case 'Connectome2'
        physical_slew_max = 598.802; physical_grad_max = 500; B0 = 2.89;
    case 'skyra'
        physical_slew_max = 180; physical_grad_max = 43; B0 = 2.89;
    case 'trio'
        physical_slew_max = 170; physical_grad_max = 38; B0 = 2.89;
    case 'CimaX'
        physical_slew_max = 200; physical_grad_max = 200; B0 = 2.89;
    case 'TerraX'
        physical_slew_max = 250; physical_grad_max = 135; B0 = 2.89;
    otherwise
        error('Undefined sys_type: %s', sys_type);
end

rfDeadTime = 100e-6;
rfRingdownTime = 100e-6;
adcDeadTime = 20e-6;
adcRasterTime = 100e-9;
rfRasterTime = 1e-6;
gradRasterTime = 10e-6;
blockDurationRaster = 10e-6;

%% Resolve paths and validate configuration

assert(ismember(Ncycles, [10 20]), ...
    'Ncycles must be one of the supported protocol options: 10 or 20.');
addSequenceToolboxPaths(pulseq_path, safe_pns_prediction_path);

assert(psfSlabThickness <= psfFovZ, ...
    'The PSF slab must not excite beyond its encoded z FOV.');
assert(all([Ncalib_y Ncalib_z Nacs_y] == ...
    round([Ncalib_y Ncalib_z Nacs_y])), ...
    'Calibration widths must be integers.');
assert(Ncalib_y > 1 && Ncalib_y <= N(2) && Ncalib_z > 1 && Nacs_y > 0, ...
    'Invalid PSF or CSM calibration width.');
assert(all([NdummyPsf NsettlePsfPerSet NdummyCsm] >= 0) && ...
    all([NdummyPsf NsettlePsfPerSet NdummyCsm] == ...
    round([NdummyPsf NsettlePsfPerSet NdummyCsm])), ...
    'Dummy and settling counts must be nonnegative integers.');
assert(mod(NdummyCsm, Nslices) == 0, ...
    'NdummyCsm must contain complete physical-slice cycles.');
assert(csmShotSpacingFloor > 0 && maximumDesiredSequenceDuration > 0, ...
    'CSM shot-spacing and sequence-duration limits must be positive.');

assert(isfile(mbPulseMatFile), ...
    'Shared MB2 RF/geometry MAT file not found: %s', mbPulseMatFile);
MB = load(mbPulseMatFile);
requiredFields = {'slice_positions_m', 'sms_groups_1based', ...
    'group_centers_m', 'interleaved_group_order_1based', ...
    'pair_separation_m', 'slice_center_spacing_m', ...
    'slice_thickness_m', 'slice_gap_m', 'physical_slice_count', ...
    'sms_factor', 'sms_group_count', 'candidate_tag', ...
    'reference_arrays_sha256', 'candidate_status'};
for iField = 1:numel(requiredFields)
    assert(isfield(MB, requiredFields{iField}), ...
        'Shared MB2 geometry MAT file is missing %s.', requiredFields{iField});
end

slicePositions = double(MB.slice_positions_m(:).');
smsGroups = double(MB.sms_groups_1based);
groupCenters = double(MB.group_centers_m(:).');
mbGroupOrder = double(MB.interleaved_group_order_1based(:).');
pairSeparation = double(MB.pair_separation_m);
candidateTag = char(string(MB.candidate_tag));
candidateChecksum = char(string(MB.reference_arrays_sha256));
candidateStatus = char(string(MB.candidate_status));
assert(numel(slicePositions) == Nslices && all(diff(slicePositions) > 0), ...
    'Expected %d strictly ascending physical slice positions.', Nslices);
assert(max(abs(diff(slicePositions)-sliceCenterSpacing)) < 1e-12, ...
    'Shared MAT slice spacing does not match the requested CSM geometry.');
assert(double(MB.physical_slice_count) == Nslices && ...
       double(MB.sms_factor) == 2 && ...
       double(MB.sms_group_count) == Nslices/2, ...
    'Shared MAT physical/SMS slice counts do not match this calibration.');
assert(abs(double(MB.slice_thickness_m)-sliceThickness) < 1e-12 && ...
       abs(double(MB.slice_gap_m)-sliceSpacing) < 1e-12, ...
    'Shared MAT slice thickness/gap does not match this calibration.');
assert(abs(mean(slicePositions)) < 1e-12, ...
    'The physical TSE slice stack must be centered.');
assert(size(smsGroups, 1) == Nslices/2 && size(smsGroups, 2) == 2, ...
    'Unexpected MB2 SMS group table.');
expectedPairSeparation = (Nslices/2)*sliceCenterSpacing;
assert(abs(pairSeparation-expectedPairSeparation) < 1e-12, ...
    'Unexpected MB2 pair separation.');
assert(strcmp(candidateTag, ...
       'ls_tbw7_d2_1e-06_beta150_phasealigned_v1') && ...
       strcmp(candidateChecksum, ...
       'f34ebbcf0c52c8638eb6dc1da6e778673776fc2465a26e71e8d6e6da5e793b22'), ...
    'Unexpected TBW7 beta150 candidate identity or array checksum.');

sliceAcquisitionOrder = [1:2:Nslices, 2:2:Nslices];
sliceByAcquisitionIndex = sliceAcquisitionOrder;
acquisitionIndexBySlice = zeros(1, Nslices);
acquisitionIndexBySlice(sliceAcquisitionOrder) = 1:Nslices;
stackExtent = Nslices*sliceThickness+(Nslices-1)*sliceSpacing;
assert(psfSlabThickness >= stackExtent && psfFovZ >= stackExtent, ...
    'The PSF slab/FOV must cover the complete physical slice stack.');

% These coordinates match the centered FFT convention used by reconstruction.
psfYCoordinates = ((0:Ncalib_y-1)-floor(Ncalib_y/2))*fov(2)/Ncalib_y;
psfZCoordinates = psfSlabCenter + ...
    ((0:Ncalib_z-1)-floor(Ncalib_z/2))*psfFovZ/Ncalib_z;

fprintf('Mixed-dimensional calibration geometry:\n');
fprintf('  PSF FOV / matrix = %.1f x %.1f x %.1f mm / %d x %d x %d\n', ...
    fov(1)*1e3, fov(2)*1e3, psfFovZ*1e3, N(1), N(2), Ncalib_z);
fprintf('  PSF slab center / thickness = %.3f / %.3f mm\n', ...
    psfSlabCenter*1e3, psfSlabThickness*1e3);
fprintf('  CSM physical slice centers [mm] = %s\n', ...
    mat2str(slicePositions*1e3, 6));

%% Setup

% The sequence-wide system is the physical scanner envelope. Only the
% active Wave samples use this slew allowance; every non-Wave gradient is
% still constructed with sys_lowPNS below.
sys_sequence = mr.opts('MaxGrad', physical_grad_max*grad_safety_margin, ...
    'GradUnit', 'mT/m', 'MaxSlew', physical_slew_max, ...
    'SlewUnit', 'T/m/s', 'rfDeadTime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, 'adcDeadTime', adcDeadTime, ...
    'adcRasterTime', adcRasterTime, 'rfRasterTime', rfRasterTime, ...
    'gradRasterTime', gradRasterTime, ...
    'blockDurationRaster', blockDurationRaster, 'B0', B0);
sys_wave = sys_sequence;
sys = sys_sequence; % timing/raster compatibility alias; not for non-Wave gradient design
sys_lowPNS = mr.opts('MaxGrad', physical_grad_max*lowPNS_grad_safety_margin, ...
    'GradUnit', 'mT/m', 'MaxSlew', ...
    physical_slew_max*lowPNS_slew_safety_margin, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, 'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime, 'adcRasterTime', adcRasterTime, ...
    'rfRasterTime', rfRasterTime, 'gradRasterTime', gradRasterTime, ...
    'blockDurationRaster', blockDurationRaster, 'B0', B0);
seq = mr.Sequence(sys_sequence);

%% RF, readout, PE, and Wave events

[rfPsf, gzPsf, gzPsfReph] = mr.makeSincPulse(alphaPsf*pi/180, ...
    'Duration', rfDurationPsf, 'SliceThickness', psfSlabThickness, ...
    'apodization', 0.42, 'timeBwProduct', 20, ...
    'system', sys_lowPNS, 'use', 'excitation');
[rfCsm, gzCsm, gzCsmReph] = mr.makeSincPulse(alphaCsm*pi/180, ...
    'Duration', rfDurationCsm, 'SliceThickness', sliceThickness, ...
    'apodization', 0.5, 'timeBwProduct', 4, ...
    'system', sys_lowPNS, 'use', 'excitation');
rfBlockDurationPsf = mr.calcDuration(rfPsf, gzPsf);
rfBlockDurationCsm = mr.calcDuration(rfCsm, gzCsm);
rfBlockEnvelopeByMode = {mr.makeDelay(rfBlockDurationPsf), ...
    mr.makeDelay(rfBlockDurationCsm)};

Nx = N(1);
Ny = N(2);
Nx_os = Nx*ro_os;
deltakX = 1/fov(1);
deltakY = 1/fov(2);
deltakZ = 1/psfFovZ;
dwell = round((ro_dur/Nx_os)/sys.adcRasterTime)*sys.adcRasterTime;
Tread = dwell*Nx_os;
assert(abs(Tread-ro_dur) <= sys.adcRasterTime, ...
    'Rasterized ADC duration does not match ro_dur.');

readoutTime = Tread+2*sys.adcDeadTime;
readoutRampTime = 250e-6;
readoutAmplitude = Nx*deltakX/readoutTime;
gxRampUp = mr.makeExtendedTrapezoid(ax.d1, ...
    'times', [0 readoutRampTime], 'amplitudes', [0 readoutAmplitude], ...
    'system', sys_lowPNS);
gxRead = mr.makeExtendedTrapezoid(ax.d1, ...
    'times', [0 readoutTime], ...
    'amplitudes', [readoutAmplitude readoutAmplitude], 'system', sys_lowPNS);
gxRampDown = mr.makeExtendedTrapezoid(ax.d1, ...
    'times', [0 readoutRampTime], 'amplitudes', [readoutAmplitude 0], ...
    'system', sys_lowPNS);
adc = mr.makeAdc(Nx_os, 'Duration', Tread, ...
    'Delay', sys.adcDeadTime, 'system', sys_sequence);

[gyWave, waveAreaY, waveTiming] = defineSineReadoutWave( ...
    ax.d2, Tread, adc.delay, readoutTime, sys_wave, Ncycles, wavePolarity, ...
    gwave_max, swave_max, physical_slew_max, adc, true, ...
    isUseWaveCalibration);
gyWave.id = seq.registerGradEvent(gyWave);
assert(abs(waveAreaY) < 1e-6, ...
    'The calibration sine must have negligible zeroth moment.');

[gzWave, waveAreaZ, cosineWaveTiming] = defineCosineReadoutWave( ...
    ax.d3, Tread, adc.delay, readoutTime, sys_wave, Ncycles, wavePolarity, ...
    gwave_max, swave_max, physical_slew_max, adc, true, ...
    isUseWaveCalibration);
gzWave.id = seq.registerGradEvent(gzWave);
[gzCosPreBridge, cosinePreBridgeTiming] = makeFixedDurationThreePointBridge( ...
    ax.d3, 0, cosineWaveTiming.edgeAmplitude, ...
    -cosineWaveTiming.areaToAdcCenter, cosineBridgeDurationBudget, sys_lowPNS);
[gzCosPostBridge, cosinePostBridgeTiming] = makeFixedDurationThreePointBridge( ...
    ax.d3, cosineWaveTiming.edgeAmplitude, 0, ...
    -(gzCosPreBridge.area+waveAreaZ), cosineBridgeDurationBudget, sys_lowPNS);
gzCosPreBridge.id = seq.registerGradEvent(gzCosPreBridge);
gzCosPostBridge.id = seq.registerGradEvent(gzCosPostBridge);
cosineBridgeDuration = max(mr.calcDuration(gzCosPreBridge), ...
    mr.calcDuration(gzCosPostBridge));
cosineCenterClosureError = gzCosPreBridge.area+ ...
    cosineWaveTiming.areaToAdcCenter;
cosineCycleClosureError = gzCosPreBridge.area+waveAreaZ+ ...
    gzCosPostBridge.area;
assert(abs(mr.calcDuration(gzCosPreBridge)-cosineBridgeDurationBudget) < ...
    sys.gradRasterTime/10 && ...
    abs(mr.calcDuration(gzCosPostBridge)-cosineBridgeDurationBudget) < ...
    sys.gradRasterTime/10, ...
    'Cosine bridges do not match the accepted TSE duration.');
assert(abs(cosineCenterClosureError) <= 1e-6 && ...
    abs(cosineCycleClosureError) <= 1e-6, ...
    'Cosine Wave moments do not close.');

gxBridgeHold = mr.makeExtendedTrapezoid(ax.d1, ...
    'times', [0 cosineBridgeDuration], ...
    'amplitudes', [readoutAmplitude readoutAmplitude], 'system', sys_lowPNS);
assert(abs(mr.calcDuration(gyWave)-readoutTime) < sys.gradRasterTime/10 && ...
    abs(mr.calcDuration(gzWave)-readoutTime) < sys.gradRasterTime/10, ...
    'Wave events do not match the complete TSE readout block.');

adcCenterWithinReadout = adc.delay+0.5*adc.numSamples*adc.dwell;
gxAreaToAdcCenter = 0.5*readoutAmplitude*mr.calcDuration(gxRampUp)+ ...
    readoutAmplitude*cosineBridgeDuration+ ...
    readoutAmplitude*adcCenterWithinReadout;
groPreArea = -gxAreaToAdcCenter;
groPreMin = mr.makeTrapezoid(ax.d1, 'Area', groPreArea, ...
    'system', sys_lowPNS);

phaseAreasY = ((0:Ny-1)-floor(Ny/2))*deltakY;
phaseAreasZ = ((0:Ncalib_z-1)-floor(Ncalib_z/2))*deltakZ;
kyPsfForTiming = centerBlockIndices(Ny, Ncalib_y);
gyPsfPreMax = mr.makeTrapezoid(ax.d2, ...
    'Area', max(abs(phaseAreasY(kyPsfForTiming))), 'system', sys_lowPNS);
kyAcsForTiming = centerBlockIndices(Ny, Nacs_y);
gyCsmPreMax = mr.makeTrapezoid(ax.d2, ...
    'Area', max(abs(phaseAreasY(kyAcsForTiming))), 'system', sys_lowPNS);
gzPsfPreMin = [
    mr.makeTrapezoid(ax.d3, 'Area', gzPsfReph.area+min(phaseAreasZ), ...
        'system', sys_lowPNS), ...
    mr.makeTrapezoid(ax.d3, 'Area', gzPsfReph.area+max(phaseAreasZ), ...
        'system', sys_lowPNS)];
preGradDurationPsf = ceil(max([mr.calcDuration(groPreMin), ...
    mr.calcDuration(gyPsfPreMax), mr.calcDuration(gzPsfPreMin(1)), ...
    mr.calcDuration(gzPsfPreMin(2))])/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
preGradDurationCsm = ceil(max([mr.calcDuration(groPreMin), ...
    mr.calcDuration(gyCsmPreMax), mr.calcDuration(gzCsmReph)])/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
groPrePsf = mr.makeTrapezoid(ax.d1, 'Area', groPreArea, ...
    'Duration', preGradDurationPsf, 'system', sys_lowPNS);
groPreCsm = mr.makeTrapezoid(ax.d1, 'Area', groPreArea, ...
    'Duration', preGradDurationCsm, 'system', sys_lowPNS);

readoutArea = readoutAmplitude*(readoutTime+readoutRampTime);
groSpoilMin = mr.makeTrapezoid(ax.d1, ...
    'Area', ro_spoil*readoutArea/2, 'system', sys_lowPNS);
gyPsfPostMin = mr.makeTrapezoid(ax.d2, ...
    'Area', max(abs(phaseAreasY(kyPsfForTiming))), 'system', sys_lowPNS);
gyCsmPostMin = mr.makeTrapezoid(ax.d2, ...
    'Area', max(abs(phaseAreasY(kyAcsForTiming))), 'system', sys_lowPNS);
psfSpoilAreaZ = 4/psfSlabThickness;
csmSpoilAreaZ = 4/sliceThickness;
gzPsfPostMin = [
    mr.makeTrapezoid(ax.d3, 'Area', psfSpoilAreaZ-min(phaseAreasZ), ...
        'system', sys_lowPNS), ...
    mr.makeTrapezoid(ax.d3, 'Area', psfSpoilAreaZ-max(phaseAreasZ), ...
        'system', sys_lowPNS)];
gzCsmPostMin = mr.makeTrapezoid(ax.d3, ...
    'Area', csmSpoilAreaZ, 'system', sys_lowPNS);
minimumPostDurationPsf = max([mr.calcDuration(groSpoilMin), ...
    mr.calcDuration(gyPsfPostMin), mr.calcDuration(gzPsfPostMin(1)), ...
    mr.calcDuration(gzPsfPostMin(2))]);
minimumPostDurationPsf = ceil(minimumPostDurationPsf/sys.gradRasterTime)* ...
    sys.gradRasterTime;
minimumPostDurationCsm = max([mr.calcDuration(groSpoilMin), ...
    mr.calcDuration(gyCsmPostMin), mr.calcDuration(gzCsmPostMin)]);
minimumPostDurationCsm = ceil(minimumPostDurationCsm/sys.gradRasterTime)* ...
    sys.gradRasterTime;

rfCenterToEndPsf = rfBlockDurationPsf-(rfPsf.delay+mr.calcRfCenter(rfPsf));
rfCenterToEndCsm = rfBlockDurationCsm-(rfCsm.delay+mr.calcRfCenter(rfCsm));
readoutLeadAfterPrephase = mr.calcDuration(gxRampUp)+ ...
    cosineBridgeDuration+adcCenterWithinReadout;

% Derive the shortest TE independently for the slab PSF and direct-2D CSM.
% Pair-matched Wave/reference SETs still use identical timing.
TEPsf = ceil((rfCenterToEndPsf+preGradDurationPsf+ ...
    readoutLeadAfterPrephase)/sys.gradRasterTime)*sys.gradRasterTime;
TECsm = ceil((rfCenterToEndCsm+preGradDurationCsm+ ...
    readoutLeadAfterPrephase)/sys.gradRasterTime)*sys.gradRasterTime;
delayTEPsf = round((TEPsf-rfCenterToEndPsf-readoutLeadAfterPrephase)/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
delayTECsm = round((TECsm-rfCenterToEndCsm-readoutLeadAfterPrephase)/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
assert(delayTEPsf >= preGradDurationPsf && ...
    delayTECsm >= preGradDurationCsm, ...
    'A derived TE is too short for its mode-specific prephaser.');
fixedDurationAfterPrephase = mr.calcDuration(gxRampUp)+ ...
    2*cosineBridgeDuration+readoutTime+mr.calcDuration(gxRampDown);

% The slab PSF uses its hardware-limited minimum TR. The slice-selective CSM
% keeps a configurable global RF-duty floor; its physical-slice revisit time
% is Nslices times the RF-to-RF shot spacing.
baseDurationPsf = rfBlockDurationPsf+delayTEPsf+fixedDurationAfterPrephase;
baseDurationCsm = rfBlockDurationCsm+delayTECsm+fixedDurationAfterPrephase;
minimumShotDurationPsf = ceil((baseDurationPsf+minimumPostDurationPsf)/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
minimumShotDurationCsm = ceil((baseDurationCsm+minimumPostDurationCsm)/ ...
    sys.gradRasterTime)*sys.gradRasterTime;
shotDurationPsf = minimumShotDurationPsf;
shotDurationCsm = max(minimumShotDurationCsm, ...
    ceil(csmShotSpacingFloor/sys.gradRasterTime)*sys.gradRasterTime);
spoilDurationPsf = shotDurationPsf-baseDurationPsf;
spoilDurationCsm = shotDurationCsm-baseDurationCsm;
assert(spoilDurationPsf >= minimumPostDurationPsf && ...
    spoilDurationCsm >= minimumPostDurationCsm, ...
    'A derived TR is too short for its mode-specific rewinders and spoilers.');
csmSliceTR = Nslices*shotDurationCsm;
fprintf(['Derived PSF TE / TR = %.6f / %.6f ms; ', ...
    'CSM TE / RF spacing / slice revisit = %.6f / %.6f / %.6f ms.\n'], ...
    TEPsf*1e3, shotDurationPsf*1e3, TECsm*1e3, ...
    shotDurationCsm*1e3, csmSliceTR*1e3);

%% Scientific acquisition parts and compact raw-data mapping

MODE_NOWAVE = 1;
MODE_WAVE_Y = 2;
MODE_WAVE_Z = 3;
RF_PSF = 1;
RF_CSM = 2;
kyCalib = centerBlockIndices(Ny, Ncalib_y);
kyCenter = centerBlockIndices(Ny, 1);
kyAcs = centerBlockIndices(Ny, Nacs_y);
kzCalib = 1:Ncalib_z;
kzCenter = centerBlockIndices(Ncalib_z, 1);

parts = struct('id', {}, 'name', {}, 'mode', {}, 'rfMode', {}, ...
    'kyList', {}, 'kzList', {}, 'zPositions', {}, 'slcLabels', {}, ...
    'nDummy', {}, 'isRef', {});
if isUseWaveCalibration
    parts(1) = makeMixedPart(0, 'nowave_y_projection_3d', ...
        MODE_NOWAVE, RF_PSF, kyCalib, kzCenter, psfSlabCenter, 0, ...
        NdummyPsf, false);
    parts(2) = makeMixedPart(1, 'wave_y_projection_3d', ...
        MODE_WAVE_Y, RF_PSF, kyCalib, kzCenter, psfSlabCenter, 0, ...
        NsettlePsfPerSet, false);
    parts(3) = makeMixedPart(2, 'nowave_z_projection_3d', ...
        MODE_NOWAVE, RF_PSF, kyCenter, kzCalib, psfSlabCenter, 0, ...
        NsettlePsfPerSet, false);
    parts(4) = makeMixedPart(3, 'wave_z_projection_3d', ...
        MODE_WAVE_Z, RF_PSF, kyCenter, kzCalib, psfSlabCenter, 0, ...
        NsettlePsfPerSet, false);
    parts(5) = makeMixedPart(4, 'nowave_slice_csm_2d', ...
        MODE_NOWAVE, RF_CSM, kyAcs, kzCenter, ...
        slicePositions(sliceAcquisitionOrder), ...
        sliceAcquisitionOrder-1, NdummyCsm, true);
else
    parts(1) = makeMixedPart(4, 'nowave_slice_csm_2d', ...
        MODE_NOWAVE, RF_CSM, kyAcs, kzCenter, ...
        slicePositions(sliceAcquisitionOrder), ...
        sliceAcquisitionOrder-1, NdummyCsm, true);
end

acqTable = struct('partArrayIdx', {}, 'setID', {}, 'mode', {}, ...
    'rfMode', {}, 'isRef', {}, 'kyPhysicalIndex', {}, ...
    'kzPhysicalIndex', {}, 'linLocal', {}, 'parLocal', {}, ...
    'slcLabel', {}, 'zPosition', {}, 'rawLin0', {});
partStart = zeros(1, numel(parts));
partStop = zeros(1, numel(parts));
for p = 1:numel(parts)
    partStart(p) = numel(acqTable)+1;
    for iKy = 1:numel(parts(p).kyList)
        for iKz = 1:numel(parts(p).kzList)
            for iPosition = 1:numel(parts(p).zPositions)
                row.partArrayIdx = p;
                row.setID = parts(p).id;
                row.mode = parts(p).mode;
                row.rfMode = parts(p).rfMode;
                row.isRef = parts(p).isRef;
                row.kyPhysicalIndex = parts(p).kyList(iKy);
                row.kzPhysicalIndex = parts(p).kzList(iKz);
                row.linLocal = iKy-1;
                row.parLocal = iKz-1;
                row.slcLabel = parts(p).slcLabels(iPosition);
                row.zPosition = parts(p).zPositions(iPosition);
                row.rawLin0 = numel(acqTable); % zero-based storage index, not physical ky
                acqTable(end+1) = row; %#ok<SAGROW>
            end
        end
    end
    partStop(p) = numel(acqTable);
end
expectedAdcCount = numel(acqTable);
storageLayoutVersion = 'compact_adc_lin_v1';
psfAdcIndices = zeros(4, max(Ncalib_y, Ncalib_z));
csmAdcIndices = zeros(Nslices, Nacs_y);
for iAdc = 1:expectedAdcCount
    row = acqTable(iAdc);
    if row.rfMode == RF_PSF
        projectionColumn = max(row.linLocal, row.parLocal)+1;
        assert(psfAdcIndices(row.setID+1, projectionColumn) == 0, ...
            'Duplicate PSF mapping at ADC %d.', iAdc);
        psfAdcIndices(row.setID+1, projectionColumn) = iAdc;
    else
        assert(csmAdcIndices(row.slcLabel+1, row.linLocal+1) == 0, ...
            'Duplicate CSM mapping at ADC %d.', iAdc);
        csmAdcIndices(row.slcLabel+1, row.linLocal+1) = iAdc;
    end
end
psfValidMask = psfAdcIndices > 0;
psfValidLengths = sum(psfValidMask, 2).';
imageAdcIndices = find(~[acqTable.isRef]);
refAdcIndices = find([acqTable.isRef]);
mappedAdcIndices = [psfAdcIndices(psfValidMask); csmAdcIndices(:)];
assert(all(csmAdcIndices(:) > 0) && ...
       isequal(sort(mappedAdcIndices), (1:expectedAdcCount).'), ...
    'PSF/CSM sidecar tables must map every ADC exactly once.');
totalDummyShots = sum([parts.nDummy]);
partAcqCounts = partStop-partStart+1;
shotDurationByMode = {shotDurationPsf, shotDurationCsm};
partShotDurations = arrayfun( ...
    @(part) shotDurationByMode{part.rfMode}, parts);
estimatedSequenceDuration = sum( ...
    (partAcqCounts+[parts.nDummy]).*partShotDurations);
assert(estimatedSequenceDuration < maximumDesiredSequenceDuration, ...
    ['Estimated calibration duration %.3f s exceeds the %.3f s design ', ...
    'limit. Review timing or acquisition counts before writing.'], ...
    estimatedSequenceDuration, maximumDesiredSequenceDuration);

lblLIN = arrayfun(@(v) mr.makeLabel('SET', 'LIN', v), ...
    0:expectedAdcCount-1, 'UniformOutput', false);
lblPAR0 = mr.makeLabel('SET', 'PAR', 0);
lblSLC0 = mr.makeLabel('SET', 'SLC', 0);
lblSET0 = mr.makeLabel('SET', 'SET', 0);
lblECO = mr.makeLabel('SET', 'ECO', 0);
lblRefOn = mr.makeLabel('SET', 'REF', true);
lblRefOff = mr.makeLabel('SET', 'REF', false);

%% RF-spoiled acquisition loop

rfByMode = {rfPsf, rfCsm};
gzByMode = {gzPsf, gzCsm};
preGradDurationByMode = {preGradDurationPsf, preGradDurationCsm};
groPreByMode = {groPrePsf, groPreCsm};
delayTEByMode = {delayTEPsf, delayTECsm};
spoilDurationByMode = {spoilDurationPsf, spoilDurationCsm};

tic;
for p = 1:numel(parts)
    partRows = acqTable(partStart(p):partStop(p));
    dummyRowIndices = mod((-parts(p).nDummy:-1), numel(partRows))+1;
    rf_phase = 0;
    rf_inc = 0;
    for iDummy = 1:numel(dummyRowIndices)
        row = partRows(dummyRowIndices(iDummy));
        addMixedDimShot(seq, row, false, [], rfByMode, gzByMode, ...
            rfBlockEnvelopeByMode, delayTEByMode, rf_phase/180*pi, ...
            gzPsfReph, gzCsmReph, ...
            phaseAreasY, phaseAreasZ, preGradDurationByMode, ...
            groPreByMode, gxRampUp, ...
            gxBridgeHold, gxRead, gxRampDown, gyWave, gzWave, ...
            gzCosPreBridge, gzCosPostBridge, spoilDurationByMode, ...
            ro_spoil*readoutArea/2, psfSpoilAreaZ, csmSpoilAreaZ, ...
            sys_lowPNS, RF_PSF, MODE_WAVE_Y, MODE_WAVE_Z);
        rf_inc = mod(rf_inc+rfSpoilingInc, 360);
        rf_phase = mod(rf_phase+rf_inc, 360);
    end

    % Reset identical Wave/reference SETs to the same RF/ADC phase evolution.
    rf_phase = 0;
    rf_inc = 0;
    for kk = partStart(p):partStop(p)
        row = acqTable(kk);
        adcThis = adc;
        adcThis.phaseOffset = rf_phase/180*pi;
        if row.isRef
            refLabel = lblRefOn;
        else
            refLabel = lblRefOff;
        end
        labels = {lblLIN{row.rawLin0+1}, lblPAR0, ...
            lblSLC0, lblSET0, lblECO, refLabel};
        addMixedDimShot(seq, row, true, {adcThis, labels}, ...
            rfByMode, gzByMode, rfBlockEnvelopeByMode, delayTEByMode, ...
            rf_phase/180*pi, ...
            gzPsfReph, gzCsmReph, phaseAreasY, phaseAreasZ, ...
            preGradDurationByMode, groPreByMode, gxRampUp, ...
            gxBridgeHold, gxRead, ...
            gxRampDown, gyWave, gzWave, gzCosPreBridge, ...
            gzCosPostBridge, spoilDurationByMode, ...
            ro_spoil*readoutArea/2, psfSpoilAreaZ, csmSpoilAreaZ, ...
            sys_lowPNS, RF_PSF, MODE_WAVE_Y, MODE_WAVE_Z);
        rf_inc = mod(rf_inc+rfSpoilingInc, 360);
        rf_phase = mod(rf_phase+rf_inc, 360);
    end
end
fprintf('Sequence blocks generated in %.3f s.\n', toc);

%% Timing, labels, trajectory, and pair validation

[ok, error_report] = seq.checkTiming;
if ~ok
    error('Pulseq timing check failed:\n%s', strjoin(error_report, '\n'));
end
if stopAfterTimingCheck
    fprintf('Timing check passed; stopAfterTimingCheck=true.\n');
    return;
end

adc_lbl = seq.evalLabels('evolution', 'adc');
requiredLabels = {'LIN', 'PAR', 'SLC', 'SET', 'ECO', 'REF'};
for iLabel = 1:numel(requiredLabels)
    assert(isfield(adc_lbl, requiredLabels{iLabel}), ...
        'Required ADC label %s was not found.', requiredLabels{iLabel});
end
assert(numel(adc_lbl.LIN) == expectedAdcCount, ...
    'Unexpected acquired ADC count.');
assert(isequal(adc_lbl.LIN(:), (0:expectedAdcCount-1).') && ...
    all(adc_lbl.PAR(:) == 0) && all(adc_lbl.SLC(:) == 0) && ...
    all(adc_lbl.SET(:) == 0) && all(adc_lbl.ECO(:) == 0) && ...
    isequal(logical(adc_lbl.REF(:)), logical([acqTable.isRef].')), ...
    'ADC label evolution does not match the compact storage contract.');
tuples = [adc_lbl.SET(:), adc_lbl.SLC(:), adc_lbl.LIN(:), adc_lbl.PAR(:)];
assert(size(unique(tuples, 'rows'), 1) == expectedAdcCount, ...
    'Duplicate [SET,SLC,LIN,PAR] tuples were found.');
storageLabelExtents = max(tuples, [], 1)+1;
assert(prod(storageLabelExtents) == expectedAdcCount, ...
    'Scanner label bounds introduce empty Cartesian-product storage slots.');
fprintf('Compact ADC storage: %d records, LIN=0...%d, PAR/SLC/SET/ECO=0.\n', ...
    expectedAdcCount, expectedAdcCount-1);

warning('off', 'mr:restoreShape');
[ktraj_adc, t_adc, ~, ~, t_excitation] = seq.calculateKspacePP();
warning('on', 'mr:restoreShape');
assert(size(ktraj_adc, 2) == expectedAdcCount*Nx_os, ...
    'Trajectory sample count does not match ADC labels.');
ktrajByReadout = reshape(ktraj_adc, 3, Nx_os, expectedAdcCount);
assert(all(squeeze(ktrajByReadout(1,end,:)) > ...
    squeeze(ktrajByReadout(1,1,:))), ...
    'Readout polarity is not consistently positive x.');

if isUseWaveCalibration
    yRefRows = partStart(1):partStop(1);
    yWaveRows = partStart(2):partStop(2);
    zRefRows = partStart(3):partStop(3);
    zWaveRows = partStart(4):partStop(4);
    assert(numel(yRefRows) == Ncalib_y && numel(zRefRows) == Ncalib_z, ...
        'Projection SET dimensions are inconsistent.');
    pairedY = ktrajByReadout(:,:,yWaveRows)-ktrajByReadout(:,:,yRefRows);
    pairedZ = ktrajByReadout(:,:,zWaveRows)-ktrajByReadout(:,:,zRefRows);
    assert(max(abs(pairedY(1,:,:)), [], 'all') < 1e-5 && ...
        max(abs(pairedY(3,:,:)), [], 'all') < 1e-5, ...
        'SET 0/1 differ outside the physical y Wave axis.');
    assert(max(abs(pairedZ(1,:,:)), [], 'all') < 1e-5 && ...
        max(abs(pairedZ(2,:,:)), [], 'all') < 1e-5, ...
        'SET 2/3 differ outside the physical z Wave axis.');
    assert(max(abs(pairedY(2,:,:)-pairedY(2,:,1)), [], 'all') < ...
        deltakY*1e-4, 'The y-Wave difference varies across ky rows.');
    assert(max(abs(pairedZ(3,:,:)-pairedZ(3,:,1)), [], 'all') < ...
        deltakZ*1e-4, 'The z-Wave difference varies across kz rows.');
end

tAdcByReadout = reshape(t_adc, Nx_os, expectedAdcCount);
adcCenters = mean(tAdcByReadout, 1);
excCenters = t_excitation(isfinite(t_excitation));
measuredTE = zeros(1, expectedAdcCount);
for iReadout = 1:expectedAdcCount
    iExc = find(excCenters < adcCenters(iReadout), 1, 'last');
    measuredTE(iReadout) = adcCenters(iReadout)-excCenters(iExc);
end
expectedTEByMode = [TEPsf TECsm];
expectedTE = arrayfun( ...
    @(row) expectedTEByMode(row.rfMode), acqTable);
assert(max(abs(measuredTE-expectedTE)) <= max(1e-9, adc.dwell/100), ...
    'Measured RF-center TE differs from its RF-mode timing contract.');
for p = 1:numel(parts)
    centersThisSet = adcCenters(partStart(p):partStop(p));
    if numel(centersThisSet) > 1
        assert(max(abs(diff(centersThisSet)-partShotDurations(p))) <= ...
            max(1e-9, adc.dwell/100), ...
            'ADC shot spacing is not constant inside SET %d.', parts(p).id);
    end
end

%% PNS/CNS and forbidden-frequency checks

% The gradient model is the system_asc_file path setting, which must match
% the sys_type entry below.
if execute_system_check_code
    sysMap = {
        'prisma',             'MP_GPA_K2309_2250V_951A_AS82.asc'
        'prisma_XA30A',       'MP_GPA_K2309_2250V_951A_AS82_XA30A_mod_twoFilesCombined.asc'
        'skyra',              'MP_GPA_K2309_2250V_793A_GC99.asc'
        'Connectome2',        'MP_GradSys_P034_c_CX600.asc'
        'C2_simulate_prisma', 'MP_GradSys_P034_c_CX600.asc'
        'CimaX',              'MP_GradSys_P034_X60.asc'
        'TerraX',             'MP_GradSys_K2298_2250V_1250A_W60_SC72CD.asc'
    };
    idx = find(strcmp(sys_type, sysMap(:,1)), 1);
    assert(~isempty(idx), 'No hardware model is configured for %s.', sys_type);
    ascName = requireSystemCheckInputs(sys_type, sysMap{idx,2}, ...
        system_asc_file, safe_pns_prediction_path);
    warning('off', 'mr:restoreShape');
    [~, tpns] = seq.calcPNS(ascName, false, 0);
    isHasCNS = strcmp(sys_type, 'CimaX');
    if isHasCNS
        [~, tcns] = seq.calcPNS(ascName, false, 1);
    else
        tcns = NaN;
    end
    warning('on', 'mr:restoreShape');
    assert(max(tpns) <= 0.95, ...
        'PNS exceeds the 0.95 acceptance threshold.');
    if isHasCNS
        assert(max(tcns) <= 0.95, ...
            'CNS exceeds the 0.95 acceptance threshold.');
    end
    forbiddenFreqCheck(seq, sys, ascName);
end

%% Definitions and write

calibrationMetadata = struct();
calibrationMetadata.StorageLayoutVersion = storageLayoutVersion;
calibrationMetadata.CalibrationRFDesign = 'Pulseq makeSincPulse';
calibrationMetadata.PsfRFDuration = rfDurationPsf;
calibrationMetadata.PsfRFTimeBandwidthProduct = 20;
calibrationMetadata.PsfRFApodization = 0.42;
calibrationMetadata.CsmRFDuration = rfDurationCsm;
calibrationMetadata.CsmRFTimeBandwidthProduct = 4;
calibrationMetadata.CsmRFApodization = 0.5;
calibrationMetadata.WaveEnabled = double(isUseWaveCalibration);
calibrationMetadata.WaveDesignedAmplitude_mTm = waveTiming.actualAmplitude_mTm;
calibrationMetadata.WaveEmittedAmplitude_mTm = ...
    double(isUseWaveCalibration)*waveTiming.actualAmplitude_mTm;
calibrationMetadata.CalibrationDimensionality = ...
    'mixed-3d-psf-projection-plus-2d-direct-slice-csm';
calibrationMetadata.FOV = [fov psfFovZ];
calibrationMetadata.Matrix = [N Ncalib_z];
calibrationMetadata.OrientationMapping = slOrientation;
calibrationMetadata.ReadoutAxis = ax.d1;
calibrationMetadata.PhaseEncodeAxis = ax.d2;
calibrationMetadata.PartitionEncodeAxis = ax.d3;
calibrationMetadata.Nx = Nx;
calibrationMetadata.Ny = Ny;
calibrationMetadata.Nx_os = Nx_os;
calibrationMetadata.ro_os = ro_os;
calibrationMetadata.ReadoutDuration = Tread;
calibrationMetadata.dwell = dwell;
calibrationMetadata.TE = [TEPsf TECsm];
calibrationMetadata.TR = [shotDurationPsf shotDurationCsm];
calibrationMetadata.PsfTE = TEPsf;
calibrationMetadata.PsfTR = shotDurationPsf;
calibrationMetadata.CsmTE = TECsm;
calibrationMetadata.CsmRfToRfSpacing = shotDurationCsm;
calibrationMetadata.CsmSliceRevisitTime = csmSliceTR;
calibrationMetadata.CsmRfToRfSpacingFloor = csmShotSpacingFloor;
calibrationMetadata.PsfFlipAngle = alphaPsf;
calibrationMetadata.CsmFlipAngle = alphaCsm;
calibrationMetadata.PsfFOV = [fov psfFovZ];
calibrationMetadata.PsfMatrix = [Nx Ny Ncalib_z];
calibrationMetadata.PsfSlabCenter = psfSlabCenter;
calibrationMetadata.PsfSlabThickness = psfSlabThickness;
calibrationMetadata.PsfYCoordinates = psfYCoordinates;
calibrationMetadata.PsfZCoordinates = psfZCoordinates;
calibrationMetadata.PsfZIsFourierEncoded = 1;
calibrationMetadata.PsfSlabAxis = ax.d3;
calibrationMetadata.Ncalib_y = Ncalib_y;
calibrationMetadata.Ncalib_z = Ncalib_z;
calibrationMetadata.CSMSetID = 4;
calibrationMetadata.CSMIsDirect2D = 1;
calibrationMetadata.Nacs_y = Nacs_y;
calibrationMetadata.NPhysicalSlices = Nslices;
calibrationMetadata.Nslices = Nslices;
calibrationMetadata.SliceThickness = sliceThickness;
calibrationMetadata.SliceGap = sliceSpacing;
calibrationMetadata.SliceCenterSpacing = sliceCenterSpacing;
calibrationMetadata.SlicePositions = slicePositions;
calibrationMetadata.CSMFOV = fov;
calibrationMetadata.CSMMatrix = N;
calibrationMetadata.CSMSlicePositions = slicePositions;
calibrationMetadata.CSMSliceThickness = sliceThickness;
calibrationMetadata.CSMSliceGap = sliceSpacing;
calibrationMetadata.PhysicalSlicePositions = slicePositions;
calibrationMetadata.TSESliceThickness = sliceThickness;
calibrationMetadata.TSESliceGap = sliceSpacing;
calibrationMetadata.TSESliceCenterSpacing = sliceCenterSpacing;
calibrationMetadata.TSEStackExtent = stackExtent;
calibrationMetadata.SliceAcquisitionOrder = sliceAcquisitionOrder;
calibrationMetadata.SliceByAcquisitionIndex = sliceByAcquisitionIndex;
calibrationMetadata.AcquisitionIndexBySlice = acquisitionIndexBySlice;
calibrationMetadata.SMSGroupsFlat = smsGroups(:).';
calibrationMetadata.SMSGroupsShape = size(smsGroups);
calibrationMetadata.SMSGroupCenters = groupCenters;
calibrationMetadata.SMSGroupOrder = mbGroupOrder;
calibrationMetadata.SMSPairSeparation = pairSeparation;
calibrationMetadata.YWaveReferenceSetID = 0;
calibrationMetadata.YWaveCalibrationSetID = 1;
calibrationMetadata.ZWaveReferenceSetID = 2;
calibrationMetadata.ZWaveCalibrationSetID = 3;
calibrationMetadata.NdummyPsf = NdummyPsf;
calibrationMetadata.NsettlePsfPerSet = NsettlePsfPerSet;
calibrationMetadata.NdummyCsm = NdummyCsm;
calibrationMetadata.NdummyCsmCycles = NdummyCsm/Nslices;
calibrationMetadata.SetDummyShots = [parts.nDummy];
calibrationMetadata.EnabledSetIDs = [parts.id];
calibrationMetadata.ExpectedADCCount = expectedAdcCount;
calibrationMetadata.ExpectedSequenceDuration = estimatedSequenceDuration;
calibrationMetadata.MaximumDesiredSequenceDuration = ...
    maximumDesiredSequenceDuration;
calibrationMetadata.WaveSinChannel = ax.d2;
calibrationMetadata.WaveCosChannel = ax.d3;
calibrationMetadata.WaveRepresentation = waveRepresentation;
calibrationMetadata.WaveAmplitude_mTm = ...
    double(isUseWaveCalibration)*waveTiming.actualAmplitude_mTm;
calibrationMetadata.WaveRequestedAmplitude_mTm = gwave_max;
calibrationMetadata.WaveSlew_Tms = swave_max;
calibrationMetadata.WaveActualSlew_Tms = double(isUseWaveCalibration)*max( ...
    waveTiming.slewPeak_T_m_s, cosineWaveTiming.slewPeak_T_m_s);
calibrationMetadata.WaveFullAmplitudeRequiredSlew_Tms = ...
    2*pi*Ncycles/Tread*gwave_max*1e-3;
calibrationMetadata.WaveCycles = Ncycles;
calibrationMetadata.WaveSupportedCycles = [10 20];
calibrationMetadata.WaveCycleParity = (-1)^Ncycles;
calibrationMetadata.WavePolarity = wavePolarity;
calibrationMetadata.WaveCosBridgeDuration = cosineBridgeDuration;
calibrationMetadata.WavePairMatchedApartFromWave = 1;
calibrationMetadata.PsfCoefficientSpatialUnits = 'rad_per_m';
calibrationMetadata.PsfGridTransfer = ...
    'fit physical y/z coordinates; evaluate on actual TSE grid; no image interpolation';
calibrationMetadata.SystemType = sys_type;
calibrationMetadata.SequenceMaxGrad_mTm = sys_sequence.maxGrad/sys_sequence.gamma*1e3;
calibrationMetadata.SequenceMaxSlew_Tms = sys_sequence.maxSlew/sys_sequence.gamma;
calibrationMetadata.WaveMaxGrad_mTm = sys_wave.maxGrad/sys_wave.gamma*1e3;
calibrationMetadata.WaveMaxSlew_Tms = sys_wave.maxSlew/sys_wave.gamma;
calibrationMetadata.LowPNSMaxGrad_mTm = sys_lowPNS.maxGrad/sys_lowPNS.gamma*1e3;
calibrationMetadata.LowPNSMaxSlew_Tms = sys_lowPNS.maxSlew/sys_lowPNS.gamma;
calibrationMetadata.MBGeometryMatFile = mbPulseMatFile;
[~, mbGeometryMatName, mbGeometryMatExt] = fileparts(mbPulseMatFile);
calibrationMetadata.MBGeometryMatFilename = ...
    [mbGeometryMatName mbGeometryMatExt];
calibrationMetadata.TargetTseRFCandidateTag = candidateTag;
calibrationMetadata.TargetTseRFReferenceArraysSHA256 = candidateChecksum;
calibrationMetadata.TargetTseRFCandidateStatus = candidateStatus;
calibrationMetadata.ExecuteSystemCheckCode = double(execute_system_check_code);
calibrationMetadata.OutputPath = out_path;

for p = 1:numel(parts)
    prefix = ['Part' num2str(parts(p).id) '_'];
    calibrationMetadata.([prefix 'Name']) = parts(p).name;
    calibrationMetadata.([prefix 'SetID']) = parts(p).id;
    calibrationMetadata.([prefix 'IsREF']) = double(parts(p).isRef);
    calibrationMetadata.([prefix 'NLinLocal']) = numel(parts(p).kyList);
    calibrationMetadata.([prefix 'NParLocal']) = numel(parts(p).kzList);
    calibrationMetadata.([prefix 'SLCLabels']) = parts(p).slcLabels;
    calibrationMetadata.([prefix 'ZPositions']) = parts(p).zPositions;
    calibrationMetadata.([prefix 'NDummy']) = parts(p).nDummy;
    calibrationMetadata.([prefix 'KyPhysIndices0']) = parts(p).kyList-1;
    calibrationMetadata.([prefix 'KzPhysIndices0']) = parts(p).kzList-1;
end

% Keep only acquisition-independent, scanner-facing definitions. Pulseq adds
% its four required rasters and TotalDuration when writing. All physical
% grids, calibration identities, diagnostics, and provenance stay in the MAT.
seq.setDefinition('FOV', [fov psfFovZ]);

if isUseWaveCalibration
    modeTag = '_waveYZ3dproj_csm2d';
else
    modeTag = '_csm2donly';
end
seqFilename = ['gre_mixed_dim_mb2_calib' modeTag '_Ncyc' num2str(Ncycles) '_' ...
    num2str(Nx) 'x' num2str(Ny) 'x' num2str(Ncalib_z) ...
    '_FOV' num2str(round(fov(1)*1e3)) '_Nacs' num2str(Nacs_y) ...
    '_os' num2str(ro_os) '_' sys_type];
seq.setDefinition('Name', seqFilename);
calibrationMetadata.GeneratedSequenceName = seqFilename;
scannerDefinitionNames = keys(seq.definitions);
calibrationMetadata.ScannerDefinitions = cell2struct( ...
    values(seq.definitions, scannerDefinitionNames), scannerDefinitionNames, 2);
sequenceFile = fullfile(out_path, [seqFilename '.seq']);
metadataFile = fullfile(out_path, [seqFilename '_geometry_labels.mat']);
reportPath = fullfile(out_path, [seqFilename '_geometry_labels.txt']);

if ~writeSequence
    fprintf('Validation-only run: waveform, labels, trajectory and metadata checked; no output written.\n');
    return;
end
assert(~isfile(sequenceFile) && ~isfile(metadataFile) && ~isfile(reportPath), ...
    ['Output already exists. Preserve the historical export and select a ', ...
     'new out_path interactively before writing compact-storage data.']);
if write_v141_format
    legacySequenceFile = fullfile(out_path, 'v141', [seqFilename '_v141.seq']);
    assert(~isfile(legacySequenceFile), ...
        'Historical legacy output already exists: %s', legacySequenceFile);
end

if ~exist(out_path, 'dir'), mkdir(out_path); end
if write_v141_format
    v141Path = fullfile(out_path, 'v141');
    if ~exist(v141Path, 'dir'), mkdir(v141Path); end
    seq.write_v141(fullfile(v141Path, [seqFilename '_v141.seq']));
end
seq.write(sequenceFile);
calibrationMetadata.SequenceSignature = seq.signatureValue;
calibrationMetadata.SequenceSignatureType = seq.signatureType;

fid = fopen(reportPath, 'w');
assert(fid >= 0, 'Could not open geometry report: %s', reportPath);
fprintf(fid, 'Name: %s\n', seqFilename);
fprintf(fid, 'StorageLayoutVersion: %s\n', storageLayoutVersion);
fprintf(fid, 'ScannerLabels: LIN=0...%d; PAR=SLC=SET=ECO=0; CSM REF=1\n', ...
    expectedAdcCount-1);
fprintf(fid, 'PSF_ValidReadoutsByRow: %s\n', mat2str(psfValidLengths));
fprintf(fid, 'PSF_ADCIndexMapShape: %s (1-based ADC records; 0=unused)\n', ...
    mat2str(size(psfAdcIndices)));
fprintf(fid, 'CSM_ADCIndexMapShape: %s (physical slice by local ky)\n', ...
    mat2str(size(csmAdcIndices)));
fprintf(fid, 'CalibrationDimensionality: mixed-3d-psf-plus-2d-csm\n');
fprintf(fid, 'PSF_FOV_mm: %.9g %.9g %.9g\n', [fov psfFovZ]*1e3);
fprintf(fid, 'PSF_Matrix: %d %d %d\n', [N Ncalib_z]);
fprintf(fid, 'PSF_SlabCenterThickness_mm: %.9g %.9g\n', ...
    psfSlabCenter*1e3, psfSlabThickness*1e3);
fprintf(fid, 'PSF_YCoordinates_mm: %s\n', mat2str(psfYCoordinates*1e3, 9));
fprintf(fid, 'PSF_ZCoordinates_mm: %s\n', mat2str(psfZCoordinates*1e3, 9));
fprintf(fid, 'CSM_SlicePositions_mm: %s\n', mat2str(slicePositions*1e3, 9));
fprintf(fid, 'TSE_StackExtent_mm: %.9g\n', stackExtent*1e3);
fprintf(fid, 'SMS_PairSeparation_mm: %.9g\n', pairSeparation*1e3);
fprintf(fid, 'MBGeometryMatFilename: %s\n', ...
    [mbGeometryMatName mbGeometryMatExt]);
fprintf(fid, 'CalibrationRFDesign: Pulseq makeSincPulse (PSF TBW20, CSM TBW4)\n');
fprintf(fid, 'TargetTseRFCandidateTag: %s\n', candidateTag);
fprintf(fid, 'TargetTseRFReferenceArraysSHA256: %s\n', candidateChecksum);
fprintf(fid, 'TargetTseRFCandidateStatus: %s\n', candidateStatus);
fprintf(fid, 'CSM_SliceAcquisitionOrder_1based: %s\n', ...
    mat2str(sliceAcquisitionOrder));
fprintf(fid, 'FlipAngles_deg_PSF_CSM: %.9g %.9g\n', alphaPsf, alphaCsm);
fprintf(fid, 'DummyShotsBySet: %s\n', mat2str([parts.nDummy]));
fprintf(fid, 'ADCCount: %d\n', expectedAdcCount);
fprintf(fid, 'PSF_TE_TR_ms_derived_minimum: %.9g %.9g\n', ...
    TEPsf*1e3, shotDurationPsf*1e3);
fprintf(fid, 'CSM_TE_RFspacing_SliceRevisit_ms: %.9g %.9g %.9g\n', ...
    TECsm*1e3, shotDurationCsm*1e3, csmSliceTR*1e3);
fprintf(fid, 'EstimatedSequenceDuration_s: %.9g\n', ...
    estimatedSequenceDuration);
fprintf(fid, 'WaveCycles: %d\n', Ncycles);
fprintf(fid, 'WaveRequestedDesignedEmittedAmplitude_mTm: %.9g %.9g %.9g\n', ...
    gwave_max, waveTiming.actualAmplitude_mTm, ...
    calibrationMetadata.WaveEmittedAmplitude_mTm);
fprintf(fid, 'WaveRequestedActualSlew_Tms: %.9g %.9g\n', ...
    swave_max, calibrationMetadata.WaveActualSlew_Tms);
fprintf(fid, 'SequenceWaveLowPNSMaxSlew_Tms: %.9g %.9g %.9g\n', ...
    sys_sequence.maxSlew/sys_sequence.gamma, ...
    sys_wave.maxSlew/sys_wave.gamma, ...
    sys_lowPNS.maxSlew/sys_lowPNS.gamma);
fprintf(fid, 'WaveCosBridgeDuration_ms: %.9g\n', cosineBridgeDuration*1e3);
for p = 1:numel(parts)
    fprintf(fid, ['ScientificPart%d: %s, localKy=%d, localKz=%d, physicalSliceLabels=%s, REF=%d, ' ...
        'dummy=%d\n'], parts(p).id, parts(p).name, ...
        numel(parts(p).kyList), numel(parts(p).kzList), ...
        mat2str(parts(p).slcLabels), parts(p).isRef, parts(p).nDummy);
end
fclose(fid);

save(metadataFile, ...
    'storageLayoutVersion', 'storageLabelExtents', 'calibrationMetadata', ...
    'psfAdcIndices', 'psfValidMask', 'psfValidLengths', 'csmAdcIndices', ...
    'imageAdcIndices', 'refAdcIndices', 'partStart', 'partStop', ...
    'acqTable', 'parts', 'slicePositions', 'sliceAcquisitionOrder', ...
    'sliceByAcquisitionIndex', 'acquisitionIndexBySlice', ...
    'psfYCoordinates', 'psfZCoordinates', 'psfFovZ', ...
    'psfSlabCenter', 'psfSlabThickness', 'phaseAreasY', 'phaseAreasZ', ...
    'smsGroups', 'groupCenters', 'mbGroupOrder', 'pairSeparation', ...
    'sliceThickness', 'sliceSpacing', 'sliceCenterSpacing', 'stackExtent', ...
    'mbPulseMatFile', 'candidateTag', 'candidateChecksum', 'candidateStatus', ...
    'Ncalib_y', 'Ncalib_z', 'Nacs_y', 'Nx', 'Ny', 'Nx_os', 'Nslices', ...
    'fov', 'N', 'ro_os', 'dwell', 'slOrientation', 'ax', ...
    'rfDurationPsf', 'rfDurationCsm', 'rfSpoilingInc', ...
    'isUseWaveCalibration', 'execute_system_check_code', 'out_path', ...
    'NdummyPsf', 'NsettlePsfPerSet', 'NdummyCsm', ...
    'alphaPsf', 'alphaCsm', 'Tread', 'TEPsf', 'TECsm', ...
    'shotDurationPsf', 'shotDurationCsm', 'csmSliceTR', ...
    'csmShotSpacingFloor', 'maximumDesiredSequenceDuration', ...
    'expectedAdcCount', 'totalDummyShots', 'estimatedSequenceDuration', ...
    'waveRepresentation', 'gwave_max', 'swave_max', 'Ncycles', ...
    'wavePolarity', 'waveAreaY', 'waveAreaZ', 'waveTiming', ...
    'cosineWaveTiming', 'cosinePreBridgeTiming', ...
    'cosinePostBridgeTiming', 'cosineCenterClosureError', ...
    'cosineCycleClosureError', 'sys_sequence', 'sys_wave', 'sys_lowPNS');

fprintf('Wrote sequence: %s\n', fullfile(out_path, [seqFilename '.seq']));
fprintf('Wrote geometry/label report: %s\n', reportPath);

return;

%% Local helpers

function part = makeMixedPart(id, name, mode, rfMode, kyList, kzList, ...
        zPositions, slcLabels, nDummy, isRef)
%MAKEMIXEDPART Create one scientific calibration part description.
% Inputs specify sidecar part identity, RF/Wave mode, physical ky/kz indices,
% direct slice positions and physical indices, dummy count, and REF routing.
% The output is one
% scalar structure. Invalid slice-position/label lengths raise an assertion.
    assert(numel(zPositions) == numel(slcLabels), ...
        'Each part must provide one physical slice index per direct RF position.');
    part.id = id;
    part.name = name;
    part.mode = mode;
    part.rfMode = rfMode;
    part.kyList = kyList;
    part.kzList = kzList;
    part.zPositions = zPositions;
    part.slcLabels = slcLabels;
    part.nDummy = nDummy;
    part.isRef = isRef;
end

function idx = centerBlockIndices(Ndim, nCenter)
%CENTERBLOCKINDICES Return centered one-based indices within a dimension.
% Inputs are a full dimension and requested width. The output is a row-vector
% index list; an oversized or nonpositive request raises an assertion.
    assert(nCenter > 0 && nCenter <= Ndim, ...
        'Requested center block does not fit its dimension.');
    centerIdx = floor(Ndim/2)+1;
    if mod(nCenter,2) == 0
        idx = (centerIdx-nCenter/2):(centerIdx+nCenter/2-1);
    else
        idx = (centerIdx-floor(nCenter/2)):(centerIdx+floor(nCenter/2));
    end
end

function rfThis = makeSliceRf(rf, gz, position, rfPhaseRad)
%MAKESLICERF Apply a physical slab/slice center and RF-spoiling phase.
% Inputs are RF and slice-gradient templates, physical center [m], and phase
% [rad]. The output is a copied RF event with consistent frequency and phase.
    rfThis = rf;
    rfThis.freqOffset = gz.amplitude*position;
    rfThis.phaseOffset = rfPhaseRad- ...
        2*pi*rfThis.freqOffset*mr.calcRfCenter(rfThis);
end

function addMixedDimShot(seq, row, acquire, adcAndLabels, rfByMode, ...
        gzByMode, rfBlockEnvelopeByMode, delayTEByMode, rfPhaseRad, ...
        gzPsfReph, gzCsmReph, ...
        phaseAreasY, phaseAreasZ, preGradDurationByMode, ...
        groPreByMode, gxRampUp, ...
        gxBridgeHold, gxRead, gxRampDown, gyWave, gzWave, ...
        gzCosPreBridge, gzCosPostBridge, spoilDurationByMode, ...
        groSpoilArea, psfSpoilAreaZ, csmSpoilAreaZ, sys, ...
        RF_PSF, MODE_WAVE_Y, MODE_WAVE_Z)
%ADDMIXEDDIMSHOT Append one dummy or acquired mixed-dimensional GRE shot.
% The row structure selects RF geometry, ky/kz encoding, physical position,
% and Wave mode. acquire controls ADC/label insertion; adcAndLabels contains
% an ADC event and label cell when acquisition is enabled. The function adds
% blocks to seq and has no output.
    rfTemplate = rfByMode{row.rfMode};
    gzRf = gzByMode{row.rfMode};
    rfBlockEnvelope = rfBlockEnvelopeByMode{row.rfMode};
    preGradDuration = preGradDurationByMode{row.rfMode};
    groPre = groPreByMode{row.rfMode};
    rfThis = makeSliceRf(rfTemplate, gzRf, row.zPosition, rfPhaseRad);
    gyArea = phaseAreasY(row.kyPhysicalIndex);
    if row.rfMode == RF_PSF
        zEncodeArea = phaseAreasZ(row.kzPhysicalIndex);
        gzPreArea = gzPsfReph.area+zEncodeArea;
        gzPostArea = psfSpoilAreaZ-zEncodeArea;
    else
        zEncodeArea = 0; %#ok<NASGU>
        gzPreArea = gzCsmReph.area;
        gzPostArea = csmSpoilAreaZ;
    end
    gyPre = mr.makeTrapezoid('y', 'Area', gyArea, ...
        'Duration', preGradDuration, 'system', sys);
    gzPre = mr.makeTrapezoid('z', 'Area', gzPreArea, ...
        'Duration', preGradDuration, 'system', sys);
    gyPost = mr.makeTrapezoid('y', 'Area', -gyArea, ...
        'Duration', spoilDurationByMode{row.rfMode}, 'system', sys);
    gzPost = mr.makeTrapezoid('z', 'Area', gzPostArea, ...
        'Duration', spoilDurationByMode{row.rfMode}, 'system', sys);
    gxPost = mr.makeTrapezoid('x', 'Area', groSpoilArea, ...
        'Duration', spoilDurationByMode{row.rfMode}, 'system', sys);

    seq.addBlock(rfThis, gzRf, rfBlockEnvelope);
    seq.addBlock(mr.align('left', ...
        mr.makeDelay(delayTEByMode{row.rfMode}), ...
        'right', groPre, gyPre, gzPre));
    seq.addBlock(gxRampUp);
    if row.mode == MODE_WAVE_Z
        seq.addBlock(gxBridgeHold, gzCosPreBridge);
    else
        seq.addBlock(gxBridgeHold);
    end

    if acquire
        adcThis = adcAndLabels{1};
        labels = adcAndLabels{2};
        if row.mode == MODE_WAVE_Y
            seq.addBlock(gxRead, gyWave, adcThis, labels{:});
        elseif row.mode == MODE_WAVE_Z
            seq.addBlock(gxRead, gzWave, adcThis, labels{:});
        else
            seq.addBlock(gxRead, adcThis, labels{:});
        end
    else
        if row.mode == MODE_WAVE_Y
            seq.addBlock(gxRead, gyWave);
        elseif row.mode == MODE_WAVE_Z
            seq.addBlock(gxRead, gzWave);
        else
            seq.addBlock(gxRead);
        end
    end

    if row.mode == MODE_WAVE_Z
        seq.addBlock(gxBridgeHold, gzCosPostBridge);
    else
        seq.addBlock(gxBridgeHold);
    end
    seq.addBlock(gxRampDown);
    seq.addBlock(gxPost, gyPost, gzPost);
end

function [gWave, area, timing] = defineSineReadoutWave(channel, Tread, adcDelay, ...
        readoutBlockDur, sys, Ncycles, wavePolarity, gwave_max, swave_max, ...
        physical_slew_max, adc, waveInfoFlag, debugFlag)
%DEFINESINEREADOUTWAVE Construct an ADC-aligned full-block sine gradient.
% Inputs describe the axis, active ADC timing, complete readout-block
% duration, polarity, scanner limits, and diagnostic flags. Outputs are the zero-delay
% arbitrary gradient, zeroth moment, and realized timing metadata. Assertions
% reject raster, ADC, block-duration, slew, or amplitude inconsistencies.
    dt = sys.gradRasterTime;
    assert(ismember(wavePolarity, [-1 1]), 'wavePolarity must be +1 or -1.');
    [G0, w, TreadRaster, actualAmplitude_mTm] = designWaveAmplitude( ...
        Tread, sys, Ncycles, gwave_max, swave_max, physical_slew_max, waveInfoFlag);
    nReadSamples = round(TreadRaster/dt);
    assert(abs(nReadSamples*dt-TreadRaster) < dt/10, ...
        'TreadRaster must contain an integer number of gradient samples.');

    nBlockSamples = round(readoutBlockDur/dt);
    readoutBlockDurRaster = nBlockSamples*dt;
    assert(abs(readoutBlockDurRaster-readoutBlockDur) < dt/10, ...
        'The readout block duration must be on the gradient raster.');

    nPadPre = round(adcDelay/dt);
    padPreDur = nPadPre*dt;
    assert(abs(padPreDur-adcDelay) < dt/10, ...
        'The ADC delay must be on the gradient raster for full-block Wave padding.');
    nPadPost = nBlockSamples-nPadPre-nReadSamples;
    assert(nPadPost >= 0, 'The active sine does not fit in the readout block.');
    padPostDur = nPadPost*dt;
    assert(nReadSamples >= 2, ...
        'At least two active samples are required for an endpoint-zero sine wave.');

    activePhase = 2*pi*Ncycles*(0:(nReadSamples-1))/(nReadSamples-1);
    activeWave = wavePolarity*G0*sin(activePhase);
    activeWave([1 end]) = 0;
    fullBlockWave = [zeros(1, nPadPre), activeWave, zeros(1, nPadPost)];
    assert(numel(fullBlockWave) == nBlockSamples, ...
        'Internal full-block calibration sine length mismatch.');

    gWave = mr.makeArbitraryGrad(channel, fullBlockWave, 'system', sys, ...
        'first', 0, 'last', 0);
    gWave.delay = 0;
    area = gWave.area;
    slewPeak_Hz_m_s = max(abs(diff(fullBlockWave)))/dt;
    waveSlewLimit = min(swave_max, physical_slew_max);
    assertArbitraryWaveformWithinSystem( ...
        fullBlockWave, 0, 0, sys, 'active sine wave');
    assert(slewPeak_Hz_m_s/sys.gamma <= waveSlewLimit*(1+1e-9), ...
        'Active sine slew exceeds the requested Wave slew limit.');

    adcCenterTime = adc.delay+0.5*adc.numSamples*adc.dwell;
    assert(adc.numSamples == 2*floor(adc.numSamples/2), ...
        'The sampled sine-center calculation requires an even ADC sample count.');
    centerSampleNumbers = [adc.numSamples/2, adc.numSamples/2+1];
    centerSampleTimes = adc.delay+(centerSampleNumbers-0.5)*adc.dwell;
    geometricCenterMoment = integrateRasterWaveformToTime( ...
        fullBlockWave, dt, 0, 0, adcCenterTime);
    centerSampleMoments = zeros(size(centerSampleTimes));
    for iCenterSample = 1:numel(centerSampleTimes)
        centerSampleMoments(iCenterSample) = integrateRasterWaveformToTime( ...
            fullBlockWave, dt, 0, 0, centerSampleTimes(iCenterSample));
    end
    sampledCenterMoment = mean(centerSampleMoments);

    timing.nReadSamples = nReadSamples;
    timing.nBlockSamples = nBlockSamples;
    timing.nPadPre = nPadPre;
    timing.nPadPost = nPadPost;
    timing.padPreDur = padPreDur;
    timing.padPostDur = padPostDur;
    timing.Twave = readoutBlockDurRaster;
    timing.Tactive = TreadRaster;
    timing.activeStartTime = padPreDur;
    timing.activeEndTime = padPreDur+TreadRaster;
    timing.G0 = G0;
    timing.omega = w;
    timing.frequency_Hz = w/(2*pi);
    timing.polarity = wavePolarity;
    timing.slewPeak_T_m_s = slewPeak_Hz_m_s/sys.gamma;
    timing.requestedFullAmplitudeSlew_T_m_s = ...
        2*pi*Ncycles/TreadRaster*gwave_max*1e-3;
    timing.actualAmplitude_mTm = actualAmplitude_mTm;
    timing.area = area;
    timing.adcCenterTime = adcCenterTime;
    timing.geometricCenterMoment = geometricCenterMoment;
    timing.centerSampleTimes = centerSampleTimes;
    timing.centerSampleMoments = centerSampleMoments;
    timing.sampledCenterMoment = sampledCenterMoment;

    adcAcqDur = adc.numSamples*adc.dwell;
    tol = max(dt/10, sys.adcRasterTime);
    assert(abs(TreadRaster-Tread) < tol, 'Wave duration does not match Tread.');
    assert(abs(adcDelay-adc.delay) < tol, 'Wave and ADC start times differ.');
    assert(abs(TreadRaster-adcAcqDur) < tol, 'Wave and ADC active durations differ.');
    assert(abs(gWave.delay) < tol, 'Full-block arbitrary Wave must have zero delay.');
    assert(abs(timing.activeStartTime-adc.delay) < tol, ...
        'Wave active start does not match the ADC start.');
    assert(abs(timing.activeEndTime-(adc.delay+adcAcqDur)) < tol, ...
        'Wave active end does not match the ADC end.');
    assert(abs(mr.calcDuration(gWave)-readoutBlockDurRaster) < tol, ...
        'Wave event duration does not match the complete readout block.');

    if debugFlag
        fprintf(['y-sine Wave: amplitude=%.6f mT/m, slew=%.6f T/m/s, ', ...
                 'cycles=%d, active/block=%.6f/%.6f ms, pre/post=%d/%d ', ...
                 'samples, center geometric/sampled=%.9g/%.9g 1/m\n'], ...
            actualAmplitude_mTm, timing.slewPeak_T_m_s, Ncycles, ...
            TreadRaster*1e3, readoutBlockDurRaster*1e3, nPadPre, nPadPost, ...
            geometricCenterMoment, sampledCenterMoment);
    end
end

function area = integrateRasterWaveformToTime( ...
        waveform, rasterTime, firstAmplitude, lastAmplitude, endTime)
%INTEGRATERASTERWAVEFORMTOTIME Integrate a Pulseq arbitrary gradient exactly.
% Inputs are the cell-centered waveform samples [Hz/m], gradient raster [s],
% explicit first/last boundary amplitudes [Hz/m], and scalar integration end
% time [s]. The output is the piecewise-linear gradient area [1/m] from zero
% through endTime. Invalid dimensions or out-of-range times raise assertions.
    waveform = waveform(:);
    duration = numel(waveform)*rasterTime;
    assert(isscalar(endTime) && endTime >= 0 && endTime <= duration, ...
        'Gradient integration time lies outside the arbitrary event.');
    knotTimes = [0; ((0:numel(waveform)-1).'+0.5)*rasterTime; duration];
    knotAmplitudes = [firstAmplitude; waveform; lastAmplitude];
    interior = knotTimes > 0 & knotTimes < endTime;
    integrationTimes = [0; knotTimes(interior); endTime];
    integrationAmplitudes = interp1( ...
        knotTimes, knotAmplitudes, integrationTimes, 'linear');
    area = trapz(integrationTimes, integrationAmplitudes);
end

function [gWave, area, timing] = defineCosineReadoutWave(channel, Tread, adcDelay, ...
        readoutBlockDur, sys, Ncycles, wavePolarity, gwave_max, swave_max, ...
        physical_slew_max, adc, waveInfoFlag, debugFlag)
%DEFINECOSINEREADOUTWAVE Construct an ADC-aligned full-block cosine gradient.
% Inputs describe the axis, ADC/readout timing, Wave design, scanner limits,
% and diagnostic flags. Outputs are the +G0-endpoint arbitrary gradient,
% its area, and realized timing/moment metadata. The cosine samples and ADC
% guard convention exactly match the separated yz-Wave TSE source.
    dt = sys.gradRasterTime;
    assert(ismember(wavePolarity, [-1 1]), 'wavePolarity must be +1 or -1.');
    [G0, w, TreadRaster, actualAmplitude_mTm] = designWaveAmplitude( ...
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
    fullBlockWave = [Gedge*ones(1, nPadPre), activeWave, ...
        Gedge*ones(1, nPadPost)];
    assert(numel(fullBlockWave) == nBlockSamples, ...
        'Internal full-block calibration cosine length mismatch.');

    gWave = mr.makeArbitraryGrad(channel, fullBlockWave, 'system', sys, ...
        'first', Gedge, 'last', Gedge);
    gWave.delay = 0;
    area = gWave.area;
    adcCenterTime = adc.delay+0.5*adc.numSamples*adc.dwell;
    nCellsToAdcCenter = round(adcCenterTime/dt);
    assert(abs(nCellsToAdcCenter*dt-adcCenterTime) < dt/10, ...
        'ADC center must lie on the gradient raster.');
    areaToAdcCenter = sum(fullBlockWave(1:nCellsToAdcCenter))*dt;

    timing.nReadSamples = nReadSamples;
    timing.nBlockSamples = nBlockSamples;
    timing.nPadPre = nPadPre;
    timing.nPadPost = nPadPost;
    timing.Twave = nBlockSamples*dt;
    timing.Tactive = TreadRaster;
    timing.activeStartTime = nPadPre*dt;
    timing.activeEndTime = nPadPre*dt+TreadRaster;
    timing.adcCenterTime = adcCenterTime;
    timing.edgeAmplitude = Gedge;
    timing.actualAmplitude_mTm = actualAmplitude_mTm;
    timing.slewPeak_T_m_s = max(abs(diff(fullBlockWave)))/dt/sys.gamma;
    waveSlewLimit = min(swave_max, physical_slew_max);
    assertArbitraryWaveformWithinSystem( ...
        fullBlockWave, Gedge, Gedge, sys, 'active cosine wave');
    assert(timing.slewPeak_T_m_s <= waveSlewLimit*(1+1e-9), ...
        'Active cosine slew exceeds the requested Wave slew limit.');
    timing.requestedFullAmplitudeSlew_T_m_s = ...
        2*pi*Ncycles/TreadRaster*gwave_max*1e-3;
    timing.area = area;
    timing.areaToAdcCenter = areaToAdcCenter;
    timing.areaAfterAdcCenter = area-areaToAdcCenter;

    if debugFlag
        fprintf(['z-cosine Wave: amplitude=%.6f mT/m, slew=%.6f T/m/s, ', ...
                 'area center/total=%.9g/%.9g 1/m, pre/post=%d/%d samples\n'], ...
            actualAmplitude_mTm, timing.slewPeak_T_m_s, areaToAdcCenter, ...
            area, nPadPre, nPadPost);
    end
end

function [gBridge, report] = makeFixedDurationThreePointBridge( ...
        channel, Gstart, Gend, areaTarget, duration, sys)
%MAKEFIXEDDURATIONTHREEPOINTBRIDGE Design a low-PNS endpoint bridge.
% The accepted bridge duration is preserved for both supported cycle cases;
% the midpoint location is chosen to minimize peak slew.
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
        slew1 = abs(Gmid-Gstart)/Tmid;
        slew2 = abs(Gend-Gmid)/(Ttotal-Tmid);
        peakSlew = max(slew1, slew2);
        if max(abs([Gstart Gmid Gend])) <= sys.maxGrad*(1+1e-12) && ...
                peakSlew <= sys.maxSlew*(1+1e-12) && peakSlew < bestSlew
            best = struct('Ttotal', Ttotal, 'Tmid', Tmid, ...
                'Gmid', Gmid, 'slew1', slew1, 'slew2', slew2, ...
                'peakSlew', peakSlew);
            bestSlew = peakSlew;
        end
    end
    if isempty(best)
        error(['No feasible separated three-point bridge was found at ', ...
               'the accepted %.3f ms duration.'], duration*1e3);
    end

    gBridge = mr.makeExtendedTrapezoid(channel, ...
        'times', [0 best.Tmid best.Ttotal], ...
        'amplitudes', [Gstart best.Gmid Gend], 'system', sys);
    areaError = gBridge.area-areaTarget;
    assert(abs(areaError) <= 1e-6, ...
        'Three-point bridge area error %.9g 1/m exceeds tolerance.', areaError);
    report = best;
    report.startAmplitude = Gstart;
    report.endAmplitude = Gend;
    report.areaTarget = areaTarget;
    report.areaActual = gBridge.area;
    report.areaError = areaError;
    report.peakAmplitude_mTm = max(abs([Gstart best.Gmid Gend]))/sys.gamma*1e3;
    report.peakSlew_T_m_s = best.peakSlew/sys.gamma;
end

function [G0_pulseq, w, TreadRaster, actualAmplitude_mTm] = designWaveAmplitude( ...
        Tread, sys, Ncycles, gwave_max, swave_max, physical_slew_max, waveInfoFlag)
%DESIGNWAVEAMPLITUDE Calculate a slew-limited sine-wave amplitude.
% Inputs specify duration, raster/system limits, cycles, and requested limits.
% Outputs are Pulseq amplitude [Hz/m], angular frequency, rasterized duration,
% and physical amplitude [mT/m]. Assertions reject invalid cycles or duration.
    assert(ismember(Ncycles, [10 20]), ...
        'Ncycles must be one of the supported protocol options: 10 or 20.');
    wavepoints = round(Tread/sys.gradRasterTime);
    TreadRaster = wavepoints*sys.gradRasterTime;
    assert(abs(TreadRaster-Tread) < sys.gradRasterTime/10, ...
        'Tread must be on the gradient raster.');
    w = 2*pi/(TreadRaster/Ncycles);

    swave_max = min(physical_slew_max, swave_max);
    activePhase = 2*pi*Ncycles*(0:(wavepoints-1))/(wavepoints-1);
    rasterSlewPerUnitAmplitude = max([ ...
        abs(diff(sin(activePhase))), abs(diff(cos(activePhase)))]) / ...
        sys.gradRasterTime;
    requestedAmplitude_T_m = gwave_max*1e-3;
    amplitude_T_m = min(requestedAmplitude_T_m, ...
        swave_max/rasterSlewPerUnitAmplitude);
    G0_pulseq = amplitude_T_m*sys.gamma;
    actualAmplitude_mTm = amplitude_T_m*1e3;

    if waveInfoFlag
        fprintf('Wave amplitude design: requested/actual %.6f/%.6f mT/m, raster slew %.6f T/m/s.\n', ...
            gwave_max, actualAmplitude_mTm, ...
            rasterSlewPerUnitAmplitude*amplitude_T_m);
    end
end

function assertArbitraryWaveformWithinSystem( ...
        waveform, firstValue, lastValue, system, segmentName)
    %ASSERTARBITRARYWAVEFORMWITHINSYSTEM Check an arbitrary gradient against limits.
    % waveform holds the raster samples [Hz/m] of a Pulseq arbitrary gradient
    % and firstValue/lastValue its block-boundary values; system supplies
    % maxGrad [Hz/m], maxSlew [Hz/m/s] and gradRasterTime, and segmentName
    % labels the error message. The slew includes the half-raster steps
    % between the boundary values and the first/last samples. No output; an
    % assertion error is raised when either limit is exceeded by more than a
    % relative 1e-10.
    waveform = waveform(:);
    slew = [2*(firstValue-waveform(1)); ...
        waveform(2:end)-waveform(1:end-1); ...
        2*(waveform(end)-lastValue)] / system.gradRasterTime;
    gradPeak = max(abs(waveform));
    slewPeak = max(abs(slew));
    tol = 1+1e-10;
    assert(gradPeak <= system.maxGrad*tol, ...
        '%s gradient %.9g exceeds assigned limit %.9g Hz/m.', ...
        segmentName, gradPeak, system.maxGrad);
    assert(slewPeak <= system.maxSlew*tol, ...
        '%s slew %.9g exceeds assigned limit %.9g Hz/m/s.', ...
        segmentName, slewPeak, system.maxSlew);
end
