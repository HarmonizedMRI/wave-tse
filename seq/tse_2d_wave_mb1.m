%% Create a matched 50-slice MB1 TSE RF-control sequence
%
% Starting point: the standard non-MB multi-slice Pulseq TSE sequence.
% This dedicated MB1 control keeps fifty independently acquired physical
% slices and provides one in-file RF design switch:
%   conventional_sinc - original 2.5 ms TSE excitation/refocusing pulses
%   sigpy_9ms         - 9 ms single-band templates saved by the MB2 generator
%
% Both modes share editable ETL/ESP/effective-TE controls and the same
% acquisition geometry as the MB2 protocol.
% The script intentionally contains no SMS groups, CAIPI, or collapsed-slice
% metadata. It must remain separate from the production MB2 sequence.
%
% Design choices:
%   - Keep the essential TSE gradient names and block structure:
%       GSex, GSref, GSspr, GSspex, GRacq, GRspr, GRspex, GRpreph,
%       GS1, GS2, GS3, GS4, GS5, GS7, GR3, GR5, GR6, GR7, GPpre, GPrew.
%   - Use GRE-style system setup with separate sys_sequence/sys_wave and
%     sys_lowPNS limits.
%   - Use sys_sequence for the sequence/ADC and sys_wave for Wave objects.
%   - Use sys_lowPNS for spoilers, crushers, PE pre/rewinders, and RF
%     slice-selective gradients.
%   - Add only gyWave. No cosine wave or wave ramp-up/down events.
%   - Represent gyWave as a zero-delay full-block arbitrary gradient with
%     zero guard samples; its active samples start and end with the ADC.
%   - Multi-slice excitation uses independent excitation/refocusing
%     slice-gradient amplitudes, interleaved slice order, and an optional
%     edge-to-edge sliceSpacing gap (default 0). Dummy trains reproduce
%     gyWave whenever wave encoding is enabled.
%   - Non-MB metadata distinguishes the regular slice sampling-grid FOV
%     (Nz times center-to-center spacing) from the smaller edge-to-edge
%     physical slice coverage when gaps are present.

clear all; close all; clc
format long

%% Path

pulseq_path='~/Code_mgh/pulseq';  % adjust path as needed
addpath([pulseq_path '/matlab']);
addpath('~/Code_mgh/sources/wave-haste/utils_seq_haste');
addpath('~/Code_mgh/sources/pulseq_basic/matlab_utils');

target_path = '~/Code_mgh/sources/pulseq_basic/matlab/';
cd(target_path);

out_path = '~/Code_mgh/scan_protocols/20261015_prisma/tse_product_control/';

% This MAT filename is shared with the existing MB2 generator. It already
% contains calibrated single-band excitation and refocusing templates.
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
mbPulseMatFile = fullfile(script_dir, 'rf_pulses', ...
    ['mb2_tse_20slices_10groups_sep60mm_sigpy_9ms_source_', ...
     'time_stretched_B1max15p000uT.mat']);

%% Parameters

% Write options:
%   write_v141_format = false  -> write only current format (v1.5.x)
%   write_v141_format = true   -> also write legacy v1.4.1 format
write_v141_format = false;
writeSequence = true;
execute_system_check_code = 1;
sys_type = 'skyra'; % prisma, skyra, Connectome2, C2_simulate_prisma, trio, prisma_XA30A, premier, CimaX

% RF control. Edit this value directly in MATLAB; do not pass it through the
% shell or environment.
%   'conventional_sinc' : original 2.5 ms TSE sinc pulses
%   'sigpy_9ms'         : exact 9 ms single-band SigPy templates
rfDesign = 'conventional_sinc';

% ------------------------- sequence timing -------------------------
alpha          = 90;                         % [deg]
ETL            = 17;                         % editable matched range: 8...17
TE1Requested   = 107e-3/7;                   % [s], ESP
TEeffRequested = 107e-3;                     % [s], snapped to nearest echo within ETL
minimumSupportedETL = 8;
maximumSupportedETL = 17;
assert(isscalar(ETL) && isfinite(ETL) && ETL == round(ETL) && ...
       ETL >= minimumSupportedETL && ETL <= maximumSupportedETL, ...
    'ETL must be an integer from %d through %d.', ...
    minimumSupportedETL, maximumSupportedETL);
assert(isscalar(TE1Requested) && isfinite(TE1Requested) && ...
       TE1Requested > 0, ...
    'TE1Requested must be one positive finite scalar in seconds.');
turboFactor    = ETL;
Nechoes        = 1;                          % ICE-visible echo dimension: TSE is single-echo/contrast
necho          = ETL;                        % keep standard TSE loop variable name
Nslices        = 50;
Ndummy         = 1;                          % original tse_pulseq.m uses one dummy train, kex=0
naverage       = 1;
rflip          = [165 150*ones(1, ETL-1)];   % first/repeated refocusing train [deg]
assert(numel(rflip) == ETL, ...
    'The refocusing flip-angle train must contain exactly ETL entries.');

sliceThickness = 3e-3;
% Additional edge-to-edge gap between adjacent slices [m]. A value of 0
% gives contiguous slices; the center-to-center spacing is
% sliceThickness + sliceSpacing.
sliceSpacing   = 0;
TE1            = TE1Requested;
TE1Baseline    = TE1;                        % legacy metadata alias
TRRequested    = 8000e-3;
TR             = TRRequested;                % updated to actual value if requested TR is infeasible
assert(isscalar(TEeffRequested) && isfinite(TEeffRequested) && ...
       TEeffRequested > 0, ...
    'TEeffRequested must be one positive finite scalar in seconds.');
effectiveEchoNumber = min(max(round(TEeffRequested/TE1), 1), ETL);
TEeff          = effectiveEchoNumber*TE1;
TEeffAdjustment = TEeff-TEeffRequested;
TEeffWasAdapted = abs(TEeffAdjustment) > eps(TEeffRequested)*8;
if TEeffWasAdapted
    warning(['Requested effective TE %.6f ms is not an available echo ', ...
        'center for ESP %.6f ms and ETL=%d. Using echo %d at %.6f ms.'], ...
        TEeffRequested*1e3, TE1*1e3, ETL, effectiveEchoNumber, TEeff*1e3);
end

% Conventional RF durations inherited from tse_pulseq.m. The custom mode
% reads its duration and dwell from mbPulseMatFile below.
tExConventional  = 2.5e-3;
tRefConventional = 3.1e-3;                  % stretched so refocusing stays below 15 uT
rfex_phase     = pi/2;
rfref_phase    = 0;

% Readout. roDuration is the ADC acquisition duration. readoutTime includes
% the ADC dead-time guard intervals before and after acquisition, matching
% the reference wave-TSE structure.
roDuration     = 6600e-6;                     % [s], ADC acquisition duration
os_factor      = 3;                          % readout oversampling factor for wave PSF-style acquisition
fspR           = 1.0;
fspS           = 0.5;
dG             = 250e-6;                     % readout/RF ramp time from tse_pulseq.m; low-PNS spoilers/crushers use auto ramps

% ------------------------- geometry -------------------------
slOrientation = 'TRA';
ax = struct;
ax.d1 = 'x';                                 % readout
ax.d2 = 'y';                                 % PE / LIN / sine wave
ax.d3 = 'z';                                 % slice select
ax.n1 = strfind('xyz', ax.d1);
ax.n2 = strfind('xyz', ax.d2);
ax.n3 = strfind('xyz', ax.d3);

fovScalar = 220e-3;
Nx = 440;
% Ny = 420;
% Nx = 210;
Ny = 220;
Nz = Nslices;

% ------------------------- acceleration / undersampling -------------------------
% GRE-style variable acceleration factor along the phase-encoding / LIN
% direction. For this 2D TSE sequence there is only one PE dimension, so we
% use Ry and no Rz. The acquired PE set is center-preserving; incomplete
% final trains use timing-only echo slots rather than dropping ky lines.
Ry = 1;                                      % acceleration along y/LIN
isUseACS = false;                            % if true and Ry>1, append central ACS ETL trains at the end
Nacs_y = 32;                                 % requested minimum number of central ACS ky lines

% ------------------------- wave -------------------------
gwave_max     = 10;                           % [mT/m]
swave_max     = 140;                         % [T/m/s]
Ncycles       = 20;                           % supported protocol options: 10 or 20
isUseWave_sin = false;
isUseWave_cos = false;
waveSinChannel = ax.d2;
waveCosChannel = '';

waveInfoFlag  = true;
waveDebugFlag = true;

%% System limits

slew_safety_magrin        = 0.7;
grad_safety_magrin        = 0.9;
lowPNS_slew_safety_margin = 0.41;
lowPNS_grad_safety_margin = grad_safety_magrin;
diff_slew_safety_margin   = 0.45;
diff_grad_safety_margin   = 0.97;

if strcmp(sys_type,'prisma') || strcmp(sys_type,'C2_simulate_prisma') || strcmp(sys_type,'prisma_XA30A')
    physical_slew_max = 200;
    physical_grad_max = 80;
    B0 = 2.89;
elseif strcmp(sys_type,'premier')
    physical_slew_max = 200;
    physical_grad_max = 70;
    B0 = 3;
elseif strcmp(sys_type,'Connectome2')
    physical_slew_max = 598.802;
    physical_grad_max = 500;
    B0 = 2.89;
elseif strcmp(sys_type,'skyra')
    physical_slew_max = 180;
    physical_grad_max = 43;
    B0 = 2.89;
elseif strcmp(sys_type,'trio')
    physical_slew_max = 170;
    physical_grad_max = 38;
    B0 = 2.89;
elseif strcmp(sys_type,'CimaX')
    physical_slew_max = 200;
    physical_grad_max = 200;
    B0 = 2.89;
elseif strcmp(sys_type,'TerraX')
    physical_slew_max = 250;
    physical_grad_max = 135;
    B0=2.89;
else
    error('Undefined sys_type: %s', sys_type)
end

isGEscanner = strcmp(sys_type,'premier');
if ~isGEscanner
    pislquant = 0; %#ok<NASGU>
end
if isGEscanner
    psd_rf_wait = 200e-6;
    rfDeadTime = 100e-6;
    rfRingdownTime = 60e-6 + psd_rf_wait;
    adcDeadTime = 20e-6;
    adcRasterTime = 2e-6;
    rfRasterTime = 2e-6;
    gradRasterTime = 4e-6;
    blockDurationRaster = 4e-6;
else
    rfDeadTime = 100e-6;
    rfRingdownTime = 100e-6;
    adcDeadTime = 20e-6;
    adcRasterTime = 100e-9;
    rfRasterTime = 1e-6;
    gradRasterTime = 10e-6;
    blockDurationRaster = 10e-6;
end

%% Resolve inputs and validate configuration

assert(ismember(Ncycles, [10 20]), ...
    'Ncycles must be one of the supported protocol options: 10 or 20.');
if Nslices < 1 || Nslices ~= round(Nslices)
    error('Nslices must be a positive integer.');
end
if sliceThickness <= 0
    error('sliceThickness must be positive.');
end
if sliceSpacing < 0
    error('sliceSpacing must be non-negative.');
end
if os_factor ~= round(os_factor) || os_factor < 1
    error('os_factor must be a positive integer.');
end
if Ry ~= round(Ry) || Ry < 1
    error('Ry must be a positive integer.');
end
if Nacs_y < 1 || Nacs_y ~= round(Nacs_y)
    error('Nacs_y must be a positive integer.');
end
if Ndummy ~= 1
    warning('This file keeps the original tse_pulseq.m dummy-loop structure, which uses one dummy train kex=0.');
end
if naverage ~= 1
    warning('This file keeps the original tse_pulseq.m averaging structure; naverage is stored only as metadata.');
end

sliceCenterSpacing = sliceThickness + sliceSpacing;

% Siemens Pulseq 2D multi-slice reconstruction uses the physical
% edge-to-edge stack extent for FOV(3):
%   max(slicePositions)-min(slicePositions)+sliceThickness.
% For equally spaced slices this is Nslices*sliceThickness plus the
% (Nslices-1) inter-slice gaps. Do not use Nslices*centerSpacing here;
% that adds one extra gap and can prevent the interpreter from recognizing
% the intended 2D slice geometry.
stackExtent = Nslices*sliceThickness + max(Nslices-1, 0)*sliceSpacing;
slicePositionsGeometric = ((1:Nslices) - 1 - (Nslices-1)/2) * sliceCenterSpacing;

% Siemens 2D reconstruction convention: SlicePositions is stored in the
% same order as acquisition, while SLC runs sequentially 0,1,...,Nslices-1.
% The interpreter then uses the positions to sort images geometrically.
sliceOrder = [1:2:Nslices, 2:2:Nslices];
assert(isequal(sort(sliceOrder), 1:Nslices), ...
    'Internal slice-order error: sliceOrder must contain each physical slice exactly once.');
slicePositions = slicePositionsGeometric(sliceOrder);

% Non-MB acquisition: raw ICE slices and physical slices are identical.
% Standard 2D multi-slice metadata use the physical stack extent, while
% Matrix(3) and SLC labels carry the independently acquired slice count.
fov = [fovScalar fovScalar stackExtent];
N = [Nx Ny Nz];
Nx_os = Nx * os_factor;

assert(Nz == Nslices && N(3) == Nslices, ...
    'Non-MB metadata requires Nz, Matrix(3), and Nslices to agree.');
assert(abs(fov(3) - (max(slicePositions)-min(slicePositions)+sliceThickness)) < 1e-12, ...
    'FOV(3) must equal max(slicePositions)-min(slicePositions)+sliceThickness.');
if Nslices > 1
    assert(max(abs(diff(slicePositionsGeometric) - sliceCenterSpacing)) < 1e-12, ...
        'Geometric slice positions must use the declared slice center spacing.');
end

fprintf(['Multi-slice geometry: Nslices=%d, thickness=%.6f mm, gap=%.6f mm, ', ...
         'center spacing=%.6f mm\n'], ...
    Nslices, sliceThickness*1e3, sliceSpacing*1e3, sliceCenterSpacing*1e3);
fprintf('  2D reconstruction FOV z [mm]       = %.6f\n', fov(3)*1e3);
fprintf('  Physical edge-to-edge coverage [mm]= %.6f\n', stackExtent*1e3);
fprintf('  Interleaved acquisition order (1-based geometric indices): %s\n', mat2str(sliceOrder));
fprintf('  Geometric slice centers [mm]: %s\n', mat2str(slicePositionsGeometric*1e3, 6));
fprintf('  Acquisition-ordered SlicePositions [mm]: %s\n', mat2str(slicePositions*1e3, 6));

%% Setup system objects

% The sequence-wide system is the physical scanner envelope. Only the
% active Wave samples use this slew allowance; every non-Wave gradient is
% still constructed with sys_lowPNS below.
sys_sequence = mr.opts('MaxGrad', physical_grad_max*grad_safety_magrin, 'GradUnit', 'mT/m', ...
    'MaxSlew', physical_slew_max, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime, ...
    'adcRasterTime', adcRasterTime, ...
    'rfRasterTime', rfRasterTime, ...
    'gradRasterTime', gradRasterTime, ...
    'blockDurationRaster', blockDurationRaster, ...
    'B0', B0);
sys_wave = sys_sequence;
sys = sys_sequence; % timing/raster compatibility alias; not for non-Wave gradient design

sys_lowPNS = mr.opts('MaxGrad', physical_grad_max*lowPNS_grad_safety_margin, 'GradUnit', 'mT/m', ...
    'MaxSlew', physical_slew_max*lowPNS_slew_safety_margin, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime, ...
    'adcRasterTime', adcRasterTime, ...
    'rfRasterTime', rfRasterTime, ...
    'gradRasterTime', gradRasterTime, ...
    'blockDurationRaster', blockDurationRaster, ...
    'B0', B0);

sys_diff = mr.opts('MaxGrad', physical_grad_max*diff_grad_safety_margin, 'GradUnit', 'mT/m', ...
    'MaxSlew', physical_slew_max*diff_slew_safety_margin, 'SlewUnit', 'T/m/s', ...
    'rfDeadTime', rfDeadTime, ...
    'rfRingdownTime', rfRingdownTime, ...
    'adcDeadTime', adcDeadTime, ...
    'adcRasterTime', adcRasterTime, ...
    'rfRasterTime', rfRasterTime, ...
    'gradRasterTime', gradRasterTime, ...
    'blockDurationRaster', blockDurationRaster, ...
    'B0', B0); %#ok<NASGU>

lims = sys_sequence; %#ok<NASGU>
seq = mr.Sequence(sys_sequence);

%% Setup RF, readout, PE, ADC, and TSE timing scalars

% Rasterize ADC duration in GRE style, while preserving the TSE variable
% names samplingTime and readoutTime.
dwell = round((roDuration / Nx_os) / sys.adcRasterTime) * sys.adcRasterTime;
Tread = dwell * Nx_os;
if abs(Tread - roDuration) > sys.adcRasterTime
    fprintf('Requested roDuration %.6f ms rasterized to Tread %.6f ms.\n', roDuration*1e3, Tread*1e3);
end
samplingTime = Tread;
readoutTime = samplingTime + 2*sys.adcDeadTime;

fprintf('Readout setup: Nx=%d, os_factor=%d, Nx_os=%d, dwell=%.3f us, Tread=%.6f ms, readoutTime=%.6f ms\n', ...
    Nx, os_factor, Nx_os, dwell*1e6, Tread*1e3, readoutTime*1e3);

%% Base gradients

%%% Slice selection
flipex = alpha*pi/180;
fliprefFirst = rflip(1)*pi/180;
fliprefRepeat = rflip(2)*pi/180;

switch rfDesign
    case 'conventional_sinc'
        tEx = tExConventional;
        tRef = tRefConventional;
        [rfex, gz] = mr.makeSincPulse(flipex, sys_lowPNS, ...
            'Duration', tEx, ...
            'SliceThickness', sliceThickness, ...
            'apodization', 0.5, ...
            'timeBwProduct', 4, ...
            'PhaseOffset', rfex_phase, ...
            'use', 'excitation');
        [rfrefFirst, gz2First] = mr.makeSincPulse(fliprefFirst, sys_lowPNS, ...
            'Duration', tRef, ...
            'SliceThickness', sliceThickness, ...
            'apodization', 0.5, ...
            'timeBwProduct', 4, ...
            'PhaseOffset', rfref_phase, ...
            'use', 'refocusing');
        [rfrefRepeat, gz2Repeat] = mr.makeSincPulse(fliprefRepeat, sys_lowPNS, ...
            'Duration', tRef, ...
            'SliceThickness', sliceThickness, ...
            'apodization', 0.5, ...
            'timeBwProduct', 4, ...
            'PhaseOffset', rfref_phase, ...
            'use', 'refocusing');
        assert(abs(gz2First.amplitude-gz2Repeat.amplitude) < ...
            max(abs(gz2First.amplitude), 1)*1e-12, ...
            'First/repeated conventional refocusing slice gradients differ.');
        rfDesignTag = sprintf( ...
            'conventionalSincEx2p5Ref3p1ms_165_150x%d', ETL-1);
        rfFilenameTag = 'convEx2p5R3p1';
        rfSourceDescription = sprintf([ ...
            'Pulseq makeSincPulse, excitation 2.5 ms / refocusing ', ...
            '3.1 ms, TBW 4, 90-165-150x%d, B1max15uT'], ETL-1);

    case 'sigpy_9ms'
        if ~isfile(mbPulseMatFile)
            error('SigPy RF template MAT file not found: %s', mbPulseMatFile);
        end
        MB = load(mbPulseMatFile);
        requiredFields = { ...
            'rf_ex_sigpy_rad_per_sample', ...
            'rf_ref_sigpy_rad_per_sample', ...
            'rf_ex_single_hz', ...
            'rf_ref_single_hz', ...
            'rf_duration_s', ...
            'rf_dwell_s', ...
            'time_bandwidth_product', ...
            'slice_thickness_m'};
        for iField = 1:numel(requiredFields)
            if ~isfield(MB, requiredFields{iField})
                error('SigPy RF template MAT file is missing field: %s', ...
                    requiredFields{iField});
            end
        end
        assert(abs(double(MB.slice_thickness_m) - sliceThickness) < 1e-12, ...
            'SigPy RF template slice thickness does not match this sequence.');

        tEx = double(MB.rf_duration_s);
        tRef = tEx;
        rfDwell = double(MB.rf_dwell_s);
        rfTBP = double(MB.time_bandwidth_product);
        assert(isscalar(tEx) && tEx >= 9e-3, ...
            'The custom MB1 control requires duration-stretched 9 ms-source templates.');
        assert(isscalar(rfDwell) && isfinite(rfDwell) && rfDwell > 0, ...
            'SigPy RF dwell must be one positive scalar.');
        assert(abs(numel(MB.rf_ex_single_hz) * rfDwell - tEx) < 1e-12, ...
            'SigPy excitation sample count, dwell, and duration disagree.');
        assert(abs(numel(MB.rf_ref_single_hz) * rfDwell - tRef) < 1e-12, ...
            'SigPy refocusing sample count, dwell, and duration disagree.');
        [rfex, gz] = mr.makeArbitraryRf( ... %#ok<ASGLU>
            MB.rf_ex_sigpy_rad_per_sample(:), flipex, ...
            'system', sys_lowPNS, ...
            'dwell', rfDwell, ...
            'timeBwProduct', rfTBP, ...
            'sliceThickness', sliceThickness, ...
            'PhaseOffset', rfex_phase, ...
            'use', 'excitation');
        [rfrefFirst, gz2First] = mr.makeArbitraryRf( ...
            MB.rf_ref_sigpy_rad_per_sample(:), fliprefFirst, ...
            'system', sys_lowPNS, ...
            'dwell', rfDwell, ...
            'timeBwProduct', rfTBP, ...
            'sliceThickness', sliceThickness, ...
            'PhaseOffset', rfref_phase, ...
            'use', 'refocusing');
        [rfrefRepeat, gz2Repeat] = mr.makeArbitraryRf( ...
            MB.rf_ref_sigpy_rad_per_sample(:), fliprefRepeat, ...
            'system', sys_lowPNS, ...
            'dwell', rfDwell, ...
            'timeBwProduct', rfTBP, ...
            'sliceThickness', sliceThickness, ...
            'PhaseOffset', rfref_phase, ...
            'use', 'refocusing');

        % Insert the already calibrated single-band signals without any
        % multiband modulation or additional normalization.
        rfex.signal = double(MB.rf_ex_single_hz(:));
        rfrefFirst.signal = double(MB.rf_ref_single_hz(:))*rflip(1)/180;
        rfrefRepeat.signal = double(MB.rf_ref_single_hz(:))*rflip(2)/180;
        assert(numel(rfex.signal) == numel(rfex.t), ...
            'SigPy excitation sample count does not match the RF event.');
        assert(numel(rfrefFirst.signal) == numel(rfrefFirst.t) && ...
               numel(rfrefRepeat.signal) == numel(rfrefRepeat.t), ...
            'SigPy refocusing sample count does not match the RF events.');
        if isfield(MB, 'slice_grad_hz_per_m')
            expectedSliceGrad = double(MB.slice_grad_hz_per_m);
            assert(abs(gz.amplitude - expectedSliceGrad) / expectedSliceGrad < 5e-3, ...
                'SigPy excitation slice gradient does not match MAT metadata.');
            assert(abs(gz2First.amplitude - expectedSliceGrad) / expectedSliceGrad < 5e-3 && ...
                   abs(gz2Repeat.amplitude - expectedSliceGrad) / expectedSliceGrad < 5e-3, ...
                'SigPy refocusing slice gradient does not match MAT metadata.');
        end
        rfDesignTag = 'sigpySingleBand9msSourceB1max15uT';
        rfFilenameTag = 'sigpy9srcB15';
        rfSourceDescription = 'Duration-stretched single-band templates from MB2 SigPy MAT';

    otherwise
        error(['Unknown rfDesign ''%s''. Use ''conventional_sinc'' or ', ...
               '''sigpy_9ms''.'], rfDesign);
end

% The slice gradients returned by makeSincPulse/makeArbitraryRf are used
% only for their amplitudes; their ramps are replaced below by GS1--GS5.
% Therefore, do not retain the RF delay that aligned the pulse to the
% discarded gradient rise. Start each RF event at its scanner dead time;
% the refocusing event is raster-centered more precisely below.
rfex.delay = rfex.deadTime;
rfrefFirst.delay = rfrefFirst.deadTime;
rfrefRepeat.delay = rfrefRepeat.deadTime;

peakRfB1_uT = max([max(abs(rfex.signal)), max(abs(rfrefFirst.signal)), ...
    max(abs(rfrefRepeat.signal))]) / 42.576;
assert(peakRfB1_uT <= 15*(1+1e-10), ...
    'MB1 RF peak %.6f uT exceeds the 15 uT protocol limit.', peakRfB1_uT);

fprintf('RF control: %s (%s), excitation/refocusing duration %.6f / %.6f ms\n', ...
    rfDesign, rfSourceDescription, tEx*1e3, tRef*1e3);

% Use the realized excitation duration for GS2. Jointly choose the
% refocusing RF delay and GS4 duration so the refocusing center lands at
% the flat-gradient block midpoint on both RF and gradient rasters. This
% prevents an alternating residual kz after successive refocusing pulses.
tExwd = ceil(mr.calcDuration(rfex) / sys.gradRasterTime) * sys.gradRasterTime;
[rfrefFirst, tRefwdFirst] = centerRefocusingRfInFlatBlock(rfrefFirst, sys);
[rfrefRepeat, tRefwdRepeat] = centerRefocusingRfInFlatBlock(rfrefRepeat, sys);
assert(abs(tRefwdFirst-tRefwdRepeat) < sys.gradRasterTime/10 && ...
       abs(rfrefFirst.delay-rfrefRepeat.delay) < sys.rfRasterTime/10, ...
    'First/repeated refocusing events do not share one centered block timing.');
tRefwd = tRefwdRepeat;
rfref = rfrefRepeat; % common timing template; sequence loops choose first/repeat RF

% Slice-select amplitudes are taken from the RF design. The actual GSex,
% GSref, GSspex, and GSspr objects are built below with one shared dGz so
% the split TSE slice-gradient chain has a consistent low-PNS ramp geometry.
GSexAmp = gz.amplitude;
GSrefAmp = gz2Repeat.amplitude;

rfExCenter = mr.calcRfCenter(rfex);
rfRefCenterFirst = mr.calcRfCenter(rfrefFirst);
rfRefCenterRepeat = mr.calcRfCenter(rfrefRepeat);
assert(abs(rfRefCenterFirst-rfRefCenterRepeat) < sys.rfRasterTime/10, ...
    'First/repeated refocusing RF temporal centers differ.');
rfRefCenter = rfRefCenterRepeat;
rfExCenterInBlock = rfex.delay + rfExCenter;
rfRefCenterInBlock = rfref.delay + rfRefCenter;
rfExCenterOffset = rfExCenterInBlock - tExwd/2;
rfRefCenterOffset = rfRefCenterInBlock - tRefwd/2;

% Compensate any remaining excitation/refocusing center offsets in emitted
% block coordinates. Symmetric conventional pulses produce zero here.
GSspexAreaCorrection = ...
    GSexAmp*rfExCenterOffset - GSrefAmp*rfRefCenterOffset;

fprintf('Slice-select RF-design amplitudes: GSex=%.9g Hz/m, GSref=%.9g Hz/m\n', ...
    GSexAmp, GSrefAmp);
fprintf('RF temporal centers: excitation/refocusing = %.6f / %.6f ms\n', ...
    rfExCenter*1e3, rfRefCenter*1e3);
fprintf('RF centers in emitted blocks              = %.6f / %.6f ms\n', ...
    rfExCenterInBlock*1e3, rfRefCenterInBlock*1e3);
fprintf('RF center offsets from block midpoint     = %+.6f / %+.6f us\n', ...
    rfExCenterOffset*1e6, rfRefCenterOffset*1e6);
fprintf('RF-center GSspex area correction           = %+.9g 1/m\n', ...
    GSspexAreaCorrection);

%%% Readout gradient
deltak = 1/fovScalar;
kWidth = Nx*deltak;

% Match the original tse_pulseq.m readout convention: GRacq is designed
% over the full constant-readout plateau, including ADC dead-time guards.
% The ADC still samples only samplingTime=Tread in the middle of this plateau.
% This keeps GRacq.amplitude, GR6, and the TSE readout prephaser moment
% consistent with the standard sequence.
GRacq = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'FlatArea', kWidth, 'FlatTime', readoutTime, 'riseTime', dG);
adc = mr.makeAdc(Nx_os, 'Duration', samplingTime, 'Delay', sys.adcDeadTime, 'system', sys_sequence);

%% Phase encoding / GRE-style Ry undersampling and optional ACS

% Build a center-preserving accelerated PE list and preserve every requested
% ky line exactly once. If its length is not divisible by ETL, timing-only
% slots in the final train omit the ADC instead of cropping outer k-space.
yCenterIdx = floor(Ny/2) + 1;
imageYIdxFull = makeAcceleratedIndexList(Ny, Ry);
imageYIdx = imageYIdxFull;
imagePESteps = imageYIdx - yCenterIdx;
nImagePE = numel(imageYIdx);
useSparseImageAdcSlots = mod(nImagePE, ETL) ~= 0;

if useSparseImageAdcSlots
    [PEindexOrderImg, LINorderImg, PEorder, phaseAreas, imageAdcMask, pe_steps, nex, k0prescr, k0curr] = makeTsePETablesWithSparseAdcSlots( ...
        imageYIdx, yCenterIdx, ETL, TEeff, TE1, deltak);
else
    [PEindexOrderImg, LINorderImg, PEorder, phaseAreas, pe_steps, nex, k0prescr, k0curr] = makeTsePETablesFromIndices( ...
        imageYIdx, yCenterIdx, ETL, TEeff, TE1, deltak);
    imageAdcMask = true(size(PEindexOrderImg));
end
echoTrainPositionImg = repmat((1:ETL).', 1, nex);
nImageEchoSlots = numel(imageAdcMask);
nImageNoAdcSlots = nImageEchoSlots-nnz(imageAdcMask);
assert(nnz(imageAdcMask) == nImagePE, ...
    'Image ADC-slot mask count does not match the acquired PE-line count.');
assert(nnz(imageAdcMask & (PEindexOrderImg == yCenterIdx)) == 1, ...
    'The image PE table must acquire the k-space center exactly once.');
assert(isequal(sort(PEindexOrderImg(imageAdcMask)).', sort(imageYIdx)), ...
    'The PE table must preserve every requested image ky line exactly once.');

% Optional ACS block appended after the image acquisition. Keep the same TSE
% ETL structure by expanding the requested central ACS width to the smallest
% unique central ky block whose size is divisible by ETL. No duplicate ACS
% ky lines are intentionally appended.
isAcquireACS = isUseACS && (Ry > 1);
if isUseACS && Ry == 1
    warning('isUseACS=true but Ry=1; ACS append is disabled to match GRE-style behavior.');
end
if isAcquireACS
    Nacs_y_actual = ceil(Nacs_y / ETL) * ETL;
    if Nacs_y_actual > Ny
        error(['Requested ACS width Nacs_y=%d rounds up to %d lines for ETL=%d, ', ...
               'which exceeds Ny=%d. Reduce Nacs_y or ETL.'], ...
            Nacs_y, Nacs_y_actual, ETL, Ny);
    end
    acsYIdx = makeCenterIndexList(Ny, Nacs_y_actual);
    acsYIdxBase = makeCenterIndexList(Ny, Nacs_y);
    acsPESteps = acsYIdx - yCenterIdx;
    nAcsPEBase = numel(acsYIdxBase);
    nAcsPE = numel(acsYIdx);
    nAcsExtraPE = nAcsPE - nAcsPEBase;
    nAcsDupPE = 0;
    [PEindexOrderACS, LINorderACS, PEorderACS, phaseAreasACS, pe_steps_acs, acsNex, k0prescrACS, k0currACS] = makeTsePETablesFromIndices( ...
        acsYIdx, yCenterIdx, ETL, TEeff, TE1, deltak);
    echoTrainPositionACS = repmat((1:ETL).', 1, acsNex);
else
    Nacs_y_actual = 0;
    acsYIdxBase = [];
    acsYIdx = [];
    acsPESteps = [];
    nAcsPEBase = 0;
    nAcsPE = 0;
    nAcsExtraPE = 0;
    nAcsDupPE = 0;
    acsNex = 0;
    PEindexOrderACS = [];
    LINorderACS = [];
    PEorderACS = [];
    pe_steps_acs = [];
    phaseAreasACS = [];
    echoTrainPositionACS = [];
    k0prescrACS = 0;
    k0currACS = 0;
end

nFullTSEPE = numel(makeAcceleratedIndexList(Ny, 1));
effectiveRy = nFullTSEPE / nImagePE;

fprintf('\nPE table / undersampling:\n');
fprintf('  Ny full matrix                     = %d\n', Ny);
fprintf('  Ry                                 = %d\n', Ry);
fprintf('  accelerated PE candidates          = %d\n', numel(imageYIdxFull));
fprintf('  acquired image PE lines            = %d\n', nImagePE);
fprintf('  image echo trains, nex             = %d\n', nex);
fprintf('  image echo slots                    = %d\n', nImageEchoSlots);
fprintf('  image no-ADC slots                  = %d\n', nImageNoAdcSlots);
fprintf('  image readouts                     = %d\n', nImagePE);
fprintf('  effective acceleration vs Ry=1 TSE = %.6f\n', effectiveRy);
if abs(effectiveRy - Ry) > 0.05
    warning(['Effective TSE acceleration %.6f differs from requested Ry=%d ', ...
             'because the centered stride grid contains %d unique lines.'], ...
             effectiveRy, Ry, nImagePE);
end
fprintf('  image PE step range                = %+d ... %+d\n', min(imagePESteps), max(imagePESteps));
fprintf('  image central echo target/native   = %d / %d\n', k0prescr, k0curr);
if isAcquireACS
    fprintf('  ACS requested central lines        = %d\n', Nacs_y);
    fprintf('  ACS actual unique central lines    = %d\n', nAcsPE);
    fprintf('  ACS extra unique lines for ETL     = %d\n', nAcsExtraPE);
    fprintf('  ACS duplicated fill lines          = %d\n', nAcsDupPE);
    fprintf('  ACS echo trains                    = %d\n', acsNex);
    fprintf('  ACS PE step range                  = %+d ... %+d\n', min(acsPESteps), max(acsPESteps));
    fprintf('  ACS central echo target/native     = %d / %d\n', k0prescrACS, k0currACS);
else
    fprintf('  ACS appended at end                = false\n');
end

%% Low-PNS spoiler/crusher/prephaser timing

% The original TSE example used a hard-coded dG=250 us for the split
% crusher/spoiler geometry. After switching these events to sys_lowPNS, the
% feasible durations have to be computed from the low-PNS limits first.
% tSpex is the block after excitation and before the first refocusing RF.
% tSp is the block on either side of each readout.
%
% Important readout convention:
%   GR5 is not the standalone GRspr trapezoid. It is a bridge waveform
%   0 -> GRspr.amplitude -> GRacq.amplitude. With low-PNS ramps, its
%   realized area differs from GRspr.area. Therefore the readout prephaser
%   GRpreph must be designed from the realized GR5 area plus the GR6 moment
%   up to the ADC center, not from GRspr.area alone.
AGRsprTarget = GRacq.area*fspR;
allPhaseAreasForTiming = phaseAreas(:);
if isAcquireACS
    allPhaseAreasForTiming = [allPhaseAreasForTiming; phaseAreasACS(:)];
end
maxPhaseArea = max(abs(allPhaseAreasForTiming));
raster = sys.gradRasterTime;

tSpTarget = 0.5*(TE1 - readoutTime - tRefwd);
tSpexTarget = 0.5*(TE1 - tExwd - tRefwd);
tSp = round(tSpTarget / raster) * raster;
tSpex = round(tSpexTarget / raster) * raster;

if tSpTarget <= 0
    error('TE1 is too short for the selected RF/readout timing: tSpTarget=%.6f ms.', tSpTarget*1e3);
end
if tSpexTarget <= 0
    error('TE1 is too short for the selected RF/readout timing: tSpexTarget=%.6f ms.', tSpexTarget*1e3);
end

% Readout spoiler and PE timing still use their shortest low-PNS durations.
% The slice-select / slice-crusher chain is handled separately below with a
% shared fixed ramp time dGz, because the split TSE z-gradient pieces must be
% geometrically consistent to avoid alternating residual kz moments.
GRsprNatural = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'area', AGRsprTarget);
GPpreMaxNatural = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
    'Area', maxPhaseArea);
tSpMinXPE = ceil(max([mr.calcDuration(GRsprNatural), mr.calcDuration(GPpreMaxNatural)]) / raster) * raster;

if tSp + raster/10 < tSpMinXPE
    TE1minRefToEcho_prelim = readoutTime + tRefwd + 2*tSpMinXPE;
    error(['TE1=%.6f ms is too short on the readout/PE side under sys_lowPNS. ', ...
           'Minimum feasible from this side is %.6f ms.'], ...
        TE1*1e3, TE1minRefToEcho_prelim*1e3);
end

% Build the selected low-PNS readout spoiler object that determines the
% realized bridge area into the readout plateau.
GRspr = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'area', AGRsprTarget, 'duration', tSp);

% Preview the actual GR5 bridge area. This bridge ends at GRacq.amplitude,
% so its moment includes the final transition into the readout plateau.
GR5timesPreview = [0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR5ampPreview = [0 GRspr.amplitude GRspr.amplitude GRacq.amplitude];
GR5areaForPreph = trapz(GR5timesPreview, GR5ampPreview);

% Moment from the start of GR6 to the ADC center. This uses the realized ADC
% delay/dwell and the constant GR6 amplitude.
GR6AreaToAdcCenter = GRacq.amplitude * (adc.delay + 0.5*adc.numSamples*adc.dwell);

% Correct readout prephaser area: after the refocusing pulse, the previous
% x moment is sign-inverted. To place kx=0 at the ADC center, the pre-180
% GRpreph moment must equal the realized post-180 moment up to ADC center.
AGRpreph = GR5areaForPreph + GR6AreaToAdcCenter;
GRprephNatural = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'Area', AGRpreph);
tSpexMinX = ceil(mr.calcDuration(GRprephNatural) / raster) * raster;

if tSpex + raster/10 < tSpexMinX
    TE1minExcToRef_prelim = tExwd + tRefwd + 2*tSpexMinX;
    error(['TE1=%.6f ms is too short on the readout-prephaser side under sys_lowPNS. ', ...
           'Minimum feasible from this side is %.6f ms.'], ...
        TE1*1e3, TE1minExcToRef_prelim*1e3);
end

% Option B for the slice-gradient chain: compute one common low-PNS-safe
% transition time dGz and use it for GSex, GSref, GSspex, and GSspr. This
% preserves the original TSE split-gradient geometry while adapting the ramp
% time to sys_lowPNS instead of allowing each standalone trapezoid to choose
% a different natural ramp time.
[GSex, GSref, GSspex, GSspr, dGz, sliceTiming] = designCommonTseSliceGradients( ...
    ax.d3, GSexAmp, GSrefAmp, tExwd, tRefwd, tSpex, tSp, fspS, ...
    GSspexAreaCorrection, dG, sys_lowPNS);
AGSex = GSex.area/2;

tSpMin = max([tSpMinXPE, sliceTiming.minTSp]);
tSpexMin = max([tSpexMinX, sliceTiming.minTSpex]);
TE1minExcToRef = tExwd + tRefwd + 2*tSpexMin;
TE1minRefToEcho = readoutTime + tRefwd + 2*tSpMin;
TE1min = max(TE1minExcToRef, TE1minRefToEcho);

fprintf('\nTSE low-PNS timing feasibility summary:\n');
fprintf('  tExwd / tRefwd RF-flat blocks       = %.6f / %.6f ms\n', tExwd*1e3, tRefwd*1e3);
fprintf('  requested TE1                       = %.6f ms\n', TE1*1e3);
fprintf('  tSpex target / low-PNS min / chosen = %.6f / %.6f / %.6f ms\n', ...
    tSpexTarget*1e3, tSpexMin*1e3, tSpex*1e3);
fprintf('  tSp target / low-PNS min / chosen   = %.6f / %.6f / %.6f ms\n', ...
    tSpTarget*1e3, tSpMin*1e3, tSp*1e3);
fprintf('  TE1 minimum from exc->ref side      = %.6f ms\n', TE1minExcToRef*1e3);
fprintf('  TE1 minimum from ref->echo side     = %.6f ms\n', TE1minRefToEcho*1e3);
fprintf('  TE1 minimum overall                 = %.6f ms\n', TE1min*1e3);

if tSp + raster/10 < tSpMin || tSpex + raster/10 < tSpexMin
    error(['TE1=%.6f ms is too short under the selected common-dGz slice timing. ', ...
           'Minimum feasible TE1 is %.6f ms (exc/ref side %.6f ms, ref/echo side %.6f ms).'], ...
        TE1*1e3, TE1min*1e3, TE1minExcToRef*1e3, TE1minRefToEcho*1e3);
end

% Now rebuild the remaining TSE readout prephaser object using the chosen
% common duration. These names intentionally match tse_pulseq.m.
GRspex = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'area', GRacq.area*(1+fspR), 'duration', tSpex); %#ok<NASGU>
GRpreph = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'Area', AGRpreph, 'duration', tSpex);
AGRspr = GRspr.area;

fprintf('\nTSE selected timing scalars:\n');
fprintf('  TE1=%.6f ms, TR=%.6f ms, tSp=%.6f ms, tSpex=%.6f ms\n', ...
    TE1*1e3, TR*1e3, tSp*1e3, tSpex*1e3);
fprintf('  common slice dGz                       = %.6f us\n', dGz*1e6);
fprintf('  slice dGz iterations / peak slew       = %d / %.6f T/m/s\n', ...
    sliceTiming.iterations, sliceTiming.peakSlew_T_m_s);
fprintf('  slice amplitudes GSex/GSspex/GSspr     = %.6f / %.6f / %.6f kHz/m\n', ...
    GSex.amplitude*1e-3, GSspex.amplitude*1e-3, GSspr.amplitude*1e-3);
fprintf('  GRspr/GPpreMax natural durations       = %.6f / %.6f ms\n', ...
    mr.calcDuration(GRsprNatural)*1e3, mr.calcDuration(GPpreMaxNatural)*1e3);
fprintf('  GRpreph natural duration               = %.6f ms\n', ...
    mr.calcDuration(GRprephNatural)*1e3);
fprintf('\nReadout moment centering diagnostics:\n');
fprintf('  kWidth                                = %.9g 1/m\n', kWidth);
fprintf('  GRacq amplitude                       = %.9g Hz/m\n', GRacq.amplitude);
fprintf('  GRacq.area used for spoiler target    = %.9g 1/m\n', GRacq.area);
fprintf('  GRspr standalone target/realized area = %.9g / %.9g 1/m\n', AGRsprTarget, GRspr.area);
fprintf('  GR5 realized bridge area preview      = %.9g 1/m\n', GR5areaForPreph);
fprintf('  GR6 area to ADC center                = %.9g 1/m\n', GR6AreaToAdcCenter);
fprintf('  GRpreph target/realized area          = %.9g / %.9g 1/m\n', AGRpreph, GRpreph.area);

%% Wave gradients: y sine

if isUseWave_sin
    [gyWave, waveAreaY, waveTiming] = defineSineReadoutWaveTSE( ...
        waveSinChannel, Tread, adc.delay, readoutTime, sys_wave, ...
        Ncycles, gwave_max, swave_max, physical_slew_max, ...
        adc, waveInfoFlag, waveDebugFlag);
    gyWave.id = seq.registerGradEvent(gyWave);
else
    gyWave = [];
    waveAreaY = 0;
    waveTiming = struct('Twave', 0, 'Tactive', 0, ...
        'peakAmplitude_mTm', 0, 'slewPeak_T_m_s', 0);
end

fprintf('\nWave carry-area summary after one TSE readout:\n');
fprintf('  y/sine area = %.9g 1/m\n', waveAreaY);
fprintf('  wave starts with ADC = %d\n', isUseWave_sin && ...
    abs(waveTiming.activeStartTime-adc.delay) < sys.gradRasterTime/10);
fprintf('  wave ends with ADC   = %d\n', isUseWave_sin && ...
    abs(waveTiming.activeEndTime-(adc.delay+adc.numSamples*adc.dwell)) < ...
    max(sys.gradRasterTime/10, sys.adcRasterTime));
if abs(waveAreaY) > 1e-6
    warning('Realized sine wave area is %.9g 1/m. PE rewinder is intentionally left unchanged in this first implementation.', waveAreaY);
end

%% Split gradients and recombine into blocks

% Slice selection and slice crushers.
GS1times = [0 GSex.riseTime];
GS1amp = [0 GSex.amplitude];
GS1 = mr.makeExtendedTrapezoid(ax.d3, 'times', GS1times, 'amplitudes', GS1amp, 'system', sys_lowPNS);

GS2times = [0 tExwd];
GS2amp = [GSex.amplitude GSex.amplitude];
GS2 = mr.makeExtendedTrapezoid(ax.d3, 'times', GS2times, 'amplitudes', GS2amp, 'system', sys_lowPNS);

GS3times = [0 GSspex.riseTime GSspex.riseTime+GSspex.flatTime GSspex.riseTime+GSspex.flatTime+GSspex.fallTime];
GS3amp = [GSex.amplitude GSspex.amplitude GSspex.amplitude GSref.amplitude];
GS3 = mr.makeExtendedTrapezoid(ax.d3, 'times', GS3times, 'amplitudes', GS3amp, 'system', sys_lowPNS);

GS4times = [0 tRefwd];
GS4amp = [GSref.amplitude GSref.amplitude];
GS4 = mr.makeExtendedTrapezoid(ax.d3, 'times', GS4times, 'amplitudes', GS4amp, 'system', sys_lowPNS);

GS5times = [0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS5amp = [GSref.amplitude GSspr.amplitude GSspr.amplitude 0];
GS5 = mr.makeExtendedTrapezoid(ax.d3, 'times', GS5times, 'amplitudes', GS5amp, 'system', sys_lowPNS);

GS7times = [0 GSspr.riseTime GSspr.riseTime+GSspr.flatTime GSspr.riseTime+GSspr.flatTime+GSspr.fallTime];
GS7amp = [0 GSspr.amplitude GSspr.amplitude GSref.amplitude];
GS7 = mr.makeExtendedTrapezoid(ax.d3, 'times', GS7times, 'amplitudes', GS7amp, 'system', sys_lowPNS);

% Readout gradient and readout spoilers/crushers. Keep TSE names.
GR3 = GRpreph;

GR5times = [0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR5amp = [0 GRspr.amplitude GRspr.amplitude GRacq.amplitude];
GR5 = mr.makeExtendedTrapezoid(ax.d1, 'times', GR5times, 'amplitudes', GR5amp, 'system', sys_lowPNS);

GR6times = [0 readoutTime];
GR6amp = [GRacq.amplitude GRacq.amplitude];
GR6 = mr.makeExtendedTrapezoid(ax.d1, 'times', GR6times, 'amplitudes', GR6amp, 'system', sys_lowPNS);

GR7times = [0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR7amp = [GRacq.amplitude GRspr.amplitude GRspr.amplitude 0];
GR7 = mr.makeExtendedTrapezoid(ax.d1, 'times', GR7times, 'amplitudes', GR7amp, 'system', sys_lowPNS);

fprintf('\nRF/block-boundary duration check for nonzero flat segments:\n');
fprintf('  GS2 duration / rfex duration   = %.6f / %.6f ms\n', mr.calcDuration(GS2)*1e3, mr.calcDuration(rfex)*1e3);
fprintf('  GS4 duration / rfref duration  = %.6f / %.6f ms\n', mr.calcDuration(GS4)*1e3, mr.calcDuration(rfref)*1e3);
fprintf('  GR6 duration / ADC start/end   = %.6f / %.6f / %.6f ms\n', ...
    mr.calcDuration(GR6)*1e3, adc.delay*1e3, (adc.delay + adc.numSamples*adc.dwell)*1e3);

%% Fill times

tex = mr.calcDuration(GS1) + mr.calcDuration(GS2) + mr.calcDuration(GS3);
tref = mr.calcDuration(GS4) + mr.calcDuration(GS5) + mr.calcDuration(GS7) + readoutTime;
tend = mr.calcDuration(GS4) + mr.calcDuration(GS5);

tETrain = tex + necho*tref + tend;
minimumTRWithoutFill = Nslices*tETrain;
requestedTRfill = TRRequested/Nslices-tETrain;
TRWasRaisedToMinimum = requestedTRfill < 0;
if TRWasRaisedToMinimum
    TRfill = 0;
else
    TRfill = sys.gradRasterTime*ceil(requestedTRfill/sys.gradRasterTime-1e-12);
end
TRActual = Nslices*(tETrain+TRfill);
TR = TRActual;
trComparisonTolerance = max(eps(TRRequested)*8, sys.gradRasterTime*1e-6);
TRWasRasterRounded = ~TRWasRaisedToMinimum && ...
    abs(TRActual-TRRequested) > trComparisonTolerance;
TRWasAdapted = TRWasRaisedToMinimum || TRWasRasterRounded;
if TRWasRaisedToMinimum
    warning(['Requested TR %.6f ms is infeasible for %d sequential MB1 ', ...
        'slices at ETL=%d. Using the exact no-fill lower bound %.6f ms.'], ...
        TRRequested*1e3, Nslices, ETL, TRActual*1e3);
end
if TRfill > 0
    delayTR = mr.makeDelay(TRfill);
else
    delayTR = [];
end

fprintf('\nTSE timing summary before sequence loop:\n');
fprintf('  tex=%.6f ms, tref=%.6f ms, tend=%.6f ms\n', tex*1e3, tref*1e3, tend*1e3);
fprintf('  tETrain=%.6f ms, TRfill=%.6f ms, TR_actual=%.6f ms\n', ...
    tETrain*1e3, TRfill*1e3, TRActual*1e3);
fprintf('  TR_requested/minimum_without_fill = %.6f / %.6f ms\n', ...
    TRRequested*1e3, minimumTRWithoutFill*1e3);
fprintf('  TE_requested/actual/echo = %.6f / %.6f ms / %d\n', ...
    TEeffRequested*1e3, TEeff*1e3, effectiveEchoNumber);
fprintf('  GR6/readout plateau duration=%.6f ms\n', mr.calcDuration(GR6)*1e3);
if isUseWave_sin
    fprintf('  gyWave object/full-block/active durations=%.6f / %.6f / %.6f ms\n', ...
        mr.calcDuration(gyWave)*1e3, waveTiming.Twave*1e3, waveTiming.Tactive*1e3);
end

%% ICE labels

% Global 0-based LIN labels. The PE index table, phase-area table, and LIN
% label table are built together so ICE labels always match the physical ky
% gradient used by that echo. ACS uses the same global LIN labels and is
% distinguished by REF=true.
lblLIN = cell(1, Ny);
for iy = 1:Ny
    lblLIN{iy} = mr.makeLabel('SET', 'LIN', iy - 1);
end

lblSLC = cell(1, Nslices);
for s = 1:Nslices
    lblSLC{s} = mr.makeLabel('SET', 'SLC', s - 1);
end

% TSE is a single-contrast/single-echo acquisition from ICE's point of
% view. The echo-train position is encoded by LIN ordering, not by ECO.
lblECO0 = mr.makeLabel('SET', 'ECO', 0);

% Siemens/Pulseq refscan flags. REF marks appended ACS readouts. IMA is kept
% false for both image and ACS, following gre_3d_wave_label_split_blocks_minprep.m.
lblSetRefScan   = mr.makeLabel('SET', 'REF', true);
lblResetRefScan = mr.makeLabel('SET', 'REF', false);
lblResetImaScan = mr.makeLabel('SET', 'IMA', false);
lblSetRefScan.id   = seq.registerLabelEvent(lblSetRefScan);
lblResetRefScan.id = seq.registerLabelEvent(lblResetRefScan);
lblResetImaScan.id = seq.registerLabelEvent(lblResetImaScan);

%% Define sequence blocks

fprintf('\nInterleaved multi-slice RF frequency table:\n');
fprintf('  acquisition position | SLC | geometric index | z [mm] | rfex offset [Hz] | rfref offset [Hz]\n');
for iSliceOrder = 1:Nslices
    geometricIndex = sliceOrder(iSliceOrder);
    slicePosition = slicePositions(iSliceOrder);
    fprintf('  %20d | %3d | %15d | %+8.3f | %+16.6f | %+17.6f\n', ...
        iSliceOrder, iSliceOrder-1, geometricIndex, slicePosition*1e3, ...
        GSex.amplitude*slicePosition, GSref.amplitude*slicePosition);
end
assert(numel(unique(round(GSex.amplitude*slicePositions, 9))) == Nslices, ...
    'Excitation RF frequency offsets are not unique across slices.');
assert(numel(unique(round(GSref.amplitude*slicePositions, 9))) == Nslices, ...
    'Refocusing RF frequency offsets are not unique across slices.');

for kex = 0:nex % kex=0 is one dummy train; kex=1:nex are Ry-undersampled image trains
    for iSliceOrder = 1:Nslices
        % Sequential SLC in acquisition order, matched to SlicePositions.
        slicePosition = slicePositions(iSliceOrder);
        rfex.freqOffset = GSex.amplitude * slicePosition;
        rfrefFirst.freqOffset = GSref.amplitude * slicePosition;
        rfrefRepeat.freqOffset = GSref.amplitude * slicePosition;
        rfex.phaseOffset = rfex_phase - 2*pi*rfex.freqOffset*mr.calcRfCenter(rfex);
        rfrefFirst.phaseOffset = rfref_phase - ...
            2*pi*rfrefFirst.freqOffset*mr.calcRfCenter(rfrefFirst);
        rfrefRepeat.phaseOffset = rfref_phase - ...
            2*pi*rfrefRepeat.freqOffset*mr.calcRfCenter(rfrefRepeat);

        seq.addBlock(GS1);
        seq.addBlock(GS2, rfex);
        seq.addBlock(GS3, GR3);

        for kech = 1:necho
            if kech == 1
                rfrefThis = rfrefFirst;
            else
                rfrefThis = rfrefRepeat;
            end
            isImageAdcSlot = kex > 0 && imageAdcMask(kech,kex);
            if isImageAdcSlot
                iy = PEindexOrderImg(kech,kex);
                linValue = LINorderImg(kech,kex); %#ok<NASGU>
                phaseArea = phaseAreas(kech,kex);
            else
                iy = yCenterIdx;
                linValue = iy - 1; %#ok<NASGU>
                phaseArea = 0;
            end
            GPpre = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                'Area', phaseArea, 'Duration', tSp);
            GPrew = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                'Area', -phaseArea, 'Duration', tSp);

            seq.addBlock(GS4, rfrefThis);
            seq.addBlock(GS5, GR5, GPpre);
            if isImageAdcSlot
                if isUseWave_sin
                    seq.addBlock(GR6, gyWave, adc, lblLIN{iy}, lblSLC{iSliceOrder}, lblECO0, lblResetRefScan, lblResetImaScan);
                else
                    seq.addBlock(GR6, adc, lblLIN{iy}, lblSLC{iSliceOrder}, lblECO0, lblResetRefScan, lblResetImaScan);
                end
            else
                % Dummy trains omit ADC/labels but reproduce the acquired
                % gradient history, including the wave gradient when enabled.
                if isUseWave_sin
                    seq.addBlock(GR6, gyWave);
                else
                    seq.addBlock(GR6);
                end
            end
            seq.addBlock(GS7, GR7, GPrew);
        end
        seq.addBlock(GS4);
        seq.addBlock(GS5);
        if TRfill > 0
            seq.addBlock(delayTR);
        end
    end
end

% Optional ACS/refscan-like block appended at the end of the image sequence.
% This keeps the same ETL train structure and uses the same readout/wave
% module as the image data. ACS readouts use global LIN labels and REF=true;
% the accompanying PEorder MAT file records the ACS order.
if isAcquireACS
    fprintf('\nAppending ACS block at end: %d trains, %d readouts.\n', acsNex, nAcsPE);
    for kexACS = 1:acsNex
        for iSliceOrder = 1:Nslices
            % Sequential SLC in acquisition order, matched to SlicePositions.
            slicePosition = slicePositions(iSliceOrder);
            rfex.freqOffset = GSex.amplitude * slicePosition;
            rfrefFirst.freqOffset = GSref.amplitude * slicePosition;
            rfrefRepeat.freqOffset = GSref.amplitude * slicePosition;
            rfex.phaseOffset = rfex_phase - 2*pi*rfex.freqOffset*mr.calcRfCenter(rfex);
            rfrefFirst.phaseOffset = rfref_phase - ...
                2*pi*rfrefFirst.freqOffset*mr.calcRfCenter(rfrefFirst);
            rfrefRepeat.phaseOffset = rfref_phase - ...
                2*pi*rfrefRepeat.freqOffset*mr.calcRfCenter(rfrefRepeat);

            seq.addBlock(GS1);
            seq.addBlock(GS2, rfex);
            seq.addBlock(GS3, GR3);

            for kech = 1:necho
                if kech == 1
                    rfrefThis = rfrefFirst;
                else
                    rfrefThis = rfrefRepeat;
                end
                iy = PEindexOrderACS(kech,kexACS);
                linValue = LINorderACS(kech,kexACS); %#ok<NASGU>
                phaseArea = phaseAreasACS(kech,kexACS);
                GPpre = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                    'Area', phaseArea, 'Duration', tSp);
                GPrew = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                    'Area', -phaseArea, 'Duration', tSp);

                seq.addBlock(GS4, rfrefThis);
                seq.addBlock(GS5, GR5, GPpre);
                if isUseWave_sin
                    seq.addBlock(GR6, gyWave, adc, lblLIN{iy}, lblSLC{iSliceOrder}, lblECO0, lblSetRefScan, lblResetImaScan);
                else
                    seq.addBlock(GR6, adc, lblLIN{iy}, lblSLC{iSliceOrder}, lblECO0, lblSetRefScan, lblResetImaScan);
                end
                seq.addBlock(GS7, GR7, GPrew);
            end
            seq.addBlock(GS4);
            seq.addBlock(GS5);
            if TRfill > 0
                seq.addBlock(delayTR);
            end
        end
    end
end

%% Timing check

[ok, error_report] = seq.checkTiming;
if ok
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% ICE label, ADC/readout, and k-space diagnostics

adc_lbl = seq.evalLabels('evolution','adc');

expectedImageReadouts = nImagePE * Nslices;
expectedAcsReadouts = acsNex * Nslices * necho;
expectedTotalReadouts = expectedImageReadouts + expectedAcsReadouts;

fprintf('\nICE label summary from ADC events:\n');
fprintf('  Expected image / ACS / total readouts = %d / %d / %d\n', ...
    expectedImageReadouts, expectedAcsReadouts, expectedTotalReadouts);
if isfield(adc_lbl, 'LIN')
    fprintf('  LIN range: %d ... %d\n', min(adc_lbl.LIN), max(adc_lbl.LIN));
end
if isfield(adc_lbl, 'SLC')
    fprintf('  SLC range: %d ... %d\n', min(adc_lbl.SLC), max(adc_lbl.SLC));
end
if isfield(adc_lbl, 'ECO')
    fprintf('  ECO range: %d ... %d\n', min(adc_lbl.ECO), max(adc_lbl.ECO));
end
if isfield(adc_lbl, 'REF')
    fprintf('  REF count: image/non-ref=%d, ref=%d\n', sum(adc_lbl.REF(:)==0), sum(adc_lbl.REF(:)~=0));
end
if isfield(adc_lbl, 'IMA')
    fprintf('  IMA count: false=%d, true=%d\n', sum(adc_lbl.IMA(:)==0), sum(adc_lbl.IMA(:)~=0));
end

assert(numel(adc_lbl.LIN) == expectedTotalReadouts, 'Unexpected total ADC label count.');
assert(all(adc_lbl.LIN >= 0) && max(adc_lbl.LIN) < Ny, 'LIN labels out of range.');
assert(all(adc_lbl.SLC >= 0) && max(adc_lbl.SLC) < Nslices, 'SLC labels out of range.');
assert(all(adc_lbl.ECO(:) == 0), 'TSE should be labeled as single-echo for ICE; ECO must stay 0.');
if isfield(adc_lbl, 'REF')
    assert(sum(adc_lbl.REF(:) ~= 0) == expectedAcsReadouts, 'REF label count does not match ACS readout count.');
    assert(all(adc_lbl.REF(1:expectedImageReadouts) == 0), 'Image readouts should have REF=false.');
    if expectedAcsReadouts > 0
        assert(all(adc_lbl.REF(expectedImageReadouts+1:end) ~= 0), 'ACS readouts should have REF=true.');
    end
end
if isfield(adc_lbl, 'IMA')
    assert(all(adc_lbl.IMA(:) == 0), 'IMA/PATRefAndIma flag should remain false for both image and ACS readouts.');
end

% Verify that every sequential acquisition-order SLC receives the same
% image and ACS LIN tables. Geometric sorting is defined by SlicePositions.
for s = 1:Nslices
    if isfield(adc_lbl, 'REF')
        imageMaskThisSlice = (adc_lbl.SLC(:) == s-1) & (adc_lbl.REF(:) == 0);
    else
        imageMaskThisSlice = (adc_lbl.SLC(:) == s-1);
    end
    assert(sum(imageMaskThisSlice) == nImagePE, ...
        'Slice %d has an unexpected number of image readouts.', s-1);
    actualImageLINThisSlice = adc_lbl.LIN(imageMaskThisSlice);
    expectedImageLINThisSlice = LINorderImg(imageAdcMask);
    assert(isequal(actualImageLINThisSlice(:), expectedImageLINThisSlice(:)), ...
        'Slice %d image LIN order does not match LINorderImg.', s-1);

    if isAcquireACS
        acsMaskThisSlice = (adc_lbl.SLC(:) == s-1) & (adc_lbl.REF(:) ~= 0);
        assert(sum(acsMaskThisSlice) == acsNex*necho, ...
            'Slice %d has an unexpected number of ACS readouts.', s-1);
        actualAcsLINThisSlice = adc_lbl.LIN(acsMaskThisSlice);
        assert(isequal(actualAcsLINThisSlice(:), LINorderACS(:)), ...
            'Slice %d ACS LIN order does not match LINorderACS.', s-1);
    end
end

warning('off', 'mr:restoreShape');
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP(); %#ok<ASGLU>
warning('on', 'mr:restoreShape');

nAdcSamples = numel(t_adc);
adcSamplesPerReadout = adc.numSamples;
if nAdcSamples > 0
    assert(mod(nAdcSamples, adcSamplesPerReadout) == 0, ...
        'ADC sample count (%d) is not divisible by ADC samples/readout (%d).', nAdcSamples, adcSamplesPerReadout);
    nReadouts = nAdcSamples / adcSamplesPerReadout;
    fprintf('\nADC/readout count summary:\n');
    fprintf('  ADC samples/readout = %d\n', adcSamplesPerReadout);
    fprintf('  acquired readouts   = %d\n', nReadouts);
    fprintf('  expected readouts   = %d\n', expectedTotalReadouts);
    assert(nReadouts == expectedTotalReadouts, 'Unexpected number of acquired readouts.');

    if size(ktraj_adc, 2) >= adcSamplesPerReadout
        firstReadoutKy = ktraj_adc(2, 1:adcSamplesPerReadout);
        fprintf('\nRepresentative first-readout ky diagnostics from calculateKspacePP:\n');
        fprintf('  ky start/center/end = %.9g / %.9g / %.9g 1/m\n', ...
            firstReadoutKy(1), firstReadoutKy(round(end/2)), firstReadoutKy(end));
        fprintf('  ky peak-to-peak during ADC = %.9g 1/m\n', max(firstReadoutKy)-min(firstReadoutKy));

        kxAdcByReadout = reshape(ktraj_adc(1,:), adcSamplesPerReadout, nReadouts);
        kyAdcByReadout = reshape(ktraj_adc(2,:), adcSamplesPerReadout, nReadouts);
        kzAdcByReadout = reshape(ktraj_adc(3,:), adcSamplesPerReadout, nReadouts);
        iCenter = round(adcSamplesPerReadout/2);
        kzCenterAll = kzAdcByReadout(iCenter,:);
        fprintf('\nSlice-axis kz diagnostics from calculateKspacePP:\n');
        fprintf('  kz ADC center range        = %+.9g ... %+.9g 1/m\n', min(kzCenterAll), max(kzCenterAll));
        kzCenterPeakToPeak = max(kzCenterAll) - min(kzCenterAll);
        fprintf('  kz ADC center peak-to-peak = %.9g 1/m\n', kzCenterPeakToPeak);
        kzEchoTolerance = 1e-2*deltak;
        assert(kzCenterPeakToPeak <= kzEchoTolerance, ...
            ['MB1 ADC kz alternates across echoes: peak-to-peak ', ...
             '%.9g 1/m exceeds %.9g 1/m.'], ...
            kzCenterPeakToPeak, kzEchoTolerance);
        nDiagReadouts = min([ETL, nReadouts, 8]);
        fprintf('\nFirst echo-train kx/ky/kz diagnostics from calculateKspacePP:\n');
        fprintf('  Echo | kx_start / kx_center / kx_end [1/m] | ky_start / ky_center / ky_end [1/m] | kz_start / kz_center / kz_end [1/m]\n');
        for iDiag = 1:nDiagReadouts
            fprintf('  %4d | %+.9g / %+.9g / %+.9g | %+.9g / %+.9g / %+.9g | %+.9g / %+.9g / %+.9g\n', ...
                iDiag, ...
                kxAdcByReadout(1,iDiag), kxAdcByReadout(iCenter,iDiag), kxAdcByReadout(end,iDiag), ...
                kyAdcByReadout(1,iDiag), kyAdcByReadout(iCenter,iDiag), kyAdcByReadout(end,iDiag), ...
                kzAdcByReadout(1,iDiag), kzAdcByReadout(iCenter,iDiag), kzAdcByReadout(end,iDiag));
        end
    end
end

%% Optional plots and report

plot_seq = false;
plot_ktraj = false;
print_test_report = false;

if plot_ktraj
    figure; plot(t_ktraj, ktraj'); title('k-space components as functions of time');
    figure; plot(ktraj(1,:), ktraj(2,:), 'b', ktraj_adc(1,:), ktraj_adc(2,:), 'r.');
    axis('equal');
    title('2D k-space');
end

if plot_seq
    seq.plot('showBlocks', true, 'stacked', true);
end

if print_test_report
    rep = seq.testReport;
    fprintf([rep{:}]);
end

%% Check PNS/CNS and forbidden frequencies

if execute_system_check_code
    ascDir = '~/Code_mgh/pulseq/matlab/idea/asc';
    sysMap = {
        'prisma',            'MP_GPA_K2309_2250V_951A_AS82.asc'
        'prisma_XA30A',      'MP_GPA_K2309_2250V_951A_AS82_XA30A_mod_twoFilesCombined.asc'
        'skyra',             'MP_GPA_K2309_2250V_793A_GC99.asc'
        'Connectome2',       'MP_GradSys_P034_c_CX600.asc'
        'C2_simulate_prisma','MP_GradSys_P034_c_CX600.asc'
        'CimaX',             'MP_GradSys_P034_X60.asc'
        'TerraX',            'MP_GradSys_K2298_2250V_1250A_W60_SC72CD.asc'
    };
    idx = strcmp(sys_type, sysMap(:,1));
    if any(idx)
        idx = find(idx, 1);
        ascName = fullfile(ascDir, sysMap{idx, 2});
    else
        if ~isGEscanner
            error('Unknown system : %s, cannot find .asc file\n', sys_type);
        end
    end

    warning('off', 'mr:restoreShape');
    isHasCNS = strcmp(sys_type, 'CimaX');
    doPlots = true;
    isCalcCNS = 0;
    [pns,tpns] = seq.calcPNS(ascName, doPlots, isCalcCNS);
    if ~isGEscanner && max(tpns) > 0.95
        warning('PNS=%.2f too high, the sequence may not run on the scanner', max(tpns));
    end
    if isHasCNS
        isCalcCNS = 1;
        [pns,tpns] = seq.calcPNS(ascName, doPlots, isCalcCNS);
        if ~isGEscanner && max(tpns) > 0.95
            warning('CNS=%.2f too high, the sequence may not run on the scanner', max(tpns));
        end
    end

    if ~isGEscanner
        tic;
        fprintf('Checking forbidden frequencies... ');
        forbiddenFreqCheck(seq, sys, ascName);
        toc;
    end
    warning('on', 'mr:restoreShape');
end

%% Sequence definitions and write

% Non-MB ICE-facing 2D multi-slice geometry. The Siemens Pulseq
% interpreter recognizes the exact standard keys SlicePositions,
% SliceThickness, and SliceGap. SliceGap is the edge-to-edge void between
% neighboring slices, not the center-to-center spacing.
seq.setDefinition('FOV', fov);
seq.setDefinition('Matrix', N);
seq.setDefinition('SliceThickness', sliceThickness);
seq.setDefinition('SliceGap', sliceSpacing);
seq.setDefinition('SliceSpacing', sliceSpacing); % custom backward-compatible alias
seq.setDefinition('SliceCenterSpacing', sliceCenterSpacing);
seq.setDefinition('SliceStackExtent', stackExtent);
seq.setDefinition('SliceCoverageExtent', stackExtent);
seq.setDefinition('SlicePositions', slicePositions);
seq.setDefinition('SliceOrder', sliceOrder);
seq.setDefinition('SliceInterleaved', double(Nslices > 1));

% Explicit aliases make clear that acquired, physical, and reconstructed
% slice counts coincide in the non-MB sequence.
seq.setDefinition('PhysicalMatrix', N);
seq.setDefinition('PhysicalSlicePositions', slicePositionsGeometric);
seq.setDefinition('PhysicalSliceStackExtent', stackExtent);
seq.setDefinition('NPhysicalSlices', Nslices);
seq.setDefinition('AcquiredSliceCount', Nslices);
seq.setDefinition('ReconstructedSlices', Nslices);
seq.setDefinition('SMSFactor', 1);
seq.setDefinition('NsliceGroups', Nslices);
seq.setDefinition('TE1', TE1);
seq.setDefinition('TE1Baseline', TE1Baseline);
seq.setDefinition('TE1Requested', TE1Requested);
seq.setDefinition('TE1MinimumFeasible', TE1min);
seq.setDefinition('TEeff', TEeff);
seq.setDefinition('TEeffRequested', TEeffRequested);
seq.setDefinition('TEeffActual', TEeff);
seq.setDefinition('TEeffEchoNumber', effectiveEchoNumber);
seq.setDefinition('TEeffAdjustment', TEeffAdjustment);
seq.setDefinition('TEeffWasAdapted', double(TEeffWasAdapted));
seq.setDefinition('TR', TR);
seq.setDefinition('TRActual', TRActual);
seq.setDefinition('TRRequested', TRRequested);
seq.setDefinition('TRMinimumWithoutFill', minimumTRWithoutFill);
seq.setDefinition('TRWasAdapted', double(TRWasAdapted));
seq.setDefinition('TRWasRaisedToMinimum', double(TRWasRaisedToMinimum));
seq.setDefinition('TRWasRasterRounded', double(TRWasRasterRounded));
seq.setDefinition('FlipAngle', alpha);
seq.setDefinition('RefocusingFlipAngle', rflip(1));
seq.setDefinition('RefocusingFlipAngleFirst', rflip(1));
seq.setDefinition('RefocusingFlipAngleRepeat', rflip(2));
seq.setDefinition('RefocusingFlipAngleTrain', rflip);
seq.setDefinition('RFDesign', rfDesign);
seq.setDefinition('RFSourceDescription', rfSourceDescription);
seq.setDefinition('RFExcitationDuration', tEx);
seq.setDefinition('RFRefocusingDuration', tRef);
seq.setDefinition('RFExcitationDelay', rfex.delay);
seq.setDefinition('RFRefocusingDelay', rfrefRepeat.delay);
seq.setDefinition('RFRefocusingFirstDelay', rfrefFirst.delay);
seq.setDefinition('RFRefocusingRepeatDelay', rfrefRepeat.delay);
seq.setDefinition('RFExcitationEventDuration', mr.calcDuration(rfex));
seq.setDefinition('RFRefocusingEventDuration', mr.calcDuration(rfrefRepeat));
seq.setDefinition('RFRefocusingFirstEventDuration', mr.calcDuration(rfrefFirst));
seq.setDefinition('RFRefocusingRepeatEventDuration', mr.calcDuration(rfrefRepeat));
seq.setDefinition('RFExcitationFlatGradientDuration', tExwd);
seq.setDefinition('RFRefocusingFlatGradientDuration', tRefwd);
seq.setDefinition('RFExcitationTemporalCenter', rfExCenter);
seq.setDefinition('RFRefocusingTemporalCenter', rfRefCenter);
seq.setDefinition('RFRefocusingFirstTemporalCenter', rfRefCenterFirst);
seq.setDefinition('RFRefocusingRepeatTemporalCenter', rfRefCenterRepeat);
seq.setDefinition('RFExcitationCenterInBlock', rfExCenterInBlock);
seq.setDefinition('RFRefocusingCenterInBlock', rfRefCenterInBlock);
seq.setDefinition('RFExcitationBlockCenterOffset', rfExCenterOffset);
seq.setDefinition('RFRefocusingBlockCenterOffset', rfRefCenterOffset);
seq.setDefinition('GSspexAreaCorrection', GSspexAreaCorrection);
if strcmp(rfDesign, 'sigpy_9ms')
    [~, rfTemplateName, rfTemplateExtension] = fileparts(mbPulseMatFile);
    seq.setDefinition('RFTemplateMatFilename', ...
        [rfTemplateName rfTemplateExtension]);
    seq.setDefinition('RFTemplateMatFile', mbPulseMatFile);
else
    seq.setDefinition('RFTemplateMatFilename', 'none');
    seq.setDefinition('RFTemplateMatFile', 'none');
end
seq.setDefinition('Nechoes', Nechoes);              % ICE-visible echo dimension
seq.setDefinition('EchoTrainLength', ETL);
seq.setDefinition('TurboFactor', turboFactor);
seq.setDefinition('ETL', ETL);
seq.setDefinition('MinimumSupportedETL', minimumSupportedETL);
seq.setDefinition('MaximumSupportedETL', maximumSupportedETL);
seq.setDefinition('necho', necho);
seq.setDefinition('Nslices', Nslices);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Ny', Ny);
seq.setDefinition('Nz', Nz);
seq.setDefinition('Nx_os', Nx_os);
seq.setDefinition('ReadoutOversamplingFactor', os_factor);
seq.setDefinition('roDuration', Tread);
seq.setDefinition('readoutTime', readoutTime);
seq.setDefinition('dwell', dwell);
seq.setDefinition('CommonSliceRampTime', dGz);
seq.setDefinition('CommonSliceRampTime_us', dGz*1e6);
seq.setDefinition('Ndummy', Ndummy);
seq.setDefinition('Averages', naverage);
seq.setDefinition('naverage', naverage);
seq.setDefinition('Ry', Ry);
seq.setDefinition('AccelerationFactor', Ry);
seq.setDefinition('AccelerationFactorPE', Ry);
seq.setDefinition('EffectiveRy', effectiveRy);
seq.setDefinition('Ny_meas', nImagePE);
seq.setDefinition('Nex', nex);
seq.setDefinition('ImageEchoSlots', nImageEchoSlots);
seq.setDefinition('ImageNoADCSlots', nImageNoAdcSlots);
seq.setDefinition('ACSAtEndAsRefscan', double(isAcquireACS));
seq.setDefinition('Nacs_y', Nacs_y);
seq.setDefinition('Nacs_y_base', nAcsPEBase);
seq.setDefinition('Nacs_y_actual', nAcsPE);
seq.setDefinition('Nacs_y_meas', nAcsPE);
seq.setDefinition('Nacs_y_extra_unique', nAcsExtraPE);
seq.setDefinition('Nacs_y_duplicated', nAcsDupPE);
seq.setDefinition('AcsNex', acsNex);
seq.setDefinition('OrientationMapping', slOrientation);
centerLineIdx = yCenterIdx;
phaseResolution = (fov(1)/Nx) / (fov(2)/Ny);
seq.setDefinition('kSpaceCenterLine', centerLineIdx - 1);
seq.setDefinition('PhaseResolution', phaseResolution);
seq.setDefinition('SliceResolution', 1);
seq.setDefinition('Slices', Nslices);
seq.setDefinition('kSpaceCenterSlice', 0);
seq.setDefinition('ReadoutPolarity', 'TSE-refocused-readout');
seq.setDefinition('WaveSinChannel', waveSinChannel);
seq.setDefinition('WaveCosChannel', waveCosChannel);
seq.setDefinition('WaveAmplitude_mTm', waveTiming.peakAmplitude_mTm);
seq.setDefinition('WaveRequestedAmplitude_mTm', gwave_max);
seq.setDefinition('WaveSlew_Tms', swave_max);
seq.setDefinition('WaveActualSlew_Tms', waveTiming.slewPeak_T_m_s);
seq.setDefinition('WaveFullAmplitudeRequiredSlew_Tms', ...
    2*pi*Ncycles/Tread*gwave_max*1e-3);
seq.setDefinition('WaveCycles', Ncycles);
seq.setDefinition('WaveSupportedCycles', [10 20]);
seq.setDefinition('WaveCycleParity', (-1)^Ncycles);
seq.setDefinition('WaveCarryY', waveAreaY);
seq.setDefinition('WaveCarryZ', 0);
seq.setDefinition('WaveActiveDur', waveTiming.Tactive);
seq.setDefinition('WaveBlockEnvelopeDur', readoutTime);
seq.setDefinition('SplitReadoutWaveBlocks', 0);
seq.setDefinition('WaveStartsWithADC', double(isUseWave_sin));
seq.setDefinition('WaveEndsWithADC', double(isUseWave_sin));
seq.setDefinition('SystemType', sys_type);
seq.setDefinition('SequenceMaxGrad_mTm', sys_sequence.maxGrad/sys_sequence.gamma*1e3);
seq.setDefinition('SequenceMaxSlew_Tms', sys_sequence.maxSlew/sys_sequence.gamma);
seq.setDefinition('WaveMaxGrad_mTm', sys_wave.maxGrad/sys_wave.gamma*1e3);
seq.setDefinition('WaveMaxSlew_Tms', sys_wave.maxSlew/sys_wave.gamma);
seq.setDefinition('LowPNSMaxGrad_mTm', sys_lowPNS.maxGrad/sys_lowPNS.gamma*1e3);
seq.setDefinition('LowPNSMaxSlew_Tms', sys_lowPNS.maxSlew/sys_lowPNS.gamma);
seq.setDefinition('ExecuteSystemCheckCode', execute_system_check_code);
seq.setDefinition('WriteV141Format', double(write_v141_format));
seq.setDefinition('MetadataRevision', 'nonMB_Siemens2D_sequentialSLC_acqOrderedPositions_v3');
seq.setDefinition('OutputPath', out_path);
seq.setDefinition('Name', 'tse_waveYsin');

if isUseWave_sin
    tag_wave = '_wY';
else
    tag_wave = '';
end
tag_cycles = ['_Ncyc' num2str(Ncycles)];
if isAcquireACS
    tag_ACS = ['_ACS' num2str(Nacs_y) 'to' num2str(nAcsPE)];
else
    tag_ACS = '';
end
teTag = strrep(sprintf('%.2f', TE1*1e3), '.', 'p');
seqFilename = ['tse_mb1_' rfFilenameTag tag_wave tag_cycles '_' slOrientation '_' ...
    num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nz) ...
    '_FOV' num2str(round(fovScalar*1e3)) ...
    '_Ry' num2str(Ry) tag_ACS ...
    '_os' num2str(os_factor) ...
    '_ETL' num2str(ETL) ...
    '_TE' teTag ...
    '_TR' num2str(round(TR*1e3)) ...
    '_' sys_type];
assert(numel([seqFilename '.seq']) <= 100, ...
    'Sequence filename exceeds 100 characters: %s.seq', seqFilename);
assert(numel([seqFilename '_PEorder.mat']) <= 100, ...
    'PE-order filename exceeds 100 characters: %s_PEorder.mat', seqFilename);
if write_v141_format
    assert(numel([seqFilename '_v141.seq']) <= 100, ...
        'Pulseq v1.4.1 filename exceeds 100 characters: %s_v141.seq', seqFilename);
end
seq.setDefinition('Name', seqFilename);

if writeSequence
    if ~exist([out_path 'v141/'], 'dir'), mkdir([out_path 'v141/']); end
    if ~exist([out_path 'v151/'], 'dir'), mkdir([out_path 'v151/']); end

    if write_v141_format
        seq.write_v141([out_path 'v141/' seqFilename '_v141.seq']);
        seq.write([out_path 'v151/' seqFilename '.seq']);
    else
        seq.write([out_path 'v151/' seqFilename '.seq']);
    end

    fprintf('\nWrote sequence name: %s\n', seqFilename);

    % Save PE ordering alongside the sequence, matching the wave-TSE reference
    % workflow where PEorder is needed during reconstruction.
    save([out_path seqFilename '_PEorder.mat'], ...
        'PEorder', 'PEindexOrderImg', 'LINorderImg', 'phaseAreas', ...
        'imageYIdx', 'imagePESteps', 'imageAdcMask', 'nImageEchoSlots', ...
        'nImageNoAdcSlots', 'echoTrainPositionImg', 'Ry', 'Nechoes', 'ETL', 'turboFactor', 'nex', ...
        'isAcquireACS', 'PEorderACS', 'PEindexOrderACS', 'LINorderACS', 'phaseAreasACS', 'echoTrainPositionACS', ...
        'acsYIdx', 'acsYIdxBase', 'acsPESteps', ...
        'Nacs_y', 'Nacs_y_actual', 'nAcsPEBase', 'nAcsPE', 'nAcsExtraPE', 'nAcsDupPE', 'acsNex', ...
        'fov', 'N', 'Nz', ...
        'sliceThickness', 'sliceSpacing', 'sliceCenterSpacing', 'stackExtent', ...
        'slicePositionsGeometric', 'slicePositions', 'sliceOrder', ...
        'rfDesign', 'rfDesignTag', 'rfFilenameTag', 'rfSourceDescription', ...
        'mbPulseMatFile', 'tEx', 'tRef', 'rflip', ...
        'minimumSupportedETL', 'maximumSupportedETL', ...
        'TE1Requested', 'TE1', 'TE1Baseline', 'TE1min', ...
        'TE1minExcToRef', 'TE1minRefToEcho', ...
        'TEeffRequested', 'TEeff', 'effectiveEchoNumber', ...
        'TEeffAdjustment', 'TEeffWasAdapted', ...
        'TRRequested', 'TRActual', 'minimumTRWithoutFill', ...
        'TRWasAdapted', 'TRWasRaisedToMinimum', 'TRWasRasterRounded', ...
        'Ncycles', 'gwave_max', 'swave_max', 'isUseWave_sin', ...
        'waveAreaY', 'waveTiming', 'sys_sequence', 'sys_wave', 'sys_lowPNS');
    fprintf('Saved PE order: %s\n', [out_path seqFilename '_PEorder.mat']);
else
    fprintf('\nValidation-only run: writeSequence=false; no sequence or PE metadata written.\n');
end

%% Local helper functions

function idx = makeAcceleratedIndexList(N, R)
    if R < 1 || R ~= round(R)
        error('Acceleration factor must be a positive integer.');
    end
    centerIdx = floor(N/2) + 1;
    idx = [];
    for ii = 1:N
        if mod(ii - centerIdx, R) == 0
            idx(end+1) = ii; %#ok<AGROW>
        end
    end
end

function idxCrop = centerCropIndexListToMultiple(idx, centerIdx, blockSize)
    if isempty(idx)
        error('Cannot crop an empty PE index list.');
    end
    if blockSize < 1 || blockSize ~= round(blockSize)
        error('blockSize must be a positive integer.');
    end
    nUse = floor(numel(idx) / blockSize) * blockSize;
    if nUse < blockSize
        error(['Only %d accelerated PE candidates are available, which is fewer than ', ...
               'one complete echo train of %d readouts. Reduce Ry or ETL.'], ...
               numel(idx), blockSize);
    end

    [~, centerPos] = min(abs(idx - centerIdx));
    startPos = centerPos - floor(nUse/2);
    startPos = max(1, min(startPos, numel(idx) - nUse + 1));
    idxCrop = idx(startPos:(startPos + nUse - 1));
end

function idx = makeCenterIndexList(N, Ncenter)
    if Ncenter < 1 || Ncenter ~= round(Ncenter)
        error('Center block size must be a positive integer.');
    end
    Ncenter = min(Ncenter, N);
    startIdx = floor(N/2) + 1 - floor(Ncenter/2);
    startIdx = max(1, min(startIdx, N - Ncenter + 1));
    idx = startIdx:(startIdx + Ncenter - 1);
end

function [PEindexOrder, LINorder, PEorder, phaseAreas, adcSlotMask, pe_steps, nex, k0prescr, k0curr] = makeTsePETablesWithSparseAdcSlots(yIdx, yCenterIdx, ETL, TEeff, TE1, deltak)
    % Preserve every requested ky index and use center-index placeholders
    % for timing-only echo slots that must not emit an ADC.
    nPE = numel(yIdx);
    if nPE < ETL
        error('PE index count (%d) must contain at least one echo train of length %d.', nPE, ETL);
    end
    nex = ceil(nPE/ETL);
    nSlots = nex*ETL;
    pe_steps = yIdx(:).'-yCenterIdx;
    [~, iPEmin] = min(abs(pe_steps));
    k0curr = floor((iPEmin-1)/nex)+1;
    k0prescr = min(max(round(TEeff/TE1), 1), ETL);

    paddedAdcMask = false(1, nSlots);
    paddedAdcMask(1:nPE) = true;
    paddedSteps = zeros(1, nSlots);
    paddedSteps(1:nPE) = pe_steps;
    peStepOrder = reshape(paddedSteps, [nex ETL]).';
    adcSlotMask = reshape(paddedAdcMask, [nex ETL]).';
    echoShift = k0prescr-k0curr;
    peStepOrder = circshift(peStepOrder, [echoShift 0]);
    adcSlotMask = circshift(adcSlotMask, [echoShift 0]);

    PEindexOrder = peStepOrder+yCenterIdx;
    LINorder = PEindexOrder-1;
    PEorder = peStepOrder;
    phaseAreas = PEorder*deltak;
    acquiredIndices = PEindexOrder(adcSlotMask);
    if numel(acquiredIndices) ~= nPE || ...
            ~isequal(sort(acquiredIndices(:)), sort(yIdx(:)))
        error('Sparse-slot PE table does not preserve every requested PE index exactly once.');
    end
    if any(PEorder(~adcSlotMask) ~= 0) || ...
            any(PEindexOrder(~adcSlotMask) ~= yCenterIdx)
        error('Sparse-slot placeholders must use zero PE area at the center index.');
    end
end


function [PEindexOrder, LINorder, PEorder, phaseAreas, pe_steps, nex, k0prescr, k0curr] = makeTsePETablesFromIndices(yIdx, yCenterIdx, ETL, TEeff, TE1, deltak)
    % Build the physical PE order and matching ICE LIN order together.
    % PEindexOrder is 1-based MATLAB/global ky index. LINorder is the
    % corresponding 0-based ICE LIN label. PEorder is centered ky step in
    % units of delta-k, matching the physical phase-area table.
    nPE = numel(yIdx);
    if nPE < ETL || mod(nPE, ETL) ~= 0
        error('PE index count (%d) must be a positive integer multiple of ETL (%d).', nPE, ETL);
    end
    nex = nPE / ETL;

    pe_steps = yIdx(:).' - yCenterIdx;
    if 0 == mod(ETL,2)
        pe_steps = circshift(pe_steps, [0, -round(nex/2)]);
    end
    [~, iPEmin] = min(abs(pe_steps));
    k0curr = floor((iPEmin-1)/nex) + 1;
    k0prescr = max(round(TEeff/TE1), 1);
    k0prescr = min(k0prescr, ETL);

    peStepOrder = circshift(reshape(pe_steps, [nex, ETL])', k0prescr-k0curr);
    PEindexOrder = peStepOrder + yCenterIdx;
    LINorder = PEindexOrder - 1;
    PEorder = peStepOrder;
    phaseAreas = PEorder * deltak;

    if any(~ismember(PEindexOrder(:), yIdx(:)))
        error('Internal PE table error: ordered PE indices are not contained in the requested PE index list.');
    end
end


function [rf, blockDuration] = centerRefocusingRfInFlatBlock(rf, sys)
    %CENTERREFOCUSINGRFINFLATBLOCK Align the RF center to the GS4 midpoint.
    % Search one gradient-raster period of RF-raster delay candidates. The
    % lower bound also guarantees that a midpoint-centered block is not
    % shorter than the RF event including ringdown.
    rfCenter = mr.calcRfCenter(rf);
    minCenteringDelay = max([rf.deadTime, ...
        rf.shape_dur + rf.ringdownTime - 2*rfCenter, 0]);
    firstDelay = ceil((minCenteringDelay - 1e-12) / ...
        sys.rfRasterTime) * sys.rfRasterTime;
    nCandidates = ceil(sys.gradRasterTime / sys.rfRasterTime) + 1;

    bestError = inf;
    bestDuration = inf;
    bestDelay = inf;
    for iCandidate = 0:nCandidates
        candidateDelay = firstDelay + iCandidate*sys.rfRasterTime;
        rfCandidate = rf;
        rfCandidate.delay = candidateDelay;
        minBlockDuration = ceil((mr.calcDuration(rfCandidate) - 1e-12) / ...
            sys.gradRasterTime) * sys.gradRasterTime;
        candidateDuration = round(2*(candidateDelay + rfCenter) / ...
            sys.gradRasterTime) * sys.gradRasterTime;
        candidateDuration = max(candidateDuration, minBlockDuration);
        candidateError = abs(candidateDelay + rfCenter - ...
            candidateDuration/2);

        isBetter = candidateError < bestError - 1e-12;
        isTieAndShorter = abs(candidateError - bestError) <= 1e-12 && ...
            candidateDuration < bestDuration - 1e-12;
        if isBetter || isTieAndShorter
            bestError = candidateError;
            bestDuration = candidateDuration;
            bestDelay = candidateDelay;
        end
    end

    rf.delay = bestDelay;
    blockDuration = bestDuration;
    assert(bestError <= sys.rfRasterTime/2 + 1e-12, ...
        ['Could not align the refocusing RF center to the flat-gradient ', ...
         'block midpoint: residual %.6f us.'], bestError*1e6);
end


function [GSex, GSref, GSspex, GSspr, dGz, report] = designCommonTseSliceGradients(channel, GSexAmp, GSrefAmp, tExwd, tRefwd, tSpex, tSp, fspS, GSspexAreaCorrection, dGmin, sys)
    %DESIGNCOMMONTSESLICEGRADIENTS Build split TSE slice gradients with one ramp time.
    % Inputs describe the gradient channel and amplitudes, RF flat times,
    % preparation durations, spoiler factor, first-refocusing bridge-area
    % correction, minimum ramp time, and scanner limits. Outputs are the
    % four slice-gradient events, shared ramp time, and validation report.
    % The original Pulseq TSE example assumes that the slice-select and
    % slice-crusher bridge waveforms share a common transition time. When
    % sys_lowPNS is used, independent natural trapezoid ramps can break this
    % geometry and leave alternating residual kz moments. This helper finds
    % one shared low-PNS-safe dGz and uses it for GSex, GSref, GSspex, GSspr.

    dt = sys.gradRasterTime;
    dGz = max(dt, ceil(dGmin / dt) * dt);
    dGz = max(dGz, ceil(abs(GSexAmp) / sys.maxSlew / dt) * dt);
    dGz = max(dGz, ceil(abs(GSrefAmp) / sys.maxSlew / dt) * dt);

    converged = false;
    maxIter = 30;
    for iter = 1:maxIter
        if 2*dGz > tSpex + dt/10
            error('Common slice ramp dGz=%.6f ms is too long for tSpex=%.6f ms.', dGz*1e3, tSpex*1e3);
        end
        if 2*dGz > tSp + dt/10
            error('Common slice ramp dGz=%.6f ms is too long for tSp=%.6f ms.', dGz*1e3, tSp*1e3);
        end

        GSex = mr.makeTrapezoid(channel, 'system', sys, ...
            'amplitude', GSexAmp, 'FlatTime', tExwd, 'riseTime', dGz);
        GSref = mr.makeTrapezoid(channel, 'system', sys, ...
            'amplitude', GSrefAmp, 'FlatTime', tRefwd, 'riseTime', dGz);
        AGSex = GSex.area/2;
        GSspexTargetArea = AGSex*fspS + GSspexAreaCorrection;
        GSspex = mr.makeTrapezoid(channel, 'system', sys, ...
            'area', GSspexTargetArea, 'duration', tSpex, 'riseTime', dGz);
        GSspr = mr.makeTrapezoid(channel, 'system', sys, ...
            'area', AGSex*(1+fspS), 'duration', tSp, 'riseTime', dGz);

        jumps = [ ...
            abs(GSex.amplitude - 0), ...
            abs(GSspex.amplitude - GSex.amplitude), ...
            abs(GSref.amplitude - GSspex.amplitude), ...
            abs(GSspr.amplitude - GSref.amplitude), ...
            abs(0 - GSspr.amplitude), ...
            abs(GSspr.amplitude - 0), ...
            abs(GSref.amplitude - GSspr.amplitude)];
        maxJump = max(jumps);
        dGzNeeded = max(dt, ceil(maxJump / sys.maxSlew / dt) * dt);
        dGzNeeded = max(dGzNeeded, ceil(dGmin / dt) * dt);

        if dGzNeeded <= dGz + dt/10
            converged = true;
            break;
        end
        dGz = dGzNeeded;
    end

    if ~converged
        error('Could not converge a common low-PNS-safe slice ramp dGz after %d iterations.', maxIter);
    end

    amps = [0, GSex.amplitude, GSspex.amplitude, GSref.amplitude, GSspr.amplitude, 0];
    peakGrad = max(abs(amps));
    peakSlew = maxJump / dGz;
    if peakGrad > sys.maxGrad * (1 + 1e-9)
        error('Common slice-gradient design exceeds maxGrad: %.6g > %.6g Hz/m.', peakGrad, sys.maxGrad);
    end
    if peakSlew > sys.maxSlew * (1 + 1e-9)
        error('Common slice-gradient design exceeds maxSlew: %.6g > %.6g Hz/m/s.', peakSlew, sys.maxSlew);
    end

    report = struct;
    report.iterations = iter;
    report.maxJump = maxJump;
    report.peakGrad = peakGrad;
    report.peakSlew = peakSlew;
    report.peakSlew_T_m_s = peakSlew / sys.gamma;
    report.minTSp = 2*dGz;
    report.minTSpex = 2*dGz;
    report.AGSex = AGSex;
    report.GSspexTargetArea = GSspexTargetArea;
    report.GSspexAreaCorrection = GSspexAreaCorrection;
end

function [gWave, area, timing] = defineSineReadoutWaveTSE( ...
        channel, Tread, adcDelay, readoutBlockDur, sys, Ncycles, ...
        gwave_max, swave_max, physical_slew_max, adc, ...
        waveInfoFlag, debugFlag)
    dt = sys.gradRasterTime;
    [G0, w, TreadRaster] = designWaveAmplitude(Tread, sys, Ncycles, ...
        gwave_max, swave_max, physical_slew_max, 'sine', waveInfoFlag);

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
        'The ADC delay must be on the gradient raster.');
    nPadPost = nBlockSamples-nPadPre-nReadSamples;
    assert(nReadSamples >= 2 && nPadPost >= 0, ...
        'The active sine does not fit in the readout block.');

    activePhase = 2*pi*Ncycles*(0:(nReadSamples-1))/(nReadSamples-1);
    activeWave = G0*sin(activePhase);
    activeWave([1 end]) = 0;
    fullBlockWave = [zeros(1,nPadPre), activeWave, zeros(1,nPadPost)];
    assert(numel(fullBlockWave) == nBlockSamples, ...
        'Internal full-block MB1 sine length mismatch.');

    gWave = mr.makeArbitraryGrad(channel, fullBlockWave, 'system', sys, ...
        'first', 0, 'last', 0);
    gWave.delay = 0;
    area = gWave.area;

    slewPeak_Hz_m_s = max(abs(diff(fullBlockWave)))/dt;
    slewPeak_T_m_s = slewPeak_Hz_m_s/sys.gamma;
    waveSlewLimit = min(swave_max, physical_slew_max);
    assertArbitraryWaveformWithinSystem( ...
        fullBlockWave, 0, 0, sys, 'active sine wave');
    assert(slewPeak_T_m_s <= waveSlewLimit*(1+1e-9), ...
        'Active sine slew %.9g T/m/s exceeds requested limit %.9g T/m/s.', ...
        slewPeak_T_m_s, waveSlewLimit);

    timing = struct;
    timing.nReadSamples = nReadSamples;
    timing.nBlockSamples = nBlockSamples;
    timing.nPadPre = nPadPre;
    timing.nPadPost = nPadPost;
    timing.Twave = readoutBlockDurRaster;
    timing.Tactive = TreadRaster;
    timing.activeStartTime = padPreDur;
    timing.activeEndTime = padPreDur+TreadRaster;
    timing.G0 = G0;
    timing.omega = w;
    timing.freq = w/(2*pi);
    timing.peakAmplitude_mTm = max(abs(activeWave))/sys.gamma*1e3;
    timing.slewPeak_T_m_s = slewPeak_T_m_s;
    timing.requestedFullAmplitudeSlew_T_m_s = ...
        2*pi*Ncycles/TreadRaster*gwave_max*1e-3;
    timing.area = area;

    tol = dt/10;
    adcAcqDur = adc.numSamples*adc.dwell;
    assert(abs(TreadRaster-Tread) < tol, ...
        'Wave duration does not match Tread.');
    assert(abs(adcDelay-adc.delay) < tol, ...
        'Wave and ADC start times differ.');
    assert(abs(TreadRaster-adcAcqDur) < max(tol,sys.adcRasterTime), ...
        'Wave and ADC active durations differ.');
    assert(abs(gWave.delay) < tol, ...
        'Full-block arbitrary Wave must have zero delay.');
    assert(abs(timing.activeStartTime-adc.delay) < tol && ...
           abs(timing.activeEndTime-(adc.delay+adcAcqDur)) < ...
           max(tol,sys.adcRasterTime), ...
        'Wave active samples do not match the ADC window.');
    assert(abs(mr.calcDuration(gWave)-readoutBlockDurRaster) < tol, ...
        'Wave event duration does not match the complete readout block.');

    if debugFlag
        fprintf('\n============= Full-block MB1 sine debug =============\n');
        fprintf('channel                          = %s\n', channel);
        fprintf('active / full-block samples       = %d / %d\n', ...
            nReadSamples, nBlockSamples);
        fprintf('pre/post zero-pad samples         = %d / %d\n', ...
            nPadPre, nPadPost);
        fprintf('active / full-block duration      = %.6f / %.6f ms\n', ...
            TreadRaster*1e3, readoutBlockDurRaster*1e3);
        fprintf('active start / end                = %.6f / %.6f ms\n', ...
            timing.activeStartTime*1e3, timing.activeEndTime*1e3);
        fprintf('actual amplitude / raster slew    = %.6f mT/m / %.6f T/m/s\n', ...
            timing.peakAmplitude_mTm, timing.slewPeak_T_m_s);
        fprintf('wave area                         = %.9g 1/m\n', area);
        fprintf('======================================================\n\n');
    end
end

function [G0_pulseq, w, TreadRaster] = designWaveAmplitude(Tread, sys, Ncycles, gwave_max, swave_max, physical_slew_max, waveName, waveInfoFlag)
    if ~ismember(Ncycles, [10 20])
        error('Ncycles must be one of the supported protocol options: 10 or 20.');
    end
    if Tread <= 0
        error('Tread must be positive.');
    end

    wavepoints = round(Tread / sys.gradRasterTime);
    TreadRaster = wavepoints * sys.gradRasterTime;
    if abs(TreadRaster - Tread) > sys.gradRasterTime/10
        error('Tread %.6f ms is not on the gradient raster %.6f us. Rasterized value would be %.6f ms.', ...
            Tread*1e3, sys.gradRasterTime*1e6, TreadRaster*1e3);
    end

    tWavePeriod = TreadRaster / Ncycles;
    w = 2*pi / tWavePeriod;

    swave_max = min(physical_slew_max, swave_max);
    activePhase = 2*pi*Ncycles*(0:(wavepoints-1))/(wavepoints-1);
    rasterSlewPerUnitAmplitude = max([ ...
        abs(diff(sin(activePhase))), abs(diff(cos(activePhase)))]) / ...
        sys.gradRasterTime;
    requestedAmplitude_T_m = gwave_max*1e-3;
    rasterSlewLimitedAmplitude_T_m = swave_max / rasterSlewPerUnitAmplitude;

    if requestedAmplitude_T_m <= rasterSlewLimitedAmplitude_T_m
        G0_T_m = requestedAmplitude_T_m;
        slewLimited = false;
    else
        G0_T_m = rasterSlewLimitedAmplitude_T_m;
        slewLimited = true;
    end

    G0_pulseq = G0_T_m * sys.gamma;
    actualSlew_T_m_s = rasterSlewPerUnitAmplitude * G0_T_m;

    if waveInfoFlag
        fprintf('\n%s wave amplitude design:\n', waveName);
        fprintf('  Tread/raster      = %.6f / %.6f ms\n', Tread*1e3, TreadRaster*1e3);
        fprintf('  Ncycles           = %d\n', Ncycles);
        fprintf('  period            = %.6f ms\n', tWavePeriod*1e3);
        fprintf('  frequency         = %.3f Hz\n', 1/tWavePeriod);
        fprintf('  omega             = %.6f rad/s\n', w);
        fprintf('  requested gmax    = %.6f mT/m\n', gwave_max);
        fprintf('  actual gmax       = %.6f mT/m\n', G0_T_m*1e3);
        fprintf('  requested slew    = %.6f T/m/s\n', swave_max);
        fprintf('  raster peak slew  = %.6f T/m/s\n', actualSlew_T_m_s);
        fprintf('  slew limited      = %d\n', slewLimited);
    end
end

function assertArbitraryWaveformWithinSystem( ...
        waveform, firstValue, lastValue, system, segmentName)
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
