%% Create an imaging TSE sequence using the 15 uT TBW7 beta150 RF candidate
%
% This scanner-facing research source is derived from the established vMB
% VERSE imaging TSE without modifying the historical production source. It
% uses the simulation-selected TBW7 beta-domain LS candidate with a
% 90-165-150xN RF train. ETL 8--17 uses a prefix of the screened ETL-17
% train. Optional Wave and CAIPI encoding remain available as explicit
% controls at the requested echo spacing below.
%
% The RF is simulation-selected but is not locked, scanner/profile approved,
% or production approved. This source records those flags in the sequence and
% requires an explicit research-verification acknowledgement below. Scanner
% timing, gradient, PNS/CNS, forbidden-frequency, RF power, and SAR checks
% remain mandatory before acquisition.
%
% Design choices:
%   - Keep the essential TSE gradient names and block structure:
%       GSspr, GSspex, GRacq, GRspr, GRspex, GRpreph,
%       GS2, GS3, GS4, GS5, GS7, GR3, GR5, GR6, GR7, GPpre, GPrew.
%   - Use GRE-style system setup with separate sys_sequence/sys_wave and
%     sys_lowPNS limits.
%   - Use sys_sequence for the sequence/ADC and sys_wave for Wave objects.
%   - Use sys_lowPNS for spoilers, crushers, PE pre/rewinders, and RF
%     slice-selective gradients.
%   - Represent gyWave and gzWave as zero-delay arbitrary gradients spanning
%     the full readout block. Gy has zero-valued guard samples. Gz remains at
%     its +G0 endpoint during the guards and is connected through separate
%     pre/post bridges, avoiding the delayed extended-trapezoid representation
%     that has failed Siemens GPA/gradient-health checks.
%   - Group-specific TBW7 MB2 excitation, first-165, and repeated-150
%     refocusing waveforms and their time-varying VERSE slice gradients are
%     loaded from the companion MAT file. Fifty contiguous physical slices
%     are acquired as 25 SMS groups:
%       [1 26], [2 27], ..., [25 50].
%   - The SLC label is sequential in acquired group order (0...24), matching
%     the acquisition-ordered SlicePositions convention required by Siemens
%     2D Pulseq reconstruction. The SLC-to-SMS-group and physical
%     slice-to-group mappings are saved in sequence definitions and PE metadata.
%   - The standard FOV reports the full extent of the fifty physical slice
%     centers plus one slice thickness. Matrix/SlicePositions/SliceGap/Nz/Nslices
%     continue to describe the 25 acquired SMS groups. Separate Physical*
%     definitions retain the fifty-slice geometry required by SMS reconstruction.
%   - CAIPI is implemented with balanced slice-axis blips around every
%     readout. This preserves a fixed CPMG refocusing phase throughout the
%     echo train. Dummy trains reproduce both Wave axes whenever enabled.

clearvars; close all; clc
format long

% This standalone entry uses the simulation-selected 15 uT TBW7 beta150
% vMB-VERSE RF candidate.
%
% Two-axis Wave timing design (2026-09-11):
%   - Gy carries the ADC-aligned sine wave used by the validated y-Wave file.
%   - Gz carries an ADC-aligned cosine wave whose readout-block endpoints
%     are +G0.
%   - The accepted slice-crusher/CAIPI topology is preserved, while its
%     durations and moments are re-derived from the VERSE RF centers. A
%     separate low-PNS three-point bridge is inserted after GS5 and before
%     the cosine readout, and its time reverse is inserted after the
%     readout and before GS7. The bridges do not overlap the crusher or
%     CAIPI waveforms.
%   - Each Skyra bridge keeps the accepted 0.37 ms duration for either
%     supported cycle count; its midpoint is solved from the realized Wave
%     amplitude and moment while remaining within sys_lowPNS.
%     The RF bundle carries system-specific minimum/recommended Wave timing;
%     do not assume that the no-Wave 14 ms validation is sufficient when
%     enabling either Wave axis.
%   - The pre-bridge cancels the realized cosine moment through the ADC
%     center. The post-bridge cancels the remaining realized cosine moment,
%     so the added Gz Wave moment closes before the next refocusing pulse.

%% Path

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
addpath(fullfile(script_dir, 'utils'));

% Machine-specific paths are entered interactively and saved beside this
% script in wave_tse_path_settings.json (ignored by Git). out_path receives
% the v151/, metadata/ and optional v141/ export folders.
pathSettings = configurePathSettings(script_dir);
pulseq_path = pathSettings.pulseq_path;
safe_pns_prediction_path = pathSettings.safe_pns_prediction_path;
out_path = ensureTrailingFilesep(pathSettings.out_path);
system_asc_file = pathSettings.system_asc_file;
addSequenceToolboxPaths(pulseq_path, safe_pns_prediction_path);

mbPulseMatFile = fullfile(script_dir, 'rf_pulses', ...
    'mb2_tse_ls_tbw7_beta150_v1verse_B1max15p000uT.mat');

%% Parameters

% Write options:
%   write_v141_format = false  -> write only current format (v1.5.x)
%   write_v141_format = true   -> also write legacy v1.4.1 format
write_v141_format = false;
writeSequence = true; % candidate remains validation-only until RF/timing review is complete
allowUnapprovedCandidateForResearchScan = true; % acknowledge only for the planned research scan
execute_system_check_code = false; % set true to run PNS/CNS and forbidden-frequency checks

sys_type      = 'skyra'; % CimaX, prisma, or skyra for automatic Wave timing

% ------------------------- sequence timing -------------------------
alpha          = 90;                         % [deg]
ETL            = 17;                         % editable, supported screened prefix: 8...17
TE1Requested   = 107e-3/7;                   % [s], ESP; [] uses MAT/system recommendation
TEeffRequested = 107e-3;                     % [s], snapped to nearest echo within ETL
allowUnvalidatedEchoSpacing = true;           % allow non-screened ESP after Pulseq timing checks
turboFactor    = ETL;
Nechoes        = 1;                          % ICE-visible echo dimension: TSE is single-echo/contrast
necho          = ETL;                        % keep standard TSE loop variable name
SMSFactor      = 2;
Nslices        = 50;                          % physical slices after SMS reconstruction
NsliceGroups   = Nslices / SMSFactor;         % acquired MB2 groups per TR
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
cosineBridgeDurationBudget = 370e-6;         % Skyra separated low-PNS bridge on each side of readout
TRRequested    = 15000e-3;
TR             = TRRequested;                % updated to actual value if requested TR is infeasible

% RF durations and effective centers are read from the vMB VERSE MAT file.
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

% ------------------------- acceleration / undersampling -------------------------
% Every Ry mode preserves its complete centered stride grid. If the number of
% unique ky lines is not divisible by ETL, timing-only slots omit the ADC.
% ------------------------- scanner verification stage -------------------------
%   Stage 1: Ry=1, isUseCaipi=false, caipiShiftFractionFullFOV=0
%   Stage 2: Ry=1, isUseCaipi=true,  caipiShiftFractionFullFOV=0.50
%   Stage 3: Ry=2, isUseCaipi=true,  caipiShiftFractionFullFOV=0.25
%   Ry=3:    Ry=3, isUseCaipi=true,  caipiShiftFractionFullFOV=1/6
%   Ry=4:    Ry=4, isUseCaipi=true,  caipiShiftFractionFullFOV=1/8
Ry = 1;                                      % acceleration along y/LIN
isUseACS = false;                            % collapsed MB ACS cannot provide independent-slice CSMs
Nacs_y = 0;                                  % no appended ACS; use five-SET FLASH calibration SET 4

% ------------------------- SMS / CAIPI -------------------------
isUseCaipi = true;                           % balanced slice-axis CAIPI blips
caipiShiftFractionFullFOV = 0.5;             % relative band-1 shift / full PE FOV
caipiPolarity = 1;                            % change sign to reverse CAIPI shift direction
caipiAdcDemodSign = -1;                       % cancel group-center common z-blip phase

% ------------------------- wave -------------------------
isUseWave_sin = false;
isUseWave_cos = false;
gwave_max     = 10;                           % [mT/m]
swave_max     = 140;                         % [T/m/s]
Ncycles       = 20;                           % supported protocol options: 10 or 20
wavePolarity  = 1;                           % +1 matches the Wave calibration
waveSinChannel = ax.d2;
waveCosChannel = ax.d3;
waveRepresentation = 'fullblock-arbitrary';
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
assert(Ndummy == 1, ...
    'This sequence loop implements exactly one dummy TSE train per SMS group.');
% ------------------------- MB2 pulse and slice geometry -------------------------
if ~isfile(mbPulseMatFile)
    error('MB2 pulse MAT file not found: %s', mbPulseMatFile);
end
MB = load(mbPulseMatFile);
requiredMbFields = { ...
    'rf_ex_single_vmb_hz', ...
    'rf_ref_first_single_vmb_hz', 'rf_ref_repeat_single_vmb_hz', ...
    'rf_ex_vmb_hz_by_group', ...
    'rf_ref_first_vmb_hz_by_group', ...
    'rf_ref_repeat_vmb_hz_by_group', ...
    'rf_ref_single_vmb_hz', 'rf_ref_vmb_hz_by_group', ...
    'rf_ex_duration_s', 'rf_ref_duration_s', ...
    'rf_ex_dwell_s', 'rf_ref_dwell_s', ...
    'rf_ex_center_s', 'rf_ref_center_s', ...
    'gz_ex_times_s', 'gz_ref_times_s', ...
    'gz_ex_vmb_hz_per_m', 'gz_ref_vmb_hz_per_m', ...
    'gz_ex_area_before_center_cycles_per_m', ...
    'gz_ex_area_after_center_cycles_per_m', ...
    'gz_ref_area_before_center_cycles_per_m', ...
    'gz_ref_area_after_center_cycles_per_m', ...
    'gradient_raster_s', 'time_bandwidth_product', ...
    'validated_echo_spacings_s', 'validated_echo_train_length', ...
    'validated_slice_spoiler_factor', ...
    'slice_positions_m', 'sms_groups_1based', ...
    'group_centers_m', 'group_center_freq_hz', ...
    'constant_gradient_equivalent_group_center_freq_hz', ...
    'absolute_position_modulation_convention_version', ...
    'absolute_position_modulation_exponent_sign', ...
    'interleaved_group_order_1based', 'pair_separation_m', ...
    'slice_thickness_m', 'slice_gap_m', 'slice_center_spacing_m', ...
    'slice_grad_hz_per_m', 'physical_slice_count', ...
    'sms_factor', 'sms_group_count', ...
    'caipi_shift_fraction_pe_fov', ...
    'peak_b1_ex_uT', 'peak_b1_ref_first_uT', ...
    'peak_b1_ref_repeat_uT', ...
    'rf_ex_energy_uT2_ms_by_group', ...
    'rf_ref_first_energy_uT2_ms_by_group', ...
    'rf_ref_repeat_energy_uT2_ms_by_group', ...
    'tse_train_rf_energy_uT2_ms_by_group', ...
    'tse_train_first_refocusing_count', ...
    'tse_train_repeat_refocusing_count', ...
    'design_max_b1_uT', 'constraint_safety_factor', ...
    'rf_frequency_offset_mode', 'candidate_tag', 'candidate_status', ...
    'reference_arrays_sha256', 'simulation_selected_candidate', ...
    'stopband_locked_candidate', 'scanner_profile_approved', ...
    'scanner_approved', 'profile_approved', 'production_approved', ...
    'rf_simulation_screened', 'b1_constrained_design', ...
    'b1_constrained_derivative', ...
    'rf_energy_is_scanner_sar_approval', ...
    'legacy_ref_alias_is_repeat_150', ...
    'target_te_s', 'target_te_echo_number'};
for iField = 1:numel(requiredMbFields)
    if ~isfield(MB, requiredMbFields{iField})
        error('MB2 pulse MAT file is missing required field: %s', requiredMbFields{iField});
    end
end
assert(double(MB.absolute_position_modulation_convention_version) == 2 && ...
       double(MB.absolute_position_modulation_exponent_sign) == 1, ...
    ['The MB2 MAT file does not use the scanner-validated absolute-z ', ...
     'modulation convention.']);
expectedCandidateTag = ...
    'ls_tbw7_d2_1e-06_beta150_phasealigned_v1';
expectedCandidateChecksum = ...
    'f34ebbcf0c52c8638eb6dc1da6e778673776fc2465a26e71e8d6e6da5e793b22';
candidateTag = char(string(MB.candidate_tag));
candidateChecksum = char(string(MB.reference_arrays_sha256));
candidateStatus = char(string(MB.candidate_status));
assert(strcmp(candidateTag, expectedCandidateTag), ...
    'Unexpected stopband candidate: %s', candidateTag);
assert(strcmp(candidateChecksum, expectedCandidateChecksum), ...
    'Unexpected TBW7 candidate-array checksum: %s', candidateChecksum);
assert(logical(MB.simulation_selected_candidate) && ...
       logical(MB.rf_simulation_screened), ...
    'The RF must be the simulation-selected, simulation-screened candidate.');
assert(~logical(MB.stopband_locked_candidate) && ...
       ~logical(MB.scanner_profile_approved) && ...
       ~logical(MB.scanner_approved) && ...
       ~logical(MB.profile_approved) && ...
       ~logical(MB.production_approved), ...
    'The candidate approval flags no longer match the unapproved RF state.');
assert(logical(MB.b1_constrained_design) && ...
       ~logical(MB.b1_constrained_derivative), ...
    'Expected an independently designed B1-constrained candidate, not a derivative.');
assert(~logical(MB.rf_energy_is_scanner_sar_approval), ...
    'RF energy must not be represented as scanner SAR approval.');
assert(logical(MB.legacy_ref_alias_is_repeat_150) && ...
       isequal(MB.rf_ref_single_vmb_hz, MB.rf_ref_repeat_single_vmb_hz) && ...
       isequal(MB.rf_ref_vmb_hz_by_group, MB.rf_ref_repeat_vmb_hz_by_group), ...
    'The legacy refocusing aliases must point exactly to repeat-150.');
if ~logical(MB.production_approved)
    warning(['The TBW7 beta150 RF is simulation-selected but is not locked, ', ...
             'scanner/profile approved, or production approved. Use this ', ...
             'export only for the planned research verification.']);
    if writeSequence
        assert(allowUnapprovedCandidateForResearchScan, ...
            ['This RF candidate is not production approved. Set ', ...
             'allowUnapprovedCandidateForResearchScan=true only for the ', ...
             'planned research verification.']);
    end
end
timingSystemKey = sys_type;
if ismember(sys_type, {'C2_simulate_prisma', 'prisma_XA30A'})
    timingSystemKey = 'prisma';
end
hasAutomaticWaveTiming = ...
    isfield(MB, 'wave_min_echo_spacing_s_by_system') && ...
    isfield(MB, 'wave_recommended_echo_spacing_s_by_system') && ...
    isfield(MB.wave_min_echo_spacing_s_by_system, timingSystemKey) && ...
    isfield(MB.wave_recommended_echo_spacing_s_by_system, timingSystemKey);
if isempty(TE1Requested) && hasAutomaticWaveTiming
    TE1 = double( ...
        MB.wave_recommended_echo_spacing_s_by_system.(timingSystemKey));
    TE1SelectionMode = 'mat_system_recommendation';
elseif isempty(TE1Requested)
    error(['TBW7 beta150 RF MAT has no automatic Wave timing estimate for ', ...
           'system %s. Regenerate it with the current Python generator.'], ...
        sys_type);
else
    assert(isscalar(TE1Requested) && isfinite(TE1Requested) && ...
        TE1Requested > 0, ...
        'TE1Requested must be empty or one positive finite scalar in seconds.');
    TE1 = TE1Requested;
    TE1SelectionMode = 'user_requested';
end
TE1RequestedResolved = TE1;
TE1RequestWasAutomatic = isempty(TE1Requested);
TE1Baseline = TE1; % legacy metadata alias; TE1Requested is the editable control
assert(isscalar(TEeffRequested) && isfinite(TEeffRequested) && ...
       TEeffRequested > 0, ...
    'TEeffRequested must be one positive finite scalar.');
effectiveEchoNumber = min(max(round(TEeffRequested/TE1), 1), ETL);
TEeff = effectiveEchoNumber*TE1;
TEeffAdjustment = TEeff-TEeffRequested;
TEeffWasAdapted = abs(TEeffAdjustment) > eps(TEeffRequested)*8;
if TEeffWasAdapted
    warning(['Requested effective TE %.6f ms is not an available echo ', ...
        'center for ESP %.6f ms and ETL=%d. Using echo %d at %.6f ms.'], ...
        TEeffRequested*1e3, TE1*1e3, ETL, effectiveEchoNumber, TEeff*1e3);
end
if hasAutomaticWaveTiming
    estimatedMinimumWaveTE1 = double( ...
        MB.wave_min_echo_spacing_s_by_system.(timingSystemKey));
    recommendedWaveTE1 = double( ...
        MB.wave_recommended_echo_spacing_s_by_system.(timingSystemKey));
else
    estimatedMinimumWaveTE1 = NaN;
    recommendedWaveTE1 = TE1;
end

tEx = double(MB.rf_ex_duration_s);
tRef = double(MB.rf_ref_duration_s);
assert(isscalar(tEx) && isfinite(tEx) && tEx > 0 && ...
       isscalar(tRef) && isfinite(tRef) && tRef > 0, ...
    'The MB2 MAT file must contain positive excitation/refocusing durations.');
rfExDurationTagUs = round(tEx*1e6);
rfRefDurationTagUs = round(tRef*1e6);
assert(abs(tEx*1e6-rfExDurationTagUs) < 1e-6 && ...
       abs(tRef*1e6-rfRefDurationTagUs) < 1e-6, ...
    'RF durations must lie on the 1 us filename-tag raster.');

smsGroups = double(MB.sms_groups_1based);
groupOrder = double(MB.interleaved_group_order_1based(:).');
slicePositions = double(MB.slice_positions_m(:).');
groupCenters = double(MB.group_centers_m(:).');
groupCenterFreqHz = double(MB.group_center_freq_hz(:).');
constantGradientEquivalentGroupCenterFreqHz = ...
    double(MB.constant_gradient_equivalent_group_center_freq_hz(:).');
pairSeparation = double(MB.pair_separation_m);
sliceCenterSpacing = double(MB.slice_center_spacing_m);

assert(SMSFactor == 2, 'This implementation is specialized for MB2.');
assert(Nslices == size(smsGroups,1)*SMSFactor, ...
    'Nslices must equal number of SMS groups times SMSFactor.');
assert(NsliceGroups == size(smsGroups,1), ...
    'NsliceGroups does not match the MB pulse group table.');
assert(numel(slicePositions) == Nslices, ...
    'MB pulse slice-position count does not match Nslices.');
assert(numel(groupCenters) == NsliceGroups, ...
    'MB pulse group-center count does not match NsliceGroups.');
assert(numel(groupCenterFreqHz) == NsliceGroups, ...
    'MB pulse group-frequency count does not match NsliceGroups.');
assert(all(abs(groupCenterFreqHz) < 1e-12), ...
    ['vMB uses group-specific kz(t) modulation. Pulseq RF frequency ', ...
     'offsets must be zero for every group.']);
assert(size(MB.rf_ex_vmb_hz_by_group, 2) == NsliceGroups && ...
       size(MB.rf_ref_first_vmb_hz_by_group, 2) == NsliceGroups && ...
       size(MB.rf_ref_repeat_vmb_hz_by_group, 2) == NsliceGroups, ...
    'The vMB RF matrices must contain one column per SMS group.');
assert(max(abs(MB.rf_ref_first_vmb_hz_by_group - ...
        (165/150)*MB.rf_ref_repeat_vmb_hz_by_group), [], 'all') < 1e-9, ...
    'Every first-165 group waveform must be 165/150 times repeat-150.');
assert(isequal(sort(groupOrder), 1:NsliceGroups), ...
    'The MB pulse group order must contain each slice group exactly once.');
assert(all(smsGroups(:) >= 1 & smsGroups(:) <= Nslices), ...
    'The MB pulse SMS group table contains an invalid physical slice index.');
assert(numel(unique(smsGroups(:))) == Nslices, ...
    'Every physical slice must occur exactly once in smsGroups.');
assert(abs(sliceThickness - 3e-3) < 1e-12, ...
    'The supplied MB2 pulse was designed for 3 mm slice thickness.');
assert(abs(sliceSpacing) < 1e-12, ...
    'The supplied MB2 pulse was designed for zero slice gap.');
assert(abs(double(MB.slice_thickness_m) - sliceThickness) < 1e-12, ...
    'MAT-file slice thickness does not match the sequence setting.');
assert(abs(double(MB.slice_gap_m) - sliceSpacing) < 1e-12, ...
    'MAT-file slice gap does not match the sequence setting.');
assert(abs(pairSeparation - 75e-3) < 1e-12, ...
    'This protocol requires exactly 75 mm separation in every MB2 pair.');
assert(double(MB.physical_slice_count) == Nslices && ...
       double(MB.sms_factor) == SMSFactor && ...
       double(MB.sms_group_count) == NsliceGroups, ...
    'MAT-file slice/group counts do not match the sequence settings.');
validatedEchoTrainLength = double(MB.validated_echo_train_length);
minimumSupportedETL = 8;
assert(isscalar(ETL) && isfinite(ETL) && ETL == round(ETL) && ...
       ETL >= minimumSupportedETL && ETL <= validatedEchoTrainLength, ...
    ['ETL must be an integer from %d through the screened maximum %d. ', ...
     'Values below %d do not retain an echo after the validated echo-7 ', ...
     'contrast, and values above %d require new RF echo-train screening.'], ...
    minimumSupportedETL, validatedEchoTrainLength, minimumSupportedETL, ...
    validatedEchoTrainLength);
assert(double(MB.tse_train_first_refocusing_count) == 1 && ...
       double(MB.tse_train_repeat_refocusing_count) >= ETL-1, ...
    ['The RF bundle does not contain enough screened repeat-150 echoes ', ...
     'for the requested ETL=%d.'], ETL);
isFullValidatedEchoTrain = ETL == validatedEchoTrainLength;
isValidatedEchoTrainPrefix = ETL < validatedEchoTrainLength;
validatedEchoSpacings = double(MB.validated_echo_spacings_s(:).');
echoSpacingMatchTolerance = 1e-9;
isGeneratorValidatedEchoSpacing = any( ...
    abs(validatedEchoSpacings - TE1) <= echoSpacingMatchTolerance) && ...
    logical(MB.rf_simulation_screened);
isUnvalidatedEchoSpacing = ~isGeneratorValidatedEchoSpacing;
assert(isGeneratorValidatedEchoSpacing || allowUnvalidatedEchoSpacing, ...
    ['ESP %.6f ms was not included in the RF echo-train simulation. ', ...
     'Set allowUnvalidatedEchoSpacing=true to permit sequence-timing ', ...
     'validation without claiming RF simulation at this ESP. Simulated ', ...
     'spacings [ms]: %s'], ...
    TE1*1e3, mat2str(validatedEchoSpacings*1e3, 9));
if isGeneratorValidatedEchoSpacing
    echoSpacingValidationMode = 'generator_validated';
else
    echoSpacingValidationMode = 'sequence_timing_only_not_rf_simulated';
    warning(['ESP %.6f ms is allowed above the actual sequence timing ', ...
        'limit, but it was not included in the RF B0/B1 echo-train ', ...
        'simulation. Simulated spacings [ms]: %s'], ...
        TE1*1e3, mat2str(validatedEchoSpacings*1e3, 9));
end
if isFullValidatedEchoTrain && isGeneratorValidatedEchoSpacing
    echoTrainValidationMode = 'generator_validated_full_train';
elseif isFullValidatedEchoTrain
    echoTrainValidationMode = 'screened_length_unvalidated_spacing';
elseif isGeneratorValidatedEchoSpacing
    echoTrainValidationMode = 'generator_validated_prefix_of_full_train';
else
    echoTrainValidationMode = 'screened_prefix_length_unvalidated_spacing';
end
currentRFTrainSimulationScreened = ...
    logical(MB.rf_simulation_screened) && ...
    isGeneratorValidatedEchoSpacing && ...
    (isFullValidatedEchoTrain || isValidatedEchoTrainPrefix);
currentTseTrainRfEnergy_uT2_ms_by_group = ...
    double(MB.rf_ex_energy_uT2_ms_by_group) + ...
    double(MB.rf_ref_first_energy_uT2_ms_by_group) + ...
    (ETL-1)*double(MB.rf_ref_repeat_energy_uT2_ms_by_group);
validatedTseTrainRfEnergy_uT2_ms_by_group = ...
    double(MB.tse_train_rf_energy_uT2_ms_by_group);
energyTolerance = 1e-10*max(1, max(validatedTseTrainRfEnergy_uT2_ms_by_group));
assert(all(currentTseTrainRfEnergy_uT2_ms_by_group <= ...
           validatedTseTrainRfEnergy_uT2_ms_by_group+energyTolerance), ...
    'Requested ETL RF energy exceeds the screened ETL-%d train.', ...
    validatedEchoTrainLength);
assert(abs(double(MB.validated_slice_spoiler_factor) - fspS) < 1e-12, ...
    'The generator echo-train validation spoiler factor does not match fspS.');
rfScreenTargetTE = double(MB.target_te_s);
rfScreenTargetEchoNumber = double(MB.target_te_echo_number);
assert(any(abs(rfScreenTargetTE - ...
        rfScreenTargetEchoNumber*validatedEchoSpacings) <= ...
        echoSpacingMatchTolerance), ...
    'The RF bundle has inconsistent reference TE/ESP simulation metadata.');
expectedSmsGroups = [(1:NsliceGroups).' ((1:NsliceGroups).' + NsliceGroups)];
assert(isequal(smsGroups, expectedSmsGroups), ...
    'Expected 50-slice pairing [1 26], [2 27], ..., [25 50].');
assert(all(abs(diff(slicePositions) - sliceCenterSpacing) < 1e-12), ...
    'The 50 physical slice centers must be uniformly spaced.');

Nz = Nslices;                                % reconstructed physical slice count
% Derive physical coverage directly from the RF geometry metadata.
stackExtent = max(slicePositions) - min(slicePositions) + sliceThickness;
fov = [fovScalar fovScalar stackExtent];      % physical fifty-slice SMS geometry
N = [Nx Ny Nz];                               % physical reconstruction matrix

% ICE receives one SLC index per acquired SMS group, not one index per
% physical slice. Siemens 2D Pulseq reconstruction expects SlicePositions in
% the same order as the sequential SLC labels. Therefore SLC=0 corresponds
% to the first acquired group in groupOrder, SLC=1 to the second acquired
% group, etc. The groupBySlc/slcByGroup maps preserve the link to the
% original geometric SMS-group table used by the RF offsets and recon.
iceNz = NsliceGroups;
groupBySlc = groupOrder(:).';                 % SLC+1 -> geometric SMS-group index
slcByGroup = zeros(1, NsliceGroups);           % geometric SMS-group index -> SLC value
slcByGroup(groupBySlc) = 0:(NsliceGroups-1);
iceSlicePositions = groupCenters(groupBySlc);  % acquisition/SLC order, required by Siemens 2D ICE
iceSliceGap = sliceSpacing;
iceStackExtent = max(iceSlicePositions) - min(iceSlicePositions) + sliceThickness;
iceFov = [fovScalar fovScalar iceStackExtent];
iceMatrix = [Nx Ny iceNz];
assert(abs(iceStackExtent - (NsliceGroups*sliceThickness + max(NsliceGroups-1,0)*sliceSpacing)) < 1e-12, ...
    'ICE stack extent must match the group centers and SliceGap.');
assert(isequal(sort(groupBySlc), 1:NsliceGroups), ...
    'groupBySlc must contain each SMS group exactly once.');

Nx_os = Nx * os_factor;

fprintf(['MB2 geometry: physical slices=%d, groups=%d, MB=%d, ', ...
         'thickness=%.6f mm, gap=%.6f mm, pair separation=%.6f mm\n'], ...
    Nslices, NsliceGroups, SMSFactor, sliceThickness*1e3, ...
    sliceSpacing*1e3, pairSeparation*1e3);
fprintf('  Physical slice centers [mm]: %s\n', mat2str(slicePositions*1e3, 6));
fprintf('  SMS groups, physical 1-based indices:\n');
disp(smsGroups);
fprintf('  Interleaved group order: %s\n', mat2str(groupOrder));
fprintf('  Group carrier centers [mm]: %s\n', mat2str(groupCenters*1e3, 6));
fprintf('  Pulseq group carrier offsets [Hz]: %s (must be zero for vMB)\n', ...
    mat2str(groupCenterFreqHz, 7));
fprintf('  Constant-gradient equivalent offsets [Hz]: %s\n', ...
    mat2str(constantGradientEquivalentGroupCenterFreqHz, 7));
fprintf(['  ICE raw-group metadata: Matrix=[%d %d %d], FOV=[%.3f %.3f %.3f] mm, ', ...
         'SliceGap=%.3f mm, SlicePositions(acq/SLC order)=%s mm\n'], ...
    iceMatrix(1), iceMatrix(2), iceMatrix(3), ...
    iceFov(1)*1e3, iceFov(2)*1e3, iceFov(3)*1e3, ...
    iceSliceGap*1e3, mat2str(iceSlicePositions*1e3, 6));
fprintf('  SLC-to-SMS-group map, groupBySlc: %s\n', mat2str(groupBySlc));
fprintf('  SMS-group-to-SLC map, slcByGroup: %s\n', mat2str(slcByGroup));
fprintf('  Candidate status: %s\n', candidateStatus);

% This TSE uses balanced Gz blips, not RF phase-state cycling, to create the
% ky-dependent relative phase. The RF MAT supplies the MB2 waveform and
% physical geometry; caipiShiftFractionFullFOV above is the applied shift.
if isUseCaipi
    caipiMode = 'blipped-z';
else
    caipiMode = 'none';
end
caipiShiftFractionAliasedFOV = caipiShiftFractionFullFOV * Ry;

if os_factor ~= round(os_factor) || os_factor < 1
    error('os_factor must be a positive integer.');
end
if ~ismember(wavePolarity, [-1 1])
    error('wavePolarity must be +1 or -1.');
end
validGradientChannels = {'x', 'y', 'z'};
assert(any(strcmp(waveSinChannel, validGradientChannels)) && ...
       any(strcmp(waveCosChannel, validGradientChannels)), ...
    'Wave channels must each be one of x, y, or z.');
if isUseWave_sin
    assert(strcmp(waveSinChannel, ax.d2), ...
        'The sine Wave channel must match the PE axis ax.d2.');
end
if isUseWave_cos
    assert(strcmp(waveCosChannel, ax.d3), ...
        ['The cosine Wave channel must match the slice axis ax.d3 because ', ...
         'its separated bridges are serialized with the slice crushers.']);
end
if isUseWave_sin && isUseWave_cos
    assert(~strcmp(waveSinChannel, waveCosChannel), ...
        'Simultaneous sine and cosine Wave events must use different channels.');
end
if Ry ~= round(Ry) || Ry < 1
    error('Ry must be a positive integer.');
end
if isUseCaipi
    assert(isscalar(caipiShiftFractionFullFOV) && ...
           isreal(caipiShiftFractionFullFOV) && ...
           isfinite(caipiShiftFractionFullFOV) && ...
           caipiShiftFractionFullFOV > 0 && ...
           caipiShiftFractionFullFOV < 1, ...
        ['Enabled CAIPI shift must be a finite real scalar strictly between ', ...
         'zero and one full PE FOV. Use caipiPolarity to reverse its direction.']);
else
    assert(caipiShiftFractionFullFOV == 0, ...
        'Set caipiShiftFractionFullFOV=0 when CAIPI is disabled.');
end
if isUseACS && (Nacs_y < 1 || Nacs_y ~= round(Nacs_y))
    error('Nacs_y must be a positive integer when ACS is enabled.');
elseif ~isUseACS && Nacs_y ~= 0
    error('Nacs_y must be zero when ACS is disabled.');
end
if Ndummy ~= 1
    warning('This file keeps the original tse_pulseq.m dummy-loop structure, which uses one dummy train kex=0.');
end
if naverage ~= 1
    warning('This file keeps the original tse_pulseq.m averaging structure; naverage is stored only as metadata.');
end

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

%%% Slice selection: group-specific TBW7 beta150 vMB-VERSE RF/Gz pulses
flipex = alpha*pi/180;
fliprefFirst = rflip(1)*pi/180;
fliprefRepeat = rflip(2)*pi/180;
rfExDwell = double(MB.rf_ex_dwell_s);
rfRefDwell = double(MB.rf_ref_dwell_s);
rfTBP = double(MB.time_bandwidth_product);

% Construct valid Pulseq RF event structures from the single-band VERSE
% templates, then insert the calibrated group-specific MB2 waveforms. The
% explicit center values are based on the realized VERSE kz trajectory.
rfex = mr.makeArbitraryRf( ...
    MB.rf_ex_single_vmb_hz(:), flipex, ...
    'system', sys_lowPNS, ...
    'dwell', rfExDwell, ...
    'center', double(MB.rf_ex_center_s), ...
    'PhaseOffset', rfex_phase, ...
    'use', 'excitation');

rfrefFirst = mr.makeArbitraryRf( ...
    MB.rf_ref_first_single_vmb_hz(:), fliprefFirst, ...
    'system', sys_lowPNS, ...
    'dwell', rfRefDwell, ...
    'center', double(MB.rf_ref_center_s), ...
    'PhaseOffset', rfref_phase, ...
    'use', 'refocusing');

rfrefRepeat = mr.makeArbitraryRf( ...
    MB.rf_ref_repeat_single_vmb_hz(:), fliprefRepeat, ...
    'system', sys_lowPNS, ...
    'dwell', rfRefDwell, ...
    'center', double(MB.rf_ref_center_s), ...
    'PhaseOffset', rfref_phase, ...
    'use', 'refocusing');

% Initialize with group 1. The sequence loops replace these signals with the
% corresponding group column before each train. CAIPI remains implemented
% by balanced z blips, not RF phase cycling.
rfex.signal = double(MB.rf_ex_vmb_hz_by_group(:, 1));
rfrefFirst.signal = double(MB.rf_ref_first_vmb_hz_by_group(:, 1));
rfrefRepeat.signal = double(MB.rf_ref_repeat_vmb_hz_by_group(:, 1));
assert(numel(rfex.signal) == numel(rfex.t), ...
    'MB2 excitation sample count does not match the Pulseq RF event.');
assert(numel(rfrefFirst.signal) == numel(rfrefFirst.t) && ...
       numel(rfrefRepeat.signal) == numel(rfrefRepeat.t), ...
    'MB2 first/repeated refocusing samples do not match their RF events.');
assert(isequal(rfrefFirst.t, rfrefRepeat.t) && ...
       abs(mr.calcRfCenter(rfrefFirst) - mr.calcRfCenter(rfrefRepeat)) < 1e-12, ...
    'First-165 and repeat-150 refocusing events must share one time base.');
refocusingScaleResidual = max(abs( ...
    rfrefFirst.signal - (165/150)*rfrefRepeat.signal));
assert(refocusingScaleResidual < 1e-9, ...
    'First-165 must equal exactly 165/150 times repeat-150 for every group.');

assert(abs(double(MB.gradient_raster_s)-sys.gradRasterTime) < 1e-12, ...
    'The generated VERSE gradient raster does not match the Pulseq system.');

% Keep the effective refocusing center near the center of its complete RF
% block. The TBW7 VERSE design has a sub-raster asymmetric center;
% a minimal pre-RF delay prevents that asymmetry from being amplified when
% one common rasterized crusher duration is used on both sides of the echo.
% This is the same alignment rule used by the matched V3 profile source.
rfRefDelayAlignment = 0;
rfRefCenterForAlignment = mr.calcRfCenter(rfrefRepeat);
bestRfRef = rfrefRepeat;
bestRfRefCenterResidual = Inf;
for iDelayRaster = 0:8
    rfRefCandidate = rfrefRepeat;
    rfRefCandidate.delay = rfrefRepeat.delay + ...
        iDelayRaster*sys.gradRasterTime;
    candidateDuration = ceil( ...
        mr.calcDuration(rfRefCandidate)/sys.gradRasterTime - 1e-12) * ...
        sys.gradRasterTime;
    candidateCenterResidual = candidateDuration - 2*( ...
        rfRefCandidate.delay + rfRefCenterForAlignment);
    if abs(candidateCenterResidual) < abs(bestRfRefCenterResidual)
        bestRfRef = rfRefCandidate;
        bestRfRefCenterResidual = candidateCenterResidual;
    end
end
rfRefDelayAlignment = bestRfRef.delay-rfrefRepeat.delay;
rfrefRepeat = bestRfRef;
rfrefFirst.delay = rfrefFirst.delay + rfRefDelayAlignment;
rfref = rfrefRepeat; % common timing template; sequence loops choose first/repeat RF

% Align each time-varying gradient with the delayed start of its RF shape.
GS2 = mr.makeExtendedTrapezoid(ax.d3, 'system', sys_lowPNS, ...
    'times', rfex.delay + double(MB.gz_ex_times_s(:).'), ...
    'amplitudes', double(MB.gz_ex_vmb_hz_per_m(:).'));
GS4 = mr.makeExtendedTrapezoid(ax.d3, 'system', sys_lowPNS, ...
    'times', rfref.delay + double(MB.gz_ref_times_s(:).'), ...
    'amplitudes', double(MB.gz_ref_vmb_hz_per_m(:).'));
assert(abs(GS2.first) < 1e-9 && abs(GS2.last) < 1e-9 && ...
       abs(GS4.first) < 1e-9 && abs(GS4.last) < 1e-9, ...
    'VERSE slice gradients must begin and end at zero.');

tExwd = ceil(max(mr.calcDuration(rfex), mr.calcDuration(GS2)) / ...
    sys.gradRasterTime - 1e-12) * sys.gradRasterTime;
tRefwd = ceil(max(mr.calcDuration(rfref), mr.calcDuration(GS4)) / ...
    sys.gradRasterTime - 1e-12) * sys.gradRasterTime;
rfExCenter = mr.calcRfCenter(rfex);
rfRefCenter = mr.calcRfCenter(rfref);
rfExCenterAbsolute = rfex.delay + rfExCenter;
rfRefCenterAbsolute = rfref.delay + rfRefCenter;
rfExCenterOffset = rfExCenter - tEx/2;
rfRefCenterOffset = rfRefCenter - tRef/2;

% Derive the standalone pre-refocusing and post-refocusing slice moments
% from the actual gradient areas on either side of the effective RF centers.
exAreaAfterCenter = double(MB.gz_ex_area_after_center_cycles_per_m);
refAreaBeforeCenter = double(MB.gz_ref_area_before_center_cycles_per_m);
refAreaAfterCenter = double(MB.gz_ref_area_after_center_cycles_per_m);
GSsprTargetArea = exAreaAfterCenter*(1+fspS);
GSspexTargetArea = refAreaAfterCenter + GSsprTargetArea - ...
    exAreaAfterCenter - refAreaBeforeCenter;
GSspexAreaCorrection = GSspexTargetArea - exAreaAfterCenter*fspS;

fprintf(['MB2 vMB RF setup: excitation/refocusing durations=%.6f / %.6f ms, ', ...
         'dwells=%.3f / %.3f us, TBW=%.3f\n'], ...
    tEx*1e3, tRef*1e3, rfExDwell*1e6, rfRefDwell*1e6, rfTBP);
fprintf('  Excitation/first-165/repeat-150 samples = %d / %d / %d\n', ...
    numel(rfex.signal), numel(rfrefFirst.signal), numel(rfrefRepeat.signal));
fprintf('  Peak excitation/first-165/repeat-150 B1 = %.3f / %.3f / %.3f uT\n', ...
    double(MB.peak_b1_ex_uT), double(MB.peak_b1_ref_first_uT), ...
    double(MB.peak_b1_ref_repeat_uT));
fprintf(['  ETL=%d current / screened-ETL-%d maximum RF train energy ', ...
         'surrogate = %.6f / %.6f uT^2 ms\n'], ...
    ETL, validatedEchoTrainLength, ...
    max(currentTseTrainRfEnergy_uT2_ms_by_group), ...
    max(validatedTseTrainRfEnergy_uT2_ms_by_group));
fprintf('  Peak excitation/refocusing Gz = %.6f / %.6f mT/m\n', ...
    max(abs(double(MB.gz_ex_vmb_hz_per_m)))/sys.gamma*1e3, ...
    max(abs(double(MB.gz_ref_vmb_hz_per_m)))/sys.gamma*1e3);
fprintf('  RF temporal centers: excitation/refocusing = %.6f / %.6f ms\n', ...
    rfExCenter*1e3, rfRefCenter*1e3);
fprintf('  RF center offsets from waveform midpoint  = %+.6f / %+.6f us\n', ...
    rfExCenterOffset*1e6, rfRefCenterOffset*1e6);
fprintf('  Added refocusing pre-RF center alignment   = %.6f us\n', ...
    rfRefDelayAlignment*1e6);
fprintf('  GSspex RF-center area correction          = %+.9g 1/m\n', ...
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

% Build a center-preserving accelerated PE list. Preserve every requested ky
% line exactly once. When its length is not divisible by ETL, fill the final
% echo table with timing-only no-ADC slots instead of cropping outer k-space.
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
nImageNoAdcSlots = nImageEchoSlots - nnz(imageAdcMask);

assert(nnz(imageAdcMask) == nImagePE, ...
    'Image ADC-slot mask count does not match the acquired PE-line count.');
assert(nnz(imageAdcMask & (PEindexOrderImg == yCenterIdx)) == 1, ...
    'The image PE table must acquire the k-space center exactly once.');
assert(isequal(sort(PEindexOrderImg(imageAdcMask)).', sort(imageYIdx)), ...
    'The PE table must preserve every requested image ky line exactly once.');
if useSparseImageAdcSlots
    assert(nex == ceil(nImagePE/ETL) && ...
           nImageNoAdcSlots == nex*ETL-nImagePE, ...
        'Sparse image ADC slots do not match the complete ky list and ETL.');
else
    assert(all(imageAdcMask(:)) && nImageNoAdcSlots == 0, ...
        'A divisible PE list must retain one ADC in every image echo slot.');
end

% Optional ACS block appended after the image acquisition. Keep the same TSE
% ETL structure by expanding the requested central ACS width to the smallest
% unique central ky block whose size is divisible by ETL. For example,
% Nacs_y=32 with ETL=20 becomes Nacs_y_actual=40. No duplicate ACS ky lines
% are intentionally appended.
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

% Keep the accepted crusher/CAIPI block separate from the cosine bridge.
% tSp is the crusher/CAIPI duration; tSpTotal includes the serial bridge.
cosineBridgeDurationPlanned = double(isUseWave_cos)*cosineBridgeDurationBudget;
adcCenterInReadoutBlock = adc.delay + 0.5*adc.numSamples*adc.dwell;
tSpexTarget = 0.5*TE1 - (tExwd-rfExCenterAbsolute) - ...
    rfRefCenterAbsolute;
tSpTargetPost = 0.5*TE1 - (tRefwd-rfRefCenterAbsolute) - ...
    cosineBridgeDurationPlanned - adcCenterInReadoutBlock;
tSpTargetPre = 0.5*TE1 - (readoutTime-adcCenterInReadoutBlock) - ...
    cosineBridgeDurationPlanned - rfRefCenterAbsolute;
tSpTarget = 0.5*(tSpTargetPost+tSpTargetPre);
tSpTotalTarget = tSpTarget + cosineBridgeDurationPlanned;
tSp = round(tSpTarget / raster) * raster;
tSpex = round(tSpexTarget / raster) * raster;
tSpTotal = tSp + cosineBridgeDurationPlanned;

if tSpTarget <= 0 || tSpTargetPost <= 0 || tSpTargetPre <= 0 || ...
        tSpTotalTarget <= 0
    error(['TE1 is too short for the separated cosine-bridge timing: ', ...
           'total half-spacing %.6f ms, crusher target %.6f ms.'], ...
        tSpTotalTarget*1e3, tSpTarget*1e3);
end
if tSpexTarget <= 0
    error('TE1 is too short for the selected RF/readout timing: tSpexTarget=%.6f ms.', tSpexTarget*1e3);
end

% Readout spoiler and PE timing still use their shortest low-PNS durations.
% The slice-select VERSE waveforms start and end at zero. Standalone
% slice-crusher trapezoids are checked separately below to avoid residual
% first-echo or alternating kz moments.
GRsprNatural = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'area', AGRsprTarget);
GPpreMaxNatural = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
    'Area', maxPhaseArea);
tSpMinXPE = ceil(max([mr.calcDuration(GRsprNatural), mr.calcDuration(GPpreMaxNatural)]) / raster) * raster;

if tSp + raster/10 < tSpMinXPE
    TE1minRefToEcho_prelim = readoutTime + tRefwd + ...
        2*(tSpMinXPE + cosineBridgeDurationPlanned);
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

% Moment from the end of GR5 to the ADC center. During the separated Gz
% pre-bridge, Gx must remain at the readout amplitude so the nonzero GR5/GR6
% boundary stays continuous. Include that serial hold in the prephaser
% correction; otherwise kx would be shifted at the ADC center.
GR6AreaToAdcCenter = GRacq.amplitude * (cosineBridgeDurationPlanned + ...
    adc.delay + 0.5*adc.numSamples*adc.dwell);

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

% The VERSE RF gradients begin and end at zero. The preparation/crusher
% gradients are therefore standalone low-PNS trapezoids whose areas are
% derived from the actual moments around the RF effective centers.
GSspexNatural = mr.makeTrapezoid(ax.d3, 'system', sys_lowPNS, ...
    'area', GSspexTargetArea);
GSsprNatural = mr.makeTrapezoid(ax.d3, 'system', sys_lowPNS, ...
    'area', GSsprTargetArea);
tSpMinZ = ceil(mr.calcDuration(GSsprNatural) / raster) * raster;
tSpexMinZ = ceil(mr.calcDuration(GSspexNatural) / raster) * raster;
tSpMin = max(tSpMinXPE, tSpMinZ);
tSpexMin = max(tSpexMinX, tSpexMinZ);
TE1minExcToRef = 2*((tExwd-rfExCenterAbsolute) + ...
    tSpexMin + rfRefCenterAbsolute);
TE1minRefToEcho = 2*((tRefwd-rfRefCenterAbsolute) + ...
    tSpMin + cosineBridgeDurationPlanned + adcCenterInReadoutBlock);
TE1minEchoToRef = 2*((readoutTime-adcCenterInReadoutBlock) + ...
    cosineBridgeDurationPlanned + tSpMin + rfRefCenterAbsolute);
TE1min = max([TE1minExcToRef, TE1minRefToEcho, TE1minEchoToRef]);

fprintf('\nTSE low-PNS timing feasibility summary:\n');
fprintf('  tExwd / tRefwd from RF events       = %.6f / %.6f ms\n', tExwd*1e3, tRefwd*1e3);
fprintf('  requested TE1                       = %.6f ms\n', TE1*1e3);
fprintf('  tSpex target / low-PNS min / chosen = %.6f / %.6f / %.6f ms\n', ...
    tSpexTarget*1e3, tSpexMin*1e3, tSpex*1e3);
fprintf('  tSp target / low-PNS min / chosen   = %.6f / %.6f / %.6f ms\n', ...
    tSpTarget*1e3, tSpMin*1e3, tSp*1e3);
fprintf('    post-ref / pre-ref targets         = %.6f / %.6f ms\n', ...
    tSpTargetPost*1e3, tSpTargetPre*1e3);
fprintf('  separated cosine bridge per side    = %.6f ms\n', ...
    cosineBridgeDurationPlanned*1e3);
fprintf('  total ref/readout half spacing       = %.6f ms\n', tSpTotal*1e3);
fprintf('  TE1 minimum from exc->ref side      = %.6f ms\n', TE1minExcToRef*1e3);
fprintf('  TE1 minimum from ref->echo side     = %.6f ms\n', TE1minRefToEcho*1e3);
fprintf('  TE1 minimum from echo->ref side     = %.6f ms\n', TE1minEchoToRef*1e3);
fprintf('  TE1 minimum overall                 = %.6f ms\n', TE1min*1e3);

if tSp + raster/10 < tSpMin || tSpex + raster/10 < tSpexMin
    error(['TE1=%.6f ms is too short for the vMB slice timing. ', ...
           'Minimum feasible TE1 is %.6f ms.'], TE1*1e3, TE1min*1e3);
end

GSspex = mr.makeTrapezoid(ax.d3, 'system', sys_lowPNS, ...
    'area', GSspexTargetArea, 'duration', tSpex);
GSspr = mr.makeTrapezoid(ax.d3, 'system', sys_lowPNS, ...
    'area', GSsprTargetArea, 'duration', tSp);
slicePrepRampTime = max(GSspex.riseTime, GSspr.riseTime);

firstEchoMomentLeft = exAreaAfterCenter + GSspex.area + refAreaBeforeCenter;
firstEchoMomentRight = refAreaAfterCenter + GSspr.area;
repeatedEchoMomentError = refAreaBeforeCenter - refAreaAfterCenter;
sliceMomentTolerance = 1e-6;
assert(abs(firstEchoMomentLeft-firstEchoMomentRight) <= sliceMomentTolerance, ...
    'The vMB first-echo slice moment does not balance: %.9g 1/m.', ...
    firstEchoMomentLeft-firstEchoMomentRight);
assert(abs(repeatedEchoMomentError) <= sliceMomentTolerance, ...
    'The vMB refocusing gradient is asymmetric by %.9g 1/m.', ...
    repeatedEchoMomentError);

excitationToRefocusingCenter = (tExwd-rfExCenterAbsolute) + ...
    tSpex + rfRefCenterAbsolute;
refocusingCenterToEcho = (tRefwd-rfRefCenterAbsolute) + tSp + ...
    cosineBridgeDurationPlanned + adcCenterInReadoutBlock;
echoToNextRefocusingCenter = (readoutTime-adcCenterInReadoutBlock) + ...
    cosineBridgeDurationPlanned + tSp + rfRefCenterAbsolute;
centerTimingTolerance = sys.gradRasterTime + sys.rfRasterTime;
assert(abs(excitationToRefocusingCenter-0.5*TE1) <= centerTimingTolerance, ...
    'Excitation-to-refocusing center timing misses TE1/2 by %.6f us.', ...
    (excitationToRefocusingCenter-0.5*TE1)*1e6);
assert(abs(refocusingCenterToEcho-0.5*TE1) <= centerTimingTolerance, ...
    'Refocusing-center-to-echo timing misses TE1/2 by %.6f us.', ...
    (refocusingCenterToEcho-0.5*TE1)*1e6);
assert(abs(echoToNextRefocusingCenter-0.5*TE1) <= centerTimingTolerance, ...
    'Echo-to-next-refocusing-center timing misses TE1/2 by %.6f us.', ...
    (echoToNextRefocusingCenter-0.5*TE1)*1e6);

% Now rebuild the remaining TSE readout prephaser object using the chosen
% common duration. These names intentionally match tse_pulseq.m.
GRspex = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'area', GRacq.area*(1+fspR), 'duration', tSpex); %#ok<NASGU>
GRpreph = mr.makeTrapezoid(ax.d1, 'system', sys_lowPNS, ...
    'Area', AGRpreph, 'duration', tSpex);
AGRspr = GRspr.area;

fprintf('\nTSE selected timing scalars:\n');
fprintf('  TE1=%.6f ms, requested TR=%.6f ms, tSpCrusher=%.6f ms, bridge=%.6f ms, tSpex=%.6f ms\n', ...
    TE1*1e3, TRRequested*1e3, tSp*1e3, cosineBridgeDurationPlanned*1e3, tSpex*1e3);
fprintf('  max standalone slice-prep ramp time    = %.6f us\n', ...
    slicePrepRampTime*1e6);
fprintf('  GSspex / GSspr areas                   = %.9g / %.9g 1/m\n', ...
    GSspex.area, GSspr.area);
fprintf('  first/repeated slice-moment errors     = %.9g / %.9g 1/m\n', ...
    firstEchoMomentLeft-firstEchoMomentRight, repeatedEchoMomentError);
fprintf('  center half-spacings exc/ref, ref/echo, echo/ref = %.6f / %.6f / %.6f ms\n', ...
    excitationToRefocusingCenter*1e3, refocusingCenterToEcho*1e3, ...
    echoToNextRefocusingCenter*1e3);
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
fprintf('  Gx bridge hold + GR6 area to ADC center = %.9g 1/m\n', GR6AreaToAdcCenter);
fprintf('  GRpreph target/realized area          = %.9g / %.9g 1/m\n', AGRpreph, GRpreph.area);

%% Wave gradients: y sine and separated z cosine

if isUseWave_sin
    [gyWave, waveAreaY, waveTiming] = defineSineReadoutWaveTSE( ...
        waveSinChannel, Tread, adc.delay, readoutTime, sys_wave, ...
        Ncycles, wavePolarity, gwave_max, swave_max, physical_slew_max, ...
        adc, waveInfoFlag, waveDebugFlag);
    gyWave.id = seq.registerGradEvent(gyWave);
else
    gyWave = [];
    waveAreaY = 0;
    waveTiming = struct();
    waveTiming.Twave = 0;
    waveTiming.Tactive = 0;
    waveTiming.peakAmplitude_mTm = 0;
    waveTiming.slewPeak_T_m_s = 0;
end

if isUseWave_cos
    [gzWave, waveAreaZ, cosineWaveTiming] = defineCosineReadoutWaveTSE( ...
        waveCosChannel, Tread, adc.delay, readoutTime, sys_wave, ...
        Ncycles, wavePolarity, gwave_max, swave_max, physical_slew_max, ...
        adc, waveInfoFlag, waveDebugFlag);
    gzWave.id = seq.registerGradEvent(gzWave);

    % These serial bridges are deliberately outside GS5/GS7. They preserve
    % the accepted crusher/CAIPI waveforms and make the additional cosine
    % moment zero both at the ADC center and before the next RF pulse.
    [gzCosPreBridge, cosinePreBridgeTiming] = makeFixedDurationThreePointBridge( ...
        waveCosChannel, 0, cosineWaveTiming.edgeAmplitude, ...
        -cosineWaveTiming.areaToAdcCenter, cosineBridgeDurationBudget, sys_lowPNS);
    cosinePostAreaTarget = -(gzCosPreBridge.area + waveAreaZ);
    [gzCosPostBridge, cosinePostBridgeTiming] = makeFixedDurationThreePointBridge( ...
        waveCosChannel, cosineWaveTiming.edgeAmplitude, 0, ...
        cosinePostAreaTarget, cosineBridgeDurationBudget, sys_lowPNS);
    gzCosPreBridge.id = seq.registerGradEvent(gzCosPreBridge);
    gzCosPostBridge.id = seq.registerGradEvent(gzCosPostBridge);

    endpointTol = 1e-9*max(1, abs(cosineWaveTiming.edgeAmplitude));
    assert(abs(gzCosPreBridge.first) <= endpointTol && ...
           abs(gzCosPreBridge.last-cosineWaveTiming.edgeAmplitude) <= endpointTol, ...
        'Separated cosine pre-bridge endpoints do not match 0 -> G0.');
    assert(abs(gzWave.first-cosineWaveTiming.edgeAmplitude) <= endpointTol && ...
           abs(gzWave.last-cosineWaveTiming.edgeAmplitude) <= endpointTol, ...
        'Cosine readout endpoints do not match G0 -> G0.');
    assert(abs(gzCosPostBridge.first-cosineWaveTiming.edgeAmplitude) <= endpointTol && ...
           abs(gzCosPostBridge.last) <= endpointTol, ...
        'Separated cosine post-bridge endpoints do not match G0 -> 0.');

    cosineCenterClosureError = gzCosPreBridge.area + ...
        cosineWaveTiming.areaToAdcCenter;
    cosineCycleClosureError = gzCosPreBridge.area + waveAreaZ + ...
        gzCosPostBridge.area;
    cosineBridgeDurationActual = max(mr.calcDuration(gzCosPreBridge), ...
        mr.calcDuration(gzCosPostBridge));

    assert(abs(mr.calcDuration(gzCosPreBridge) - cosineBridgeDurationBudget) < ...
            sys.gradRasterTime/10, ...
        ['Separated cosine pre-bridge duration %.6f ms differs from the ', ...
         'accepted %.6f ms timing budget. Update TE1 only after reviewing the design.'], ...
        mr.calcDuration(gzCosPreBridge)*1e3, cosineBridgeDurationBudget*1e3);
    assert(abs(mr.calcDuration(gzCosPostBridge) - cosineBridgeDurationBudget) < ...
            sys.gradRasterTime/10, ...
        ['Separated cosine post-bridge duration %.6f ms differs from the ', ...
         'accepted %.6f ms timing budget. Update TE1 only after reviewing the design.'], ...
        mr.calcDuration(gzCosPostBridge)*1e3, cosineBridgeDurationBudget*1e3);
    assert(abs(cosineCenterClosureError) <= 1e-6, ...
        'Cosine moment at the ADC center is not closed: %.9g 1/m.', ...
        cosineCenterClosureError);
    assert(abs(cosineCycleClosureError) <= 1e-6, ...
        'Cosine moment over the complete separated module is not closed: %.9g 1/m.', ...
        cosineCycleClosureError);
    fprintf('\nSeparated cosine-bridge design (no crusher/CAIPI overlap):\n');
    fprintf('  pre duration / midpoint time / midpoint Gz = %.6f / %.6f ms / %+.6f mT/m\n', ...
        cosinePreBridgeTiming.Ttotal*1e3, cosinePreBridgeTiming.Tmid*1e3, ...
        cosinePreBridgeTiming.Gmid/sys.gamma*1e3);
    fprintf('  post duration / midpoint time / midpoint Gz = %.6f / %.6f ms / %+.6f mT/m\n', ...
        cosinePostBridgeTiming.Ttotal*1e3, cosinePostBridgeTiming.Tmid*1e3, ...
        cosinePostBridgeTiming.Gmid/sys.gamma*1e3);
    fprintf('  pre/post peak slew = %.6f / %.6f T/m/s\n', ...
        cosinePreBridgeTiming.peakSlew_T_m_s, ...
        cosinePostBridgeTiming.peakSlew_T_m_s);
else
    gzWave = [];
    gzCosPreBridge = [];
    gzCosPostBridge = [];
    waveAreaZ = 0;
    cosineCenterClosureError = 0;
    cosineCycleClosureError = 0;
    cosineBridgeDurationActual = 0;
    cosineWaveTiming = struct('Twave', 0, 'Tactive', 0, ...
        'peakAmplitude_mTm', 0, 'slewPeak_T_m_s', 0, ...
        'activeStartTime', 0, 'activeEndTime', 0, ...
        'edgeAmplitude', 0, 'areaToAdcCenter', 0);
    cosinePreBridgeTiming = struct();
    cosinePostBridgeTiming = struct();
end

waveStartsWithADC = isUseWave_sin && ...
    abs(waveTiming.activeStartTime - adc.delay) < sys.gradRasterTime/10;
waveEndsWithADC = isUseWave_sin && ...
    abs(waveTiming.activeEndTime - (adc.delay + adc.numSamples*adc.dwell)) < ...
        max(sys.gradRasterTime/10, sys.adcRasterTime);
cosineStartsWithADC = isUseWave_cos && ...
    abs(cosineWaveTiming.activeStartTime - adc.delay) < sys.gradRasterTime/10;
cosineEndsWithADC = isUseWave_cos && ...
    abs(cosineWaveTiming.activeEndTime - (adc.delay + adc.numSamples*adc.dwell)) < ...
        max(sys.gradRasterTime/10, sys.adcRasterTime);
fprintf('\nWave carry-area summary after one TSE readout:\n');
fprintf('  y/sine area = %.9g 1/m\n', waveAreaY);
fprintf('  z/cosine readout area = %.9g 1/m\n', waveAreaZ);
fprintf('  z/cosine area through ADC center = %.9g 1/m\n', ...
    cosineWaveTiming.areaToAdcCenter);
fprintf('  z/cosine center closure error = %.9g 1/m\n', cosineCenterClosureError);
fprintf('  z/cosine full-cycle closure error = %.9g 1/m\n', cosineCycleClosureError);
fprintf('  wave representation = %s\n', waveRepresentation);
fprintf('  wave starts with ADC = %d\n', waveStartsWithADC);
fprintf('  wave ends with ADC   = %d\n', waveEndsWithADC);
fprintf('  cosine starts with ADC = %d\n', cosineStartsWithADC);
fprintf('  cosine ends with ADC   = %d\n', cosineEndsWithADC);
if isUseWave_sin
    assert(waveStartsWithADC && waveEndsWithADC, ...
        'The active sine waveform must start and end with the ADC window.');
    assert(abs(mr.calcDuration(gyWave) - readoutTime) < sys.gradRasterTime/10, ...
        'The full-block sine event must span the complete GR6 readout block.');
    assert(abs(waveAreaY) <= 1e-6, ...
        ['The realized sine area %.9g 1/m is not sufficiently balanced. ', ...
         'The PE rewinder does not compensate Wave carry area.'], waveAreaY);
end
if isUseWave_cos
    assert(cosineStartsWithADC && cosineEndsWithADC, ...
        'The active cosine waveform must start and end with the ADC window.');
    assert(abs(mr.calcDuration(gzWave) - readoutTime) < sys.gradRasterTime/10, ...
        'The full-block cosine event must span the complete GR6 readout block.');
    assert(abs(cosineBridgeDurationActual - cosineBridgeDurationPlanned) < ...
            sys.gradRasterTime/10, ...
        'Actual and planned separated cosine bridge durations differ.');
end

%% Assemble the standalone vMB slice gradients and crushers

% GS2 and GS4 are the explicit time-varying VERSE gradients created with
% the RF events. Since both end at zero, GS3/GS5/GS7 can use standalone
% trapezoids without hidden boundary areas.
GS3 = GSspex;
GS5 = GSspr;
GS7 = GSspr;
if isUseWave_cos
    zBoundaryTol = 1e-9*max(1, abs(cosineWaveTiming.edgeAmplitude));
    assert(abs(GS5.last-gzCosPreBridge.first) <= zBoundaryTol && ...
           abs(gzCosPostBridge.last-GS7.first) <= zBoundaryTol, ...
        'Gz is not continuous at the crusher/separated-bridge boundaries.');
end

% Readout gradient and readout spoilers/crushers. Keep TSE names.
GR3 = GRpreph;

GR5times = [0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR5amp = [0 GRspr.amplitude GRspr.amplitude GRacq.amplitude];
GR5 = mr.makeExtendedTrapezoid(ax.d1, 'times', GR5times, 'amplitudes', GR5amp, 'system', sys_lowPNS);

GR6times = [0 readoutTime];
GR6amp = [GRacq.amplitude GRacq.amplitude];
GR6 = mr.makeExtendedTrapezoid(ax.d1, 'times', GR6times, 'amplitudes', GR6amp, 'system', sys_lowPNS);

if isUseWave_cos
    % Hold Gx at the readout plateau while the separated Gz bridge runs.
    % This preserves continuity from GR5 to GR6 and from GR6 to GR7.
    GRWaveBridgeHold = mr.makeExtendedTrapezoid(ax.d1, ...
        'times', [0 cosineBridgeDurationActual], ...
        'amplitudes', [GRacq.amplitude GRacq.amplitude], ...
        'system', sys_lowPNS);
else
    GRWaveBridgeHold = [];
end

GR7times = [0 GRspr.riseTime GRspr.riseTime+GRspr.flatTime GRspr.riseTime+GRspr.flatTime+GRspr.fallTime];
GR7amp = [GRacq.amplitude GRspr.amplitude GRspr.amplitude 0];
GR7 = mr.makeExtendedTrapezoid(ax.d1, 'times', GR7times, 'amplitudes', GR7amp, 'system', sys_lowPNS);

if isUseWave_cos
    xEndpointTol = 1e-9*max(1, abs(GRacq.amplitude));
    assert(abs(GR5.last-GRWaveBridgeHold.first) <= xEndpointTol && ...
           abs(GRWaveBridgeHold.last-GR6.first) <= xEndpointTol && ...
           abs(GR6.last-GRWaveBridgeHold.first) <= xEndpointTol && ...
           abs(GRWaveBridgeHold.last-GR7.first) <= xEndpointTol, ...
        'Gx is not continuous across GR5/bridge/GR6/bridge/GR7.');
end

fprintf('\nRF/block-boundary duration check for nonzero flat segments:\n');
fprintf('  GS2 duration / rfex duration   = %.6f / %.6f ms\n', mr.calcDuration(GS2)*1e3, mr.calcDuration(rfex)*1e3);
fprintf('  GS4 duration / rfref duration  = %.6f / %.6f ms\n', mr.calcDuration(GS4)*1e3, mr.calcDuration(rfref)*1e3);
fprintf('  GR6 duration / ADC start/end   = %.6f / %.6f / %.6f ms\n', ...
    mr.calcDuration(GR6)*1e3, adc.delay*1e3, (adc.delay + adc.numSamples*adc.dwell)*1e3);
if isUseWave_cos
    fprintf('  separated Gz bridge / Gx hold  = %.6f / %.6f ms\n', ...
        cosineBridgeDurationActual*1e3, mr.calcDuration(GRWaveBridgeHold)*1e3);
end

%% Fill times

tex = mr.calcDuration(GS2, rfex) + mr.calcDuration(GS3);
tref = mr.calcDuration(GS4) + mr.calcDuration(GS5) + mr.calcDuration(GS7) + ...
    readoutTime + 2*cosineBridgeDurationActual;
tend = mr.calcDuration(GS4) + mr.calcDuration(GS5);

tETrain = tex + necho*tref + tend;
minimumTRWithoutFill = NsliceGroups*tETrain;
requestedTRfill = TRRequested/NsliceGroups - tETrain;
TRWasRaisedToMinimum = requestedTRfill < 0;
if TRWasRaisedToMinimum
    TRfill = 0;
else
    TRfill = sys.gradRasterTime * ceil( ...
        requestedTRfill / sys.gradRasterTime - 1e-12);
end
TRActual = NsliceGroups*(tETrain+TRfill);
TR = TRActual;
trComparisonTolerance = max(eps(TRRequested)*8, sys.gradRasterTime*1e-6);
TRWasRasterRounded = ~TRWasRaisedToMinimum && ...
    abs(TRActual-TRRequested) > trComparisonTolerance;
TRWasAdapted = TRWasRaisedToMinimum || TRWasRasterRounded;
if TRWasRaisedToMinimum
    warning(['Requested TR %.6f ms is infeasible for %d sequential MB2 ', ...
        'groups at ETL=%d. The sequence uses the exact no-fill lower-bound ', ...
        'TR %.6f ms.'], ...
        TRRequested*1e3, NsliceGroups, ETL, TRActual*1e3);
elseif TRWasRasterRounded
    fprintf(['Requested TR %.6f ms was rounded up to %.6f ms on the ', ...
        'per-group gradient raster.\n'], TRRequested*1e3, TRActual*1e3);
else
    fprintf('TR fill per group = %.6f ms\n', TRfill*1e3);
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
    fprintf(['  gyWave object/full-block/active durations=', ...
             '%.6f / %.6f / %.6f ms\n'], ...
        mr.calcDuration(gyWave)*1e3, waveTiming.Twave*1e3, ...
        waveTiming.Tactive*1e3);
end

%% ICE labels

% Global 0-based LIN labels.
lblLIN = cell(1, Ny);
for iy = 1:Ny
    lblLIN{iy} = mr.makeLabel('SET', 'LIN', iy - 1);
end

% In SMS, one ADC contains both physical slices in an acquired group. Use
% sequential SLC labels in acquisition order, matching iceSlicePositions.
lblSLC = cell(1, NsliceGroups);
for iSlc = 1:NsliceGroups
    lblSLC{iSlc} = mr.makeLabel('SET', 'SLC', iSlc - 1);
end

lblECO0 = mr.makeLabel('SET', 'ECO', 0);
lblSetRefScan   = mr.makeLabel('SET', 'REF', true);
lblResetRefScan = mr.makeLabel('SET', 'REF', false);
lblResetImaScan = mr.makeLabel('SET', 'IMA', false);
lblSetRefScan.id   = seq.registerLabelEvent(lblSetRefScan);
lblResetRefScan.id = seq.registerLabelEvent(lblResetRefScan);
lblResetImaScan.id = seq.registerLabelEvent(lblResetImaScan);

%% Precompute balanced z-blipped CAIPI objects

caipiPhaseByLin = zeros(1, Ny);
caipiAreaByLin = zeros(1, Ny);
GS5CaipiByLin = cell(1, Ny);
GS7CaipiByLin = cell(1, Ny);
for iy = 1:Ny
    if isUseCaipi
        caipiPhaseUnwrapped = 2*pi*(iy - yCenterIdx)*caipiShiftFractionFullFOV;
        caipiPhaseByLin(iy) = mod(caipiPhaseUnwrapped + pi, 2*pi) - pi;
        caipiAreaByLin(iy) = caipiPolarity * caipiPhaseByLin(iy) / ...
            (2*pi*pairSeparation);
    end
    [GS5CaipiByLin{iy}, GS7CaipiByLin{iy}] = ...
        makeBalancedBlippedCaipiGradients( ...
            GS5, GS7, GSspr, caipiAreaByLin(iy), ax.d3, tSp, sys_lowPNS);
end

% Verify the complete corrected crusher objects, not only the standalone
% triangular CAIPI corrections. This catches duration, sign, or moment
% regressions before the objects enter dummy, imaging, or ACS trains.
caipiMomentTolerance = 1e-6;
for iy = 1:Ny
    caipiPreArea = GS5CaipiByLin{iy}.area - GS5.area;
    caipiRewArea = GS7CaipiByLin{iy}.area - GS7.area;
    assert(abs(caipiPreArea-caipiAreaByLin(iy)) <= caipiMomentTolerance, ...
        'LIN %d CAIPI pre-readout area mismatch: %.9g 1/m.', ...
        iy-1, caipiPreArea-caipiAreaByLin(iy));
    assert(abs(caipiRewArea+caipiAreaByLin(iy)) <= caipiMomentTolerance, ...
        'LIN %d CAIPI rewinder area mismatch: %.9g 1/m.', ...
        iy-1, caipiRewArea+caipiAreaByLin(iy));
    assert(abs(caipiPreArea+caipiRewArea) <= caipiMomentTolerance, ...
        'LIN %d CAIPI pre/rewinder moments do not close.', iy-1);
    assert(abs(mr.calcDuration(GS5CaipiByLin{iy})-tSp) < ...
            sys.gradRasterTime/10 && ...
           abs(mr.calcDuration(GS7CaipiByLin{iy})-tSp) < ...
            sys.gradRasterTime/10, ...
        'LIN %d CAIPI crusher duration does not match tSp.', iy-1);
end

fprintf('\nMB2 CAIPI encoding:\n');
fprintf('  enabled                           = %d\n', isUseCaipi);
fprintf('  mode                              = %s\n', caipiMode);
fprintf('  nominal MAT shift                 = %.6f of full PE FOV\n', ...
    double(MB.caipi_shift_fraction_pe_fov));
fprintf('  applied shift                     = %.6f of full PE FOV\n', ...
    caipiShiftFractionFullFOV);
fprintf('  applied shift                     = %.6f of accelerated alias FOV\n', ...
    caipiShiftFractionAliasedFOV);
fprintf('  pair separation                   = %.6f mm\n', pairSeparation*1e3);
fprintf('  ADC common-phase demodulation     = %d, sign %+d\n', isUseCaipi, caipiAdcDemodSign);
fprintf('  CAIPI z area range                = %+.9g ... %+.9g 1/m\n', ...
    min(caipiAreaByLin), max(caipiAreaByLin));
imageCaipiPhases = unique(round(caipiPhaseByLin(imageYIdx), 12));
if isUseCaipi
    assert(numel(imageCaipiPhases) >= 2, ...
        ['CAIPI is enabled, but the acquired Ry=%d ky grid sees only one ', ...
         'relative phase state at shift %.9g of the full FOV. Select a ', ...
         'shift compatible with the accelerated ky stride.'], ...
        Ry, caipiShiftFractionFullFOV);
else
    assert(isequal(imageCaipiPhases, 0), ...
        'Disabled CAIPI must produce one zero relative-phase state.');
end
fprintf('  image-line relative phases [rad]  = %s\n', mat2str(imageCaipiPhases, 6));
if isAcquireACS
    acsCaipiPhases = unique(round(caipiPhaseByLin(acsYIdx), 12));
    fprintf('  ACS relative phases [rad]         = %s\n', mat2str(acsCaipiPhases, 6));
end

%% Define sequence blocks

fprintf('\nInterleaved MB2 group-specific vMB RF table:\n');
fprintf('  acquisition position | SLC | SMS group | physical slices | centers [mm] | group center [mm] | Pulseq offset [Hz]\n');
for iGroupOrder = 1:NsliceGroups
    g = groupBySlc(iGroupOrder);
    physical = smsGroups(g,:);
    fprintf('  %20d | %3d | %9d | %s | %s | %+14.6f | %+12.6f\n', ...
        iGroupOrder, iGroupOrder-1, g, mat2str(physical), ...
        mat2str(slicePositions(physical)*1e3, 6), iceSlicePositions(iGroupOrder)*1e3, groupCenterFreqHz(g));
end

readoutWaveEvents = {GR6};
readoutWaveChannels = {ax.d1};
if isUseWave_sin
    readoutWaveEvents{end+1} = gyWave;
    readoutWaveChannels{end+1} = waveSinChannel;
end
if isUseWave_cos
    readoutWaveEvents{end+1} = gzWave;
    readoutWaveChannels{end+1} = waveCosChannel;
end
assert(numel(unique(readoutWaveChannels)) == numel(readoutWaveChannels), ...
    ['GR6 and enabled Wave events must occupy distinct physical gradient ', ...
     'channels within the common readout block.']);
nExcitationBlocks = 0;
nFirstRefocusingBlocks = 0;
nRepeatRefocusingBlocks = 0;

for kex = 0:nex % kex=0 is one dummy train; kex=1:nex are image trains
    for iGroupOrder = 1:NsliceGroups
        g = groupBySlc(iGroupOrder);
        rfex.signal = double(MB.rf_ex_vmb_hz_by_group(:, g));
        rfrefFirst.signal = double(MB.rf_ref_first_vmb_hz_by_group(:, g));
        rfrefRepeat.signal = double(MB.rf_ref_repeat_vmb_hz_by_group(:, g));
        rfex.freqOffset = 0;
        rfrefFirst.freqOffset = 0;
        rfrefRepeat.freqOffset = 0;
        rfex.phaseOffset = rfex_phase;
        rfrefFirst.phaseOffset = rfref_phase;
        rfrefRepeat.phaseOffset = rfref_phase;

        seq.addBlock(GS2, rfex);
        nExcitationBlocks = nExcitationBlocks + 1;
        seq.addBlock(GS3, GR3);

        for kech = 1:necho
            if kech == 1
                rfrefThis = rfrefFirst;
                nFirstRefocusingBlocks = nFirstRefocusingBlocks + 1;
            else
                rfrefThis = rfrefRepeat;
                nRepeatRefocusingBlocks = nRepeatRefocusingBlocks + 1;
            end
            isImageAdcSlot = kex > 0 && imageAdcMask(kech,kex);
            if isImageAdcSlot
                iy = PEindexOrderImg(kech,kex);
                phaseArea = phaseAreas(kech,kex);
            else
                % Dummy and sparse Ry3 slots keep identical echo timing but
                % use zero PE/CAIPI moment and do not acquire an ADC.
                iy = yCenterIdx;
                phaseArea = 0;
            end
            GPpre = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                'Area', phaseArea, 'Duration', tSp);
            GPrew = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                'Area', -phaseArea, 'Duration', tSp);
            GS5This = GS5CaipiByLin{iy};
            GS7This = GS7CaipiByLin{iy};
            adcThis = adc;
            commonCaipiPhase = 2*pi*caipiAreaByLin(iy)*groupCenters(g);
            adcThis.phaseOffset = angle(exp(1i*caipiAdcDemodSign*commonCaipiPhase));

            seq.addBlock(GS4, rfrefThis);
            seq.addBlock(GS5This, GR5, GPpre);
            if isUseWave_cos
                seq.addBlock(gzCosPreBridge, GRWaveBridgeHold);
            end
            if isImageAdcSlot
                seq.addBlock(readoutWaveEvents{:}, adcThis, lblLIN{iy}, lblSLC{iGroupOrder}, ...
                    lblECO0, lblResetRefScan, lblResetImaScan);
            else
                seq.addBlock(readoutWaveEvents{:});
            end
            if isUseWave_cos
                seq.addBlock(gzCosPostBridge, GRWaveBridgeHold);
            end
            seq.addBlock(GS7This, GR7, GPrew);
        end
        seq.addBlock(GS4);
        seq.addBlock(GS5);
        if TRfill > 0
            seq.addBlock(delayTR);
        end
    end
end

% ACS/refscan block, using the same MB2 groups and continuous CAIPI phase
% ramp evaluated at every fully sampled ACS ky location.
if isAcquireACS
    fprintf('\nAppending MB2 ACS block at end: %d trains/group, %d lines/group.\n', ...
        acsNex, nAcsPE);
    for kexACS = 1:acsNex
        for iGroupOrder = 1:NsliceGroups
            g = groupBySlc(iGroupOrder);
            rfex.signal = double(MB.rf_ex_vmb_hz_by_group(:, g));
            rfrefFirst.signal = double(MB.rf_ref_first_vmb_hz_by_group(:, g));
            rfrefRepeat.signal = double(MB.rf_ref_repeat_vmb_hz_by_group(:, g));
            rfex.freqOffset = 0;
            rfrefFirst.freqOffset = 0;
            rfrefRepeat.freqOffset = 0;
            rfex.phaseOffset = rfex_phase;
            rfrefFirst.phaseOffset = rfref_phase;
            rfrefRepeat.phaseOffset = rfref_phase;

            seq.addBlock(GS2, rfex);
            nExcitationBlocks = nExcitationBlocks + 1;
            seq.addBlock(GS3, GR3);

            for kech = 1:necho
                if kech == 1
                    rfrefThis = rfrefFirst;
                    nFirstRefocusingBlocks = nFirstRefocusingBlocks + 1;
                else
                    rfrefThis = rfrefRepeat;
                    nRepeatRefocusingBlocks = nRepeatRefocusingBlocks + 1;
                end
                iy = PEindexOrderACS(kech,kexACS);
                phaseArea = phaseAreasACS(kech,kexACS);
                GPpre = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                    'Area', phaseArea, 'Duration', tSp);
                GPrew = mr.makeTrapezoid(ax.d2, 'system', sys_lowPNS, ...
                    'Area', -phaseArea, 'Duration', tSp);
                GS5This = GS5CaipiByLin{iy};
                GS7This = GS7CaipiByLin{iy};
                adcThis = adc;
                commonCaipiPhase = 2*pi*caipiAreaByLin(iy)*groupCenters(g);
                adcThis.phaseOffset = angle(exp(1i*caipiAdcDemodSign*commonCaipiPhase));

                seq.addBlock(GS4, rfrefThis);
                seq.addBlock(GS5This, GR5, GPpre);
                if isUseWave_cos
                    seq.addBlock(gzCosPreBridge, GRWaveBridgeHold);
                end
                seq.addBlock(readoutWaveEvents{:}, adcThis, lblLIN{iy}, lblSLC{iGroupOrder}, ...
                    lblECO0, lblSetRefScan, lblResetImaScan);
                if isUseWave_cos
                    seq.addBlock(gzCosPostBridge, GRWaveBridgeHold);
                end
                seq.addBlock(GS7This, GR7, GPrew);
            end
            seq.addBlock(GS4);
            seq.addBlock(GS5);
            if TRfill > 0
                seq.addBlock(delayTR);
            end
        end
    end
end

expectedTseTrainCount = (nex + Ndummy + acsNex)*NsliceGroups;
assert(nExcitationBlocks == expectedTseTrainCount && ...
       nFirstRefocusingBlocks == expectedTseTrainCount && ...
       nRepeatRefocusingBlocks == expectedTseTrainCount*(ETL-1), ...
    'Constructed RF event counts do not match 90-165-150x%d per TSE train.', ...
    ETL-1);
fprintf(['Constructed RF train counts (90-165-150x%d): ', ...
         'trains/excitations/first-165/repeat-150 = %d / %d / %d / %d\n'], ...
    ETL-1, ...
    expectedTseTrainCount, nExcitationBlocks, ...
    nFirstRefocusingBlocks, nRepeatRefocusingBlocks);

%% Timing check

[ok, error_report] = seq.checkTiming;
if ok
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
    error('Sequence timing check failed; no sequence files were written.');
end

%% ICE label, ADC/readout, and k-space diagnostics

adc_lbl = seq.evalLabels('evolution','adc');

expectedImageReadouts = nImagePE * NsliceGroups;
expectedAcsReadouts = acsNex * NsliceGroups * necho;
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
assert(all(adc_lbl.SLC >= 0) && max(adc_lbl.SLC) < NsliceGroups, 'SLC labels out of range.');
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

% Verify the LIN ordering independently for every acquired SLC. SLC is
% sequential in acquisition order; groupBySlc maps SLC+1 back to the
% geometric SMS group table used by RF offsets and recon.
for iSlc = 1:NsliceGroups
    g = groupBySlc(iSlc);
    if isfield(adc_lbl, 'REF')
        imageMaskThisGroup = (adc_lbl.SLC(:) == iSlc-1) & (adc_lbl.REF(:) == 0);
    else
        imageMaskThisGroup = (adc_lbl.SLC(:) == iSlc-1);
    end
    assert(sum(imageMaskThisGroup) == nImagePE, ...
        'SLC %d / SMS group %d has an unexpected number of image readouts.', iSlc-1, g);
    actualImageLINThisGroup = adc_lbl.LIN(imageMaskThisGroup);
    expectedImageLINThisGroup = LINorderImg(imageAdcMask);
    assert(isequal(actualImageLINThisGroup(:), expectedImageLINThisGroup(:)), ...
        'SLC %d / SMS group %d image LIN order does not match LINorderImg.', iSlc-1, g);

    if isAcquireACS
        acsMaskThisGroup = (adc_lbl.SLC(:) == iSlc-1) & (adc_lbl.REF(:) ~= 0);
        assert(sum(acsMaskThisGroup) == acsNex*necho, ...
            'SLC %d / SMS group %d has an unexpected number of ACS readouts.', iSlc-1, g);
        actualAcsLINThisGroup = adc_lbl.LIN(acsMaskThisGroup);
        assert(isequal(actualAcsLINThisGroup(:), LINorderACS(:)), ...
            'SLC %d / SMS group %d ACS LIN order does not match LINorderACS.', iSlc-1, g);
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
        adcLinOneBased = double(adc_lbl.LIN(:)) + 1;
        assert(numel(adcLinOneBased) == nReadouts, ...
            'ADC LIN label count does not match the trajectory readout count.');
        assert(all(adcLinOneBased >= 1 & adcLinOneBased <= Ny), ...
            'ADC LIN labels fall outside the encoded PE matrix.');
        expectedCaipiKz = caipiAreaByLin(adcLinOneBased);
        kzCenterResidual = kzCenterAll(:) - expectedCaipiKz(:);
        kzResidualPeakToPeak = max(kzCenterResidual) - min(kzCenterResidual);
        kzEchoTolerance = 1e-2*deltak;
        fprintf('\nSlice-axis kz diagnostics from calculateKspacePP:\n');
        fprintf('  kz ADC center range        = %+.9g ... %+.9g 1/m\n', min(kzCenterAll), max(kzCenterAll));
        fprintf('  kz ADC center peak-to-peak = %.9g 1/m\n', max(kzCenterAll) - min(kzCenterAll));
        fprintf('  expected CAIPI kz range    = %+.9g ... %+.9g 1/m\n', ...
            min(expectedCaipiKz), max(expectedCaipiKz));
        fprintf('  CAIPI-adjusted kz residual = %+.9g ... %+.9g 1/m\n', ...
            min(kzCenterResidual), max(kzCenterResidual));
        fprintf('  residual peak-to-peak      = %.9g 1/m (limit %.9g)\n', ...
            kzResidualPeakToPeak, kzEchoTolerance);
        assert(kzResidualPeakToPeak <= kzEchoTolerance, ...
            ['MB2 ADC kz contains an echo-dependent residual after removing ', ...
             'the designed CAIPI state: %.9g 1/m exceeds %.9g 1/m.'], ...
            kzResidualPeakToPeak, kzEchoTolerance);
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

% Optional scanner-specific PNS/CNS and forbidden-frequency checks. The
% gradient model is the system_asc_file path setting, which must match the
% sys_type entry below, and the sequence definition records whether these
% checks were enabled for this export.
if execute_system_check_code
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
        ascName = requireSystemCheckInputs(sys_type, sysMap{idx, 2}, ...
            system_asc_file, safe_pns_prediction_path);
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

% Report the full physical coverage through the standard FOV for the scanner
% UI. The remaining standard/ICE-facing geometry describes the 25 acquired
% SMS groups.
% Siemens 2D Pulseq reconstruction recognizes the exact keys
% SlicePositions, SliceThickness, and SliceGap. SlicePositions must be in
% sequential SLC/acquisition order, not geometric group order.
seq.setDefinition('FOV', [iceFov(1:2) fov(3)]);
seq.setDefinition('Matrix', iceMatrix);
seq.setDefinition('SliceThickness', sliceThickness);
seq.setDefinition('SliceGap', iceSliceGap);
seq.setDefinition('SliceSpacing', sliceSpacing); % custom backward-compatible alias
seq.setDefinition('SliceCenterSpacing', sliceCenterSpacing);
seq.setDefinition('SliceStackExtent', iceStackExtent);
seq.setDefinition('SliceCoverageExtent', iceStackExtent);
seq.setDefinition('SlicePositions', iceSlicePositions);
seq.setDefinition('SliceOrder', groupBySlc);

% Preserve the complete physical fifty-slice geometry under explicit custom
% names for SMS-aware reconstruction and downstream bookkeeping.
seq.setDefinition('PhysicalFOV', fov);
seq.setDefinition('PhysicalMatrix', N);
seq.setDefinition('PhysicalSliceThickness', sliceThickness);
seq.setDefinition('PhysicalSliceSpacing', sliceSpacing);
seq.setDefinition('PhysicalSliceCenterSpacing', sliceCenterSpacing);
seq.setDefinition('PhysicalSliceStackExtent', stackExtent);
seq.setDefinition('PhysicalSlicePositions', slicePositions);
seq.setDefinition('GroupFOV', iceFov);
seq.setDefinition('GroupMatrix', iceMatrix);
seq.setDefinition('GroupSliceGap', iceSliceGap);
seq.setDefinition('GroupSliceStackExtent', iceStackExtent);
seq.setDefinition('GroupSlicePositions', iceSlicePositions);
seq.setDefinition('GroupSlicePositionsGeometric', groupCenters);
seq.setDefinition('SMSGroupsFlat', smsGroups(:).');
seq.setDefinition('SMSGroupsShape', size(smsGroups));
seq.setDefinition('SMSGroupOrder', groupOrder);
seq.setDefinition('SMSGroupBySLC', groupBySlc);
seq.setDefinition('SMSSLCByGroup', slcByGroup);
seq.setDefinition('SMSGroupCenters', groupCenters);
seq.setDefinition('SMSGroupCentersAcqOrder', iceSlicePositions);
seq.setDefinition('SMSGroupCenterFreqHz', groupCenterFreqHz);
seq.setDefinition('SMSConstantGradientEquivalentCenterFreqHz', ...
    constantGradientEquivalentGroupCenterFreqHz);
seq.setDefinition('VMBGroupSpecificRF', 1);
seq.setDefinition('VMBPulseqFrequencyOffsetHz', 0);
seq.setDefinition('SMSPairSeparation', pairSeparation);
seq.setDefinition('SliceInterleaved', double(NsliceGroups > 1));
seq.setDefinition('TE1', TE1);
seq.setDefinition('TE1Baseline', TE1Baseline);
seq.setDefinition('TE1Requested', TE1RequestedResolved);
seq.setDefinition('TE1RequestWasAutomatic', double(TE1RequestWasAutomatic));
seq.setDefinition('TE1SelectionMode', TE1SelectionMode);
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
seq.setDefinition('TSETrainCount', expectedTseTrainCount);
seq.setDefinition('ExcitationRFEventCount', nExcitationBlocks);
seq.setDefinition('FirstRefocusingRFEventCount', nFirstRefocusingBlocks);
seq.setDefinition('RepeatRefocusingRFEventCount', nRepeatRefocusingBlocks);
seq.setDefinition('Nechoes', Nechoes);              % ICE-visible echo dimension
seq.setDefinition('EchoTrainLength', ETL);
seq.setDefinition('TurboFactor', turboFactor);
seq.setDefinition('ETL', ETL);
seq.setDefinition('MinimumSupportedETL', minimumSupportedETL);
seq.setDefinition('ValidatedEchoTrainLength', validatedEchoTrainLength);
seq.setDefinition('EchoTrainValidationMode', echoTrainValidationMode);
seq.setDefinition('ValidatedEchoTrainPrefix', ...
    double(isValidatedEchoTrainPrefix));
seq.setDefinition('CurrentRFTrainSimulationScreened', ...
    double(currentRFTrainSimulationScreened));
seq.setDefinition('necho', necho);
% Generic slice-count fields are raw acquisition dimensions and therefore
% use the number of SMS groups. Physical/reconstructed counts remain explicit.
seq.setDefinition('Nslices', NsliceGroups);
seq.setDefinition('NPhysicalSlices', Nslices);
seq.setDefinition('NsliceGroups', NsliceGroups);
seq.setDefinition('SMSFactor', SMSFactor);
seq.setDefinition('Nx', Nx);
seq.setDefinition('Ny', Ny);
seq.setDefinition('Nz', iceNz);
seq.setDefinition('PhysicalNz', Nz);
seq.setDefinition('Nx_os', Nx_os);
seq.setDefinition('ReadoutOversamplingFactor', os_factor);
seq.setDefinition('roDuration', Tread);
seq.setDefinition('readoutTime', readoutTime);
seq.setDefinition('dwell', dwell);
seq.setDefinition('SlicePrepRampTime', slicePrepRampTime);
seq.setDefinition('SlicePrepRampTime_us', slicePrepRampTime*1e6);
seq.setDefinition('Ndummy', Ndummy);
seq.setDefinition('Averages', naverage);
seq.setDefinition('naverage', naverage);
seq.setDefinition('Ry', Ry);
seq.setDefinition('AccelerationFactor', Ry);
seq.setDefinition('AccelerationFactorPE', Ry);
seq.setDefinition('EffectiveRy', effectiveRy);
seq.setDefinition('Ny_meas', nImagePE);
seq.setDefinition('Nex', nex);
seq.setDefinition('PartialFourierEnabled', 0);
seq.setDefinition('PartialFourierFactor', 1);
seq.setDefinition('ImageUniqueKyCount', nImagePE);
seq.setDefinition('ImageFullKyGridPreserved', ...
    double(Ry == 1 && nImagePE == Ny));
if useSparseImageAdcSlots
    seq.setDefinition('ImageEchoSlotCount', nImageEchoSlots);
    seq.setDefinition('ImageADCCount', nImagePE);
    seq.setDefinition('ImageNoADCSlotCount', nImageNoAdcSlots);
    seq.setDefinition('ImageADCSlotMask', double(imageAdcMask(:).'));
end
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
phaseResolution = (iceFov(1)/Nx) / (iceFov(2)/Ny);
seq.setDefinition('kSpaceCenterLine', centerLineIdx - 1);
seq.setDefinition('PhaseResolution', phaseResolution);
seq.setDefinition('SliceResolution', 1);
seq.setDefinition('Slices', NsliceGroups);
seq.setDefinition('ReconstructedSlices', Nslices);
seq.setDefinition('kSpaceCenterSlice', 0);
seq.setDefinition('ReadoutPolarity', 'TSE-refocused-readout');
seq.setDefinition('CAIPIMode', caipiMode);
seq.setDefinition('CAIPIEnabled', double(isUseCaipi));
seq.setDefinition('CAIPIPolarity', caipiPolarity);
seq.setDefinition('CAIPIAdcDemodSign', caipiAdcDemodSign);
seq.setDefinition('CAIPICommonPhaseDemodulation', double(isUseCaipi));
seq.setDefinition('CAIPIShiftFractionFullFOV', caipiShiftFractionFullFOV);
seq.setDefinition('CAIPIShiftFractionAliasedFOV', caipiShiftFractionAliasedFOV);
seq.setDefinition('CAIPIImagePhaseStates', imageCaipiPhases);
seq.setDefinition('CAIPIAreaRange', [min(caipiAreaByLin) max(caipiAreaByLin)]);
seq.setDefinition('CAIPIMomentValidationComplete', 1);
seq.setDefinition('CAIPIMomentTolerance', caipiMomentTolerance);
seq.setDefinition('MBPulseMatFile', mbPulseMatFile);
seq.setDefinition('RFDesign', 'ls_tbw7_beta150_phasealigned_v1verse_B1max15uT');
seq.setDefinition('RFCandidateTag', candidateTag);
seq.setDefinition('RFCandidateStatus', candidateStatus);
seq.setDefinition('RFReferenceArraysSha256', candidateChecksum);
seq.setDefinition('StopbandCandidateTag', candidateTag);
seq.setDefinition('StopbandLockedCandidate', ...
    double(logical(MB.stopband_locked_candidate)));
seq.setDefinition('SimulationSelectedCandidate', ...
    double(logical(MB.simulation_selected_candidate)));
seq.setDefinition('RFBundleSimulationScreened', ...
    double(logical(MB.rf_simulation_screened)));
seq.setDefinition('RFSimulationScreened', ...
    double(currentRFTrainSimulationScreened));
seq.setDefinition('B1ConstrainedDesign', ...
    double(logical(MB.b1_constrained_design)));
seq.setDefinition('B1ConstrainedDerivative', ...
    double(logical(MB.b1_constrained_derivative)));
seq.setDefinition('RFOnlyEchoSpacingValidationComplete', ...
    double(isGeneratorValidatedEchoSpacing));
seq.setDefinition('StopbandScannerProfileApproved', ...
    double(logical(MB.scanner_profile_approved)));
seq.setDefinition('ScannerApproved', double(logical(MB.scanner_approved)));
seq.setDefinition('ProfileApproved', double(logical(MB.profile_approved)));
seq.setDefinition('StopbandProductionApproved', ...
    double(logical(MB.production_approved)));
seq.setDefinition('ResearchVerificationOnly', ...
    double(~logical(MB.production_approved)));
seq.setDefinition('RFEnergyIsScannerSARApproval', ...
    double(logical(MB.rf_energy_is_scanner_sar_approval)));
seq.setDefinition('ValidatedEchoSpacings', validatedEchoSpacings);
seq.setDefinition('GeneratorValidatedEchoSpacing', ...
    double(isGeneratorValidatedEchoSpacing));
seq.setDefinition('UnvalidatedEchoSpacingAllowed', ...
    double(allowUnvalidatedEchoSpacing));
seq.setDefinition('CurrentEchoSpacingRFSimulated', ...
    double(~isUnvalidatedEchoSpacing));
seq.setDefinition('EchoSpacingValidationMode', echoSpacingValidationMode);
seq.setDefinition('RFScreenTargetTE', rfScreenTargetTE);
seq.setDefinition('RFScreenTargetEchoNumber', rfScreenTargetEchoNumber);
seq.setDefinition('VMBAbsolutePositionConventionVersion', ...
    double(MB.absolute_position_modulation_convention_version));
seq.setDefinition('VMBAbsolutePositionExponentSign', ...
    double(MB.absolute_position_modulation_exponent_sign));
seq.setDefinition('MBPeakB1Ex_uT', double(MB.peak_b1_ex_uT));
seq.setDefinition('MBPeakB1Ref_uT', double(MB.peak_b1_ref_first_uT));
seq.setDefinition('MBPeakB1RefFirst_uT', double(MB.peak_b1_ref_first_uT));
seq.setDefinition('MBPeakB1RefRepeat_uT', double(MB.peak_b1_ref_repeat_uT));
seq.setDefinition('MBDesignMaxB1_uT', double(MB.design_max_b1_uT));
if isfinite(estimatedMinimumWaveTE1)
    seq.setDefinition('WaveEstimatedMinimumTE1', estimatedMinimumWaveTE1);
end
seq.setDefinition('WaveRecommendedTE1', recommendedWaveTE1);
seq.setDefinition('MBRFEnergyEx_uT2ms', ...
    max(double(MB.rf_ex_energy_uT2_ms_by_group)));
seq.setDefinition('MBRFEnergyRefFirst_uT2ms', ...
    max(double(MB.rf_ref_first_energy_uT2_ms_by_group)));
seq.setDefinition('MBRFEnergyRefRepeat_uT2ms', ...
    max(double(MB.rf_ref_repeat_energy_uT2_ms_by_group)));
seq.setDefinition('MBRFTrainEnergy_uT2ms', ...
    max(currentTseTrainRfEnergy_uT2_ms_by_group));
seq.setDefinition('MBRFValidatedFullTrainEnergy_uT2ms', ...
    max(validatedTseTrainRfEnergy_uT2_ms_by_group));
seq.setDefinition('WaveSinChannel', waveSinChannel);
seq.setDefinition('WaveCosChannel', waveCosChannel);
seq.setDefinition('WaveEnabled', double(isUseWave_sin || isUseWave_cos));
seq.setDefinition('WaveSinEnabled', double(isUseWave_sin));
seq.setDefinition('WaveCosEnabled', double(isUseWave_cos));
seq.setDefinition('WaveCombinedChannelValidationComplete', 1);
seq.setDefinition('WaveRepresentation', waveRepresentation);
seq.setDefinition('WaveAmplitude_mTm', max( ...
    waveTiming.peakAmplitude_mTm, cosineWaveTiming.peakAmplitude_mTm));
seq.setDefinition('WaveRequestedAmplitude_mTm', gwave_max);
seq.setDefinition('WaveActualAmplitude_mTm', max( ...
    waveTiming.peakAmplitude_mTm, cosineWaveTiming.peakAmplitude_mTm));
seq.setDefinition('WaveSlew_Tms', swave_max);
seq.setDefinition('WaveActualSlew_Tms', max( ...
    waveTiming.slewPeak_T_m_s, cosineWaveTiming.slewPeak_T_m_s));
seq.setDefinition('WaveFullAmplitudeRequiredSlew_Tms', ...
    2*pi*Ncycles/Tread*gwave_max*1e-3);
seq.setDefinition('WaveCycles', Ncycles);
seq.setDefinition('WaveSupportedCycles', [10 20]);
seq.setDefinition('WaveCycleParity', (-1)^Ncycles);
seq.setDefinition('WavePolarity', wavePolarity);
seq.setDefinition('WaveCarryY', waveAreaY);
seq.setDefinition('WaveReadoutCarryZ', waveAreaZ);
seq.setDefinition('WaveCarryZ', cosineCycleClosureError);
seq.setDefinition('WaveCosAreaToAdcCenter', cosineWaveTiming.areaToAdcCenter);
seq.setDefinition('WaveCosCenterClosureError', cosineCenterClosureError);
seq.setDefinition('WaveCosCycleClosureError', cosineCycleClosureError);
seq.setDefinition('WaveCosBridgeDuration', cosineBridgeDurationActual);
seq.setDefinition('WaveCosBridgeDesign', 'fixed-duration-lowPNS-three-point');
seq.setDefinition('WaveCosBridgeOverlapsCrusherCaipi', 0);
seq.setDefinition('CrusherCaipiDuration', tSp);
seq.setDefinition('RefReadoutHalfSpacing', tSpTotal);
seq.setDefinition('WaveActiveDur', max( ...
    waveTiming.Tactive, cosineWaveTiming.Tactive));
seq.setDefinition('WaveBlockEnvelopeDur', readoutTime);
seq.setDefinition('SplitReadoutWaveBlocks', double(isUseWave_cos));
seq.setDefinition('WaveStartsWithADC', double(waveStartsWithADC));
seq.setDefinition('WaveEndsWithADC', double(waveEndsWithADC));
seq.setDefinition('WaveCosStartsWithADC', double(cosineStartsWithADC));
seq.setDefinition('WaveCosEndsWithADC', double(cosineEndsWithADC));
seq.setDefinition('SystemType', sys_type);
seq.setDefinition('SequenceMaxGrad_mTm', sys_sequence.maxGrad/sys_sequence.gamma*1e3);
seq.setDefinition('SequenceMaxSlew_Tms', sys_sequence.maxSlew/sys_sequence.gamma);
seq.setDefinition('WaveMaxGrad_mTm', sys_wave.maxGrad/sys_wave.gamma*1e3);
seq.setDefinition('WaveMaxSlew_Tms', sys_wave.maxSlew/sys_wave.gamma);
seq.setDefinition('LowPNSMaxGrad_mTm', sys_lowPNS.maxGrad/sys_lowPNS.gamma*1e3);
seq.setDefinition('LowPNSMaxSlew_Tms', sys_lowPNS.maxSlew/sys_lowPNS.gamma);
seq.setDefinition('ExecuteSystemCheckCode', execute_system_check_code);
seq.setDefinition('WriteV141Format', double(write_v141_format));
seq.setDefinition('SimulationOnly', 0);
seq.setDefinition('OutputPath', out_path);
seq.setDefinition('RFDuration', tRef); % backward-compatible refocusing duration
seq.setDefinition('RFExcitationDuration', tEx);
seq.setDefinition('RFRefocusingDuration', tRef);
seq.setDefinition('RFExcitationTemporalCenter', rfExCenter);
seq.setDefinition('RFRefocusingTemporalCenter', rfRefCenter);
seq.setDefinition('GSspexAreaCorrection', GSspexAreaCorrection);
seq.setDefinition('VMBFirstEchoSliceMomentError', ...
    firstEchoMomentLeft-firstEchoMomentRight);
seq.setDefinition('VMBRepeatedEchoSliceMomentError', repeatedEchoMomentError);
seq.setDefinition('VMBExcitationToRefocusingCenter', ...
    excitationToRefocusingCenter);
seq.setDefinition('VMBRefocusingCenterToEcho', refocusingCenterToEcho);
seq.setDefinition('VMBEchoToNextRefocusingCenter', ...
    echoToNextRefocusingCenter);

if isUseWave_sin && isUseWave_cos
    tag_wave = '_wYZ';
elseif isUseWave_sin
    tag_wave = '_wY';
elseif isUseWave_cos
    tag_wave = '_wZ';
else
    tag_wave = '';
end
if isUseCaipi
    caipiFovDenominator = round(1/caipiShiftFractionFullFOV);
    isReciprocalIntegerShift = caipiFovDenominator >= 1 && ...
        abs(caipiShiftFractionFullFOV - 1/caipiFovDenominator) < 1e-12;
    if isReciprocalIntegerShift
        tag_caipi = ['_cF' num2str(caipiFovDenominator)];
    else
        caipiShiftTag = strrep(sprintf('%.8g', caipiShiftFractionFullFOV), '.', 'p');
        caipiShiftTag = strrep(caipiShiftTag, '-', 'm');
        caipiShiftTag = strrep(caipiShiftTag, '+', 'p');
        tag_caipi = ['_cS' caipiShiftTag];
    end
else
    tag_caipi = '_noC';
end
tag_cycles = ['_N' num2str(Ncycles)];
if isAcquireACS
    tag_ACS = ['_ACS' num2str(Nacs_y) 'to' num2str(nAcsPE)];
else
    tag_ACS = '';
end
if Ry == 3
    % Distinguish the complete 73-line grid from historical 66-line Ry3
    % exports while keeping the scanner-facing sequence filename compact.
    tag_acceleration = '_R3p73';
    tag_etl = ['_E' num2str(ETL)];
else
    tag_acceleration = ['_R' num2str(Ry)];
    tag_etl = ['_E' num2str(ETL)];
end
b1SequenceValueTag = sprintf('%.0f', double(MB.design_max_b1_uT));
tag_b1 = ['_B' b1SequenceValueTag];
espTag = strrep(sprintf('%.2f', TE1*1e3), '.', 'p');
effectiveTeTag = strrep(sprintf('%.0f', TEeff*1e3), '.', 'p');
seqFilename = ['tse_m' num2str(SMSFactor) 'b7v1' ...
    tag_b1 ...
    tag_caipi tag_wave tag_cycles '_' slOrientation '_' ...
    num2str(Nx) 'x' num2str(Ny) 'x' num2str(Nslices) ...
    '_F' num2str(round(fovScalar*1e3)) ...
    tag_acceleration tag_ACS ...
    '_os' num2str(os_factor) ...
    tag_etl ...
    '_ES' espTag ...
    '_TE' effectiveTeTag ...
    '_TR' num2str(round(TR*1e3)) ...
    '_' sys_type];
assert(numel([seqFilename '.seq']) <= 100, ...
    'Sequence filename exceeds 100 characters: %s.seq', seqFilename);
seq.setDefinition('Name', seqFilename);

if writeSequence
    v151Path = fullfile(out_path, 'v151');
    metadataPath = fullfile(out_path, 'metadata');
    if ~exist(v151Path, 'dir'), mkdir(v151Path); end
    if ~exist(metadataPath, 'dir'), mkdir(metadataPath); end
    if write_v141_format
        v141Path = fullfile(out_path, 'v141');
        if ~exist(v141Path, 'dir'), mkdir(v141Path); end
        seq.write_v141(fullfile(v141Path, [seqFilename '_v141.seq']));
    end
    sequenceFile = fullfile(v151Path, [seqFilename '.seq']);
    metadataFile = fullfile(metadataPath, [seqFilename '_PEorder.mat']);
    seq.write(sequenceFile);

    % Save PE ordering and candidate provenance for reconstruction and audit.
    save(metadataFile, ...
        'PEorder', 'PEindexOrderImg', 'LINorderImg', 'phaseAreas', ...
        'imageYIdx', 'imagePESteps', 'echoTrainPositionImg', ...
        'useSparseImageAdcSlots', 'imageAdcMask', ...
        'nImageEchoSlots', 'nImageNoAdcSlots', ...
        'Ry', 'Nechoes', 'ETL', 'turboFactor', 'nex', ...
        'isAcquireACS', 'PEorderACS', 'PEindexOrderACS', 'LINorderACS', 'phaseAreasACS', 'echoTrainPositionACS', ...
        'acsYIdx', 'acsYIdxBase', 'acsPESteps', ...
        'Nacs_y', 'Nacs_y_actual', 'nAcsPEBase', 'nAcsPE', 'nAcsExtraPE', 'nAcsDupPE', 'acsNex', ...
        'sliceThickness', 'sliceSpacing', 'sliceCenterSpacing', 'stackExtent', 'slicePositions', ...
        'iceNz', 'iceStackExtent', 'iceFov', 'iceMatrix', 'iceSliceGap', 'iceSlicePositions', ...
        'SMSFactor', 'NsliceGroups', 'smsGroups', 'groupOrder', 'groupBySlc', 'slcByGroup', ...
        'groupCenters', 'groupCenterFreqHz', 'constantGradientEquivalentGroupCenterFreqHz', ...
        'pairSeparation', 'isUseCaipi', 'caipiMode', 'caipiPolarity', 'caipiAdcDemodSign', ...
        'caipiShiftFractionFullFOV', 'caipiShiftFractionAliasedFOV', ...
        'caipiPhaseByLin', 'caipiAreaByLin', 'caipiMomentTolerance', ...
        'mbPulseMatFile', ...
        'candidateTag', 'candidateChecksum', 'candidateStatus', ...
        'validatedEchoSpacings', 'validatedEchoTrainLength', ...
        'minimumSupportedETL', 'isFullValidatedEchoTrain', ...
        'isValidatedEchoTrainPrefix', 'echoTrainValidationMode', ...
        'currentRFTrainSimulationScreened', ...
        'currentTseTrainRfEnergy_uT2_ms_by_group', ...
        'validatedTseTrainRfEnergy_uT2_ms_by_group', 'MB', 'rflip', ...
        'expectedTseTrainCount', 'nExcitationBlocks', ...
        'nFirstRefocusingBlocks', 'nRepeatRefocusingBlocks', ...
        'allowUnapprovedCandidateForResearchScan', ...
        'allowUnvalidatedEchoSpacing', ...
        'isGeneratorValidatedEchoSpacing', 'isUnvalidatedEchoSpacing', ...
        'echoSpacingValidationMode', ...
        'rfScreenTargetTE', 'rfScreenTargetEchoNumber', ...
        'tEx', 'tRef', 'rfExCenter', 'rfRefCenter', 'GSspexAreaCorrection', ...
        'firstEchoMomentLeft', 'firstEchoMomentRight', 'repeatedEchoMomentError', ...
        'TE1min', 'TE1minExcToRef', 'TE1minRefToEcho', ...
        'TE1minEchoToRef', ...
        'excitationToRefocusingCenter', 'refocusingCenterToEcho', ...
        'echoToNextRefocusingCenter', ...
        'isUseWave_sin', 'isUseWave_cos', 'readoutWaveChannels', ...
        'waveRepresentation', 'wavePolarity', 'gwave_max', 'swave_max', ...
        'Ncycles', 'waveAreaY', 'waveAreaZ', 'waveTiming', 'cosineWaveTiming', ...
        'cosinePreBridgeTiming', 'cosinePostBridgeTiming', ...
        'cosineCenterClosureError', 'cosineCycleClosureError', ...
        'cosineBridgeDurationActual', 'cosineBridgeDurationBudget', ...
        'sys_sequence', 'sys_wave', 'sys_lowPNS', 'TE1Requested', ...
        'TE1RequestedResolved', 'TE1RequestWasAutomatic', ...
        'TE1SelectionMode', 'TE1', 'TE1Baseline', ...
        'TEeffRequested', 'TEeff', 'effectiveEchoNumber', ...
        'TEeffAdjustment', 'TEeffWasAdapted', ...
        'TRRequested', 'TRActual', 'minimumTRWithoutFill', ...
        'TRWasAdapted', 'TRWasRaisedToMinimum', 'TRWasRasterRounded');
    fprintf('\nWrote imaging sequence: %s\n', sequenceFile);
    fprintf('Saved PE order and RF provenance: %s\n', metadataFile);
else
    fprintf(['\nSequence checks completed without writing output because ', ...
             'writeSequence=false.\n']);
    fprintf('After confirming all checks and plots, set writeSequence=true to write %s.\n', ...
        [seqFilename '.seq']);
end

%% Local helper functions

function [gPre, gRew] = makeBalancedBlippedCaipiGradients(GS5, GS7, GSspr, area, channel, duration, sys)
    %MAKEBALANCEDBLIPPEDCAIPIGRADIENTS Add balanced CAIPI crusher corrections.
    % Inputs are the base pre/post crusher events, crusher template, signed
    % CAIPI area, physical channel, common duration, and Pulseq system.
    % Outputs are the corrected pre-readout and post-readout gradients.
    %
    % Do not add a full-duration trapezoid on top of GS5/GS7. The slice
    % crusher bridge already uses most of the allowed slew during its
    % endpoint transitions, and a simultaneously ramping CAIPI trapezoid can
    % push the summed waveform slightly above the slew limit. Instead, place
    % a small triangular CAIPI correction only inside the flat GSspr portion
    % of the bridge. Then the added CAIPI slew occurs while the base crusher
    % is flat, so the base transition slew is not increased.
    if abs(area) < 1e-12
        gPre = GS5;
        gRew = GS7;
        return;
    end

    dt = sys.gradRasterTime;
    if GSspr.flatTime <= 2*dt
        error('CAIPI blip cannot fit: GSspr.flatTime %.6f us is too short.', GSspr.flatTime*1e6);
    end

    % Stay at least one raster point away from the base-gradient ramp
    % boundaries when possible. This avoids numerical/raster merging of the
    % CAIPI ramp with the slice-crusher ramp in mr.addGradients().
    guard = dt;
    if GSspr.flatTime - 2*guard < 2*dt
        guard = 0;
    end

    n0 = round((GSspr.riseTime + guard) / dt);
    n2 = round((GSspr.riseTime + GSspr.flatTime - guard) / dt);
    if n2 - n0 < 2
        error('CAIPI blip cannot fit: only %d raster intervals available inside GSspr flat top.', n2-n0);
    end
    nMid = floor((n0 + n2) / 2);
    nMid = max(n0 + 1, min(nMid, n2 - 1));

    t0 = n0 * dt;
    tMid = nMid * dt;
    t2 = n2 * dt;
    durationRaster = round(duration / dt) * dt;
    if t2 > durationRaster + dt/10
        error('Internal CAIPI timing error: correction end %.6f ms exceeds duration %.6f ms.', t2*1e3, durationRaster*1e3);
    end

    % Area of a triangular lobe is 0.5*peak*(t2-t0). The peak can be
    % positive or negative to realize the requested signed area.
    peak = 2 * area / (t2 - t0);
    corrSlew = max(abs(peak)/(tMid - t0), abs(peak)/(t2 - tMid));
    if corrSlew > sys.maxSlew * (1 + 1e-9)
        error(['CAIPI triangular correction exceeds maxSlew: %.6g > %.6g Hz/m/s. ', ...
               'Increase tSp/TE1 or reduce CAIPI shift.'], corrSlew, sys.maxSlew);
    end
    if abs(peak) > sys.maxGrad * (1 + 1e-9)
        error(['CAIPI triangular correction exceeds maxGrad: %.6g > %.6g Hz/m. ', ...
               'Increase tSp/TE1 or reduce CAIPI shift.'], abs(peak), sys.maxGrad);
    end

    times = [0, t0, tMid, t2, durationRaster];
    ampsPre = [0, 0, peak, 0, 0];
    ampsRew = [0, 0, -peak, 0, 0];

    gCaipiPre = mr.makeExtendedTrapezoid(channel, 'times', times, 'amplitudes', ampsPre, 'system', sys);
    gCaipiRew = mr.makeExtendedTrapezoid(channel, 'times', times, 'amplitudes', ampsRew, 'system', sys);

    gPre = mr.addGradients({GS5, gCaipiPre}, sys);
    gRew = mr.addGradients({GS7, gCaipiRew}, sys);
end

function idx = makeAcceleratedIndexList(N, R)
    %MAKEACCELERATEDINDEXLIST Select center-preserving accelerated PE indices.
    % Inputs are the full PE matrix size N and integer acceleration R. The
    % output is a row vector of one-based PE indices containing the center.
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
    %CENTERCROPINDEXLISTTOMULTIPLE Crop PE indices to complete echo trains.
    % Inputs are an ordered index list, its center index, and the required
    % block size. The output is the centered subset with divisible length.
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
    %MAKECENTERINDEXLIST Return a centered contiguous one-based index list.
    % Inputs are the full length N and requested central length Ncenter. The
    % output is the clipped central index vector.
    if Ncenter < 1 || Ncenter ~= round(Ncenter)
        error('Center block size must be a positive integer.');
    end
    Ncenter = min(Ncenter, N);
    startIdx = floor(N/2) + 1 - floor(Ncenter/2);
    startIdx = max(1, min(startIdx, N - Ncenter + 1));
    idx = startIdx:(startIdx + Ncenter - 1);
end


function [PEindexOrder, LINorder, PEorder, phaseAreas, adcSlotMask, pe_steps, nex, k0prescr, k0curr] = makeTsePETablesWithSparseAdcSlots(yIdx, yCenterIdx, ETL, TEeff, TE1, deltak)
    %MAKETSEPETABLESWITHSPARSEADCSLOTS Build a padded TSE PE table.
    % Inputs specify the complete PE index list, center, ETL, effective/echo
    % spacing, and delta-k. Outputs are ETL-by-train PE/label/area tables, a
    % logical ADC-slot mask, the acquired PE steps, train count, and the
    % requested/native center-echo positions. Safe placeholder slots use
    % zero PE area and the center index; callers must omit ADC in those slots.
    nPE = numel(yIdx);
    if nPE < ETL
        error('PE index count (%d) must contain at least one echo train of length %d.', nPE, ETL);
    end
    nex = ceil(nPE / ETL);
    nSlots = nex * ETL;

    pe_steps = yIdx(:).' - yCenterIdx;
    [~, iPEmin] = min(abs(pe_steps));
    k0curr = floor((iPEmin-1)/nex) + 1;
    k0prescr = max(round(TEeff/TE1), 1);
    k0prescr = min(k0prescr, ETL);

    nMissing = nSlots - nPE;
    % Append timing-only placeholders before reshaping. The subsequent echo
    % shift moves both PE steps and ADC flags together, so no acquired ky line
    % is duplicated or dropped and k-space center remains on the target echo.
    paddedAdcMask = false(1, nSlots);
    paddedAdcMask(1:nPE) = true;
    paddedSteps = zeros(1, nSlots);
    paddedSteps(1:nPE) = pe_steps;

    peStepOrder = reshape(paddedSteps, [nex, ETL])';
    adcSlotMask = reshape(paddedAdcMask, [nex, ETL])';
    echoShift = k0prescr-k0curr;
    peStepOrder = circshift(peStepOrder, [echoShift, 0]);
    adcSlotMask = circshift(adcSlotMask, [echoShift, 0]);

    PEindexOrder = peStepOrder + yCenterIdx;
    LINorder = PEindexOrder - 1;
    PEorder = peStepOrder;
    phaseAreas = PEorder * deltak;

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
    %MAKETSEPETABLESFROMINDICES Build physical and ICE TSE PE tables.
    % Inputs specify candidate PE indices, center, ETL, effective/echo
    % spacing, and delta-k. Outputs provide ordered indices and labels,
    % centered PE steps/areas, train count, and requested/native k0 echoes.
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


function [gWave, area, timing] = defineSineReadoutWaveTSE( ...
        channel, Tread, adcDelay, readoutBlockDur, sys, Ncycles, ...
        wavePolarity, gwave_max, swave_max, physical_slew_max, adc, ...
        waveInfoFlag, debugFlag)
    %DEFINESINEREADOUTWAVETSE Build an ADC-aligned full-block sine gradient.
    % Inputs specify the physical gradient channel, active ADC duration and
    % delay, complete readout-block duration, Pulseq system limits, sine
    % cycles/polarity/amplitude/slew, ADC event, and diagnostic flags. Outputs
    % are the zero-delay arbitrary gradient, its zeroth moment, and realized timing
    % and hardware metrics. Validation errors are raised for raster, ADC,
    % block-duration, gradient-amplitude, or slew inconsistencies.
    dt = sys.gradRasterTime;
    if ~ismember(wavePolarity, [-1 1])
        error('wavePolarity must be +1 or -1.');
    end
    [G0, w, TreadRaster] = designWaveAmplitude(Tread, sys, Ncycles, ...
        gwave_max, swave_max, physical_slew_max, 'sine-full-block-arbitrary', waveInfoFlag);

    nReadSamples = round(TreadRaster / dt);
    if abs(nReadSamples*dt - TreadRaster) > dt/10
        error('TreadRaster is not an integer number of gradient raster samples.');
    end

    nBlockSamples = round(readoutBlockDur / dt);
    readoutBlockDurRaster = nBlockSamples * dt;
    if abs(readoutBlockDurRaster - readoutBlockDur) > dt/10
        error(['readoutBlockDur %.6f ms is not on the gradient raster %.6f us. ', ...
               'Rasterized value would be %.6f ms.'], ...
            readoutBlockDur*1e3, dt*1e6, readoutBlockDurRaster*1e3);
    end

    nPadPre = round(adcDelay / dt);
    padPreDur = nPadPre * dt;
    if abs(padPreDur - adcDelay) > dt/10
        error(['adcDelay %.6f ms is not on the gradient raster %.6f us. ', ...
               'Rasterized value would be %.6f ms.'], ...
            adcDelay*1e3, dt*1e6, padPreDur*1e3);
    end

    nPadPost = nBlockSamples - nPadPre - nReadSamples;
    if nPadPost < 0
        error(['Full-block sine cannot fit in the readout block: ', ...
               'nBlockSamples=%d, nPadPre=%d, nReadSamples=%d.'], ...
            nBlockSamples, nPadPre, nReadSamples);
    end
    padPostDur = nPadPost * dt;

    if nReadSamples < 2
        error('At least two active samples are required for an endpoint-zero sine wave.');
    end

    % makeArbitraryGrad stores one amplitude per gradient-raster cell. Use an
    % integer-cycle active shape with explicit zero endpoints, then pad it to
    % the complete GR6 block so the event has zero delay and boundary values.
    activePhase = 2*pi*Ncycles * (0:(nReadSamples-1)) / (nReadSamples-1);
    activeWave = wavePolarity * G0 * sin(activePhase);
    activeWave(1) = 0;
    activeWave(end) = 0;
    fullBlockWave = [zeros(1, nPadPre), activeWave, zeros(1, nPadPost)];
    assert(numel(fullBlockWave) == nBlockSamples, ...
        'Internal full-block sine waveform length mismatch.');

    gWave = mr.makeArbitraryGrad(channel, fullBlockWave, 'system', sys, ...
        'first', 0, 'last', 0);
    gWave.delay = 0;
    area = gWave.area;

    slewPeak_Hz_m_s = max(abs(diff(fullBlockWave))) / dt;
    slewPeak_T_m_s = slewPeak_Hz_m_s / sys.gamma;
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
    timing.padPreDur = padPreDur;
    timing.padPostDur = padPostDur;
    timing.Twave = readoutBlockDurRaster;
    timing.Tactive = TreadRaster;
    timing.activeStartTime = padPreDur;
    timing.activeEndTime = padPreDur + TreadRaster;
    timing.G0 = G0;
    timing.omega = w;
    timing.freq = w/(2*pi);
    timing.polarity = wavePolarity;
    timing.peakAmplitude_mTm = max(abs(activeWave)) / sys.gamma * 1e3;
    timing.slewPeak_T_m_s = slewPeak_T_m_s;
    timing.requestedFullAmplitudeSlew_T_m_s = ...
        2*pi*Ncycles/TreadRaster*gwave_max*1e-3;
    timing.area = area;

    tol = dt/10;
    adcAcqDur = adc.numSamples * adc.dwell;
    assert(abs(TreadRaster - Tread) < tol, ...
        'Wave Tread %.6f ms does not match requested Tread %.6f ms.', TreadRaster*1e3, Tread*1e3);
    assert(abs(adcDelay - adc.delay) < tol, ...
        'Requested sine active-start delay and adc.delay do not match.');
    assert(abs(TreadRaster - adcAcqDur) < max(tol, sys.adcRasterTime), ...
        'Active sine duration and ADC acquisition duration do not match.');
    assert(abs(gWave.delay) < tol, ...
        'Full-block arbitrary sine must have zero event delay.');
    assert(abs(timing.activeStartTime - adc.delay) < tol, ...
        'Active sine start %.6f ms does not match ADC start %.6f ms.', ...
        timing.activeStartTime*1e3, adc.delay*1e3);
    assert(abs(timing.activeEndTime - (adc.delay + adcAcqDur)) < ...
            max(tol, sys.adcRasterTime), ...
        'Active sine end %.6f ms does not match ADC end %.6f ms.', ...
        timing.activeEndTime*1e3, (adc.delay + adcAcqDur)*1e3);
    assert(abs(mr.calcDuration(gWave) - readoutBlockDurRaster) < tol, ...
        'Full-block arbitrary sine must span the complete readout block.');

    if debugFlag
        fprintf('\n============= Full-block sine readout wave debug: TSE =============\n');
        fprintf('channel                          = %s\n', channel);
        fprintf('active / full-block samples       = %d / %d\n', nReadSamples, nBlockSamples);
        fprintf('pre/post zero-pad samples         = %d / %d\n', nPadPre, nPadPost);
        fprintf('Tread input/raster               = %.6f / %.6f ms\n', Tread*1e3, TreadRaster*1e3);
        fprintf('readout block duration            = %.6f ms\n', readoutBlockDurRaster*1e3);
        fprintf('gyWave.delay                     = %.6f ms\n', gWave.delay*1e3);
        fprintf('adc.delay                        = %.6f ms\n', adc.delay*1e3);
        fprintf('active start / active end         = %.6f / %.6f ms\n', ...
            timing.activeStartTime*1e3, timing.activeEndTime*1e3);
        fprintf('ADC acquisition duration          = %.6f ms\n', adcAcqDur*1e3);
        fprintf('active start - ADC start          = %.6f us\n', ...
            (timing.activeStartTime - adc.delay)*1e6);
        fprintf('active end - ADC end              = %.6f us\n', ...
            (timing.activeEndTime - (adc.delay + adcAcqDur))*1e6);
        fprintf('max amplitude                     = %.6f mT/m\n', timing.peakAmplitude_mTm);
        fprintf('peak slew from raster waveform    = %.6f T/m/s\n', slewPeak_T_m_s);
        fprintf('wave object duration incl delay   = %.6f ms\n', mr.calcDuration(gWave)*1e3);
        fprintf('wave area from object             = %.9g 1/m\n', area);
        fprintf('===================================================================\n\n');
    end
end

function [gWave, area, timing] = defineCosineReadoutWaveTSE( ...
        channel, Tread, adcDelay, readoutBlockDur, sys, Ncycles, ...
        wavePolarity, gwave_max, swave_max, physical_slew_max, adc, ...
        waveInfoFlag, debugFlag)
    %DEFINECOSINEREADOUTWAVETSE Build an ADC-aligned full-block cosine gradient.
    % Inputs specify the gradient channel, ADC/readout timing, Pulseq system,
    % Wave cycles/polarity/limits, ADC event, and diagnostic flags. Outputs
    % are the nonzero-endpoint arbitrary gradient, its readout-block area,
    % and realized timing/moment metrics. The waveform holds +G0 during ADC
    % guard intervals and uses an endpoint-inclusive integer-cycle cosine
    % during acquisition. Validation errors are raised for raster, hardware,
    % ADC-alignment, or moment-index inconsistencies.
    dt = sys.gradRasterTime;
    if ~ismember(wavePolarity, [-1 1])
        error('wavePolarity must be +1 or -1.');
    end

    [G0, w, TreadRaster] = designWaveAmplitude(Tread, sys, Ncycles, ...
        gwave_max, swave_max, physical_slew_max, ...
        'cosine-full-block-arbitrary', waveInfoFlag);
    Gedge = wavePolarity * G0;

    nReadSamples = round(TreadRaster / dt);
    readoutBlockDurRaster = round(readoutBlockDur / dt) * dt;
    nBlockSamples = round(readoutBlockDurRaster / dt);
    nPadPre = round(adcDelay / dt);
    nPadPost = nBlockSamples - nPadPre - nReadSamples;
    if nPadPost < 0
        error(['Full-block cosine cannot fit in the readout block: ', ...
               'nBlockSamples=%d, nPadPre=%d, nReadSamples=%d.'], ...
            nBlockSamples, nPadPre, nReadSamples);
    end
    if nReadSamples < 2
        error('At least two active samples are required for the cosine wave.');
    end

    activePhase = 2*pi*Ncycles * (0:(nReadSamples-1)) / (nReadSamples-1);
    activeWave = wavePolarity * G0 * cos(activePhase);
    activeWave(1) = Gedge;
    activeWave(end) = Gedge;
    fullBlockWave = [Gedge*ones(1, nPadPre), activeWave, ...
        Gedge*ones(1, nPadPost)];
    assert(numel(fullBlockWave) == nBlockSamples, ...
        'Internal full-block cosine waveform length mismatch.');

    gWave = mr.makeArbitraryGrad(channel, fullBlockWave, 'system', sys, ...
        'first', Gedge, 'last', Gedge);
    gWave.delay = 0;
    area = gWave.area;

    adcCenterTime = adc.delay + 0.5*adc.numSamples*adc.dwell;
    nCellsToAdcCenter = round(adcCenterTime / dt);
    assert(abs(nCellsToAdcCenter*dt - adcCenterTime) < dt/10, ...
        'ADC center %.6f ms is not on the gradient raster.', adcCenterTime*1e3);
    assert(nCellsToAdcCenter >= 1 && nCellsToAdcCenter <= nBlockSamples, ...
        'ADC-center gradient-cell index is outside the cosine waveform.');
    areaToAdcCenter = sum(fullBlockWave(1:nCellsToAdcCenter))*dt;

    slewPeak_Hz_m_s = max(abs(diff(fullBlockWave))) / dt;
    slewPeak_T_m_s = slewPeak_Hz_m_s / sys.gamma;
    waveSlewLimit = min(swave_max, physical_slew_max);
    assertArbitraryWaveformWithinSystem( ...
        fullBlockWave, Gedge, Gedge, sys, 'active cosine wave');
    assert(slewPeak_T_m_s <= waveSlewLimit*(1+1e-9), ...
        'Active cosine slew %.9g T/m/s exceeds requested limit %.9g T/m/s.', ...
        slewPeak_T_m_s, waveSlewLimit);
    adcAcqDur = adc.numSamples * adc.dwell;

    timing = struct;
    timing.nReadSamples = nReadSamples;
    timing.nBlockSamples = nBlockSamples;
    timing.nPadPre = nPadPre;
    timing.nPadPost = nPadPost;
    timing.Twave = readoutBlockDurRaster;
    timing.Tactive = TreadRaster;
    timing.activeStartTime = nPadPre*dt;
    timing.activeEndTime = nPadPre*dt + TreadRaster;
    timing.adcCenterTime = adcCenterTime;
    timing.nCellsToAdcCenter = nCellsToAdcCenter;
    timing.edgeAmplitude = Gedge;
    timing.peakAmplitude_mTm = max(abs(activeWave)) / sys.gamma * 1e3;
    timing.slewPeak_T_m_s = slewPeak_T_m_s;
    timing.requestedFullAmplitudeSlew_T_m_s = ...
        2*pi*Ncycles/TreadRaster*gwave_max*1e-3;
    timing.area = area;
    timing.areaToAdcCenter = areaToAdcCenter;
    timing.areaAfterAdcCenter = area - areaToAdcCenter;

    tol = dt/10;
    assert(abs(TreadRaster - adcAcqDur) < max(tol, sys.adcRasterTime), ...
        'Active cosine duration and ADC acquisition duration do not match.');
    assert(abs(timing.activeStartTime - adc.delay) < tol, ...
        'Active cosine start does not match the ADC start.');
    assert(abs(timing.activeEndTime - (adc.delay + adcAcqDur)) < ...
            max(tol, sys.adcRasterTime), ...
        'Active cosine end does not match the ADC end.');
    assert(abs(mr.calcDuration(gWave) - readoutBlockDurRaster) < tol, ...
        'Full-block arbitrary cosine must span the complete readout block.');

    if debugFlag
        fprintf('\n=========== Full-block cosine readout wave debug: TSE ===========\n');
        fprintf('channel                              = %s\n', channel);
        fprintf('active / full-block samples           = %d / %d\n', ...
            nReadSamples, nBlockSamples);
        fprintf('pre/post +G0 guard samples            = %d / %d\n', nPadPre, nPadPost);
        fprintf('edge / peak amplitude                 = %.6f / %.6f mT/m\n', ...
            Gedge/sys.gamma*1e3, timing.peakAmplitude_mTm);
        fprintf('peak slew from raster waveform        = %.6f T/m/s\n', slewPeak_T_m_s);
        fprintf('readout area                          = %.9g 1/m\n', area);
        fprintf('area through ADC center               = %.9g 1/m\n', areaToAdcCenter);
        fprintf('area after ADC center                 = %.9g 1/m\n', ...
            timing.areaAfterAdcCenter);
        fprintf('=================================================================\n\n');
    end
end

function [gBridge, report] = makeFixedDurationThreePointBridge( ...
        channel, Gstart, Gend, areaTarget, duration, sys)
    %MAKEFIXEDDURATIONTHREEPOINTBRIDGE Design a low-PNS endpoint bridge.
    % The accepted bridge duration is preserved for both supported cycle
    % cases; the midpoint location is chosen to minimize peak slew.
    dt = sys.gradRasterTime;
    nIntervals = round(duration/dt);
    Ttotal = nIntervals*dt;
    assert(nIntervals >= 2 && abs(Ttotal-duration) <= dt/10, ...
        'Cosine bridge duration must contain at least two raster intervals.');
    best = [];
    bestSlew = inf;
    for iMid = 1:(nIntervals-1)
        Tmid = iMid*dt;
        Gmid = (2*areaTarget - Tmid*Gstart - ...
            (Ttotal-Tmid)*Gend) / Ttotal;
        slew1 = abs(Gmid-Gstart) / Tmid;
        slew2 = abs(Gend-Gmid) / (Ttotal-Tmid);
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
               'the accepted %.3f ms duration for area %.9g 1/m.'], ...
              duration*1e3, areaTarget);
    end

    gBridge = mr.makeExtendedTrapezoid(channel, ...
        'times', [0 best.Tmid best.Ttotal], ...
        'amplitudes', [Gstart best.Gmid Gend], ...
        'system', sys);
    areaError = gBridge.area - areaTarget;
    assert(abs(areaError) <= 1e-6, ...
        'Three-point bridge area error %.9g 1/m exceeds tolerance.', areaError);

    report = best;
    report.startAmplitude = Gstart;
    report.endAmplitude = Gend;
    report.areaTarget = areaTarget;
    report.areaActual = gBridge.area;
    report.areaError = areaError;
    report.peakAmplitude_mTm = max(abs([Gstart best.Gmid Gend])) / ...
        sys.gamma * 1e3;
    report.peakSlew_T_m_s = best.peakSlew / sys.gamma;
end

function [G0_pulseq, w, TreadRaster] = designWaveAmplitude(Tread, sys, Ncycles, gwave_max, swave_max, physical_slew_max, waveName, waveInfoFlag)
    %DESIGNWAVEAMPLITUDE Select a rasterized, hardware-limited sine amplitude.
    % Inputs specify the active duration, Pulseq system, integer cycle count,
    % requested gradient/slew limits, physical slew ceiling, diagnostic name,
    % and print flag. Outputs are amplitude in Pulseq units, angular
    % frequency, and rasterized duration. Invalid timing or limits raise an
    % error; the returned amplitude respects the requested rasterized slew.
    if ~ismember(Ncycles, [10 20])
        error('Ncycles must be one of the supported protocol options: 10 or 20.');
    end
    if Tread <= 0
        error('Tread must be positive.');
    end
    if gwave_max < 0
        error('gwave_max must be non-negative.');
    end
    if swave_max <= 0 || physical_slew_max <= 0
        error('swave_max and physical_slew_max must be positive.');
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
