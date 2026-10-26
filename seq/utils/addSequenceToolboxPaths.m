function pulseqMatlabPath = addSequenceToolboxPaths(pulseq_path, safe_pns_prediction_path)
%ADDSEQUENCETOOLBOXPATHS Add Pulseq and optional safe_pns_prediction to the MATLAB path.
%
%   pulseqMatlabPath = addSequenceToolboxPaths(pulseq_path, safe_pns_prediction_path)
%
% Inputs:
%   pulseq_path               Pulseq repository root or its matlab/ folder.
%   safe_pns_prediction_path  safe_pns_prediction folder (for example the
%                             external/safe_pns_prediction submodule), or ''
%                             when the PNS/CNS checks are not run. seq.calcPNS
%                             needs it when execute_system_check_code is true.
%
% Output:
%   pulseqMatlabPath          Folder that directly contains the +mr package.
%
% Side effects:
%   Prepends the Pulseq and safe_pns_prediction folders to the MATLAB path.
%
% Errors:
%   Raises an error when +mr is found neither in pulseq_path/matlab nor in
%   pulseq_path, or when a non-empty safe_pns_prediction_path is not a folder.

    candidates = {fullfile(pulseq_path, 'matlab'), pulseq_path};
    hasMrPackage = cellfun(@(folder) isfolder(fullfile(folder, '+mr')), candidates);
    if ~any(hasMrPackage)
        error('Pulseq +mr package not found in %s or its matlab/ folder.', pulseq_path);
    end
    pulseqMatlabPath = candidates{find(hasMrPackage, 1)};
    addpath(pulseqMatlabPath);

    if ~isempty(safe_pns_prediction_path)
        if ~isfolder(safe_pns_prediction_path)
            error('safe_pns_prediction folder not found: %s', safe_pns_prediction_path);
        end
        addpath(safe_pns_prediction_path);
    end
end
