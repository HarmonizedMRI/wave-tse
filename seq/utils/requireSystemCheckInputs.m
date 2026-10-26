function ascName = requireSystemCheckInputs(sys_type, expectedAscName, system_asc_file, safe_pns_prediction_path)
%REQUIRESYSTEMCHECKINPUTS Validate the inputs of the PNS/CNS and forbidden-frequency checks.
%
%   ascName = requireSystemCheckInputs(sys_type, expectedAscName, ...
%       system_asc_file, safe_pns_prediction_path)
%
% Inputs:
%   sys_type                  Scanner key of the generator, e.g. 'skyra'.
%   expectedAscName           Gradient-system model file name that the
%                             generator's table binds to sys_type.
%   system_asc_file           system_asc_file path setting.
%   safe_pns_prediction_path  safe_pns_prediction_path path setting.
%
% Output:
%   ascName                   system_asc_file, validated for sys_type.
%
% Errors:
%   Raises an error when system_asc_file is missing, when its file name
%   differs from expectedAscName (PNS/CNS would then be evaluated with
%   another scanner's gradient model), or when safe_pns_prediction_path is
%   empty, because seq.calcPNS calls safe_gwf_to_pns and safe_plot.

    if isempty(system_asc_file) || ~isfile(system_asc_file)
        error(['System checks for sys_type=%s need system_asc_file (%s). ', ...
               'Rerun the generator and set it in the path settings.'], ...
            sys_type, expectedAscName);
    end
    [~, ascBase, ascExt] = fileparts(system_asc_file);
    if ~strcmp([ascBase ascExt], expectedAscName)
        error('system_asc_file %s does not match sys_type=%s, which expects %s.', ...
            system_asc_file, sys_type, expectedAscName);
    end
    if isempty(safe_pns_prediction_path)
        error(['System checks need safe_pns_prediction_path for seq.calcPNS, ', ...
               'for example the external/safe_pns_prediction submodule.']);
    end
    ascName = system_asc_file;
end
