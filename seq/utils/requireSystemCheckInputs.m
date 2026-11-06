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
%   ascName                   Gradient model of sys_type: system_asc_file,
%                             or expectedAscName in the folder of
%                             system_asc_file when the saved file belongs to
%                             another system.
%
% Side effects:
%   Prints a note when it substitutes the model of the selected system.
%
% Errors:
%   Raises an error when system_asc_file is missing, when it names another
%   system's model and expectedAscName is not in its folder (PNS/CNS would
%   otherwise be evaluated with another scanner's gradient model), or when
%   safe_pns_prediction_path is empty, because seq.calcPNS calls
%   safe_gwf_to_pns and safe_plot.

    if isempty(system_asc_file) || ~isfile(system_asc_file)
        error(['System checks for sys_type=%s need system_asc_file (%s). ', ...
               'Rerun the generator and set it in the path settings.'], ...
            sys_type, expectedAscName);
    end
    [ascFolder, ascBase, ascExt] = fileparts(system_asc_file);
    if ~strcmp([ascBase ascExt], expectedAscName)
        % The system is selected on every run, but one .asc file is saved;
        % use the model of the selected system from the same folder.
        candidate = fullfile(ascFolder, expectedAscName);
        if ~isfile(candidate)
            error(['system_asc_file %s does not match sys_type=%s, which expects %s, ', ...
                   'and that file is not in the same folder. Select the system of the ', ...
                   'saved file at the prompt, or set system_asc_file to %s in the path settings.'], ...
                system_asc_file, sys_type, expectedAscName, expectedAscName);
        end
        fprintf('system_asc_file %s belongs to another system; using %s for sys_type=%s.\n', ...
            system_asc_file, candidate, sys_type);
        system_asc_file = candidate;
    end
    if isempty(safe_pns_prediction_path)
        error(['System checks need safe_pns_prediction_path for seq.calcPNS, ', ...
               'for example the external/safe_pns_prediction submodule.']);
    end
    ascName = system_asc_file;
end
