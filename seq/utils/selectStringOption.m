function value = selectStringOption(varName, promptText, options, defaultValue)
%SELECTSTRINGOPTION Ask for one entry of a list of names in the Command Window.
%
%   value = selectStringOption(varName, promptText, options, defaultValue)
%
% Same prompt as the Wave-MPRAGE and Wave-GRE generators. The options are
% printed as a numbered menu; the reply may be a menu number or an option
% name, and pressing Enter keeps the default.
%
% Inputs:
%   varName       Name of the base-workspace variable that receives the
%                 answer, e.g. 'sys_type'. A valid char/string value already
%                 held by that variable becomes the default; an invalid one
%                 is ignored with a warning.
%   promptText    Heading printed above the menu.
%   options       Cell array of allowed char values.
%   defaultValue  Default used when the workspace holds no valid value
%                 (optional; default options{1}).
%
% Output:
%   value         Selected option (char).
%
% Side effects:
%   Prompts in the Command Window and assigns value to varName in the base
%   workspace.
%
% Errors:
%   Raises an error when the reply is neither a valid menu number nor one of
%   the option names.

    if nargin < 4 || isempty(defaultValue)
        defaultValue = options{1};
    end

    value = defaultValue;
    if evalin('base', sprintf('exist(''%s'', ''var'')', varName))
        workspaceValue = evalin('base', varName);
        if isstring(workspaceValue), workspaceValue = char(workspaceValue); end
        if ischar(workspaceValue) && any(strcmp(workspaceValue, options))
            value = workspaceValue;
        else
            warning('Ignoring invalid workspace value for %s.', varName);
        end
    end

    defaultIdx = find(strcmp(value, options), 1);
    if isempty(defaultIdx), defaultIdx = 1; end

    fprintf('\n%s:\n', promptText);
    for ii = 1:numel(options)
        fprintf('  %2d) %s\n', ii, options{ii});
    end
    reply = input(sprintf('Select 1-%d [default %d: %s]: ', numel(options), defaultIdx, options{defaultIdx}), 's');

    if ~isempty(strtrim(reply))
        numericChoice = str2double(reply);
        if ~isnan(numericChoice) && numericChoice == round(numericChoice) && numericChoice >= 1 && numericChoice <= numel(options)
            value = options{numericChoice};
        else
            trimmedReply = strtrim(reply);
            idx = find(strcmp(trimmedReply, options), 1);
            if isempty(idx)
                error('Invalid %s selection: %s', varName, reply);
            end
            value = options{idx};
        end
    else
        value = options{defaultIdx};
    end

    assignin('base', varName, value);
end
