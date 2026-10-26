function pathValue = normalizeUserPath(pathValue)
%NORMALIZEUSERPATH Trim a user-entered path and expand a leading '~'.
%
%   pathValue = normalizeUserPath(pathValue)
%
% Input:
%   pathValue  char or string path; empty input is returned unchanged.
%
% Output:
%   pathValue  char path without surrounding whitespace, with '~' or '~/...'
%              replaced by the home directory (HOME, or USERPROFILE on
%              Windows).
    if isempty(pathValue)
        return;
    end
    if isstring(pathValue), pathValue = char(pathValue); end
    pathValue = strtrim(pathValue);
    if startsWith(pathValue, ['~' filesep]) || strcmp(pathValue, '~')
        homeDir = getenv('HOME');
        if isempty(homeDir) && ispc
            homeDir = getenv('USERPROFILE');
        end
        if strcmp(pathValue, '~')
            pathValue = homeDir;
        else
            pathValue = fullfile(homeDir, pathValue(3:end));
        end
    end
end
