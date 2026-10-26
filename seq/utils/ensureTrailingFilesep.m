function pathValue = ensureTrailingFilesep(pathValue)
%ENSURETRAILINGFILESEP Append the platform file separator to a folder path.
%
%   pathValue = ensureTrailingFilesep(pathValue)
%
% Input:
%   pathValue  char folder path; empty input is returned unchanged.
%
% Output:
%   pathValue  The path ending in filesep, so that [pathValue name] forms a
%              file path inside the folder.
    if isempty(pathValue)
        return;
    end
    if pathValue(end) ~= filesep
        pathValue = [pathValue filesep];
    end
end
