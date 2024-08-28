%
% Extract file/directory/object name from path.
%
%
% ARGUMENTS
% =========
% path
%       NOTE: Will interpret anything after last filesep as filename, including
%       ".." or ".", or empty string for paths that end with filesep.
%       NOTE: Path is not required to exist, since it is only used for a string
%       operation.
%       IMPORTANT NOTE: If path ends with slash, then the name is an empty
%       string.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created circa 2020-04-29.
%
function [name, parentPath] = get_name(path)
% PROPOSAL: Generalize to accepting cell array of paths.
%   PRO: More convenient.
%   CON: Easy to use cellfun.
%       Ex: cellfun(@irf.fs.get_name, {DsmdArray1.path}', 'UniformOutput', false)
%   CON: Provides no extra value, e.g. speed (except convenience).
%   CON: Output would have to change to cell arrays for cell array input
%        (special/separate case).
%
% PROPOSAL: Remove trailing slash.

[parentPath, basename, suffix] = fileparts(path);
name = [basename, suffix];
end
