%
% Write empty file, including parent directories recursively if not already
% pre-existing. Also return path to newly create file.
%
% This function is useful for creating automated tests and debugging.
%
%
% ARGUMENTS
% =========
% pathPartsCa
%       Path to file as separate parts (strings) which are merged into one path
%       with fullfile().
%
%
% RETURN VALUE
% ============
% filePath
%       Path to file just created, based on argument.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function filePath = write_empty_file(pathPartsCa)
% IMPLEMENTATION NOTE: Having one argument instead of using varargin in order to
%                      make it possible to add optional arguments in the future.

% PROPOSAL: Error if fails to create file.
% PROPOSAL: Policy arguments
%   (1) Whether to permit path not available:
%       Permit overwrite
%       Assert no overwrite
%   (2) What happens when createing, writing file:
%       Assert write success (does not assert no overwrite)
%       Permit write failure (does not assert no overwrite), return ~error code/boolean

assert(iscell(pathPartsCa), 'Argument pathPartsCa is not a cell array.')
filePath = fullfile(pathPartsCa{:});

[parentDir, ~, ~] = fileparts(filePath);
if ~exist(parentDir, 'dir')
  mkdir(parentDir)
end

fclose(fopen(filePath, 'w'));
end
