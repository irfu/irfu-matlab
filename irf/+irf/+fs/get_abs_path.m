%
% Convert (relative/absolute) path to an absolute path (not necessarily
% canonical without symlinks).
%
%
% NOTES
% =====
% NOTE: MATLAB does indeed seem to NOT have a function for getting the absolute
%       path!
% NOTE: Also converts "~" to the home directory.
% NOTE: The resulting path will NOT end with slash/backslash unless it is the
%       system root directory on Linux ("/").
% NOTE: Only works with filesep = "/".
%
%
% ARGUMENT
% ========
% path
%       NOTE: Only path to object which parent directory exists. The object
%       itself does not have to exist.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-06-09.
%
function path = get_abs_path(path)
% PROPOSAL: Uses Linux's "readlink -f".
%   CON: Platform-dependent.
%
% PROPOSAL: Use "what".
%   NOTE: Does not work on files.
% PROPOSAL: Use "fileparts" to make it work for files in existing directories.
%
% NOTE: Different use cases, operations.
%   PROPOSAL: Convert (absolute or relative) to absolute path: add current directory to path.
%       NOTE: Does not require existing object.
%   PROPOSAL: Find canonical path: Replace symlinks with non-links.
%       NOTE: Requires existing object (or at least up to last link).
%   PROPOSAL: Rationalize away .. and .
%       NOTE: Does not require existing object.
%   PROPOSAL: Replace ~ with home dir.
%       NOTE: Does not require existing object.
%
% PROPOSAL: Separate function for replacing "~" with home dir. Cf python.
%
% PROPOSAL: Test code.

try
  if strcmp(path, '~') || (length(path)>=2 && strcmp(path(1:2), ['~', filesep]))
    % Not entirely sure if may have trailing slash.
    homeDir = getenv('HOME');
    path = fullfile(homeDir, path(2:end));
  end

  if ~exist(path, 'dir')
    % NOTE: If ends with slash, then everything is assigned to dirPath!
    [dirPath, basename, suffix] = fileparts(path);

    % Uses MATLAB trick to convert path to absolute path.
    %dirPath = cd(cd(dirPath));
    whatInfo = what(dirPath);
    dirPath = whatInfo.path;
    path = fullfile(dirPath, [basename, suffix]);
  else
    %path = cd(cd(path));
    whatInfo = what(path);
    path = whatInfo.path;
  end

  % Remove trailing slashes, except for system root.
  path = regexp(path, '^(/.*[^/]|/)', 'match');
  path = path{1};

catch Exc
  error('get_abs_path:FailedToConvertPath', ...
    ['Failed to convert path "%s" to an absolute path.\n', ...
    'Exc.message = "%s"'], ...
    path, Exc.message)
end

end
