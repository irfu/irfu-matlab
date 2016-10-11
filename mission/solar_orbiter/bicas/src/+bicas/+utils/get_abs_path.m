% path = get_abs_path(path)   Get absolute directory path.
%
% Convert (relative/absolute) path to an absolute (canonical) path.
% (MATLAB does indeed seem to NOT have a function for doing this(!).)
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
% NOTE: Will also convert ~ (home directory).
% NOTE: Only works for valid paths (paths which exist in the file system).
% NOTE: Only works for directory paths (not regular files).
% NOTE: Will return string that does not end with slash, except for the ~Linux system root "/" (Windows "C:\"?).
%
function path = get_abs_path(path)

try
    % Uses MATLAB trick to convert path to absolute path.
    path = cd(cd(path));
catch e
    error('get_abs_path:PathNotFound', 'Failed to convert path "%s" to absolute path.\nException message: "%s"', path, e.message)
end

end
