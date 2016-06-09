% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
% Convert (relative/absolute) path to an absolute (canonical) path.
%
% NOTE: Will also convert ~ (home directory).
% (MATLAB does indeed seem to NOT have a function for doing this(!).)
%
function path = get_abs_path(path)
%
% PROPOSAL: Rethrow exception via errorp with amended message somehow?

global ERROR_CODES

try
    % Uses MATLAB trick to convert path to absolute path.
    path = cd(cd(path));
catch e
    errorp(ERROR_CODES.PATH_NOT_FOUND, 'Failed to convert path "%s" to absolute path.\nException message: "%s"', path, e.message)
end

end
