%
% Create empty file. This is useful for debugging purposes sometimes, e.g.
% if batch code selects the correct files to be created.
%
%
% ARGUMENTS
% =========
% path
%       Path to file.
%       NOTE: Parent directory has to pre-exist.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function create_empty_file(path)
    % PROPOSAL: Error if fails to create file.
    % PROPOSAL: Policy arguments
    %   (1) Whether to permit path not available:
    %       Permit overwrite
    %       Assert no overwrite
    %   (2) What happens when createing, writing file:
    %       Assert write success (does not assert no overwrite)
    %       Permit write failure (does not assert no overwrite), return ~error code/boolean

    fclose(fopen(path, 'w'));
end
