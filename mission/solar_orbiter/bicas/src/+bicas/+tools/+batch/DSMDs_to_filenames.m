%
% Convert multiple DSMDs to filenames (not paths).
%
%
% ARGUMENTS
% =========
% DsmdArray
%       Array of solo.adm.DSMD.
%
%
% RETURN VALUES
% =============
% filenamesCa
%       Cell array of same size as DsmdArray. Filenames (not paths) in
%       DsmdArray.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function filenamesCa = DSMDs_to_filenames(DsmdArray)
% PROPOSAL: Automatic test code.
% PROPOSAL: Move to DSMD class.

assert(isa(DsmdArray, 'solo.adm.DSMD'))

filenamesCa = cellfun(...
  @irf.fs.get_name, {DsmdArray.path}, ...
  'UniformOutput', false);
end
