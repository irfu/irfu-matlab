%
% Get the filename for the master CDF file for a given dataset ID and version.
%
% This function DECIDES the filename of the master CDF file that BICAS should
% use.
%
%
% ARGUMENTS
% =========
% dsi
%       DSI
% skeletonVersionStr
%       Two-character two-digit string specifying version of the master CDF
%       file.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-07-28
%
function masterCdfFilename = get_master_CDF_filename(dsi, skeletonVersionStr)

bicas.assert_BICAS_DSI(dsi)
bicas.assert_skeleton_version(skeletonVersionStr)

masterCdfFilename = [dsi, '_V', skeletonVersionStr, '.cdf'];
end
