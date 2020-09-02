% Get the filename for the master CDF file for a given dataset ID and version.
%
% This function _decides_ the filename of the master CDF file that BICAS should use.
%
%
% ARGUMENTS
% =========
% datasetId          : DATASET_ID
% skeletonVersionStr : Two-character two-digit string specifying version of the master CDF file.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-07-28
%
function masterCdfFilename = get_master_CDF_filename(datasetId, skeletonVersionStr)
% PROPOSAL: Move function into ~constants.
% PROPOSAL: Standard assertion for skeletonVersionStr.
% PROPOSAL: Update BICAS so that both lowercase and uppercase .cdf is accepted.
%           ==> Has to be done in other code which finds set of matching files.

    bicas.assert_DATASET_ID(datasetId)
    bicas.assert_skeleton_version(skeletonVersionStr)
    
    masterCdfFilename = [datasetId, '_V', skeletonVersionStr, '.cdf'];
end
