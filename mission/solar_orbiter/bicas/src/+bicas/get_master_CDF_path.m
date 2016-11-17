% get_master_CDF_path Get the path to the master CDF file for a given dataset ID and version.
%
% [master_CDF_path, master_filename] = get_master_CDF_path(dataset_ID, skeleton_version_str)
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-28
%
% This function _decides_ the path and filename of the master CDF file that BICAS should use.
% NOTE: The S/W descriptor needs the filename only. The execution of a S/W mode needs the entire path. Hence the return
% values.
%
% ARGUMENTS AND RETURN VALUES
% skeleton_version_str : Two-character two-digit string specifying version of the master CDF file.
%
function [master_CDF_path, master_filename] = get_master_CDF_path(dataset_ID, skeleton_version_str)
% QUESTION: Should define the master CDFs name (construct it with code), or look it up among the constants?
% PROPOSAL: Move function into constants.
% PROPOSAL: Retrieve value from constants somehow (or rather, verify with constants).
%
    global CONSTANTS
    
    master_filename = [dataset_ID, '_V', skeleton_version_str, '.cdf'];
    master_CDF_path = fullfile(CONSTANTS.C.BICAS_ROOT_PATH, CONSTANTS.C.MASTER_CDFS_RELATIVE_DIR, master_filename);
end
