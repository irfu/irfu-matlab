%
% Iterate over files. Try to parse filenames as if they were datasets using
% standard filenaming conventions
% (solo.adm.parse_dataset_filename()). Saves the result plus path.
% Ignores unparsable filenames without error.
%
% NOTE: Function accepts FILENAMES and PATHS, not just paths (as the argument
% name implies), and not just filenames (like
% solo.adm.parse_dataset_filename()). It is therefore not entirely
% wrong that the function name mentions filenames, not paths.
%
%
% ARGUMENT
% ========
% filePathCa
%       Column cell array of paths to files. Can be both datasets and not.
%
%
% RETURN VALUE
% ============
% fiCa
%       Nx1 cell array of structs. FI = File Info
%       {iDataset} : Struct. Fields from
%           solo.adm.parse_dataset_filename() plus extra field
%           below:
%               .path : Path in filePathList{iFile}.
% bIsDatasetArray
%       Logical column array. Same size as argument. True iff the corresponding
%       input path was interpreted as a dataset (was translated into a DSMD).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-25.
%
function [fiCa, bIsDatasetArray] = parse_dataset_filename_many(filePathCa)
% PROPOSAL: Change name
%   PROPOSAL: parse_dataset_filenames_many  ("FILENAMES" in plural)
%   PROPOSAL: parse_dataset_filename_many_paths
%   PROPOSAL: parse_dataset_paths_many
%   PROPOSAL: Change to accepting many filenames, not paths.
% PROPOSAL: Drop "_many" from function name and merge with
%        solo.adm.parse_dataset_filename().
%   PRO: Works for singular case with only basename (not path) too.
%   CON: Function adds functionality of filtering/ignoring non-datasets.
%
% PROPOSAL: Policy argument for reacting to unparsable filenames.
%   TODO-DEC: Different policy for suffix .cdf and not?
% PROPOSAL: Assertion for parsable *.cdf filenames. Ignore filenames without
%           suffix ".cdf".
%
% NOTE: Has no separate test code. Is indirectly tested by
%       solo.adm.paths_to_DSMD_array___UTEST.

assert(iscell(filePathCa),   'filePathCa is not a cell array.')
assert(iscolumn(filePathCa), 'filePathCa is not a column array.')

fiCa            = cell(0, 1);
bIsDatasetArray = false(numel(filePathCa), 1);
for iFile = 1:numel(filePathCa)

  filename = irf.fs.get_name(filePathCa{iFile});

  R = solo.adm.parse_dataset_filename(filename);

  if ~isempty(R)
    % CASE: File can be identified as a dataset.

    % NOTE: Ignores .versionStr
    Fi      = R;
    Fi.path = filePathCa{iFile};

    fiCa{           end+1, 1} = Fi;
    bIsDatasetArray(iFile, 1) = true;
  else
    % CASE: Can not identify as dataset filename.

    % Do nothing. (Silently ignore files that can not be identified as
    % datasets.)
    bIsDatasetArray(iFile, 1) = false;
  end
end
end
