%
% Iterate over files. Try to parse filenames as if they were datasets using
% standard filenaming conventions (solo.adm.dsfn.DatasetFilename). Ignores
% unparsable filenames without raising error.
%
% NOTE: Function accepts FILENAMES and PATHS, not just paths (as the argument
% name implies), and not just filenames. It is therefore not entirely wrong that
% the function name mentions filenames, not paths.
%
%
% ARGUMENT
% ========
% filePathCa
%       Column cell array of paths to files. Can be both datasets and not.
%
%
% RETURN VALUES
% =============
% DfCa
%       Nx1 cell array of solo.adm.dsfn.DatasetFilename.
% bIsDatasetArray
%       Logical column array. Same size as argument. True iff the corresponding
%       input path was interpreted as a dataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-25.
%
function [DfCa, bIsDatasetArray] = parse_dataset_filename_many(filePathCa)
% PROPOSAL: Change name
%   PROPOSAL: parse_dataset_filenames_many  ("FILENAMES" in plural)
%   PROPOSAL: parse_dataset_filename_many_paths
%   PROPOSAL: parse_dataset_paths_many
%   PROPOSAL: Change to accepting many filenames, not paths.
%
% PROPOSAL: Policy argument for reacting to unparsable filenames.
% PROPOSAL: Convert into method.
%   CON: Not generalizable to having multiple filenaming conventions in separate
%        classes.
%
% NOTE: Has no separate test code. Is indirectly tested by
%       solo.adm.paths_to_DSMD_array___UTEST.

assert(iscell(filePathCa),   'filePathCa is not a cell array.')
assert(iscolumn(filePathCa), 'filePathCa is not a column array.')

DfCa            = cell(0, 1);
bIsDatasetArray = false(numel(filePathCa), 1);
for iFile = 1:numel(filePathCa)

  filename = irf.fs.get_name(filePathCa{iFile});

  Df = solo.adm.dsfn.DatasetFilename.parse_filename(filename);

  if ~isempty(Df)
    % CASE: File can be identified as a dataset.

    DfCa{           end+1, 1} = Df;
    bIsDatasetArray(iFile, 1) = true;
  else
    % CASE: File can *NOT* be identified as dataset.

    % Do nothing. (Silently ignore files that can not be identified as
    % datasets.)
    bIsDatasetArray(iFile, 1) = false;
  end
end

end
