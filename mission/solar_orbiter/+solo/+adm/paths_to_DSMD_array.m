%
% Convert paths into array of DSMDs, using filenames only, assuming filenaming
% conventions.
% * Filenames without filenaming convention ==> Ignore file
% * Filenames with    filenaming convention
%   * Filenames without time ==> Ignore file (important!)
%   * Filenames with time:   ==> Use filename to derive time interval.
%       Note: Filenames with only date are assumed to cover the entire day.
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
% DsmdArray
%       Column array of DSMD objects for those input filenames which follow
%       naming conventions.
% bIsDatasetArray
%       Logical column array. Same size as argument. True iff the corresponding
%       input path was interpreted as a dataset (was translated into a DSMD).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-08.
%
function [DsmdArray, bIsDatasetArray] = paths_to_DSMD_array(filePathCa)
% PROPOSAL: Rename
%   PROPOSAL: DSMDs_from_paths()
%   PROPOSAL: paths_to_DSMDs().
%       CON: Information mostly from filename, not path.
%
% PROPOSAL: Refactor to static method for DSMD class.
%   PRO: Is like secondary constructor.
%       CON: Can handle many paths.
%   CON: Function's/method's identifier path becomes long if replaces DSMD-->DatasetMetadata.
%       Ex: solo.adm.DatasetMetadata.DSMD_from_paths().
%       Ex: solo.adm.DSMD.array_from_paths().
%       Ex: solo.adm.DSMD_from_paths().
%   CON: Ties class to filenaming convention.
%       CON: Weak argument (but valid).
%
% PROPOSAL: Preallocate DsmdArray. Remove excess size at the end.
%
% PROPOSAL: Policy for how to handle not being able to recognize filenaming
%           convention for a .cdf file (as could be expected).
%   CON: Do not want to recognize RCTs if applying to entire ROC data/ dir.
%   NOTE: Needs support in parse_dataset_filename(_many).

% FI = File Info
[DfCa, bIsDatasetArray] = solo.adm.dsfn.parse_dataset_filename_many(filePathCa);
datasetPathCa           = filePathCa(bIsDatasetArray);

DsmdArray = solo.adm.DSMD.empty(0, 1);

for i = 1:numel(DfCa)
  Df = DfCa{i};

  Dsmd = solo.adm.DSMD(...
    datasetPathCa{i}, Df.datasetId, Df.versionNbr, Df.isCdag, Df.Dt1, Df.Dt2);

  DsmdArray(end+1, 1) = Dsmd;
end

end
