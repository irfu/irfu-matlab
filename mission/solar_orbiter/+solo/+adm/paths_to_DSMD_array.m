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
[fiCa, bIsDatasetArray] = solo.adm.parse_dataset_filename_many(filePathCa);

DsmdArray = solo.adm.DSMD.empty(0, 1);

for i = 1:numel(fiCa)
  Fi = fiCa{i};

  %===========================================
  % Interpret file info from parsing filename
  % (mostly begin & end time)
  %===========================================
  hasDv12 = isfield(Fi, 'dateVec1') && isfield(Fi, 'dateVec2');
  hasDv   = isfield(Fi, 'dateVec');

  if hasDv12 && ~hasDv
    dv1Len = numel(Fi.dateVec1);
    dv2Len = numel(Fi.dateVec2);
    if (dv1Len == 6) && (dv2Len == 6)
      dv1 = Fi.dateVec1;
      dv2 = Fi.dateVec2;
    elseif (dv1Len == 3) && (dv2Len == 3)
      dv1 = [Fi.dateVec1, 0, 0, 0];
      % NOTE: Adding one day since length-3 dateVec2 specifies
      % midnight, not just day.
      dv2 = datevec(datenum(Fi.dateVec2) + 1);
    else
      error(...
        ['Can not interpret parsed filename for "s".', ...
        ' Illegal date vector lengths.'], ...
        filePathCa)
    end

  elseif ~hasDv12 && hasDv
    assert(numel(Fi.dateVec) == 3)

    dv1 = [Fi.dateVec, 0, 0, 0];
    dv2 = datevec(datenum(Fi.dateVec)+1);
  else
    % Skip file since can not derive any time from it.
    continue
  end

  assert(numel(dv1) == 6)
  assert(numel(dv2) == 6)

  Dsmd = solo.adm.DSMD(...
    Fi.path, Fi.datasetId, str2double(Fi.versionStr), Fi.isCdag, ...
    datetime(dv1, 'TimeZone', 'UTCLeapSeconds'), ...
    datetime(dv2, 'TimeZone', 'UTCLeapSeconds') ...
    );

  DsmdArray(end+1, 1) = Dsmd;
end
end
