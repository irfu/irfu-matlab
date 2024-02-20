%
% Find all groups of datasets (DSMDs) such that for every group
% * there is exactly one dataset of each specified DATASET_ID, and
% * all datasets overlap in time (there is some time interval which is
%   represented in all datasets).
%   NOTE: Overlap must have non-zero length.
%
% NOTE: Algorithm only works under the assumption that datasets with the same
%       DATASET_ID do not overlap (in time) (assertion). Can thus not handle
%       SOLO_L2_RPW-LFR-SBM1/2-CWF-E.
%
%
% ARGUMENTS
% =========
% DsmdArray
% datasetIdList
%       Cell array of unique DATASET_IDs.
%
%
% RETURN VALUE
% ============
% dsmdGroupsCa
%       Nx1 cell array.
%       dsmdGroupsCa{iGroup} = DSMD array of datasets overlapping in time.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-13.
%
function dsmdGroupsCa = find_overlapping_DSMD_groups(DsmdArray, datasetIdList)
% PROPOSAL: Better name
%   PROPOSAL: find_time_overlapping_groups, find_time_groups
%   PROPOSAL: Name find_overlapping_DSMD_groups
%   PROPOSAL: ~DSMD
%
% PROPOSAL: Convert to DSMD method.
%   PROPOSAL: Instance method?
%   PROPOSAL: Static method?
%
% PROPOSAL: Assertion for identifying dataset groups with only unique DATASET_ID.
% PROPOSAL: Assertion for that datasets with the same DATASET_ID do not overlap in time.
%
% NEED: Be able to find overlapping datasets of same DATASET_ID.
%   Ex: Compare versions of datasets.
%       Ex: Sweeps.
%   PROBLEM: Not well defined problem. A+B, B+C may overlap while A+C d not overlap.

assert(isa(DsmdArray, 'solo.adm.DSMD'))
solo.adm.assert_no_time_overlap(DsmdArray);

DsmdArray = solo.adm.filter_DSMD_DATASET_ID(DsmdArray, datasetIdList);
DsmdArray = DsmdArray(:);

if isempty(DsmdArray)
  dsmdGroupsCa = cell(0, 1);
else
  DT0 = datetime(2000, 1, 1, 0, 0, 0, 'TimeZone', 'UTCLeapSeconds');

  % IMPLEMENTATION NOTE: Does not work for empty DsmdArray in which case
  % vertcat(DsmdArray.dt1) etc. does not produce a datetime array.
  t1Array = seconds(vertcat(DsmdArray.dt1) - DT0);
  t2Array = seconds(vertcat(DsmdArray.dt2) - DT0);

  [setsCa, nArray, oiT1Array, oiT2Array] = irf.utils.find_interval_overlaps(...
    t1Array, t2Array);

  % Only keep overlaps with (1) the exact number of datasets, and (2) non-zero
  % length overlap.
  bKeep = (nArray == numel(datasetIdList)) & (oiT1Array ~= oiT2Array);
  setCa = setsCa(bKeep);

  nGroups = numel(setCa);

  dsmdGroupsCa = cell(nGroups, 1);
  for i = 1:nGroups
    dsmdGroupsCa{i,1} = DsmdArray(setCa{i});
  end
end
end
