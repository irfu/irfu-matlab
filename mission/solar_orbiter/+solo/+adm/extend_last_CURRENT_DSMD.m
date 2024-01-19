%
% Given an array of DSMDs, find the last CURRENT datasets (last=highest
% DSMD.dt2; "CURRENT" = SOLO_L1_RPW-BIA-CURRENT), and extend its time coverage
% as stated in the DSMD by a specified amount of time.
%
% This is a hack to make it easier to keep using the latest CURRENT dataset
% longer than officially specified. It is useful for e.g. automatic processing.
%
%
% ARGUMENTS
% =========
% DsmdArray
%       Array of DSMDs. May contain any DATASET_IDs. Must not contain any time
%       overlap for the "CURRENT" datasets (SOLO_L1_RPW-BIA-CURRENT). If there
%       is no CURRENT dataset, then the array is returned without modification
%       and without error.
% timeExtensionDays
%       Number of calendar days to add to the end of the time coverage of the
%       last CURRENT dataset.
%
%
% RETURN VALUES
% =============
% DsmdArray
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-06-23.
%
function DsmdArray = extend_last_CURRENT_DSMD(DsmdArray, timeExtensionDays)
% TODO-DEC: How handle the case that there is no CURRENT dataset?
%   PROPOSAL: Do nothing.
%       PRO: bicas_batch:
%           Should be possible to not have any datasets.
%           Should be possible (in principle) to only specify datasets that define modes that do not use CURRENT.

CURRENT_DSI = 'SOLO_L1_RPW-BIA-CURRENT';

irf.assert.vector(DsmdArray)
assert(isscalar( timeExtensionDays))
assert(isnumeric(timeExtensionDays))
assert(timeExtensionDays >= 0)

iCurArray = find(strcmp({DsmdArray.datasetId}, CURRENT_DSI));

if isempty(iCurArray)
  return
end

CurDsmdArray = DsmdArray(iCurArray);
solo.adm.assert_no_time_overlap(CurDsmdArray)   % Overkill?

[~, iiLast] = max(vertcat(CurDsmdArray.dt2));
iLast = iCurArray(iiLast);

Dsmd1 = DsmdArray(iLast);

% NOTE: daysadd() adds calender days, i.e. 86400+-1 secmonds, depending on
%       leap seconds (for UTCLeapSeconds).
dt2 = daysadd(Dsmd1.dt2, timeExtensionDays);

Dsmd2 = solo.adm.DSMD(...
  Dsmd1.path, Dsmd1.datasetId, Dsmd1.versionNbr, ...
  Dsmd1.isCdag, Dsmd1.dt1, dt2);

DsmdArray(iLast) = Dsmd2;
end
