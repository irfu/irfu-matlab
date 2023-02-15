%
% This function groups CURRENT DSMDs by calender month. For each month, only keeps the DSMDs with the greatest time
% coverage within each month. This may be multiple datasets per calendar month.
%
% NOTE: solo.adm.group_sort_DSMD_versions calls this function.
% NOTE: Mixes CDAG and non-CDAG (ignores it). Maybe it ideally should, but it likely does not matter.
%
% NOTE: Files on different naming conventions convert times into time intervals (DSMD) differently, which may be
% somewhat unintuitive (but consistent). Below files have the same DSMD.
%   solo_L1_rpw-bia-current-cdag_20200211T000000-20200301T000000_V01.cdf
%   solo_L1_rpw-bia-current-cdag_20200211-20200229_V01.cdf
%
%
% RATIONALE
% =========
% CURRENT datasets have different time coverage in filename for same month. Version number is incremented only when
% using the same time coverage. Therefore, the version number does not represent the true version for CURRENT datasets.
%   Ex: solo_L1_rpw-bia-current-cdag_20200401T000000-20200421T000000_V01.cdf
%       solo_L1_rpw-bia-current-cdag_20200401T000000-20200501T000000_V01.cdf
% The official naming convention will change to eliminate this problem.   /XB e-mail 2020-05-27
% This function exists to make it possible to work around this problem as a "better hack".
% solo.adm.group_sort_DSMD_versions should do the remaining identification of latest versions. This
% function does not need to distinguish between version numbers.
%
%
% ARGUMENTS
% =========
% DsmdArray1
%       DSMD 1D array containing both CURRENT and non-CURRENT datasets.
%
%
% RETURN VALUES
% =============
% DsmdArray2
%       DSMD 1D array. Same as DsmdArray1 (but column vector), where CURRENT
%       datasets with non-largest time overage have been removed.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-06-10.
%
function DsmdArray2 = filter_DSMD_CURRENT_largest_time_coverage(DsmdArray1)
    % NOTE: 2020-07-16: ROC generated
    %   solo_L1_rpw-bia-current-cdag_20200501-20200531_V01.cdf
    % after
    %   former_versions/solo_L1_rpw-bia-current-cdag_20200501T000000-20200601T000000_V02.cdf

    CURRENT_DSI = 'SOLO_L1_RPW-BIA-CURRENT';

    DsmdArray1 = DsmdArray1(:);   % Algorithm requires(?) 1D array.

    bCurrent = ismember({DsmdArray1.datasetId}, CURRENT_DSI);
    bKeep    = ~bCurrent;

    iCurrentArray   = find(bCurrent);
    iLvArray        = find_largest_time_coverages(DsmdArray1(iCurrentArray));
    iLvArray        = iCurrentArray(iLvArray);
    bKeep(iLvArray) = true;

    DsmdArray2 = DsmdArray1(bKeep);
end



% Find exactly one DSMD per calender month represented.
% ASSUMPTION: Only CURRENT DSMDs.
function iLvArray = find_largest_time_coverages(DsmdArray)

    % For each DSMD, find
    % (1) length of time coverage
    % (2) calender month
    DurArray           = duration.empty(0,1);
    monthBeginDtArray = datetime.empty(0,1);
    monthBeginDtArray.TimeZone = 'UTCLeapSeconds';
    for iCurDsmd = 1:numel(DsmdArray)
        Dt1 = vertcat(DsmdArray(iCurDsmd).dt1);
        Dt2 = vertcat(DsmdArray(iCurDsmd).dt2);

        MonthBeginDt = dateshift(Dt1, 'start', 'month');

        % ASSERTION: DSMD does not cover/overlap with more than one CALENDER month.
        assert(Dt2 <= MonthBeginDt + calmonths(1), ...
            'Dataset covers more than one calender month. Can not handle this case. "%s"', ...
            DsmdArray(iCurDsmd).path);

        % Use timestamp for beginning of month as unique ID for calender month.
        DurArray(          iCurDsmd, 1) = Dt2 - Dt1;
        monthBeginDtArray(iCurDsmd, 1)  = MonthBeginDt;
    end



    % For every calender month, find DSMD with greatest time coverage.
    uniqueMonthBeginDtArray = unique(monthBeginDtArray);
    nSets = numel(uniqueMonthBeginDtArray);
    iLvArray = zeros(0, 1);
    % Iterate over sets (months) of datasets. Find one latest version dataset in
    % each set.
    for iSet = 1:nSets
        iSetDsmd = find(monthBeginDtArray == uniqueMonthBeginDtArray(iSet));

        iSetLvArray = find_largest_time_coverage(DsmdArray(iSetDsmd), DurArray(iSetDsmd));
        iLvArray = [iLvArray; iSetDsmd(iSetLvArray)];
    end

end



% Find all DSMDs with the greatest time coverage.
%
% ASSUMPTION: DsmdArray only contains CURRENT DSMDs for one and same calender
% month.
function iLvArray = find_largest_time_coverage(DsmdArray, DurArray)
    % ASSERTION
    assert(numel(DsmdArray) == numel(DurArray))

    MaxDur   = max(DurArray);
    iLvArray = find(DurArray == MaxDur);

    % ASSERTION: Nore more than one maximum duration dataset.
%     if numel(iLvArray) ~= 1
%         errorMsg = sprintf(...
%             ['Found multiple comparable datasets with the equivalent greatest duration (time coverage) %s', ...
%             ' (substitute for version number for CURRET datasets).".', char(MaxDur)]);
%
%         for iLv = 1:numel(iLvArray)
%             path = DsmdArray(iLvArray(iLv)).path;
%             errorMsg = sprintf('%s\n    %s', errorMsg, path);
%         end
%         error('filter_DSMD_CURRENT_largest_time_coverage:MultipleHighestVersionDatasets:Assertion', errorMsg)
%     end
%
%     iLv = iLvArray;
end
