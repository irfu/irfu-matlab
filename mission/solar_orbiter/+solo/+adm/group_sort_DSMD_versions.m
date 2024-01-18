%
% Group DSMDs into "comparable datasets" and sort them according to "version"
% but optionally use more than "version number" to work around certain practical
% problems.
%
%
% NOTE: Might be slow (~minutes) if applying to all ROC data files, including
% former_versions/. This is likely (?) due to irf.utils.find_equalities().
%
% NOTE: Includes "better hacks" to resolve uncertainties, ambiguities and
% complexities in the official naming & versioning scheme. See ALGORITHM.
%
% NOTE: Can not find comparable DSMDs among DMSDs that vary in time coverage
% (except for the CURRENT hack).
%
%
% ALGORITM
% ========
% NOTE: What is the latest version is ambiguous (2020-06-10).
%   EJ's RCS telecon notes 2020-06-04:
%       """"XB: CDAG files will have CDAG removed and start over from V01. Might
%       continue after creation of non-CDAG to again update/reprocess'' CDAG
%       file (and non-CDAG file in parallel). Not clear algorithm/procedure.""""
% NOTE: Empirically, comparing version numbers for both CDAG and non-CDAG (as of
% 2020-06-10) makes no sense.
% --
% NOTE: Only datasets with identical
% (1) DATASET_ID, and
% (2) beginning and end time
% are compared with each other.
% NOTE: Uses solo.adm.filter_DSMD_CURRENT_largest_time_coverage()
% as a workaround for CURRENT datasets which have varying time coverage.
% --
% Sorts in order (might change):
% (1) CDAG vs non-CDAG
% (2) version number (DSMD.versionNbr)
% (3) if enabled: whether datasets is under a former_versions/ directory
%     according to solo.adm.DSMD.path.
%
%
% LOCAL ABBREVIATIONS
% ===================
% FVD = former_versions/ directory.
%
%
% ARGUMENTS
% =========
% DsmdArray1
% mode
%       String constant.
%       latest' : Result = DSMD 1D array. Only latest versions.
%       'all'   : Result = 1D cell array of DSMD 1D arrays. {iSet}(iVersion).
%       NOTE: Latest version = max(iVersion).
% varargin
%       Settings as interpreted by
%       irf.utils.interpret_settings_args(). See implementation.
%       "cdagAlgorithm"
%           String constant determining which algorithm/definition of highest
%           version to use.
%
%
% RETURN VALUE
% ============
% Result
%       See argument "mode".
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-08.
%
function Result = group_sort_DSMD_versions(DsmdArray1, mode, varargin)
    %
    % =================
    % SORTING ALGORITHM
    % =================
    % PROBLEM: How handle DSMDs with time intervals that change somewhat from version to version?
	%   Ex: Sweeps (potentially)
    %   Ex: DSMD time intervals based on content (not filenames).
    %   PROPOSAL: Use some form of clustering of the time intervals. All files within same cluster are regarded as being
    %       equal w.r.t. time.
    %       CON: Inexact.
    %           CON: Yes, but there is no rigorous criterion for determining that two datasets are the same but of
    %                different versions, not even in theory.
    %       PROPOSAL: Use clustering of midpoints.
    %           PRO: Easier to implement clustering algorithm in 1D.
    %   PROPOSAL: Assume that datasets that overlap in time and share DATASET_ID should be compared.
    %       ~CON: Overlap in time is not a transitive relation ==> Not well defined groups of overlapping datasets.
    %           CON-PROPOSAL: Compare in sequence. V01 overlaps with V02, V02 overlaps with V03, ...
    %               CON: Could have missing versions.
    %
    % PROBLEM: How handle DSMDs with exactly the same time interval, same version number, but different naming
    %          conventions?
    %   Ex: CURRENT changing naming convention:
    %       former_versions/solo_L1_rpw-bia-current-cdag_20200211T000000-20200301T000000_V01.cdf
    %       solo_L1_rpw-bia-current-cdag_20200211-20200229_V01.cdf
    %   PROPOSAL: Avoid problem by never calling such sets of files.
    %       PROPOSAL: Filter DSMDs and remove those with file under former_versions/ in path.
    %       PROPOSAL: Never search former_versions directories for files to turn into DSMDs.
    %   PROPOSAL: Use former_versions/ in path when comparing files with same version number.
    %       CON: Should not really use path. Not part of the "real" DSMD. No equivalent when there is no path.
    %       PROPOSAL: Setting for algorithm hack.
    %   PROPOSAL: Add and use file modification date.
    %   PROPOSAL: Add and use glob.attr. Generation_date.
    %   PROPOSAL: Prefer the later filenamning convention over the former.
    %       CON: That information is (presently) only found in path. Would need to parse.
    %       CON: Should be independent of filenaming.
    %   PROPOSAL: Use external DSMD filter hack that removes former_versions/.
    %
    % PROPOSAL: Use filename strings to identify identical datasets (except
    %           version): DATASET_ID + Datetime
    %
    % PROPOSAL: Use glob. attr. Generation_date.
    %   CON: Bad for automatic testing since has to read file.
    %   CON: Slower, since has to read file.
    %   CON: Not necessarily the latest version in a general set of files (same DATASET_ID). Could have re-generated
    %       an old file for debugging.
    %
    % PROPOSAL: Use alternate/modified version of solo.adm.filter_DSMD_CURRENT_largest_time_coverage that
    %           produces rank wrt. CURRENT largest time coverage. Use for analogously to other criteria.
    %   NOTE: Only for CURRENT; no other DATASET_IDs.
    %   PROBLEM: Can not use DSMD fields to identify groups of comparable CURRENT datasets.
    %       NOTE: solo.adm.filter_DSMD_CURRENT_largest_time_coverage uses that CURRENT files should stay
    %       within the same calendar month.
    %
    % =====
    % OTHER
    % =====
    % PROPOSAL: Turn helper functions into separate, generic functions.
    %   NOTE: Somewhat non-standard for a filter_DSMD_* functions in that they returns indices, not DsmdArray.
    %
    % PROPOSAL: sort_DSMDs_wrt_algorithm as separate public function.
    %   PRO: Easier to test.
    %   TODO-DEC: Names?
    %       PROPOSAL: sort_DSMD_algorithmic_version
    %   PROBLEM: How handle different assertions for "latest" vs "all"?
    %       PROPOSAL: Return fhArray
    %       PROPOSAL: Argument n. Check n first in list if they are sorted.
    %
    % PROPOSAL: Separate function for grouping into "comparable datasets".
    %
    % PROPOSAL: For all non-version number sorting criteria: Setting for algorithm
    %   alt 1: Use specified sorting algorithm/order
    %   alt 2: Use specified value. Ignore other values.
    %       Ex: Use only cdag
    %       Ex: Use only non-former_versions/.
    %
    % PROBLEM: Order of DSMDs (when returning 'all' version) in result is confusing. Reverse?
    %
    % PROPOSAL: Not convert objects to structs. Extract arrays for the relevant properties instead.
    %   NOTE: irf.utils.find_equalities works with multiple arrays.
    %
    % PROPOSAL: Abolish cdagAlgorithm, sortWrtFormerVersionsDir.
    % PROPOSAL: Abolish currentHack?

    SETTINGS.cdagAlgorithm            = 'CDAG then non-CDAG';
    SETTINGS.sortWrtFormerVersionsDir = false;
    SETTINGS.currentHack              = true;
    Settings = irf.utils.interpret_settings_args(SETTINGS, varargin);
    irf.assert.struct(Settings, fieldnames(SETTINGS), {})
    assert(islogical(SETTINGS.sortWrtFormerVersionsDir))
    assert(islogical(SETTINGS.currentHack))



    assert(isa(DsmdArray1, 'solo.adm.DSMD'))
    DsmdArray1 = DsmdArray1(:);

    if Settings.currentHack
        DsmdArray1 = solo.adm.filter_DSMD_CURRENT_largest_time_coverage(DsmdArray1);
    end

    DsmdArray1 = DsmdArray1(:);    % Algorithm requires vector.

    %======================================================================
    % Find comparable datasets
    % ------------------------
    % Use subset of fields, including time, to identify sets of comparable
    % datasets.
    %======================================================================
    % IMPORTANT NOTE: irf.utils.find_equalities() can be very slow
    % for large numbers of DSMDs.
    % IMPLEMENTATION NOTE: Can not use vertcat for strings.

    % Arbitrary "epoch" which value should not affect the result.
    DT_EPOCH = datetime('2020-01-01T00:00:00.000Z', 'TimeZone', 'UTCLeapSeconds');
    if isempty(DsmdArray1)
        t1 = zeros(0, 1);
        t2 = zeros(0, 1);
    else
        % IMPLEMENTATION NOTE: vertcat() for empty array returns double (not
        %                      datetime).
        t1 = seconds(vertcat(DsmdArray1.dt1) - DT_EPOCH);
        t2 = seconds(vertcat(DsmdArray1.dt2) - DT_EPOCH);
    end
    datasetIdCa = {DsmdArray1.datasetId};
    datasetIdCa = datasetIdCa(:);

    % IMPORTANT: The call to irf.utils.find_equalities() is critical for
    % performance. The arguments and their order are chosen to speed it up.
    % Fastest arrays first. time_double < time_strings << time_datetime. In
    % particular, avoiding datetime.
    % t = tic();
    fhArray = irf.utils.find_equalities(Inf, t1, t2, datasetIdCa);

    fhUniques = unique(fhArray);
    nSets     = numel(fhUniques);

    %======================================
    % Sort each set of comparable datasets
    %======================================
    iKeepArray = [];
    sortedDsmdArraysCa = cell(nSets, 1);
    % Iterate over sets of datasets. Find one latest version dataset in each set.
    for iSet = 1:nSets
        iFhArray = find(fhArray == fhUniques(iSet));

        iSortArray = sort_DSMDs_wrt_algorithm(DsmdArray1(iFhArray), mode, Settings);

        %iKeepArray(end+1)        = iSortArray(end);   % BUGG
        iKeepArray(end+1)        = iFhArray(iSortArray(end));
        sortedDsmdArraysCa{iSet} = DsmdArray1(iFhArray(iSortArray));
    end



    switch(mode)
        case 'latest'
            Result = DsmdArray1(iKeepArray, 1);
            assert(numel(unique({Result.path})) == nSets)
        case 'all'
            Result = sortedDsmdArraysCa;

        otherwise
            error('group_sort_DSMD_versions:IllegalArgument', ...
                'Illegal argument mode="%s".', mode)
    end

    assert(numel(Result) == nSets)
end



% Sort datasets internally, assuming they are all comparable.
%
% DsmdArray(iSortArray) is the sorted array. Last = Latest version.
%
function iSortArray = sort_DSMDs_wrt_algorithm(DsmdArray, mode, Settings)

    %===================
    % Sorting algorithm
    %===================
    rankCdagArray   = rank_DSMDs_wrt_CDAG(DsmdArray, Settings.cdagAlgorithm);
    rankVerNbrArray = rank_DSMDs_wrt_versionNbr(DsmdArray);

    rankMatrix = [rankCdagArray(:), rankVerNbrArray(:)];

    if Settings.sortWrtFormerVersionsDir
        rankFvdArray = rank_DSMDs_wrt_former_versions_dir(DsmdArray);
        rankMatrix   = [rankMatrix, rankFvdArray(:)];
    end

    [sortedRankMatrix, iSortArray] = sortrows(rankMatrix);



    %=====================================================
    % ASSERTION: Not more than one latest version dataset
    %=====================================================
    nDsmd = numel(DsmdArray);
    switch(mode)
        case 'latest'

            if (nDsmd >= 2) && (size(unique(sortedRankMatrix(end-1:end, :), 'rows'), 1) == 1)
                % CASE: Latest-version, and second-latest-version have identical rank.
                errorMsg = ...
                    ['Algorithm found a set of datasets that should be', ...
                    ' comparable, and that should include exactly ONE', ...
                    ' "latest version" (not just "V02" etc).', ...
                    ' The algorithm failed to find ONE such dataset within this set.'];
                throw_error_can_not_rank(errorMsg, DsmdArray(iSortArray))
            end

        case 'all'

            if (nDsmd >= 2) && (size(unique(sortedRankMatrix, 'rows'), 1) < nDsmd)
                errorMsg = [...
                    'Algorithm found a set of datasets that should be', ...
                    ' comparable, and that should be entirely sortable', ...
                    ' w.r.t. "latest version" (not just "V02" etc).', ...
                    ' The algorithm failed to sort this set..'];
                throw_error_can_not_rank(errorMsg, DsmdArray(iSortArray))
            end

        otherwise
            error('group_sort_DSMD_versions:IllegalArgument', ...
                'Illegal argument mode="%s".', mode)
    end

end



% NOTE: Prefereably submit DSMDs in ranking order.
function throw_error_can_not_rank(errorMsg, DsmdArray)
    for i = 1:numel(DsmdArray)
        path     = DsmdArray(i).path;
        errorMsg = sprintf('%s\n    %s', errorMsg, path);
    end

    error('group_sort_DSMD_versions:MultipleLatestVersionDatasets:Assertion', errorMsg)
end



% Find the one DSMD with highest version number (DsmdArray.versionNbr) among all
% submitted DSMDs.
%
function rankArray = rank_DSMDs_wrt_versionNbr(DsmdArray)
    % PROPOSAL: Remove function since simple.

    rankArray = [DsmdArray.versionNbr];
end



function rankArray = rank_DSMDs_wrt_CDAG(DsmdArray, cdagAlgorithm)
    isCdagArray = [DsmdArray.isCdag];

    switch(cdagAlgorithm)
        case 'CDAG then non-CDAG'
            rankArray = isCdagArray;
        case 'non-CDAG then CDAG'
            rankArray = ~isCdagArray;
        otherwise
            error('group_sort_DSMD_versions:IllegalArgument:Assertion', ...
                'Illegal argument cdagAlgorithm="%s"', cdagAlgorithm)
    end

    rankArray = double(rankArray);
end



function rankArray = rank_DSMDs_wrt_former_versions_dir(DsmdArray)
    FVD_NAME = 'former_versions';

    pathArray = {DsmdArray.path};
    rankArray = zeros(size(pathArray));
    for i = 1:numel(pathArray)
        pathPartsArray = strsplit(pathArray{i}, filesep);
        rankArray(i)   = double(~any(ismember(FVD_NAME, pathPartsArray)));
    end
end
