%
% Assert that there is no overlap in time between datasets for every unique
% DATASET_ID separately.
%
% NOTE: Algorithm does not care about dataset version.
% NOTE: Overlap must have non-zero length.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-15.
%
function assert_no_time_overlap(DsmdArray)
    
    assert(isa(DsmdArray, 'solo.adm.DSMD'))

    datasetIdCa = unique({DsmdArray.datasetId});    % Unique DATASET_IDs.

    for iDsi = 1:numel(datasetIdCa)
        DsiDsmdArray = solo.adm.filter_DSMD_DATASET_ID(...
            DsmdArray, {datasetIdCa{iDsi}});

        DsiDsmdArray = DsiDsmdArray(:);

        assert_no_time_overlap_DSI(DsiDsmdArray, datasetIdCa{iDsi})
    end
end



% Assert not time overlaps, assuming all DSMDs have the same DATASET_ID.
function assert_no_time_overlap_DSI(DsmdArray, datasetId)
    % PROPOSAL: Refactor to "validate" arbitrary list of DSMDs and return list
    %           of lists of overlapping DSMDs. Not raise error. Let parent
    %           function handle errors.
    %   PRO: Cleaner
    %   PRO: Can extend code to generate error message including all overlapping
    %        DSMDs (for all DATASET_IDs).

    % Convert timestamps to scalar with arbitrary epoch and unit.
    DT0 = datetime('2000-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
    t1Array = seconds(vertcat(DsmdArray.dt1) - DT0);
    t2Array = seconds(vertcat(DsmdArray.dt2) - DT0);

    [setsCa, nArray, oiT1Array, oiT2Array] = ...
        irf.utils.find_interval_overlaps(t1Array, t2Array);

    % Remove zero-length output intervals.
    bKeep = (oiT1Array ~= oiT2Array);
    clear oiT1Array oiT2Array
    setsCa = setsCa(bKeep);
    nArray = nArray(bKeep);

    %===========
    % ASSERTION
    %===========
    iOverlapsArray = find(nArray > 1);
    if numel(iOverlapsArray) >= 1

        % NOTE: The same dataset may be present in multiple sets/groups of
        % datasets that overlap.
        errorMsg = sprintf(...
            'Files with DATASET_ID="%s" overlap.', datasetId);
        for iOverlap = iOverlapsArray(:)'
            errorMsg = sprintf(...
                '%s\n    Datasets that overlap in a given time interval.', ...
                errorMsg);

            set = setsCa{iOverlap};
            for iDataset = 1:numel(set)
                errorMsg = sprintf('%s\n        %s', ...
                    errorMsg, DsmdArray(set(iDataset)).path);
            end

        end
        error('assert_no_time_overlap:Assertion', errorMsg);

    end

end
