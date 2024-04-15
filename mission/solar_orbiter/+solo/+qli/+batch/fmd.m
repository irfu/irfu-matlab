%
% Code for obtaining days to reprocess QLIs based on file modification dates of
% datasets and pre-existing QLIs.
%
%
% NOTE: There are two types of timestamps in the implementation: (1) timestamps
% referring to when data was measured (UTC; present in filenames and quicklook
% filenames), and (2) timestamps referring to FMDs (~local time; not UTC).
%
%
% ADDITIONAL NAMING CONVENTIONS
% =============================
% DFMDD
%   Day-to-FMD Dictionary.
%   Dictionary with
%   * keys   = datetime (UTC; midnight) representing specific days of measured
%              data
%   * values = datetime (no timezone) representing relevant FMDs (e.g. most
%              recent FMDs for all input datasets on the day specified in the
%              key).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef fmd
  % NOTE: Uses BICAS code: bicas.tools.batch.get_file_paths().
  %   PROPOSAL: Refactor it to a generic function.
  %
  % PROPOSAL: Argument for specifying some kind of upper limit to the number of
  %             days returned.
  %   PROPOSAL: Specify interval of dates. -- IMPLEMENTED
  %   PROPOSAL: Max number of dates. -- IMPLEMENTED
  %   PROPOSAL: Separate function in solo.qli.batch.utils for filtering lists of
  %             dates. -- IMPLEMENTED
  %     PRO: Easily testable.
  %     PRO: If functionality is filtering, then it should be separate from
  %          solo.qli.batch.fmd.
  %
  % PROPOSAL: Function for generating list of all dates for which quicklook FMDs
  %           are in a specific time range. Expose that functionality to
  %           user.
  %   PRO: Can regenerate QLIs which were generated by a specific code version
  %        (e.g. one that is affected by a known bug).
  %   PRO: Can regenerate only the oldest QLIs, e.g. for regenerating all QLIs
  %        in multiple steps (separate batch runs).
  %   PROBLEM: Can not be implemented with external filter function for list of
  %            days.
  %     PROPOSAL: Return list of days, AND corresponding QLI FMDs.
  %       PRO: Can create filtering function which considers FMDs.
  %
  % PROPOSAL: Use solo.db_list_files() for obtaining paths to datasets.
  %   TODO-NI: Unclear if performance is reasonable.
  %   NOTE: Needs to rewrite tests.
  %   PRO: Avoids dependence on explicitly specifying dataset directories.
  %        Uses exactly the CDF directories which SolO DB would use, and as
  %        function of DSI.
  %   PRO: Resolves bug ("Implementetion checks for *all* DSIs in *all*
  %        specified dataset directories.")
  %   --
  %   CON: Preliminary tests (on local mount) indicate that it takes much more
  %        time than using solo.adm on concrete paths (current implementation
  %        2024-04-12).
  %     Ex:
  %         TINT = irf.tint('2020-02-12T00:00:00/2021-04-13T00:00:00');
  %         FsoiArray = solo.db_list_files('solo_L2_mag-rtn-normal', TINT);
  %         ==>  274.6 s
  %   CON: Does not iterate over DSIs, but "dataset prefixes" (w/wo. cdag;
  %        always uppercase archiving level), since that is the
  %        solo.db_list_files() argument. Can not (elegantly) convert DSI to
  %        dataset prefix.
  %   --
  %   PROPOSAL: Call it once for the entire (mission) time range in one call (per DSI).
  %             Still needs to use solo.adm+DSMDs.
  %     PRO: Likely faster than calling once per day.
  %   PROPOSAL: Call it once for each day (per DSI).
  %             Can avoid(?) using solo.adm+DSMDs.
  %     CON: Likely slower.
  %     NOTE: Can yield multiple day files due to midnight overlap.
  %
  % PROPOSAL: Merge
  %   solo.qli.batch.fmd.construct_QLI_DFMDD()
  %   solo.qli.batch.fmd.get_QLI_DFMDD()
  % PROPOSAL: Merge in more of the related functionality from the one call into
  %   solo.qli.batch.fmdget_dataset_DFMDD_for_all_DSIs(),
  %
  % PROPOSAL: Class for DFMDDs.
  %
  %
  %
  % ~BUG: Bad thinking: Implementetion checks for *all* DSIs in *all* specified
  %       dataset directories.
  %   If using both IRFU and LESIA dataset directories, then some datasets are
  %   found in both. Ideally, the algorithm should only check for the right
  %   datasets in the right directory.
  %   --
  %   PROPOSAL: For every directory, store a set of DSIs.
  %     NOTE: If using subdirectories (for speeding up), e.g.
  %           /data/solo/remote/data/L3/lfr_efield/ for EFIELD,
  %           then these directories are separate from the space of the
  %           directory IDs (SOAR, LESIA, IRFU) which specify logs.
  %   PROPOSAL: Caller should specify a list of dataset directories which are so
  %             granular that there is no overlap in datasets.
  %     Ex: /data/solo/data_irfu/latest/rpw/L3/ and
  %         /data/solo/remote/data/L2/
  %     CON: Risk of manual mistakes.
  %   --
  %   PROPOSITION: ~It is a harmless bug.
  %     PRO: Only leads to false positives (never false negatives), i.e. will
  %          only lead to generating more dates, not fewer.
  %     PRO: Not yet using IRFU and LESIA simultaneously.
  %     PRO: If using IRFU and LESIA simultaneously, then the overlapping DSIs
  %          (currently BIAS L3) are more likely to be more recent in IRFU
  %          directory, than in LESIA directory.
  %       CON: Not just after a delivery of L3 datasets.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function QliDfmdd = get_days_from_QLI_FMD_interval(qliDir, startInclFmdDt, stopExclFmdDt, Fsr)
      assert(startInclFmdDt <= stopExclFmdDt)

      QliDfmdd = solo.qli.batch.fmd.get_QLI_DFMDD(qliDir, Fsr);
      QliDfmdd = solo.qli.batch.fmd.filter_DFMDD_by_FMD(QliDfmdd, startInclFmdDt, stopExclFmdDt);
    end



    % For a given a DFMDD, filter it to only keep the FMDs in a certain range.
    function QliDfmdd2 = filter_DFMDD_by_FMD(QliDfmdd1, startInclFmdDt, stopExclFmdDt)
      assert(strcmp(startInclFmdDt.TimeZone, ''))
      assert(strcmp(stopExclFmdDt.TimeZone,  ''))

      QliDfmdd2 = solo.qli.batch.DayDayDictionary();

      DayDtArray = QliDfmdd1.DaysDtArray();
      FmdDtArray = QliDfmdd1.FmdDtArray();
      for i = 1:QliDfmdd1.n
        DayDt = DayDtArray(i);
        FmdDt = FmdDtArray(i);

        if (startInclFmdDt <= FmdDt) && (FmdDt < stopExclFmdDt)
          QliDfmdd2(DayDt) = FmdDt;
        end
      end
    end



    % Derive array of dates from the DMRQ algorithm and reading file system
    % data.
    function DaysDtArray = get_days_from_DMRQ_and_FS(datasetDirsCa, qliDir, dsiCa, Fsr)
      assert(iscell(datasetDirsCa) && iscolumn(datasetDirsCa))
      assert(ischar(qliDir))

      %=========================================
      % Obtain information from the file system
      %=========================================
      % Datasets
      irf.log('n', 'Collecting paths and FMDs for datasets.')
      [datasetPathsCa, DatasetFmdDtArray] = Fsr.get_file_paths_FMDs(datasetDirsCa);
      irf.log('n', 'Obtaining DSMDs from paths.')
      [DsmdArray, bIsDatasetArray]         = solo.adm.paths_to_DSMD_array(datasetPathsCa);
      DatasetFmdDtArray = DatasetFmdDtArray(bIsDatasetArray);
      DatasetFmdDtArray = DatasetFmdDtArray(:);
      DatasetsDfmdd = solo.qli.batch.fmd.get_dataset_DFMDD_for_all_DSIs(...
        DsmdArray, DatasetFmdDtArray, dsiCa);

      % QLIs
      irf.log('n', 'Collecting paths and FMDs for QLIs.')
      QliDfmdd = solo.qli.batch.fmd.get_QLI_DFMDD(qliDir, Fsr);

      %==============
      % Derive dates
      %==============
      irf.log('n', 'Determining days for which QLIs could/should be updated (DMRQ).')
      DaysDtArray = solo.qli.batch.fmd.get_days_from_DMRQ_algorithm(...
        DatasetsDfmdd, QliDfmdd);

      solo.qli.utils.assert_UTC_midnight_datetime(DaysDtArray)
    end



    % Derive array of dates from the DMRQ algorithm and arguments containing
    % file system data.
    function ChangedDatasetsDtArray = get_days_from_DMRQ_algorithm(...
        DatasetsDfmdd, QliDfmdd)

      % IMPLEMENTATION NOTE: An empty dictionary can not specify timezone in
      % keys/values. Must therefore always normalize to UTC first.
      AllDatasetsDtArray = intersect(...
        datetime(DatasetsDfmdd.DaysDtArray(), 'TimeZone', 'UTCLeapSeconds'), ...
        datetime(QliDfmdd.DaysDtArray(),      'TimeZone', 'UTCLeapSeconds'));

      % Preallocate.
      ChangedDatasetsDtArray = NaT(...
        [numel(AllDatasetsDtArray), 1], 'TimeZone', 'UTCLeapSeconds');

      nChangedDatasets = 0;
      for iDatasetDt = 1:numel(AllDatasetsDtArray)
        DatasetDt = AllDatasetsDtArray(iDatasetDt);

        if DatasetsDfmdd(DatasetDt) >= QliDfmdd(DatasetDt)
          nChangedDatasets = nChangedDatasets + 1;
          ChangedDatasetsDtArray(nChangedDatasets, 1) = DatasetDt;
        end
      end

      ChangedDatasetsDtArray = ChangedDatasetsDtArray(1:nChangedDatasets, 1);
    end



    % Given DSMDs and corresponding FMDs, get DFMDD for the most recent dataset
    % FMDs for multiple DSIs.
    function Dfmdd = get_dataset_DFMDD_for_all_DSIs(DsmdArray, FmdDtArray, dsiCa)
      DfmddCa = cell(0, 1);

      for iDsi = 1:numel(dsiCa)
        DfmddCa{iDsi, 1} = solo.qli.batch.fmd.get_dataset_DFMDD_for_one_DSI(...
          DsmdArray, FmdDtArray, dsiCa{iDsi});
      end

      Dfmdd = solo.qli.batch.DayDayDictionary.merge_max(DfmddCa);
    end



    % Given DSMDs and corresponding FMDs, get DFMDD for the most recent dataset
    % FMDs for one specific DSI.
    %
    % ARGUMENTS
    % =========
    % FmdDtArray
    %       Column array of FMDs for every DSMD.
    %
    function Dfmdd = get_dataset_DFMDD_for_one_DSI(DsmdArray, FmdDtArray, dsi)
      assert(isa(DsmdArray,  'solo.adm.DSMD'))
      assert(isa(FmdDtArray, 'datetime')     )
      assert(ischar(dsi))
      irf.assert.sizes(...
        DsmdArray,  [-1], ...
        FmdDtArray, [-1])

      bKeep      = strcmp({DsmdArray.datasetId}, dsi);
      DsmdArray  = DsmdArray(bKeep);
      FmdDtArray = FmdDtArray(bKeep);

      Dfmdd      = solo.qli.batch.DayDayDictionary();

      for iDsmd = 1:numel(DsmdArray)
        % IMPLEMENTATION NOTE: Handle datasets which cover an arbitrary length
        % of time, though that should not be needed yet, since all datasets used
        % for QLI are day-long datasets (midnight to midnight; as of
        % 2024-04-04).
        DatasetDt1 = dateshift(DsmdArray(iDsmd).dt1,                 'start', 'days');
        DatasetDt2 = dateshift(DsmdArray(iDsmd).dt2-milliseconds(1), 'start', 'days');
        DatasetDt2 = max(DatasetDt1, DatasetDt2);
        DatasetDtArray = DatasetDt1:caldays(1):DatasetDt2;

        for iDatasetDt = 1:numel(DatasetDtArray)
          Dfmdd = Dfmdd.set_if_greater(DatasetDtArray(iDatasetDt), FmdDtArray(iDsmd));
        end
      end
    end



    % Given FMDs for paths to potential QLI files, get DFMDD for the most recent
    % QLI FMDs.
    %
    function Dfmdd = construct_QLI_DFMDD(qliPathsCa, qliFmdDtArray)
      assert(iscell(qliPathsCa))
      assert(isa(qliFmdDtArray, 'datetime'))
      irf.assert.sizes(...
        qliPathsCa,    [-1], ...
        qliFmdDtArray, [-1])

      %Dfmdd = dictionary(datetime.empty, datetime.empty);
      Dfmdd = solo.qli.batch.DayDayDictionary();

      for iFile = 1:numel(qliPathsCa)
        [FilenameDt1, FilenameDt2] = solo.qli.utils.parse_quicklook_filename(...
          irf.fs.get_name(qliPathsCa{iFile}));

        if isempty(FilenameDt1)
          continue
        end

        FilenameDt1     = dateshift(FilenameDt1,                 'start', 'day');
        FilenameDt2     = dateshift(FilenameDt2-milliseconds(1), 'start', 'day');
        FilenameDt2     = max(FilenameDt1, FilenameDt2);
        FilenameDtArray = FilenameDt1:caldays(1):FilenameDt2;

        % IMPLEMENTATION NOTE: May have multiple FMDs for the same date due to
        % having
        % (1) filenames for overlapping time intervals (2h+6h+24h+1 week),
        % (2) the same filename in multiple locations (should not happen), and
        % (3) the same path multiple times (should not happen).
        for iDay = 1:numel(FilenameDtArray)
          Dfmdd = Dfmdd.set_if_smaller(FilenameDtArray(iDay), qliFmdDtArray(iFile));
        end
      end
    end



    function QliDfmdd = get_QLI_DFMDD(qliDir, Fsr)
      irf.log('n', 'Collecting paths and FMDs for QLIs.')
      [qliPathsCa, qliFmdDtArray] = Fsr.get_file_paths_FMDs({qliDir});
      qliFmdDtArray = qliFmdDtArray(:);    % Normalize to column vector.
      QliDfmdd      = solo.qli.batch.fmd.construct_QLI_DFMDD(...
        qliPathsCa, qliFmdDtArray);
    end



  end    % methods(Static)



end
