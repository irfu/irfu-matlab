%
% Code for obtaining days to reprocess QLIs based on file modification dates
% (FMD) of datasets and pre-existing QLIs.
%
%
% NOTE: There are two types of timestamps in the implementation: (1) timestamps
% referring to when data was measured (UTC; present in filenames and quicklook
% filenames), and (2) timestamps referring to FMDs (~local time; not UTC).
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
  %   solo.qli.batch.fmd.construct_QLI_UFD()
  %   solo.qli.batch.fmd.get_QLI_UFD()
  % PROPOSAL: Merge in more of the related functionality from the one call into
  %   solo.qli.batch.fmdget_dataset_UFD_for_all_DSIs(),
  %
  %
  %
  % ~BUG: Bad thinking: Implementation checks for *all* DSIs in *all* specified
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



    function QliUfd = get_days_from_QLI_FMD_interval(qliDir, startInclFmdDt, stopExclFmdDt, Fsr)
      assert(startInclFmdDt <= stopExclFmdDt)

      QliUfd = solo.qli.batch.fmd.get_QLI_UFD(qliDir, Fsr);
      QliUfd = solo.qli.batch.fmd.filter_UFD_by_FMD(QliUfd, startInclFmdDt, stopExclFmdDt);
    end



    % For a given a UFD, filter it to only keep the FMDs in a certain range.
    function QliUfd2 = filter_UFD_by_FMD(QliUfd1, startInclFmdDt, stopExclFmdDt)
      assert(strcmp(startInclFmdDt.TimeZone, ''))
      assert(strcmp(stopExclFmdDt.TimeZone,  ''))

      QliUfd2 = solo.qli.batch.UmdFmdDictionary();

      UmdDtArray = QliUfd1.UmdDtArray();
      FmdDtArray = QliUfd1.FmdDtArray();
      for i = 1:QliUfd1.n
        UmdDt = UmdDtArray(i);
        FmdDt = FmdDtArray(i);

        if (startInclFmdDt <= FmdDt) && (FmdDt < stopExclFmdDt)
          QliUfd2(UmdDt) = FmdDt;
        end
      end
    end



    % Derive array of dates from the DMRQ and by reading file system data
    % itself.
    function UmdDtArray = get_days_from_DMRQ_and_FS(datasetDirsCa, qliDir, dsiCa, Fsr)
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
      DatasetsUfd = solo.qli.batch.fmd.get_dataset_UFD_for_all_DSIs(...
        DsmdArray, DatasetFmdDtArray, dsiCa);

      % QLIs
      irf.log('n', 'Collecting paths and FMDs for QLIs.')
      QliUfd = solo.qli.batch.fmd.get_QLI_UFD(qliDir, Fsr);

      %==============
      % Derive dates
      %==============
      irf.log('n', 'Determining days for which QLIs could/should be updated (DMRQ).')
      UmdDtArray = solo.qli.batch.fmd.get_days_from_DMRQ_algorithm(...
        DatasetsUfd, QliUfd);

      solo.qli.utils.assert_UMD_DT(UmdDtArray)
    end



    % Derive array of dates from the DMRQ and arguments containing file system
    % data (i.e. do not read from file system itself).
    function ChangedDatasetsDtArray = get_days_from_DMRQ_algorithm(...
        DatasetsUfd, QliUfd)

      AllDatasetsDtArray = intersect(...
        DatasetsUfd.UmdDtArray(), QliUfd.UmdDtArray());

      % Preallocate.
      ChangedDatasetsDtArray = NaT(...
        size((AllDatasetsDtArray)), 'TimeZone', 'UTCLeapSeconds');

      nChangedDatasets = 0;
      for iDatasetDt = 1:numel(AllDatasetsDtArray)
        DatasetDt = AllDatasetsDtArray(iDatasetDt);

        if DatasetsUfd(DatasetDt) >= QliUfd(DatasetDt)
          nChangedDatasets = nChangedDatasets + 1;
          ChangedDatasetsDtArray(nChangedDatasets, 1) = DatasetDt;
        end
      end

      ChangedDatasetsDtArray = ChangedDatasetsDtArray(1:nChangedDatasets, 1);
    end



    % Given DSMDs and corresponding FMDs, get UFD for the most recent dataset
    % FMDs for multiple DSIs.
    function Ufd = get_dataset_UFD_for_all_DSIs(DsmdArray, FmdDtArray, dsiCa)
      UfdCa = cell(0, 1);

      for iDsi = 1:numel(dsiCa)
        UfdCa{iDsi, 1} = solo.qli.batch.fmd.get_dataset_UFD_for_one_DSI(...
          DsmdArray, FmdDtArray, dsiCa{iDsi});
      end

      Ufd = solo.qli.batch.UmdFmdDictionary.merge_max(UfdCa);
    end



    % Given DSMDs and corresponding FMDs, get UFD for the most recent dataset
    % FMDs for one specific DSI.
    %
    % ARGUMENTS
    % =========
    % FmdDtArray
    %       Column array of FMDs for every DSMD.
    %
    function Ufd = get_dataset_UFD_for_one_DSI(DsmdArray, FmdDtArray, dsi)
      assert(isa(DsmdArray,  'solo.adm.DSMD'))
      assert(isa(FmdDtArray, 'datetime')     )
      assert(ischar(dsi))
      irf.assert.sizes(...
        DsmdArray,  [-1], ...
        FmdDtArray, [-1])

      bKeep      = strcmp({DsmdArray.datasetId}, dsi);
      DsmdArray  = DsmdArray(bKeep);
      FmdDtArray = FmdDtArray(bKeep);

      Ufd        = solo.qli.batch.UmdFmdDictionary();

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
          Ufd = Ufd.set_if_greater(DatasetDtArray(iDatasetDt), FmdDtArray(iDsmd));
        end
      end
    end



    % Given FMDs for paths to potential QLI files, get UFD for the most recent
    % QLI FMDs.
    %
    function Ufd = construct_QLI_UFD(qliPathsCa, qliFmdDtArray)
      assert(iscell(qliPathsCa))
      assert(isa(qliFmdDtArray, 'datetime'))
      irf.assert.sizes(...
        qliPathsCa,    [-1], ...
        qliFmdDtArray, [-1])

      Ufd = solo.qli.batch.UmdFmdDictionary();

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
          Ufd = Ufd.set_if_smaller(FilenameDtArray(iDay), qliFmdDtArray(iFile));
        end
      end
    end



    function QliUfd = get_QLI_UFD(qliDir, Fsr)
      irf.log('n', 'Collecting paths and FMDs for QLIs.')
      [qliPathsCa, qliFmdDtArray] = Fsr.get_file_paths_FMDs({qliDir});
      qliFmdDtArray = qliFmdDtArray(:);    % Normalize to column vector.
      QliUfd        = solo.qli.batch.fmd.construct_QLI_UFD(...
        qliPathsCa, qliFmdDtArray);
    end



  end    % methods(Static)



end
