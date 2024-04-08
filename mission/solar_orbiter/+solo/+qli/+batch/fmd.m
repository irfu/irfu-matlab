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
%   * keys   = datetime (UTC; midnight) representing days of measured data
%   * values = datetime (no timezone) representing relevant FMD (e.g. most
%     recent for datasets on the day specified in the key).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef fmd
  % NOTE: Uses BICAS code: bicas.tools.batch.get_file_paths().
  %   PROPOSAL: Refactor it as generic function.
  %
  % PROPOSAL: Argument for specifying some kind of upper limit to the number of
  %             days returned.
  %   PROPOSAL: Specify interval of dates.
  %   PROPOSAL: Max number of dates.
  %   PROPOSAL: Separate function in solo.qli.batch.utils for filtering lists of
  %             dates.
  %     PRO: Easily testable.
  %     PRO: If functionality is filtering, then it should be separate from
  %          solo.qli.batch.fmd.
  %
  % ~BUG: Bad thinking: Implementetion checks for all DSIs in all specified
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
  %   --
  %   PROPOSITION: ~Harmless bug.
  %     PRO: Only leads to false positives (never false negatives). Will only
  %          lead to generating more dates, not fewer.
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



    function DaysDtArray = get_days_from_FMDs(datasetDirsCa, qliDir, dsiCa)
      assert(iscell(datasetDirsCa) && iscolumn(datasetDirsCa))
      assert(ischar(qliDir))

      %=========================================
      % Obtain information from the file system
      %=========================================
      % Datasets
      irf.log('n', 'Collecting paths to datasets.')
      [datasetPathsCa, DatasetFsoiArray] = bicas.tools.batch.get_file_paths(datasetDirsCa);
      irf.log('n', 'Obtaining DSMDs from paths.')
      [DsmdArray, bIsDatasetArray]       = solo.adm.paths_to_DSMD_array(datasetPathsCa);
      DatasetFsoiArray   = DatasetFsoiArray(bIsDatasetArray);
      DatasetFmdSdnArray = [DatasetFsoiArray.datenum];
      DatasetFmdSdnArray = DatasetFmdSdnArray(:);
      DatasetFmdDtArray  = datetime(DatasetFmdSdnArray, 'ConvertFrom', 'datenum');

      % QLIs
      irf.log('n', 'Collecting paths to quicklooks.')
      [QliPathsCa, QliFsoiArray] = bicas.tools.batch.get_file_paths({qliDir});
      QliFmdSdnArray = [QliFsoiArray.datenum];
      QliFmdSdnArray = QliFmdSdnArray(:);    % Normalize to column vector.
      QliFmdDtArray  = datetime(QliFmdSdnArray, 'ConvertFrom', 'datenum');

      %==============
      % Derive dates
      %==============
      irf.log('n', 'Determining days for which quicklooks could/should be updated.')
      DaysDtArray = solo.qli.batch.fmd.get_days_from_FMDs_from_file_info(...
        DsmdArray, DatasetFmdDtArray, dsiCa, QliPathsCa, QliFmdDtArray);
    end



    % Get array of days for which to generate QLIs from file system information
    % which is only specified in arguments.
    %
    % IMPLEMENTATION NOTE: This function separates the algorithms/logic from the
    % file system reading so as to have a function that is nice for automated
    % testing.
    %
    function DaysDtArray = get_days_from_FMDs_from_file_info(...
        DsmdArray, DatasetFmdDtArray, dsiCa, QliPathsCa, QliFmdDtArray)

      DatasetsFmdDict = solo.qli.batch.fmd.get_dataset_DFMDD_all(...
        DsmdArray, DatasetFmdDtArray, dsiCa);
      QliFmdDict      = solo.qli.batch.fmd.get_QLI_DFMDD(...
        QliPathsCa, QliFmdDtArray);

      DaysDtArray = solo.qli.batch.fmd.get_days_from_FMDs_algorithm(DatasetsFmdDict, QliFmdDict);

      solo.qli.utils.assert_UTC_midnight_datetime(DaysDtArray)
    end



    function ChangedDatasetsDtArray = get_days_from_FMDs_algorithm(...
        DatasetsFmdDict, QliFmdDict)

      % IMPLEMENTATION NOTE: An empty dictionary can not specify timezone in
      % keys/values. Must therefore always normalize to UTC first.
      AllDatasetsDtArray = intersect(...
        datetime(DatasetsFmdDict.keys, 'TimeZone', 'UTCLeapSeconds'), ...
        datetime(QliFmdDict.keys,      'TimeZone', 'UTCLeapSeconds'));

      % Preallocate.
      ChangedDatasetsDtArray = NaT(...
        [numel(AllDatasetsDtArray), 1], 'TimeZone', 'UTCLeapSeconds');

      nChangedDatasets = 0;
      for iDatasetDt = 1:numel(AllDatasetsDtArray)
        DatasetDt = AllDatasetsDtArray(iDatasetDt);

        if DatasetsFmdDict(DatasetDt) >= QliFmdDict(DatasetDt)
          nChangedDatasets = nChangedDatasets + 1;
          ChangedDatasetsDtArray(nChangedDatasets, 1) = DatasetDt;
        end
      end

      ChangedDatasetsDtArray = ChangedDatasetsDtArray(1:nChangedDatasets, 1);
    end



    % Same concept as solo.qli.batch.fmd.get_dataset_DFMDD_DSI(),
    % except that it covers multiple DSIs and only keeps the latest FMD for data
    % timestamp (dictionary key) collisions.
    function DayFmdDict = get_dataset_DFMDD_all(DsmdArray, FmdDtArray, dsiCa)
      DayFmdDictCa = cell(0, 1);

      for iDsi = 1:numel(dsiCa)
        dsi = dsiCa{iDsi};
        DayFmdDictCa{iDsi, 1} = solo.qli.batch.fmd.get_dataset_DFMDD_DSI(DsmdArray, FmdDtArray, dsi);
      end

      DayFmdDict = solo.qli.batch.utils.merge_dictionaries_max(DayFmdDictCa, datetime.empty, datetime.empty);
    end



    % Given DSMDs and FMDs, get dictionary of the latest FMDs for every day for
    % a specified DSI.
    %
    % ARGUMENTS
    % =========
    % FmdDtArray
    %       Column array of FMDs for every DSMD.
    %
    % RETURN VALUE
    % ============
    % DayFmdDict
    %       Dictionary day-->FMD
    %
    function DayFmdDict = get_dataset_DFMDD_DSI(DsmdArray, FmdDtArray, dsi)
      assert(isa(DsmdArray,  'solo.adm.DSMD'))
      assert(isa(FmdDtArray, 'datetime')     )
      assert(ischar(dsi))
      irf.assert.sizes(...
        DsmdArray,  [-1], ...
        FmdDtArray, [-1])

      bKeep      = strcmp({DsmdArray.datasetId}, dsi);
      DsmdArray  = DsmdArray(bKeep);
      FmdDtArray = FmdDtArray(bKeep);

      DayFmdDict = dictionary(datetime.empty, datetime.empty);

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
          DayFmdDict = solo.qli.batch.utils.dictionary_set_value_max(...
            DayFmdDict, DatasetDtArray(iDatasetDt), FmdDtArray(iDsmd));
        end
      end
    end



    % Given information returned from dir() for potential QLI files, get
    % dictionary of the latest FMDs for every day.
    %
    function DayFmdDict = get_QLI_DFMDD(QliPathsCa, QliFmdDtArray)
      assert(isa(QliFmdDtArray, 'datetime'))
      irf.assert.sizes(...
        QliPathsCa,    [-1], ...
        QliFmdDtArray, [-1])

      DayFmdDict = dictionary(datetime.empty, datetime.empty);

      for iFile = 1:numel(QliPathsCa)
        [FilenameDt1, FilenameDt2] = solo.qli.utils.parse_quicklook_filename(...
          irf.fs.get_name(QliPathsCa{iFile}));

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
          DayFmdDict = solo.qli.batch.utils.dictionary_set_value_max(...
            DayFmdDict, FilenameDtArray(iDay), QliFmdDtArray(iFile));
        end
      end
    end



  end    % methods(Static)



end
