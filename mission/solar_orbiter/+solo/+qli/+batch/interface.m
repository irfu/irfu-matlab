%
% Code associated with the "interface".
%
% Many functions have string arguments suitable for being passed on from the
% user in e.g. bash scripts calling MATLAB code. Code should thus have
% human-readable error messages for bad arguments.
%
% Therefore also giving better error messages for illegal arguments (and number
% of arguments) supplied from the top-level caller (e.g. bash wrapper).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef interface
  % PROPOSAL: Better automatic tests.
  % PROPOSAL: Merge TIME_INTERVAL and QLI_FMD_INTERVAL.
  %   maxNDays fmdMinDate fmdMaxDate dataMinDate dataMaxDate



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function DaysDtArray = get_days_from_selected_algorithm(...
        Settings, fmdQliDir, dateSelectionAlgorithmId, algorithmArgumentsCa)
      % PROPOSAL: Replace Settings-->Separate arguments.

      % NOTE: "Settings" is deliberately not passed on to any function, so as to
      % not hide its usage (and implicitly meaning).

      assert(iscell(algorithmArgumentsCa))

      switch(dateSelectionAlgorithmId)
        case 'TIME_INTERVAL'
          DaysDtArray = solo.qli.batch.interface.get_days_from_time_interval(...
            algorithmArgumentsCa);

        case 'LOGS'
          DaysDtArray = solo.qli.batch.interface.get_days_from_logs(...
            Settings.LogFileDirPatternDict, ...
            algorithmArgumentsCa);

        case 'DMRQ'
          DaysDtArray = solo.qli.batch.interface.get_days_from_DMRQ(...
            Settings.datasetDirsCa, fmdQliDir, Settings.Fsr, ...
            algorithmArgumentsCa);

        case 'QLI_FMD_INTERVAL'
          DaysDtArray = solo.qli.batch.interface.get_days_from_QLI_FMD_interval( ...
            fmdQliDir, Settings.Fsr, algorithmArgumentsCa);

        otherwise
          error('Illegal argument dateSelectionAlgorithmId="%s"', dateSelectionAlgorithmId)
      end
    end



    % Interpret argument for main function interface.
    %
    function value = interpret_boolean_flag(arg)
      % NOTE: num2str() converts string/number-->string.
      assert(isscalar(arg), 'Flag argument "%s" is not scalar.',   num2str(arg))
      assert(ischar(arg),   'Flag argument "%s" is not a string.', num2str(arg))

      if     ischar(arg) && arg=='0'
        value = false;
      elseif ischar(arg) && arg=='1'
        value = true;
      else
        error('Can not interpret argument flag="%s". Illegal value.', arg)
      end
    end



    % Function for checking whether string-valued argument is a date on form
    % YYYY-MM-DD. Function is made to make it easy to generate more
    % easy-to-understand error messages.
    function check_interface_date_str(dateStr)
      DATE_RE = '^[0-9][0-9][0-9][0-9]-[0-1][0-9]-[0-3][0-9]$';

      assert(ischar(dateStr))

      i = regexp(dateStr, DATE_RE, 'once');
      if isempty(i)
        error('String "%s" is not on the form YYYY-MM-DD.', dateStr)
      end
    end



    function DaysDtArray = filter_days_array(...
        DaysDtArray, maxNDaysStr, beginDayUtcInclStr, endDayUtcExclStr)

      solo.qli.utils.assert_UTC_midnight_datetime(DaysDtArray)

      solo.qli.batch.interface.check_interface_date_str(beginDayUtcInclStr)
      solo.qli.batch.interface.check_interface_date_str(endDayUtcExclStr)

      Dt1      = solo.qli.utils.umdt(beginDayUtcInclStr);
      Dt2      = solo.qli.utils.umdt(endDayUtcExclStr);
      maxNDays = str2double(maxNDaysStr);

      if ~isnumeric(maxNDays)
        error('Argument "%s" could not be interpreted as a number.', maxNDaysStr)
      end

      % IMPLEMENTATION NOTE: Filtering by time interval BEFORE filtering by max
      % number of dates (since that makes sense).

      % Remove date outside specified time interval.
      bKeep       = (Dt1 <= DaysDtArray) & (DaysDtArray < Dt2);
      DaysDtArray = DaysDtArray(bKeep);

      % Remove last dates, if exceeding number of dates.
      bKeep       = [1:min(maxNDays, numel(DaysDtArray))]';
      DaysDtArray = DaysDtArray(bKeep);
    end



    % Generate array of dates for all dates within a specified time interval.
    %
    %
    % NOTE: Arguments are designed to partly handle arguments deriving from
    %       the bash/the OS and are therefore on a string format.
    %
    %
    % ARGUMENTS
    % =========
    % algorithmArgumentsCa{1} = beginDayUtcInclStr
    % algorithmArgumentsCa{2} = endDayUtcExclStr
    %       Strings on format YYYY-MM-DD.
    %       Beginning and end of time interval for which quicklooks should be
    %       generated.
    %       beginDayUtcInclStr is INCLUSIVE. endDayUtcExclStr is EXCLUSIVE, i.e.
    %       the day AFTER the last day of the time interval.
    %
    %
    % RETURN VALUES
    % =============
    % (None)
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    % First created 2022-08-30.
    %
    function DaysDtArray = get_days_from_time_interval(...
      algorithmArgumentsCa)

      assert(numel(algorithmArgumentsCa) == 2, 'Illegal number of algorithm arguments.')
      beginDayUtcInclStr = algorithmArgumentsCa{1};
      endDayUtcExclStr   = algorithmArgumentsCa{2};

      solo.qli.batch.interface.check_interface_date_str(beginDayUtcInclStr)
      solo.qli.batch.interface.check_interface_date_str(endDayUtcExclStr)

      BeginDayInclDt = solo.qli.utils.umdt(beginDayUtcInclStr);
      EndDayExclDt   = solo.qli.utils.umdt(endDayUtcExclStr);

      % NOTE: Indirectly assertion on the string timestamps.
      solo.qli.utils.assert_UTC_midnight_datetime(BeginDayInclDt)
      solo.qli.utils.assert_UTC_midnight_datetime(EndDayExclDt)



      % IMPLEMENTATION NOTE: Subtracting one day from argument for end timestamp to
      % ensure that it is in accordance with the definition of the corresponding
      % argument.

      % IMPLEMENTATION NOTE: Needs to use "caldays()" not "days()" for handling leap
      %                      seconds.
      EndDayInclDt = EndDayExclDt - caldays(1);



      % Construct array of timestamps, where every timestamp represents one day
      % beginning on that timestamp.
      % -----------------------------------------------------------------------
      % IMPLEMENTATION NOTE: Needs to use "caldays()" not "days()" for handling leap
      %                      seconds.
      DaysDtArray = [BeginDayInclDt:caldays(1):EndDayInclDt]';
    end



    %
    % Generate array of dates derived from implicitly specified log files, which
    % should indicate updated datasets.
    %
    %
    % ARGUMENTS
    % =========
    % Settings
    %       Struct with fields for relevant settings.
    % algorithmArgumentsCa
    %       String IDs representing different dataset source directories. The
    %       function will use the logs for the specified directories to derive dates
    %       for which QLIs should be generated.
    %
    %
    % RETURN VALUES
    % =============
    % DaysDtArray
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    %
    function DaysDtArray = get_days_from_logs(LogFileDirPatternDict, algorithmArgumentsCa)
      assert(isequal(...
        sort(LogFileDirPatternDict.keys), ...
        sort(solo.qli.batch.const.SOURCE_DSI_DICT.keys)), ...
        'Settings.LogFileDirPatternDict defines the wrong set of keys.')

      DaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;

      for i = 1:numel(algorithmArgumentsCa)
        datasetsSourceId = algorithmArgumentsCa{i};

        assert(ischar(datasetsSourceId), 'datasetsSourceId{%i} is not a string.', i)
        if ~solo.qli.batch.const.SOURCE_DSI_DICT.isKey(datasetsSourceId)
          error('Illegal datasetsSourceId{%i}="%s"', i, datasetsSourceId)
        end

        dsiCaCa           = solo.qli.batch.const.SOURCE_DSI_DICT(datasetsSourceId);
        dsiCa             = dsiCaCa{1};
        logFileDirPattern = LogFileDirPatternDict(datasetsSourceId);

        SourceDaysDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
          logFileDirPattern, dsiCa);

        DaysDtArray = [DaysDtArray; SourceDaysDtArray];
      end

      DaysDtArray = unique(DaysDtArray);
    end



    % Generate array of dates derived from file modification dates, which should
    % indicate datasets which are newer than the corresponding quicklooks.
    %
    %
    % ARGUMENTS
    % =========
    % Settings
    %       Struct with fields for relevant settings.
    %
    %
    % RETURN VALUES
    % =============
    % DaysDtArray
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    %
    function DaysDtArray = get_days_from_DMRQ(...
        datasetDirsCa, qliDir, Fsr, algorithmArgumentsCa)

      assert(numel(algorithmArgumentsCa) == 3, 'Illegal number of algorithm arguments.')
      maxNDaysStr        = algorithmArgumentsCa{1};
      beginDayUtcInclStr = algorithmArgumentsCa{2};
      endDayUtcExclStr   = algorithmArgumentsCa{3};

      solo.qli.batch.interface.check_interface_date_str(beginDayUtcInclStr)
      solo.qli.batch.interface.check_interface_date_str(endDayUtcExclStr)

      dsiCa = [solo.qli.batch.const.SOURCE_DSI_DICT.values{:}]';

      % Get raw list of DMRQ days.
      DaysDtArray = solo.qli.batch.fmd.get_days_from_DMRQ_and_FS(...
        datasetDirsCa, qliDir, dsiCa, Fsr);

      % Filter list of days.
      DaysDtArray = solo.qli.batch.interface.filter_days_array(...
        DaysDtArray, maxNDaysStr, beginDayUtcInclStr, endDayUtcExclStr);
    end



    function DaysDtArray = get_days_from_QLI_FMD_interval(...
        qliDir, Fsr, algorithmArgumentsCa)

      assert(numel(algorithmArgumentsCa) == 3, 'Illegal number of algorithm arguments.')

      maxNDaysStr    = algorithmArgumentsCa{1};
      startInclFmdDt = datetime(algorithmArgumentsCa{2});
      stopExclFmdDt  = datetime(algorithmArgumentsCa{3});

      QliDfmdd = solo.qli.batch.fmd.get_days_from_QLI_FMD_interval( ...
        qliDir, startInclFmdDt, stopExclFmdDt, Fsr);
      if QliDfmdd.numEntries == 0
        DaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;
      else
        DaysDtArray = QliDfmdd.keys;
      end

      DaysDtArray = solo.qli.batch.interface.filter_days_array(...
        DaysDtArray, maxNDaysStr, '0000-01-01', '9999-12-31');
    end



  end    % methods(Static)



end
