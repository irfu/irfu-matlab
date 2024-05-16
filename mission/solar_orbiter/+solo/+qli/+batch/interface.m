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
  % PROBLEM: Name "interface" could imply that users should look for functions
  %          to call in this file (which is wrong).
  %          Cf. JUICE/RPWI GS TM-to-L1a (python).
  %   PROPOSAL: Better name.
  %     ~interface
  %     ~utils
  %     ~DASA
  %
  % PROPOSAL: Merge TIME_INTERVAL and QLI_FMD_INTERVAL.
  %   maxNDays fmdMinDate fmdMaxDate dataMinDate dataMaxDate
  %
  % PROPOSAL: Include generate_quicklooks_syntax:list_operation



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    %###############
    %###############
    % DAY SELECTION
    %###############
    %###############



    % ARGUMENTS
    % =========
    % Fsr
    %     solo.qli.batch.FileSystemReader* object. Should be an instance of
    %     solo.qli.batch.FileSystemReaderImplementation object in the nominal
    %     case. Arugment exists to facilitate automated tests.
    function DaysDtArray = get_days_from_DASA(...
        datasetDirsCa, LogFileDirPatternDict, Fsr, ...
        fmdQliDir, dasaid, dasaArgumentsCa)

      assert(iscell(datasetDirsCa) & iscolumn(datasetDirsCa))
      assert(isa(LogFileDirPatternDict, 'dictionary'))
      assert(isa(Fsr, 'solo.qli.batch.FileSystemReaderAbstract'))
      assert(iscell(datasetDirsCa) & iscolumn(datasetDirsCa))
      assert(ischar(dasaid))
      assert(iscell(dasaArgumentsCa) & iscolumn(dasaArgumentsCa))

      switch(dasaid)
        case 'TIME_INTERVAL'
          DaysDtArray = solo.qli.batch.interface.get_days_from_time_interval(...
            dasaArgumentsCa);

        case 'LOGS'
          DaysDtArray = solo.qli.batch.interface.get_days_from_logs(...
            Settings.LogFileDirPatternDict, ...
            dasaArgumentsCa);

        case 'DMRQ'
          DaysDtArray = solo.qli.batch.interface.get_days_from_DMRQ(...
            Settings.datasetDirsCa, fmdQliDir, Settings.Fsr, ...
            dasaArgumentsCa);

        case 'QLI_FMD_INTERVAL'
          DaysDtArray = solo.qli.batch.interface.get_days_from_QLI_FMD_interval( ...
            fmdQliDir, Settings.Fsr, dasaArgumentsCa);

        otherwise
          error('Illegal argument dasaid="%s"', dasaid)
      end
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
    % dasaArgumentsCa{1} = beginDayUtcInclStr
    % dasaArgumentsCa{2} = endDayUtcExclStr
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
      dasaArgumentsCa)

      solo.qli.batch.interface.check_nbr_of_DASA_arguments(dasaArgumentsCa, 2)

      beginDayUtcInclStr = dasaArgumentsCa{1};
      endDayUtcExclStr   = dasaArgumentsCa{2};

      solo.qli.batch.interface.check_interface_date_str(beginDayUtcInclStr)
      solo.qli.batch.interface.check_interface_date_str(endDayUtcExclStr)

      BeginDayInclDt = solo.qli.utils.umddt(beginDayUtcInclStr);
      EndDayExclDt   = solo.qli.utils.umddt(endDayUtcExclStr);

      % NOTE: Indirectly assertion on the string timestamps.
      solo.qli.utils.assert_UMD_DT(BeginDayInclDt)
      solo.qli.utils.assert_UMD_DT(EndDayExclDt)



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
    % dasaArgumentsCa
    %       String IDs representing different dataset source directories. The
    %       function will use the logs for the specified directories to derive
    %       dates for which QLIs should be generated.
    %
    %
    % RETURN VALUES
    % =============
    % DaysDtArray
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    %
    function DaysDtArray = get_days_from_logs(LogFileDirPatternDict, dasaArgumentsCa)
      assert(isequal(...
        sort(LogFileDirPatternDict.keys), ...
        sort(solo.qli.batch.const.SOURCE_DSI_DICT.keys)), ...
        'LogFileDirPatternDict defines the wrong set of keys.')



      DaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;

      for i = 1:numel(dasaArgumentsCa)
        datasetsSourceId = dasaArgumentsCa{i};

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
    % RETURN VALUES
    % =============
    % DaysDtArray
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    %
    function DaysDtArray = get_days_from_DMRQ(...
        datasetDirsCa, qliDir, Fsr, dasaArgumentsCa)

      solo.qli.batch.interface.check_nbr_of_DASA_arguments(dasaArgumentsCa, 3)

      maxNDaysStr        = dasaArgumentsCa{1};
      beginDayUtcInclStr = dasaArgumentsCa{2};
      endDayUtcExclStr   = dasaArgumentsCa{3};

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
        qliDir, Fsr, dasaArgumentsCa)

      solo.qli.batch.interface.check_nbr_of_DASA_arguments(dasaArgumentsCa, 3)

      maxNDaysStr    = dasaArgumentsCa{1};
      startInclFmdDt = datetime(dasaArgumentsCa{2});
      stopExclFmdDt  = datetime(dasaArgumentsCa{3});



      QliUfd = solo.qli.batch.fmd.get_days_from_QLI_FMD_interval( ...
        qliDir, startInclFmdDt, stopExclFmdDt, Fsr);

      if QliUfd.n == 0
        DaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;
      else
        DaysDtArray = QliUfd.DaysDtArray;
      end

      DaysDtArray = solo.qli.batch.interface.filter_days_array(...
        DaysDtArray, maxNDaysStr, '0000-01-01', '9999-12-31');
    end



    %###############
    %###############
    % MISCELLANEOUS
    %###############
    %###############



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
      assert(ischar(dateStr))

      i = regexp(dateStr, solo.qli.batch.const.DATE_RE, 'once');
      if isempty(i)
        error('String "%s" is not on the form YYYY-MM-DD.', dateStr)
      end
    end



    function check_nbr_of_DASA_arguments(dasaArgumentsCa, expNArgs)
      if numel(dasaArgumentsCa) ~= expNArgs
        error('Illegal number of algorithm (DASA) arguments. Expected %i arguments', ...
          expNArgs)
      end
    end



    % Utility function for filtering an array of days.
    function DaysDtArray = filter_days_array(...
        DaysDtArray, maxNDaysStr, beginDayUtcInclStr, endDayUtcExclStr)

      solo.qli.utils.assert_UMD_DT(DaysDtArray)

      solo.qli.batch.interface.check_interface_date_str(beginDayUtcInclStr)
      solo.qli.batch.interface.check_interface_date_str(endDayUtcExclStr)

      Dt1      = solo.qli.utils.umddt(beginDayUtcInclStr);
      Dt2      = solo.qli.utils.umddt(endDayUtcExclStr);
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



  end    % methods(Static)



end
