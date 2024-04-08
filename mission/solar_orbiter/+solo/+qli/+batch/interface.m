%
% Code associated with the "interface".
%
% Many functions have string arguments suitable for being passed on from the
% user in e.g. bash scripts calling MATLAB code.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef interface
  % PROPOSAL: Automatic test code.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function DaysDtArray = get_days_from_selected_algorithm(...
        Settings, dateSelectionAlgorithmId, algorithmArgumentsCa)
      % PROPOSAL: Replace Settings-->Separate arguments.

      assert(iscell(algorithmArgumentsCa))

      switch(dateSelectionAlgorithmId)
        case 'TIME_INTERVAL'
          DaysDtArray = solo.qli.batch.interface.get_days_from_time_interval(algorithmArgumentsCa{:});

        case 'LOGS'
          DaysDtArray = solo.qli.batch.interface.get_days_from_logs(Settings, algorithmArgumentsCa{:});

        case 'FMDS'
          DaysDtArray = solo.qli.batch.interface.get_days_from_FMDs(Settings, algorithmArgumentsCa{:});

        otherwise
          error('Illegal dateSelectionAlgorithmId="%s"', dateSelectionAlgorithmId)
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



    % Generate array of dates for all dates within a specified time interval.
    %
    %
    % NOTE: Arguments are designed to partly handle arguments deriving from
    %       the bash/the OS and are therefore on a string format.
    %
    %
    % ARGUMENTS
    % =========
    % beginDayUtcInclStr, endDayUtcExclStr
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
      beginDayUtcInclStr, endDayUtcExclStr)
      
      assert(ischar(beginDayUtcInclStr))
      assert(ischar(endDayUtcExclStr))
      
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
    % varargin
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
    function DaysDtArray = get_days_from_logs(Settings, varargin)      
      assert(isequal(...
        sort(Settings.LogFileDirPatternDict.keys), ...
        sort(solo.qli.batch.const.SOURCE_DSI_DICT.keys)), ...
        'Settings.LogFileDirPatternDict defines the wrong set of keys.')
      
      
      
      DaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;
      
      for i = 1:numel(varargin)
        datasetsSourceId = varargin{i};
        assert(ischar(datasetsSourceId), 'logFilesId %i is not a string.', i)
      
        if ~solo.qli.batch.const.SOURCE_DSI_DICT.isKey(datasetsSourceId)
          error('Illegal datasetsSourceId="%s"', datasetsSourceId)
        end
      
        dsiCaCa           = solo.qli.batch.const.SOURCE_DSI_DICT(datasetsSourceId);
        dsiCa             = dsiCaCa{1};
        logFileDirPattern = Settings.LogFileDirPatternDict(datasetsSourceId);
      
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
    function DaysDtArray = get_days_from_FMDs(Settings, varargin)    
      assert(numel(varargin) == 0, 'Illegal number of additional arguments.')
      
      dsiCa = [solo.qli.batch.const.SOURCE_DSI_DICT.values{:}]';
      
      DaysDtArray = solo.qli.batch.fmd.get_days_from_FMDs(...
        Settings.datasetDirsCa, ...
        Settings.qliDir, ...
        dsiCa);
    end



  end    % methods(Static)



end
