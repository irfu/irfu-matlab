%
% Generate array of dates derived from implicitly specified log files, which
% should indicate updated datasets.
%
%
% NOTE: Several arguments are designed to partly handle arguments deriving from
%       the bash/the OS and are therefore on a string format.
%
% NOTE: The dataset IDs referred to in the implementation must be consistent
%       with the use of datasets in
%       solo.qli.generate_quicklooks_24h_6h_2h_using_DB_SPICE() and
%       solo.qli.generate_quicklook_7days().
%
%
% ARGUMENTS
% =========
% outputDir
%       Path to output directory.
% varargin
%       String IDs representing different dataset source directories. The
%       function will use the logs for the specified directories to derive dates
%       for which QLI should be generated.
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

DaysDtArray = solo.qli.const.EMPTY_DT_ARRAY;

for i = 1:numel(varargin)
  logFilesId = varargin{i};
  assert(ischar(logFilesId), 'logFilesId %i is not a string.', i)

  switch(logFilesId)

%     case 'IRFU'
%       % '/home/erjo/logs/so_bicas_batch_cron.*.log'
%       logFileDirPattern = Settings.irfuLogFileDirPattern;
%       dsiCa             = solo.qli.batch.const.IRFU_LOGS_DSI_CA;

    case 'LESIA'
      logFileDirPattern = Settings.lesiaLogFileDirPattern;
      dsiCa             = solo.qli.batch.const.LESIA_LOGS_DSI_CA;

    case 'SOAR'
      logFileDirPattern = Settings.soarLogFileDirPattern;
      dsiCa             = solo.qli.batch.const.SOAR_LOGS_DSI_CA;

    otherwise
      error('Illegal datasetsSourceId="%s"', logFilesId)
  end

  SourceDaysDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
    logFileDirPattern, dsiCa);


  DaysDtArray = [DaysDtArray; SourceDaysDtArray];
end

DaysDtArray = unique(DaysDtArray);

end
