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
% Settings
%       Struct with fields for relevant settings.
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
