function mms_sdc_sdp_log_init(procStr, fileStr, runTime)
% MMS_SDC_SDP_LOG_INIT Create log files for MMS processing at SDC.
% Setup log writing for MMS processing at SDC using the designated path and
% name, following naming convention of the files being processed in "procStr".
%
% Creates a log file at destination depending on which process is running
% and what data it is processing. Log files follow the pattern:
% $LOG_PATH_ROOT / mmsX / edp / {tmmode}/ {level} / {descriptor} / {yyyy}
%  / {mm} / mmsX_edp_tmmode_level_descriptor_yyyymmddhhmmss_runtime.log
% (with an extra "/{dd}" added for burst data).
% Where X is one of 1, 2, 3, 4 and (determined by fileStr)
% tmmode is one of fast, slow, brst, comm (determined by fileStr)
% level is one of "ql", "l2a", "l2pre", "l2" (determined by procStr)
% descriptor is one of "dce2d", "scpot" (determined by procStr)
%
% Input:
%  - procStr = Process name, should be one of ql, l2a, l2pre, scpot
%  - fileStr = File name of data being processed.
%  - runTime = Start time of processing, in format yyyymmddTHHMMSS.
%
% Example:
%   MMS_SDC_SDP_LOG_INIT('ql', 'mms2_edp_fast_l1b_dce_20170102_v1.2.3.cdf', '20170105T082536');
% would cause logs to be written as
% "$LOG_PATH_ROOT/mms2/edp/fast/ql/dce2d/2017/01/mms2_edp_fast_ql_dce2d_20170102_20170105T082536.log".
% where "$LOG_PATH_ROOT" is determined by environment variable.
%
% If this fails, logs will be written to
%  $LOG_PATH_ROOT/mmsX/edp/runtime_IRFU.log.
%
% See also: MMS_SDC_SDP_INIT, MMS_CONSTANTS.

narginchk(3,3);
if ~regexp(runTime,'\d{8,8}T\d{6,6}'), runTime = char(datetime("now","Format","uuuuMMdd'T'HHmmss")); end
procStr = lower(procStr);

% Get environment variables
global ENVIR
if isempty(ENVIR), ENVIR = mms_sdc_sdp_init; end

% Identify file, date and time beging processed.
fileInfo = regexp(fileStr, 'mms(?<SCid>[1-4])_edp_(?<tmmode>fast|slow|comm|brst)_(?<dataLevIn>\w+)_(?<descIn>\w+)_(?<dateStr>\d{8,})_\w+','names');

if isempty(fileInfo)
  % Failed to identify file, fall back to old path
  if isempty(ENVIR.LOG_PATH_ROOT)
    irf.log('warning','Environment var LOG_PATH_ROOT not set: logging to screen')
    irf.log('log_out','screen')
  elseif ~exist(ENVIR.LOG_PATH_ROOT, 'dir')
    irf.log('warning',['Logging directory LOG_PATH_ROOT (' ...
      ENVIR.LOG_PATH_ROOT ') doest not exist: logging to screen'])
    irf.log('log_out','screen')
  else
    logDir = [ENVIR.LOG_PATH_ROOT, filesep, 'mms', fileStr(4), filesep, 'edp'];
    if ~exist(logDir, 'dir')
      [s,m] = mkdir(logDir);
      if ~s, error(['Cannot create log directory ', logDir, ': ', m]); end
    end
    irf.log('log_out', [logDir, filesep, runTime, '_IRFU.log']);
  end
else
  % New path, as agreed upon 2017/11/16.
  % $ENVIR.LOG_PATH_ROOT/ mms$SCid/ edp/ $tmmode/ $procIdStr/ $descriptor/
  % $dateStr[1-4]/$dateStr[5-6]/
  switch procStr
    case {'ql', 'l2a', 'l2pre'}
      descr = 'dce2d';
    case 'scpot'
      procStr = 'l2';
      descr = 'scpot';
    otherwise
      error(['Unexpected process input: ', procStr]);
  end
  if strcmp(fileInfo.tmmode, 'brst')
    day = [fileInfo.dateStr(7:8), filesep];
  else
    day = '';
  end
  logDir = [ENVIR.LOG_PATH_ROOT, filesep, 'mms', fileInfo.SCid, filesep, ...
    'edp', filesep, fileInfo.tmmode, filesep, procStr, filesep, descr, ...
    filesep, fileInfo.dateStr(1:4), filesep, fileInfo.dateStr(5:6), filesep, ...
    day];
  if ~exist(logDir, 'dir')
    [s,m] = mkdir(logDir);
    if ~s, error(['Cannot create log directory ', logDir, ': ', m]); end
  end
  logFile = ['mms', fileInfo.SCid, '_edp_', fileInfo.tmmode, '_', ...
    procStr, '_', descr, '_', fileInfo.dateStr, '_', runTime, '.log'];
  irf.log('log_out', [logDir, filesep, logFile]);
end


irf.log('notice'); % Set log level to notice

% Store information in the log file about which version of Matlab is used
% and which version of IRFU-MATLAB.
irf.log('notice', ['Matlab version used is: ' version]);
irf.log('notice', ['irfu-matlab version used is:', irf('version')]);
% Store information in log about the leapsecond table used.
leapsStatus = evalc('spdfcdfleapsecondsinfo');
irf.log('notice', leapsStatus(1:end-1)); % Skip the "\n".

end