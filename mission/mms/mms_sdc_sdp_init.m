function ENVIR = mms_sdc_sdp_init(scNumberStr)
% MMS_SDC_SDP_INIT  initialize environment for MMS FIELDS processing
%
% 	[ENVIR, MMS_CONST] = MMS_SDC_SDP_INIT(scNumberStr) returns environment
%   variables and constants useful for MMS processing. Input argument 
%   should be the sc number (as a string), i.e. '1' for mms1 and '2' for 
%   mms2 etc. It also configures logging to "$LOG_PATH_ROOT/ mmsX/ edp/
%   date_IRFU.log".
%
%       The struct ENVIR will contain the following:
%	  .DATA_PATH_ROOT       - Root dir of data files
%	  .CDF_BASE             - Root dir of CDF tools
%	  .DROPBOX_ROOT         - Root dir of our output files (temporary location)
%	  .LOG_PATH_ROOT        - Root dir of log files
% 	  .CAL_PATH_ROOT        - Root dir of calibration files
%
%	Example:
%		ENVIR = MMS_SDC_SDP_INIT('1');
%
% 	See also MMS_CONSTANTS.

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

narginchk(1,3); % SC number to ensure log is put in right place.

ENVIR.CDF_BASE = getenv('CDF_BASE'); % Get path to CDF tools.
ENVIR.DATA_PATH_ROOT = getenv('DATA_PATH_ROOT'); % The final path of data.
ENVIR.LOG_PATH_ROOT = getenv('LOG_PATH_ROOT'); % Get path to logs.
ENVIR.DROPBOX_ROOT = getenv('DROPBOX_ROOT'); % Get path to output location,
% DROPBOX_ROOT is only the temporary location, file are then to be moved by
% other scripts to their final location when processing is done.
ENVIR.CAL_PATH_ROOT = getenv('CAL_PATH_ROOT'); % Get path to cal.

% Setup logging.
% Create a logfile at $LOG_PATH_ROOT / mmsX / edp /
% named after current run day yyyymmdd and _IRFU.log. If this fails
% create it at $LOG_PATH_ROOT and include full date with seconds.
if ~ischar(scNumberStr) || ...
    isempty(intersect(str2double(scNumberStr), MMS_CONST.MMSids))
  
  irf.log('log_out', [ENVIR.LOG_PATH_ROOT, filesep, ...
        datestr(now,'yyyymmddTHHMMSS'), '_IRFU.log']);
    errStr = ['invalid input: scNumber idetermined as: ',scNumberStr];
    irf.log('critical', errStr);
    error('Matlab:MMS_SDC_SDP_INIT',['MMS_SDC_SDP_INIT recieved an ', ...
        'unexpected sc number string. Input to processing should be ', ...
        'fullpath/filename.cdf according to MMS standard, mmsX_whatever',...
        ' where X = 1, 2, 3 or 4.']);
end

% Check to verify that output dir exists for logging. If not, create
% it.
if isempty(ENVIR.LOG_PATH_ROOT)
  irf.log('warning','Environment var LOG_PATH_ROOT not set: logging to screen')
  irf.log('log_out','screen')
elseif ~exist(ENVIR.LOG_PATH_ROOT,'dir')
  irf.log('warning',['Logging directory LOG_PATH_ROOT (' ...
    ENVIR.LOG_PATH_ROOT ') doest not exist: logging to screen'])
  irf.log('log_out','screen')
else
  logDir = [ENVIR.LOG_PATH_ROOT, filesep, 'mms', scNumberStr, filesep, 'edp'];
  if ~exist(logDir, 'dir')
    [s,m] = mkdir(logDir);
    if ~s
      error('Matlab:MMS_SDC_SDP_INIT:io',...
        ['Cannot create log directory ' logDir ': ' m])
    end
  end
  irf.log('log_out',...
    [logDir, filesep, datestr(now,'yyyymmdd'), '_IRFU.log']);
end
irf.log('notice'); % XXX: Set log level to notice

% Store information in the log file about which version of Matlab is used
% and which version of IRFU-MATLAB.
irf.log('notice', ['Matlab version used is: ' version]);
irf.log('notice', ['irfu-matlab version used is:', irf('version')]);
