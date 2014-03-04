function [ENVIR, MMS_CONST] = mms_sdc_sdp_init(scNumber)
% MMS_SDC_SDP_INIT reads initial environment and constants for MMS FIELDS processing
% 	[ENVIR, MMS_CONST] = MMS_SDC_SDP_INIT(scNumber) returns environment variables 
%       and constants useful for MMS processing. Input argument should be
%       the sc number (as a string), i.e. '1' for mms1 and '2' for mms2 etc.
%
%       The struct ENVIR will contain the following:
%	  .DATA_PATH_ROOT	  - Root dir of data files
%	  .CDF_BASE		  - Root dir of CDF tools
%	  .DROPBOX_ROOT		  - Root dir of our output files (temporary location)
%	  .LOG_PATH_ROOT	  - Root dir of log files
% 	  .CAL_PATH_ROOT	  - Root dir of calibration files
%
%	The struct MMS_CONST will contain the following:
%	  .Version.X		  - Major Software version used.
%		  .Y		  - Major Calibration version used.
%		  .Z		  - File version (should perhaps be removed).
%	  .Bitmask.OnlyDCE = 0x01 - Only DCE was found at these points in time. 
%
%	Example:
%		[ENVIR, MMS_CONST] = MMS_SDC_SDP_INIT('1');
%

narginchk(1,1); % SC number to ensure log is put in right place.

ENVIR = [];
MMS_CONST = [];

% Version numbering, start with X, Y, Z = 0, 0, 0. When releasing new
% software update values here and subsequent output files created will have
% these numbers. 
% When simply re-running a dataset, the Z value should be increased by one.

MMS_CONST.Version.X = 0; % Major new Software version
MMS_CONST.Version.Y = 0; % New Calibration version
MMS_CONST.Version.Z = 0; % File revision, increased by 1 for each re-run.

% Bitmask constant values
MMS_CONST.Bitmask.OnlyDCE = 1; % Bits 0x01.

ENVIR.CDF_BASE = getenv('CDF_BASE'); % get environment variable.
ENVIR.DATA_PATH_ROOT = getenv('DATA_PATH_ROOT'); % Get path to data.
ENVIR.LOG_PATH_ROOT = getenv('LOG_PATH_ROOT'); % Get path to logs.
ENVIR.DROPBOX_ROOT = getenv('DROPBOX_ROOT'); % Get path to output location, (temporary location, other scripts then move it once fully written and our script is done). 
ENVIR.CAL_PATH_ROOT = getenv('CAL_PATH_ROOT'); % Get path to cal.

% Setup logging.
% Create a logfile at $LOG_PATH_ROOT / mmsX / sdp /
% named after current run day yyyymmdd and _IRFU.log. If this fails
% create it at $LOG_PATH_ROOT and include full date with seconds.
if(str2double(scNumber)>1||str2double(scNumber)<4)
    % Check to verify that output dir exists for logging. If not, create
    % it.
    if(~exist([ENVIR.LOG_PATH_ROOT, '/mms',scNumber,'/sdp'],'dir'))
        mkdir([ENVIR.LOG_PATH_ROOT, '/mms',scNumber],'sdp');
    end
    irf.log('log_out',strcat(ENVIR.LOG_PATH_ROOT,'/mms',scNumber,'/sdp/',datestr(now,'yyyymmdd'),'_IRFU.log'));
    % Set log level to debug initially.
    irf.log('debug');
else
    irf.log('log_out',strcat(ENVIR.LOG_PATH_ROOT,'/',datestr(now,'yyyymmddTHHMMSS'),'_IRFU.log'));
    irf.log('debug');
    irf.log('critical',['Matlab:MMS_SDC_SDP_INIT:InputArg scNumber incorrectly determined as: ',scNumber]);
    error('Matlab:MMS_SDC_SDP_INIT','MMS_SDC_SDP_INIT recieved an unexpected sc number string. Input to processing should be fullpath/filename.cdf according to MMS standard, mmsX_whatever where X = 1, 2, 3 or 4.');
end
