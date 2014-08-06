function [ENVIR, MMS_CONST] = mms_sdc_sdp_init(scNumber)
% MMS_SDC_SDP_INIT reads initial environment and constants for MMS FIELDS processing
% 	[ENVIR, MMS_CONST] = MMS_SDC_SDP_INIT(scNumber) returns environment
%   variables and constants useful for MMS processing. Input argument 
%   should be the sc number (as a string), i.e. '1' for mms1 and '2' for 
%   mms2 etc. It also configures logging to "$LOG_PATH_ROOT/ mmsX/ sdp/
%   date_IRFU.log".
%
%       The struct ENVIR will contain the following:
%	  .DATA_PATH_ROOT       - Root dir of data files
%	  .CDF_BASE             - Root dir of CDF tools
%	  .DROPBOX_ROOT         - Root dir of our output files (temporary location)
%	  .LOG_PATH_ROOT        - Root dir of log files
% 	  .CAL_PATH_ROOT        - Root dir of calibration files
%
%	The struct MMS_CONST will contain the following:
%	  .Version.X		- Major Software version used.
%		  .Y		    - Major Calibration version used.
%		  .Z		    - File version (should perhaps be removed).
%	  .Bitmask.OnlyDCE  - Only DCE was found at these points in time. 
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

% Spin rate max and min, nominally 3.0 rpm +/-0.2 rpm.
MMS_CONST.Spinrate.max = 3.2; % Rev per Minute.
MMS_CONST.Spinrate.min = 2.8; % Rev per Minute.

% Bitmask values; 2^(bit_number - 1):
MMS_CONST.Bitmask.SIGNAL_OFF               =  1;       % Bit 1
MMS_CONST.Bitmask.BAD_BIAS                 =  2;       % Bit 2
MMS_CONST.Bitmask.PROBE_SATURATION         =  4;       % Bit 3
MMS_CONST.Bitmask.LOW_DENSITY_SATURATION   =  8;       % Bit 4
MMS_CONST.Bitmask.SWEEP_DATA               =  16;      % Bit 5

% % DC V source bitmasks
% %for 16 ks/s channels, up to 6 channels at the same time:
% MMS_CONST.Source.SCM1 = 1;      % Bit 0x00 = SCM1  enable/disable
% MMS_CONST.Source.SCM1 = 2;      % Bit 0x01 = SCM2  enable/disable
% MMS_CONST.Source.SCM3 = 4;      % Bit 0x02 = SCM3  enable/disable
% MMS_CONST.Source.V1 = 8;        % Bit 0x03 = V1    enable/disable
% MMS_CONST.Source.V2 = 16;       % Bit 0x04 = V2    enable/disable
% MMS_CONST.Source.V3 = 32;       % Bit 0x05 = V3    enable/disable
% MMS_CONST.Source.V4 = 64;       % Bit 0x06 = V4    enable/disable
% MMS_CONST.Source.V5 = 128;      % Bit 0x07 = V5    enable/disable
% MMS_CONST.Source.V6 = 256;      % Bit 0x08 = V6    enable/disable
% MMS_CONST.Source.E12DC = 512;   % Bit 0x09 = E12DC enable/disable
% MMS_CONST.Source.E34DC = 1024;  % Bit 0x10 = E34DC enable/disable
% MMS_CONST.Source.E56DC = 2048;  % Bit 0x11 = E56DC enable/disable
% 
% % DC E source bitmasks
% %for 256 ks/s channels (ACE and High Speed Burst), up to 3 channels at 
% %the same time:
% MMS_CONST.Source.E12_AC = 1;    % Bit 0x00 = E12_AC enable/disable
% MMS_CONST.Source.E34_AC = 2;    % Bit 0x01 = E34_AC enable/disable
% MMS_CONST.Source.E56_AC = 4;    % Bit 0x02 = E56_AC enable/disable
% MMS_CONST.Source.V1_AC = 8;     % Bit 0x03 = V1_AC  enable/disable
% MMS_CONST.Source.V2_AC = 16:    % Bit 0x04 = V2_AC  enable/disable

ENVIR.CDF_BASE = getenv('CDF_BASE'); % Get path to CDF tools.
ENVIR.DATA_PATH_ROOT = getenv('DATA_PATH_ROOT'); % The final path of data.
ENVIR.LOG_PATH_ROOT = getenv('LOG_PATH_ROOT'); % Get path to logs.
ENVIR.DROPBOX_ROOT = getenv('DROPBOX_ROOT'); % Get path to output location,
% DROPBOX_ROOT is only the temporary location, file are then to be moved by
% other scripts to their final location when processing is done.
ENVIR.CAL_PATH_ROOT = getenv('CAL_PATH_ROOT'); % Get path to cal.

% Setup logging.
% Create a logfile at $LOG_PATH_ROOT / mmsX / sdp /
% named after current run day yyyymmdd and _IRFU.log. If this fails
% create it at $LOG_PATH_ROOT and include full date with seconds.
if( str2double(scNumber)>1 || str2double(scNumber)<4 )
    % Check to verify that output dir exists for logging. If not, create
    % it.
    if(~exist([ENVIR.LOG_PATH_ROOT, filesep, 'mms', scNumber, filesep, ...
            'sdp'], 'dir'))
        mkdir([ENVIR.LOG_PATH_ROOT, filesep, 'mms', scNumber], 'sdp');
    end
    irf.log('log_out', [ENVIR.LOG_PATH_ROOT, filesep, 'mms', ...
        scNumber, filesep, 'sdp', filesep, datestr(now,'yyyymmdd'),...
        '_IRFU.log']);
    mms_sdc_sdp_datamanager('init',str2double(scNumber))
    % Set log level
    irf.log('notice');
else
    irf.log('log_out', [ENVIR.LOG_PATH_ROOT, filesep, ...
        datestr(now,'yyyymmddTHHMMSS'), '_IRFU.log']);
    irf.log('debug');
    err_str = ['Matlab:MMS_SDC_SDP_INIT:InputArg scNumber incorrectly ',...
        'determined as: ',scNumber];
    irf.log('critical', err_str);
    error('Matlab:MMS_SDC_SDP_INIT',['MMS_SDC_SDP_INIT recieved an ', ...
        'unexpected sc number string. Input to processing should be ', ...
        'fullpath/filename.cdf according to MMS standard, mmsX_whatever',...
        ' where X = 1, 2, 3 or 4.']);
end

% Store information in the log file about which version of Matlab is used
% and which version of IRFU-MATLAB.
irf.log('notice', ['Matlab version used is: ' version]);
irf.log('notice', ['irfu-matlab version used is:', irf('version')]);
