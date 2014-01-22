function [ENVIR, MMS_CONST] = mms_init(remoteSDC)
% Function to return environment variables in a struct:
% 
% ENVIR
%      .DATA_PATH_ROOT  - Root dir of data files
%      .CDF_BASE        - Root dir of CDF tools
%      .DROPBOX_ROOT    - Root dir of our output files (temporary location)
%      .LOG_PATH_ROOT   - Root dir of log files
%      .CAL_PATH_ROOT   - Root dir of calibration files
% MMS_CONST
%      .Version.X       - Major Software version used.
%              .Y       - Major Calibration version used.
%              .Z       - File version (should perhaps be removed).
%      .Bitmask.OnlyDCE = 0x01 - Only DCE was found at these points in time. 
%
%        ? - Something more perhaps..
%
% Input: remoteSDC == 1, if run at SDC, == 0 if run locally.
%
% Date of latest modification: 2014/01/16
%


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

% FIXME: SET PATH TO IRFU CODE (add path ../irfu-matlab/...)

if(remoteSDC==1)
    ENVIR.CDF_BASE = getenv('CDF_BASE'); % get environment variable.
    ENVIR.DATA_PATH_ROOT = getenv('DATA_PATH_ROOT'); % Get path to data.
    ENVIR.LOG_PATH_ROOT = getenv('LOG_PATH_ROOT'); % Get path to logs.
    ENVIR.DROPBOX_ROOT = getenv('DROPBOX_ROOT'); % Get path to output location, (temporary location, other scripts then move it once fully written and our script is done). 
    ENVIR.CAL_PATH_ROOT = getenv('CAL_PATH_ROOT'); % Get path to cal.

    % Setup logging.
    % Create a logfile at log_path_root, named after current run date and IRFU.log
    irf.log('log_out',strcat(ENVIR.LOG_PATH_ROOT,'/',date,'_IRFU.log'));
    % Set log level to debug initially.
    irf.log('loggingLevel',4);

elseif(remoteSDC==0)
    % For local code running on IRFU computers
    ENVIR.CDF_BASE = '/home/thoni/MMS/cdf35_0-dist';
    ENVIR.DATA_PATH_ROOT = '/home/thoni/MMS/MMS_cdf';
    ENVIR.LOG_PATH_ROOT = '/home/thoni/MMS/log';
    ENVIR.DROPBOX_ROOT = '/home/thoni/MMS/fileoutput';
    %ENVIR.CAL_PATH_ROOT = getenv('CAL_PATH_ROOT');
    addpath(genpath('/home/thoni/MMS/irfu-matlab'));

    % Setup logging.
    % Create a logfile at log_path_root, named after current run date and IRFU.log
    irf.log('log_out',strcat(ENVIR.LOG_PATH_ROOT,'/',date,'_IRFU.log'));
    % Set log level to debug initially.
    irf.log('loggingLevel',4);

else
    error('Could not setup environment or log');
    %irf.log('warning','mms_init recived input. Was remoteSDC set?');
end