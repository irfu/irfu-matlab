function ENVIR = mms_sdc_sdp_init
% MMS_SDC_SDP_INIT  initialize environment for MMS FIELDS processing
%
% 	ENVIR = MMS_SDC_SDP_INIT returns environment variables and constants
%           useful for MMS processing.
%
%       The struct ENVIR will contain the following:
%	  .DATA_PATH_ROOT       - Root dir of data files
%	  .CDF_BASE             - Root dir of CDF tools
%	  .DROPBOX_ROOT         - Root dir of our output files (temporary location)
%	  .LOG_PATH_ROOT        - Root dir of log files
%	  .CAL_PATH_ROOT        - Root dir of calibration files
%
%	Example:
%		ENVIR = MMS_SDC_SDP_INIT;
%
% 	See also MMS_CONSTANTS, MMS_SDC_SDP_LOG_INIT.

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

ENVIR.CDF_BASE = getenv('CDF_BASE'); % Get path to CDF tools.
ENVIR.DATA_PATH_ROOT = getenv('DATA_PATH_ROOT'); % The final path of data.
ENVIR.LOG_PATH_ROOT = getenv('LOG_PATH_ROOT'); % Get path to logs.
ENVIR.DROPBOX_ROOT = getenv('DROPBOX_ROOT'); % Get path to output location,
% DROPBOX_ROOT is only the temporary location, file are then to be moved by
% other scripts to their final location when processing is done.
ENVIR.CAL_PATH_ROOT = getenv('CAL_PATH_ROOT'); % Get path to cal.

end