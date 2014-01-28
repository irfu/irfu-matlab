function filename_output = mms_usc(filename_dcv_source_file)
% This is the main file to be run at SDC, from this file all other substeps
% are performed.
% Input(s):
%    filename_dcv_source_file = filename of DCV data cdf file. Not required
%    for basic processing but if not included it is assumed this file does
%    not yet exist and only partial processing can be done.
% Output(s):
%    cdf files created containing Usc potential DCV data.
%

% narginchk - Min 1 (dcv), max 1 (dcv)
narginchk(1,1);

% Store runTime when script was called.
runTime = datestr(now,'yyyymmddHHMMSS');

% FIXME: Set to 0 if running locally at IRFU, set to 1 if running at SDC.
remoteRun =1;

global ENVIR;
global MMS_CONST;

% ENVIR & MMS_CONST structs created by init script.
[ENVIR, MMS_CONST] = mms_init(remoteRun);

    irf.log('debug',['mms_usc trying mms_cdf_in_process on input file: ', filename_dcv_source_file]);
    [dcv_source, dcv_source_fileData] = mms_cdf_in_process(filename_dcv_source_file,'sci');
    
    % Set bitmask for all times in dce_source.
    bitmask = mms_bitmasking(dcv_source);
    
    % FIXME: Do some proper quality assesment.
    quality = bitmask; 
    
    % Get sunpulse data from the same interval, using dce_source_fileData
    % for start time, MMS S/C id, etc.
    %hk_sunpulse = mms_cdf_in_process(hk_fileSOMETHING,'ancillary');
    
    %FIXME Do some processing... Actually De-filter, De-spin, etc.
    
    %
    HeaderInfo = [];
    HeaderInfo.calledBy = 'usc'; % or 'sitl_dce' if running sitl instead of 'ql'.
    HeaderInfo.scId = dcv_source_fileData.scId;
    HeaderInfo.instrumentId = dcv_source_fileData.instrumentId;
    HeaderInfo.dataMode = dcv_source_fileData.dataMode;
    HeaderInfo.dataLevel = dcv_source_fileData.dataLevel;
    HeaderInfo.startTime = dcv_source_fileData.startTime;
    HeaderInfo.vXYZ = dcv_source_fileData.vXYZ;
    HeaderInfo.numberOfSources = 1;
    HeaderInfo.parents_1 = filename_dcv_source_file;
   
    irf.log('debug', 'mms_usc trying mms_cdf_write');
    filename_output = mms_cdf_writing(dcv_source, bitmask(:,2), HeaderInfo, quality(:,2));
    
    
    % Write out filename as an empty logfile so it can be easily found by
    % SDC scripts.  scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z_runTime.log
    
    unix(['touch',' ', ENVIR.LOG_PATH_ROOT,'/',filename_output,'_',runTime,'.log']);
    
