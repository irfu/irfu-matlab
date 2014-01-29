function mms_ql_dce(filename_dce_source_file, filename_dcv_source_file)
% This is the main file to be run at SDC, from this file all other substeps
% are performed.
% Input(s):
%    filename_dce_source_file = filename of DCE data cdf file, required.
%    filename_dcv_source_file = filename of DCV data cdf file. Not required
%    for basic processing but if not included it is assumed this file does
%    not yet exist and only partial processing can be done.
% Output(s):
%    cdf files created containing QuickLook DCE data.
%

% narginchk - Min 1 (dce), max 2 (dcv)
narginchk(1,2);

% Store runTime when script was called.
runTime = datestr(now,'yyyymmddHHMMSS');

% FIXME: Set to 0 if running locally at IRFU, set to 1 if running at SDC.
remoteRun = 1;

global ENVIR;
global MMS_CONST;

% ENVIR & MMS_CONST structs created by init script.
[ENVIR, MMS_CONST] = mms_init(remoteRun);



% If only one is found we cannot do all data processing
if(nargin==1)
    % Log message so we know we only got one input.
    irf.log('warning','mms_ql_dce recieved only one input argument. Can perform some but not all processing.');
    
    irf.log('debug',['mms_ql_dce trying mms_cdf_in_process on input file :', filename_dce_source_file]);

    [dce_source, dce_source_fileData] = mms_cdf_in_process(filename_dce_source_file,'sci');
    
    % Set bitmask for all times in dce_source.
    bitmask = mms_bitmasking(dce_source);
    
    % FIXME: Do some proper quality assesment.
    quality = bitmask; 
    
    % Get sunpulse data from the same interval, using dce_source_fileData
    % for start time, MMS S/C id, etc.
    %hk_sunpulse = mms_cdf_in_process(hk_fileSOMETHING,'ancillary');
    
    %FIXME Do some processing... Actually De-filter, De-spin, etc.
    
    %
    HeaderInfo = [];
    HeaderInfo.calledBy = 'ql'; % or 'sitl_dce' if running sitl instead of 'ql'.
    HeaderInfo.scId = dce_source_fileData.scId;
    HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
    HeaderInfo.dataMode = dce_source_fileData.dataMode;
    HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
    HeaderInfo.startTime = dce_source_fileData.startTime;
    HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
    HeaderInfo.numberOfSources = 1;
    HeaderInfo.parents_1 = filename_dce_source_file;
   
    irf.log('debug', 'mms_ql_dce trying mms_cdf_write');

    filename_output = mms_cdf_writing(dce_source, bitmask(:,2), HeaderInfo, quality(:,2));
    
    % Write out filename as an empty logfile so it can be easily found by
    % SDC scripts.  scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z_runTime.log
    
    unix(['touch',' ', ENVIR.LOG_PATH_ROOT,'/',filename_output,'_',runTime,'.log']);
    
elseif(nargin==2)
    % Log message so we know we got both.
    irf.log('notice','mms_ql_dce recieved two input arguments. Can perform full processing.');
    
    % First get dce data
    irf.log('debug',['mms_ql_dce trying mms_cdf_in_process on input file :', filename_dce_source_file]);

    [dce_source, dce_source_fileData] = mms_cdf_in_process(filename_dce_source_file,'sci');
    
    % Then get dcv data
    irf.log('debug',['mms_ql_dce trying mms_cdf_in_process on input file :', filename_dcv_source_file]);

    [dcv_source, dcv_source_fileData] = mms_cdf_in_process(filename_dcv_source_file,'sci');
    
    % Set bitmask for all times that dcv_source do not match dce_source,
    % priority goes to dce_source.
    bitmask = mms_bitmasking(dce_source, dcv_source);
    
    % FIXME: Do some proper quality assesment.
    quality = bitmask; 
    % Get sunpulse data from the same interval, using dce_source_fileData
    % for start time, MMS S/C id, etc.
    %hk_sunpulse = mms_cdf_in_process(hk_fileSOMETHING,'ancillary');
    
    %FIXME Do some processing... Actually De-filter, De-spin, etc.
    
    
    HeaderInfo = [];
    HeaderInfo.calledBy = 'ql';
    HeaderInfo.scId = dce_source_fileData.scId;
    HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
    HeaderInfo.dataMode = dce_source_fileData.dataMode;
    HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
    HeaderInfo.startTime = dce_source_fileData.startTime;
    HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
    HeaderInfo.numberOfSources = 2;
    HeaderInfo.parents_1 = filename_dce_source_file;
    HeaderInfo.parents_2 = filename_dcv_source_file;
   
    irf.log('debug', 'mms_ql_dce trying mms_cdf_write');
    filename_output = mms_cdf_writing(dce_source, bitmask(:,2), HeaderInfo, quality(:,2));
    
    % Write out filename as empty logfile so it can be easily found by SDC
    % scripts.
    unix(['touch',' ', ENVIR.LOG_PATH_ROOT,'/',filename_output,'_',runTime,'.log']);
    
elseif(nargin>2)
    % Log message so we know it went wrong... Should not happen as
    % narginchk(1,2) has check it. But if we add more input arguments
    % later..
    irf.log('warning','mms_ql_dce received more then two input. What is what?');
end
