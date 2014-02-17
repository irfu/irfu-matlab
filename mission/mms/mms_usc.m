function mms_usc(filename_dcv_source_file)
% MMS_USC startpoint and main function for MMS SDC USC processing. 
%	MMS_USC(filename_dcv_source_file) takes input fullpath filename of DCV file and runs processing to determine SC potential for MMS.
%
% 	MMS_USC(filename_dcv_source_file) using various subroutines the input file is read, processed and then the final output is written
%	to a corresponding CDF file in accordiance with MMS CDF Format Guide and MMS SDC Developer Guide.
%
%	Example:
%		mms_usc('/full/path/to/source_dcv_file.cdf');
%
% 	See also MMS_INIT, MMS_CDF_IN_PROCESS. MMS_CDF_WRITE.


% narginchk - Min 1 (dcv), max 1 (dcv)
narginchk(1,1);

% Store runTime when script was called.
runTime = datestr(now,'yyyymmddHHMMSS');

global ENVIR;
global MMS_CONST;

% ENVIR & MMS_CONST structs created by init script.
[ENVIR, MMS_CONST] = mms_init();

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
    HeaderInfo.parents_1 = dcv_source_fileData.filename;
   
    irf.log('debug', 'mms_usc trying mms_cdf_write');
    filename_output = mms_cdf_writing(dcv_source, bitmask(:,2), HeaderInfo, quality(:,2));
    
    
    % Write out filename as an empty logfile so it can be easily found by
    % SDC scripts.  scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z_runTime.log
    
    unix(['touch',' ', ENVIR.LOG_PATH_ROOT,'/',filename_output,'_',runTime,'.log']);
    
