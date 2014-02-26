function mms_ql_dce(filename_dce_source_file, filename_dcv_source_file)
% MMS_QL_DCE start point and main function for MMS SDC QuickLook DCE processing.
%	MMS_QL_DCE(filename_dce_source_file, filename_dcv_source_file) takes input fullpath filenames 
%	of one DCE and one DCV file for corresponding time period and runs processing to determine 
%	electric field for MMS SITL.
%
% 	MMS_QL_DCE(filename_dce_source_file, filename_dcv_source_file) using various subroutines 
%	the input files are read, processed and then the final output is written
%	to a corresponding CDF file in accordiance with MMS CDF Format Guide and MMS SDC Developer Guide.
%	Bitmask information is added to the output file to mark if the two datasources do not fully overlap
%	in time. And compared with mms_sitl_dce a quality flag is also added.
%	MMS_QL_DCE can perform some of the processsing if only provided with the DCE file but not all.
%
%	Example:
%		mms_ql_dce('/full/path/to/source_dce_file.cdf');
%		mms_ql_dce('/full/path/to/source_dce_file.cdf', '/full/path/to/source_dcv_file.cdf');
%
% 	See also MMS_INIT, MMS_CDF_IN_PROCESS. MMS_CDF_WRITING, MMS_SITL_DCE.

% narginchk - Min 1 (dce), max 2 (dcv)
narginchk(1,2);

% Store runTime when script was called.
runTime = datestr(now,'yyyymmddHHMMSS');

global ENVIR;
global MMS_CONST;

% Quickly try to find which spacecraft, in order to place log correctly.
[pathFile, nameFile, extensionFile] = fileparts(filename_dce_source_file);
HeaderInfo = [];
HeaderInfo.numberStr = nameFile(4);
clear pathFile nameFile extensionFile

% ENVIR & MMS_CONST structs created by init script.
[ENVIR, MMS_CONST] = mms_init(HeaderInfo.numberStr);


% If only one is found we cannot do all data processing
if(nargin==1)
    % Log message so we know we only got one input.
    irf.log('warning','mms_ql_dce received only one input argument. Can perform some but not all processing.');
    
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
    
    HeaderInfo.calledBy = 'ql'; % or 'sitl_dce' if running sitl instead of 'ql'.
    HeaderInfo.scId = dce_source_fileData.scId;
    HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
    HeaderInfo.dataMode = dce_source_fileData.dataMode;
    HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
    HeaderInfo.startTime = dce_source_fileData.startTime;
    HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
    HeaderInfo.numberOfSources = 1;
    HeaderInfo.parents_1 = dce_source_fileData.filename;
   
    irf.log('debug', 'mms_ql_dce trying mms_cdf_write');

    filename_output = mms_cdf_writing(dce_source, bitmask(:,2), HeaderInfo, quality(:,2));
    
elseif(nargin==2)
    % Log message so we know we got both.
    irf.log('notice','mms_ql_dce received two input arguments. Can perform full processing.');
    
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
    
    
    HeaderInfo.calledBy = 'ql';
    HeaderInfo.scId = dce_source_fileData.scId;
    HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
    HeaderInfo.dataMode = dce_source_fileData.dataMode;
    HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
    HeaderInfo.startTime = dce_source_fileData.startTime;
    HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
    HeaderInfo.numberOfSources = 2;
    HeaderInfo.parents_1 = dce_source_fileData.filename;
    HeaderInfo.parents_2 = dcv_source_fileData.filename;
   
    irf.log('debug', 'mms_ql_dce trying mms_cdf_write');
    filename_output = mms_cdf_writing(dce_source, bitmask(:,2), HeaderInfo, quality(:,2));
    
elseif(nargin>2)
    % Log message so we know it went wrong... Should not happen as
    % narginchk(1,2) has check it. But if we add more input arguments
    % later..
    irf.log('warning','mms_ql_dce received more then two input. What is what?');
end


% Write out filename as empty logfile so it can be easily found by SDC
% scripts.
unix(['touch',' ', ENVIR.LOG_PATH_ROOT,'/mms',HeaderInfo.numberStr,'/sdp/',filename_output,'_',runTime,'.log']);
