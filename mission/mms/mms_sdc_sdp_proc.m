function mms_sdc_sdp_proc( varargin )
% MMS_SDC_SDP_PROC main starting point for MMS SDC processing.
%	MMS_SDC_SDP_PROC('processingType', '/pathTo/input/file1.cdf', ...
%   '/pathTo/input/file2.cdf', ...); will start the MMS processing for
%   processing types "ql", "sitl" and "usc" which is to be performed at
%   SDC. If not all needed cdf files are provided, but all required cdf 
%   required cdf files are only some parts of the processing can occur. 
%
%	Example:
%   Full quicklook ("ql") processing, with all needed cdf files as input,
%		mms_sdc_sdp_proc('ql',...
%         '/path/mms2_sdp_fast_dce_20150410_v0.0.1.cdf, ...
%         '/path/mms2_sdp_fast_dcv_20150410_v0.0.0.cdf', ...
%         '/path/mm2_fields_hk_101_20150410_v0.0.2.cdf');
%	Only partial scientist in the loop ("sitl") processing, with only the
%	required cdf files as input, (note for full sitl include DCV).
%       mms_sdc_sdp_proc('sitl',...
%         '/path/mms2_sdp_fast_dce_20150410_v0.0.1.cdf, ...
%         '/path/mm2_fields_hk_101_20150410_v0.0.2.cdf');
%
% 	See also MMS_SDC_SDP_INIT, MMS_SDC_SDP_BITMASKING.


% Store runTime when script was called.
runTime = datestr(now,'yyyymmddHHMMSS'); % Use this time to ease for SDC to
% find output file created. An empty file is to be created with its
% filename being the same as that of the dataproduct created appended with
% _runTime.txt and placed in ENVIR.LOG_PATH_ROOT.

global ENVIR;
global MMS_CONST;

% This is a new test for combined data processing.
narginchk(3,4); % FIXME: correct number of required input arguments.

% First argument is and should always be which mode to run.
% QuickLook, Usc or SITL. (Is static in the bash script that starts Matlab)
if( ischar( varargin{1} ) )
    HeaderInfo = [];
    HeaderInfo.calledBy = lower(varargin{1});
    if( any( [strcmp(HeaderInfo.calledBy, 'usc'), strcmp(HeaderInfo.calledBy, 'ql'), strcmp(HeaderInfo.calledBy, 'sitl')] ) )
        % Ok
    else
        error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC valid values for first argument is: "ql", "sitl" or "usc".');
    end
else
    error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC first argument must be a string identifying which process to run. Valid values are: "ql", "sitl" or "usc"');
end

for i=2:nargin
    % Go through each input argument and find if it is a HK_SunpulseFile,
    % DCV_File or DCE_File.
    if( ischar( varargin{i} ) )
        % Ok, it is a string. Lets check which string it is.
        
        [pathIn, fileIn, extIn] = fileparts(varargin{i});
        if( any( [isempty(pathIn), isempty(fileIn), isempty(extIn)] ) )
            error('Matlab:MMS_SDC_SDP_PROC',['MMS_SDC_SDP_PROC input arguments except first is supposed to be a string for cdf files with full path names. Was read as: ', varargin{i}]);
        elseif(i==2)
            % Setup log and environment.
            HeaderInfo.numberStr = fileIn(4); % Store S/C Number X from mmsX_... as a string in HeaderInfo
            % If fileIn(4) is not as expected {1, 2, 3, or 4} log file will
            % be placed in LOG_PATH_ROOT and an error will be issued.
            [ENVIR, MMS_CONST] = mms_sdc_sdp_init(HeaderInfo.numberStr);
            irf.log('notice',['MMS_SDC_SDP_PROC was called with first argument HeaderInfo.calledBy as: ', HeaderInfo.calledBy, '.']);
        elseif(i>2)
            if(~strcmp(HeaderInfo.numberStr, fileIn(4)))
                irf.log('critical',['MMS_SDC_SDP_PROC was called with one file from SC number ', HeaderInfo.numberStr, ' and one other file from SC number ', fileIn(4),'.' ]);
                irf.log('warning',['MMS_SDC_SDP_PROC file argument ', varargin{i},' does not match expected sc numbering. Aborting with error.']);
                error('Matlab:MMS_SDC_SDP_PROC:Input',['MMS_SDC_SDP_PROC was called with one file from SC number ', HeaderInfo.numberStr, ' and one other file from SC number ', fileIn(4),'.']);
            end
        end

        % Is it HK_Sunpulse cdf file
        % Files are named as: mmsX_fields_hk_l1b_101_20150410_v0.0.1.cdf
        sunHKFilePat = '_101_';
        % Files are named as: mmsX_..... v0.0.1.cdf
        dcvFilePat = '_dcv_';
        % Filess are named as: mmsX_.... v0.1.2.cdf
        dceFilePat = '_dce_';
        
        if( ~isempty( regexpi(fileIn, sunHKFilePat) ) )
            if(exist('HK_SunpulseFile','var'))
                irf.log('critical',['Several input argument files identified as HK_Sunpulse files. This one as well: ', varargin{i}]);
                error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC received two or more HK_Sunpulse files.');
            else
                % It is the HK_Sunpulse file
                HK_SunpulseFile = varargin{i};
                irf.log('notice',['HK_Sunpulse file identified as: ', HK_SunpulseFile]);
            end
        elseif( ~isempty( regexpi(fileIn, dcvFilePat) ) )
            if(exist('DCV_File','var'))
                irf.log('critical',['Several input argument files identified as DCV files. This one as well: ', varargin{i}]);
                error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC received two or more DCV files.');
            else
                % It is the DCV file
                DCV_File = varargin{i};
                irf.log('notice',['DCV file identified as: ', DCV_File]);
            end
        elseif( ~isempty( regexpi(fileIn, dceFilePat) ) )
            if(exist('DCE_File','var'))
                irf.log('critical',['Several input argument files identified as DCE files. This one as well: ', varargin{i}]);
                error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC received two or more DCE files.');
            else
                % It is the DCE file
                DCE_File = varargin{i};
                irf.log('notice',['DCE file identified as: ', DCE_File]);
            end
            
        else
            % Unidentified input argument
            irf.log('warning',['MMS_SDC_SDP_PROC received one input it could not identify: ', varargin{i}]);
        end
    else
        error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC input arguments must be strings.');
    end
end

% All input arguments read. All files required identified correct?
if( ~all( [exist('HK_SunpulseFile', 'var'), exist('DCE_File','var'), exist('DCV_File','var')] ) )
    irf.log('warning','MMS_SDC_SDP_PROC not all input variables were provided and identified.');
    for i=1:nargin
        irf.log('warning',['MMS_SDC_SDP_PROC received input argument: ', varargin{i}]);
    end
end


% Begin actual processing for Usc or QL or SITL.
switch(HeaderInfo.calledBy)
    
    case('usc')
        
        if( all( [exist('DCV_File', 'var'), exist('HK_SunpulseFile', 'var')] ) )
            
            irf.log('debug',['MMS_SDC_SDP_PROC Usc using mms_sdc_sdp_cdf_in_process on input file: ', DCV_File]);
            [dcv_source, dcv_source_fileData] = mms_sdc_sdp_cdf_in_process(DCV_File,'sci');

            irf.log('debug',['MMS_SDC_SDP_PROC Usc using dataobj on input file: ', HK_SunpulseFile]);
            hk_sunpulsecdf = dataobj(HK_SunpulseFile,'tint',0,'true');

            % Set bitmask for all times in dce_source.
            bitmask = mms_sdc_sdp_bitmasking(dcv_source);

            % FIXME: Do some proper quality assesment.
            quality = bitmask; 

            % Get sunpulse data from the same interval, using dce_source_fileData
            % for start time, MMS S/C id, etc.
            %hk_sunpulse = mms_sdc_sdp_cdf_in_process(hk_fileSOMETHING,'ancillary');

            %FIXME Do some processing... Actually De-filter, De-spin, etc.

            HeaderInfo.scId = dcv_source_fileData.scId;
            HeaderInfo.instrumentId = dcv_source_fileData.instrumentId;
            HeaderInfo.dataMode = dcv_source_fileData.dataMode;
            HeaderInfo.dataLevel = dcv_source_fileData.dataLevel;
            HeaderInfo.startTime = dcv_source_fileData.startTime;
            HeaderInfo.vXYZ = dcv_source_fileData.vXYZ;
            HeaderInfo.numberOfSources = 1;
            HeaderInfo.parents_1 = dcv_source_fileData.filename;

            irf.log('debug', 'MMS_SDC_SDP_PROC Usc using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(dcv_source, bitmask(:,2), HeaderInfo, quality(:,2));
        
        else
            
            irf.log('critical','MMS_SDC_SDP_PROC USC received some unclear input arguments. DCV and HK_Sunpulse CDF files are required,');
            irf.log('warning',['Variable HK_SunpulseFile exists: ', num2str(exist('HK_SunpulseFile','var')), ...
                ', variable DCE_File exists: ', num2str(exist('DCE_File','var')), ...
                ', variable DCV_File exists: ', num2str(exist('DCV_File','var'))]);
            error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC USC did not recieve all required input file arguments or was unsuccesful in identifying them. Please see log file.');
        
        end
    case('sitl')
   
        % Check if all required and needed files are sent as input and have
        % been identified properly.
        if( all( [exist('DCE_File', 'var'), exist('HK_SunpulseFile', 'var'), ~exist('DCV_File', 'var')] ) )
            
            % Log message so we know we are missing one input.
            irf.log('warning','MMS_SDC_SDP_PROC SITL received DCE and HK_Sunpulse but no DCV file argument. Can perform some but not all processing.');

            irf.log('debug',['MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_in_process on input file: ', DCE_File]);
            [dce_source, dce_source_fileData] = mms_sdc_sdp_cdf_in_process(DCE_File,'sci');

            irf.log('debug',['MMS_SDC_SDP_PROC SITL using dataobj on input file: ', HK_SunpulseFile]);
            hk_sunpulsecdf = dataobj(HK_SunpulseFile,'tint',0,'true');
            
            % Set bitmask for all times in dce_source.
            bitmask = mms_sdc_sdp_bitmasking(dce_source);
            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.

            HeaderInfo.scId = dce_source_fileData.scId;
            HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
            HeaderInfo.dataMode = dce_source_fileData.dataMode;
            HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
            HeaderInfo.startTime = dce_source_fileData.startTime;
            HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
            HeaderInfo.numberOfSources = 1;
            HeaderInfo.parents_1 = dce_source_fileData.filename;

            irf.log('debug', 'MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(dce_source, bitmask(:,2), HeaderInfo);

        elseif( all( [exist('HK_SunpulseFile','var'), exist('DCE_File','var'), exist('DCV_File','var')] ) )
            
            % Log message so we know we got both.
            irf.log('notice','MMS_SDC_SDP_PROC SITL received all expected input arguments, DCE, DCV and HK_Sunpulse file arguments. Can perform full processing.');

            % First get dce data
            irf.log('debug',['MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_in_process on input file: ', DCE_File]);
            [dce_source, dce_source_fileData] = mms_sdc_sdp_cdf_in_process(DCE_File,'sci');

            % Then get dcv data
            irf.log('debug',['MMS_SDC_SDP_PROC SITL trying mms_sdc_sdp_cdf_in_process on input file: ', DCV_File]);
            [dcv_source, dcv_source_fileData] = mms_sdc_sdp_cdf_in_process(DCV_File,'sci');

            irf.log('debug',['MMS_SDC_SDP_PROC SITL using dataobj on input file: ', HK_SunpulseFile]);
            hk_sunpulsecdf = dataobj(HK_SunpulseFile,'tint',0,'true');
            
            % Set bitmask for all times that dcv_source do not match dce_source,
            % priority goes to dce_source.
            bitmask = mms_sdc_sdp_bitmasking(dce_source, dcv_source);

            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.
          
            
            HeaderInfo.scId = dce_source_fileData.scId;
            HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
            HeaderInfo.dataMode = dce_source_fileData.dataMode;
            HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
            HeaderInfo.startTime = dce_source_fileData.startTime;
            HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
            HeaderInfo.numberOfSources = 2;
            HeaderInfo.parents_1 = dce_source_fileData.filename;
            HeaderInfo.parents_2 = dcv_source_fileData.filename;

            irf.log('debug', 'MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(dce_source, bitmask(:,2), HeaderInfo);

        else
            
            irf.log('critical','MMS_SDC_SDP_PROC SITL received some unclear input arguments. DCE and HK_Sunpulse CDF files are required, DCV needed for some but not all processing.');
            irf.log('warning',['Variable HK_SunpulseFile exists: ', num2str(exist('HK_SunpulseFile','var')), ...
                ', variable DCE_File exists: ', num2str(exist('DCE_File','var')), ...
                ', variable DCV_File exists: ', num2str(exist('DCV_File','var'))]);
            error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC SITL did not recieve all required input file arguments or was unsuccesful in identifying them. Please see log file.');
        
        end
        
        
    case('ql')

        % Check if all required and needed files are sent as input and have
        % been identified properly.
        if( all( [exist('DCE_File', 'var'), exist('HK_SunpulseFile', 'var'), ~exist('DCV_File', 'var')] ) )
            
            % Log message so we know we are missing one input.
            irf.log('warning','MMS_SDC_SDP_PROC QL received DCE and HK_Sunpulse but no DCV file argument. Can perform some but not all processing.');

            irf.log('debug',['MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_in_process on input file: ', DCE_File]);
            [dce_source, dce_source_fileData] = mms_sdc_sdp_cdf_in_process(DCE_File,'sci');

            irf.log('debug',['MMS_SDC_SDP_PROC QL using dataobj on input file: ', HK_SunpulseFile]);
            hk_sunpulsecdf = dataobj(HK_SunpulseFile,'tint',0,'true');
            
            % Set bitmask for all times in dce_source.
            bitmask = mms_sdc_sdp_bitmasking(dce_source);
            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.

            % FIXME: Do some proper quality assesment.
            quality = bitmask; 

            HeaderInfo.scId = dce_source_fileData.scId;
            HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
            HeaderInfo.dataMode = dce_source_fileData.dataMode;
            HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
            HeaderInfo.startTime = dce_source_fileData.startTime;
            HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
            HeaderInfo.numberOfSources = 1;
            HeaderInfo.parents_1 = dce_source_fileData.filename;

            irf.log('debug', 'MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(dce_source, bitmask(:,2), HeaderInfo, quality(:,2));

        elseif( all( [exist('HK_SunpulseFile','var'), exist('DCE_File','var'), exist('DCV_File','var')] ) )
            
            % Log message so we know we got both.
            irf.log('notice','MMS_SDC_SDP_PROC QL received all expected input arguments, DCE, DCV and HK_Sunpulse file arguments. Can perform full processing.');

            % First get dce data
            irf.log('debug',['MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_in_process on input file: ', DCE_File]);
            [dce_source, dce_source_fileData] = mms_sdc_sdp_cdf_in_process(DCE_File,'sci');

            % Then get dcv data
            irf.log('debug',['MMS_SDC_SDP_PROC QL trying mms_sdc_sdp_cdf_in_process on input file: ', DCV_File]);
            [dcv_source, dcv_source_fileData] = mms_sdc_sdp_cdf_in_process(DCV_File,'sci');

            irf.log('debug',['MMS_SDC_SDP_PROC QL using dataobj on input file: ', HK_SunpulseFile]);
            hk_sunpulsecdf = dataobj(HK_SunpulseFile,'tint',0,'true');
            
            % Set bitmask for all times that dcv_source do not match dce_source,
            % priority goes to dce_source.
            bitmask = mms_sdc_sdp_bitmasking(dce_source, dcv_source);

            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.

            % FIXME: Do some proper quality assesment.
            quality = bitmask; 
            
            
            HeaderInfo.scId = dce_source_fileData.scId;
            HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
            HeaderInfo.dataMode = dce_source_fileData.dataMode;
            HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
            HeaderInfo.startTime = dce_source_fileData.startTime;
            HeaderInfo.vXYZ = dce_source_fileData.vXYZ;
            HeaderInfo.numberOfSources = 2;
            HeaderInfo.parents_1 = dce_source_fileData.filename;
            HeaderInfo.parents_2 = dcv_source_fileData.filename;

            irf.log('debug', 'MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(dce_source, bitmask(:,2), HeaderInfo, quality(:,2));

        else
            
            irf.log('critical','MMS_SDC_SDP_PROC QL received some unclear input arguments. DCE and HK_Sunpulse CDF files are required, DCV needed for some but not all processing.');
            irf.log('warning',['Variable HK_SunpulseFile exists: ', num2str(exist('HK_SunpulseFile','var')), ...
                ', variable DCE_File exists: ', num2str(exist('DCE_File','var')), ...
                ', variable DCV_File exists: ', num2str(exist('DCV_File','var'))]);
            error('Matlab:MMS_SDC_SDP_PROC:Input','MMS_SDC_SDP_PROC QL did not recieve all required input file arguments or was unsuccesful in identifying them. Please see log file.');
        
        end

        
    otherwise
        
        % This should not happen as we have checked for HeaderInfo.calledBy before.
        % But nevertheless.
        irf.log('critical',['MMS_SDC_SDP_PROC recieved some unknown HeaderInfo.calledBy: ', HeaderInfo.calledBy]);
        error('Matlab:MMS_SDC_SDP_PROC:Input',['MMS_SDC_SDP_PROC recieved some unknown HeaderInfo.calledBy: ', HeaderInfo.calledBy]);

end


% Write out filename as empty logfile so it can be easily found by SDC
% scripts.
unix(['touch',' ', ENVIR.LOG_PATH_ROOT,'/mms',HeaderInfo.numberStr,'/sdp/',filename_output,'_',runTime,'.log']);

end