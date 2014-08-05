function mms_sdc_sdp_proc( procName, varargin)
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

narginchk(4,4);

init_matlab_path()

HK_101_File = '';
DCV_File = '';
DCE_File = '';

% First argument is and should always be which mode to run.
% QuickLook, Usc or SITL. (Is static in the bash script that starts Matlab)
if ~ischar(procName)
    error('Matlab:MMS_SDC_SDP_PROC:Input', ...
    'MMS_SDC_SDP_PROC first argument must be a string');
end
procName = lower(procName);
HeaderInfo = []; HeaderInfo.calledBy = procName;
if isempty(intersect(procName,{'usc','ql','sitl'}))
    error('Matlab:MMS_SDC_SDP_PROC:Input', ...
    'MMS_SDC_SDP_PROC first argument must be one of: "ql", "sitl" or "usc"');
end

irf.log('notice', ['MMS_SDC_SDP_PROC process name: ', procName]);
    
for i=1:nargin-1
    if isempty(varargin{i}), continue, end
    
    % Go through each input argument and find if it is a HK_101_File,
    % DCV_File or DCE_File.
    if ~ischar(varargin{i})
      error('Matlab:MMS_SDC_SDP_PROC:Input', ...
          'MMS_SDC_SDP_PROC input arguments must be strings.');
    end
    
    [pathIn, fileIn, extIn] = fileparts(varargin{i});
    if any(isempty([pathIn, fileIn, extIn]))
        error('Matlab:MMS_SDC_SDP_PROC', ...
            ['MMS_SDC_SDP_PROC expecting cdf file (full path), got: ', ...
            varargin{i}]);
    end
      
    if i==1
        % Setup log and environment.
        numberStr = fileIn(4);
        HeaderInfo.numberStr = numberStr; % Store S/C Number X from mmsX_... as a string in HeaderInfo
        % If fileIn(4) is not as expected {1, 2, 3, or 4} log file will
        % be placed in LOG_PATH_ROOT and an error will be issued.
        [ENVIR, MMS_CONST] = mms_sdc_sdp_init(HeaderInfo.numberStr);
    else
        if ~strcmp(HeaderInfo.numberStr, fileIn(4))
            err_str = ['MMS_SDC_SDP_PROC was called with one file from SC number ',...
                HeaderInfo.numberStr, ' and one other file from SC number ', ...
                fileIn(4), '.'];
            irf.log('critical', err_str);
            irf.log('critical', ...
                ['MMS_SDC_SDP_PROC file argument ', varargin{i}, ...
                ' does not match expected sc numbering. Aborting with error.']);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
    end

    if regexpi(fileIn, '_101_') % 101, mmsX_fields_hk_l1b_101_20150410_v0.0.1.cdf
        if ~isempty(HK_101_File)
            err_str = ['Received multiple HK_101 files in input (',...
                HK_101_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        % It is the HK_101 file
        HK_101_File = varargin{i};
        irf.log('notice', ['HK_101 file identified as: ', ...
            HK_101_File]);
        
    elseif regexpi(fileIn, '_dcv_') %DCV
        if ~isempty(DCV_File)
            err_str = ['Received multiple DC V files in input (',...
                DCV_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        DCV_File = varargin{i};
        irf.log('notice', ['DCV file identified as: ', DCV_File]);
        
    elseif regexpi(fileIn, '_dce_') % DCE
        if ~isempty(DCE_File)
            err_str = ['Received multiple DC E files in input (',...
                DCE_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        DCE_File = varargin{i};
        irf.log('notice', ['DCE file identified as: ', DCE_File]);
        
    else
        % Unidentified input argument
        err_str = ['MMS_SDC_SDP_PROC unrecognized input file: ',...
            varargin{i}];
        irf.log('critical', err_str);
        error('Matlab:MMS_SDC_SDP_PROC:Input', err_str)
    end
end

% All input arguments read. All files required identified correct?
if any(isempty([HK_101_File, DCE_File, DCV_File]))
    irf.log('warning', 'MMS_SDC_SDP_PROC missing some input.');
    for i=1:nargin-1
        irf.log('warning',...
            ['MMS_SDC_SDP_PROC received input argument: ', varargin{i}]);
    end
end

if isempty(HK_101_File)
    errMsg = 'MMS_SDC_SDP_PROC missing required input: HK_101_File';
    irf.log('critical', errMsg);
    error('Matlab:MMS_SDC_SDP_PROC:Input', errMsg);
end

% Begin actual processing for Usc or QL or SITL.
switch(procName)
    case('usc')
        if isempty(DCV_File)
            errMsg = 'MMS_SDC_SDP_PROC missing required input for USC: DCV';
            irf.log('critical', errMsg);
            error('Matlab:MMS_SDC_SDP_PROC:Input', errMsg);
        end

        irf.log('debug', ['MMS_SDC_SDP_PROC input DCV file: ', DCV_File]);
        dcv_source_fileData = mms_sdc_sdp_cdf_in_process(DCV_File, 'sci', 'dcv');

        irf.log('debug',['MMS_SDC_SDP_PROC input HK_101 file: ',...
            HK_101_File]);
        mms_sdc_sdp_cdf_in_process(HK_101_File, 'sci', 'hk_101');

        
        % Set bitmask
        bitmask = mms_sdc_sdp_bitmasking(mms_sdc_sdp_datamanager('dcvtime'));

        % FIXME: Do some proper quality assesment.
        quality = bitmask;


        %FIXME Do some processing... Actually De-filter, De-spin, etc.

        %copy_header(1)
        HeaderInfo.scId = dcv_source_fileData.scId;
        HeaderInfo.instrumentId = dcv_source_fileData.instrumentId;
        HeaderInfo.dataMode = dcv_source_fileData.dataMode;
        HeaderInfo.dataLevel = dcv_source_fileData.dataLevel;
        HeaderInfo.startTime = dcv_source_fileData.startTime;
        HeaderInfo.numberOfSources = 1;
        HeaderInfo.parents_1 = dcv_source_fileData.filename;
    
        irf.log('debug', ...
            'MMS_SDC_SDP_PROC Usc using mms_sdc_sdp_cdf_writing');
        filename_output = mms_sdc_sdp_cdf_writing(bitmask(:,2), ...
          HeaderInfo, quality(:,2));
  
    case('sitl')

        % Check if all required and needed files are sent as input and have
        % been identified properly.
        if( all( [~isempty(DCE_File), ~isempty(HK_101_File), ...
                isempty(DCV_File)] ) )
            
            % Log message so we know we are missing one input.
            irf.log('warning', ...
                'MMS_SDC_SDP_PROC SITL received DCE and HK_101 but no DCV file argument. Can perform some but not all processing.');

            irf.log('debug',...
                ['MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_in_process on input file: ',...
                DCE_File]);
            dce_source_fileData = mms_sdc_sdp_cdf_in_process(DCE_File, ...
                'sci', 'dce');
            
            % Then the sunpulse file
            irf.log('debug', ...
                ['MMS_SDC_SDP_PROC SITL trying mms_sdc_sdp_cdf_in_process on input file: ',...
                HK_101_File]);
            mms_sdc_sdp_cdf_in_process(HK_101_File, 'sci', 'hk_101');
            
            irf.log('debug', 'MMS_SDC_SDP_PROC Request the phase');
            dcephase = mms_sdc_sdp_datamanager('dcephase');
            
            % Set bitmask
            bitmask = mms_sdc_sdp_bitmasking(mms_sdc_sdp_datamanager('dcetime'));
            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.

            copy_header(1)
            irf.log('debug', ...
                'MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(bitmask(:,2), HeaderInfo);

        elseif( all( [~isempty(HK_101_File), ...
                ~isempty(DCE_File), ~isempty(DCV_File)] ) )
            
            % Log message so we know we got both.
            irf.log('notice', 'MMS_SDC_SDP_PROC SITL received all expected input arguments, DCE, DCV and HK_101 file arguments. Can perform full processing.');

            % First get dce data
            irf.log('debug', ...
                ['MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_in_process on input file: ',...
                DCE_File]);
            dce_source_fileData = mms_sdc_sdp_cdf_in_process(DCE_File, 'sci', 'dce');

            % Then get dcv data
            irf.log('debug', ...
                ['MMS_SDC_SDP_PROC SITL trying mms_sdc_sdp_cdf_in_process on input file: ',...
                DCV_File]);
            dcv_source_fileData = mms_sdc_sdp_cdf_in_process(DCV_File, 'sci', 'dcv');

            % Then the sunpulse file
            irf.log('debug', ...
                ['MMS_SDC_SDP_PROC SITL trying mms_sdc_sdp_cdf_in_process on input file: ',...
                HK_101_File]);
            mms_sdc_sdp_cdf_in_process(HK_101_File, 'sci', 'hk_101');
            irf.log('debug', 'MMS_SDC_SDP_PROC Request the phase');
            dcephase = mms_sdc_sdp_datamanager('dcephase');
            
            % Set bitmask for overlap or not
            bitmask = mms_sdc_sdp_bitmasking(mms_sdc_sdp_datamanager('dcetime'), mms_sdc_sdp_datamanager('dcvtime'));

            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.
          
            
            copy_header(2)
       
            
            irf.log('debug', ...
                'MMS_SDC_SDP_PROC SITL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(bitmask(:,2), HeaderInfo);

        else
            
            irf.log('critical',...
                'MMS_SDC_SDP_PROC SITL received some unclear input arguments. DCE and HK_101 CDF files are required, DCV needed for some but not all processing.');
            irf.log('warning',...
                ['Variable HK_101_File exists: ', ...
                num2str(exist('HK_101_File','var')), ...
                ', variable DCE_File exists: ', ...
                num2str(exist('DCE_File','var')), ...
                ', variable DCV_File exists: ', ...
                num2str(exist('DCV_File','var'))]);
            error('Matlab:MMS_SDC_SDP_PROC:Input',...
                'MMS_SDC_SDP_PROC SITL did not recieve all required input file arguments or was unsuccesful in identifying them. Please see log file.');
        
        end
        
        
    case('ql')

        % Check if all required and needed files are sent as input and have
        % been identified properly.
        if( all( [~isempty(DCE_File), ~isempty(HK_101_File), ...
                isempty(DCV_File)] ) )
            
            % Log message so we know we are missing one input.
            irf.log('warning',...
                'MMS_SDC_SDP_PROC QL received DCE and HK_101 but no DCV file argument. Can perform some but not all processing.');

            irf.log('debug',...
                ['MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_in_process on input file: ',...
                DCE_File]);
            dce_source_fileData = mms_sdc_sdp_cdf_in_process(DCE_File, 'sci', 'dce');

            irf.log('debug',...
                ['MMS_SDC_SDP_PROC QL trying mms_sdc_sdp_cdf_in_process on input file: ',...
                HK_101_File]);
            mms_sdc_sdp_cdf_in_process(HK_101_File, 'sci', 'hk_101');
            irf.log('debug', 'MMS_SDC_SDP_PROC Request the phase');
            dcephase = mms_sdc_sdp_datamanager('dcephase');
            
            % Set bitmask for all times in dcetime.
            bitmask = mms_sdc_sdp_bitmasking(mms_sdc_sdp_datamanager('dcetime'));
            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.

            % FIXME: Do some proper quality assesment.
            quality = bitmask; 

            copy_header(1)
            irf.log('debug', ...
                'MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(bitmask(:,2), ...
                HeaderInfo, quality(:,2));

            
        elseif( all( [~isempty(HK_101_File), ~isempty(DCE_File), ...
                ~isempty(DCV_File)] ) )
            
            % Log message so we know we got both.
            irf.log('notice','MMS_SDC_SDP_PROC QL received all expected input arguments, DCE, DCV and HK_101 file arguments. Can perform full processing.');

            % First get dce data
            irf.log('debug', ...
                ['MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_in_process on input file: ', ...
                DCE_File]);
            dce_source_fileData = mms_sdc_sdp_cdf_in_process(DCE_File, 'sci', 'dce');

            % Then get dcv data
            irf.log('debug', ...
                ['MMS_SDC_SDP_PROC QL trying mms_sdc_sdp_cdf_in_process on input file: ', ...
                DCV_File]);
            dcv_source_fileData = mms_sdc_sdp_cdf_in_process(DCV_File, 'sci', 'dcv');

            irf.log('debug', ...
                ['MMS_SDC_SDP_PROC QL trying mms_sdc_sdp_cdf_in_process on input file: ', ...
                HK_101_File]);
            mms_sdc_sdp_cdf_in_process(HK_101_File, 'sci', 'hk_101');
            irf.log('debug', 'MMS_SDC_SDP_PROC Request the phase');
            dcephase = mms_sdc_sdp_datamanager('dcephase');
            
            % Set bitmask for overlap or not.
            bitmask = mms_sdc_sdp_bitmasking(mms_sdc_sdp_datamanager('dcetime'), mms_sdc_sdp_datamanager('dcvtime'));

            
            % FIXME Do some processing... Actually De-filter, De-spin, etc.

            % FIXME: Do some proper quality assesment.
            quality = bitmask; 
            
            
            copy_header(2)
            irf.log('debug', ...
                'MMS_SDC_SDP_PROC QL using mms_sdc_sdp_cdf_writing');
            filename_output = mms_sdc_sdp_cdf_writing(bitmask(:,2), ...
                HeaderInfo, quality(:,2));

        else
            
            irf.log('critical',...
                'MMS_SDC_SDP_PROC QL received some unclear input arguments. DCE and HK_101 CDF files are required, DCV needed for some but not all processing.');
            irf.log('warning',...
                ['Variable HK_101_File exists: ', ...
                num2str(~isempty(HK_101_File)), ...
                ', variable DCE_File exists: ', ...
                num2str(~isempty(DCE_File)), ...
                ', variable DCV_File exists: ', ...
                num2str(~isempty(DCV_File))]);
            error('Matlab:MMS_SDC_SDP_PROC:Input',...
                'MMS_SDC_SDP_PROC QL did not recieve all required input file arguments or was unsuccesful in identifying them. Please see log file.');
        
        end

        
    otherwise
        % Should not be here
        error('Matlab:MMS_SDC_SDP_PROC:Input', 'wrong process name');
end


% Write out filename as empty logfile so it can be easily found by SDC
% scripts.
unix(['touch', ' ', ENVIR.LOG_PATH_ROOT, filesep 'mms', ...
    HeaderInfo.numberStr, filesep, 'sdp', filesep, filename_output, ...
    '_',runTime,'.log']);

  function init_matlab_path()
    irfPath = [irf('path') filesep];
    irfDirectories = {'irf',...
      ['mission' filesep 'mms'],...
      ['mission' filesep 'cluster'],...
      ['contrib' filesep 'nasa_cdf_patch'],...
      };
    for iPath = 1:numel(irfDirectories)
      pathToAdd = [irfPath irfDirectories{iPath}];
      addpath(pathToAdd);
      disp(['Added to path: ' pathToAdd]);
    end
  end

  function copy_header(nSrc)
    if nSrc>2 || nSrc<1, error('nSrc must be 1 or 2'), end
    HeaderInfo.scId = dce_source_fileData.scId;
    HeaderInfo.instrumentId = dce_source_fileData.instrumentId;
    HeaderInfo.dataMode = dce_source_fileData.dataMode;
    HeaderInfo.dataLevel = dce_source_fileData.dataLevel;
    HeaderInfo.startTime = dce_source_fileData.startTime;
    HeaderInfo.numberOfSources = nSrc;
    HeaderInfo.parents_1 = dce_source_fileData.filename;
    if nSrc==2
      HeaderInfo.parents_2 = dcv_source_fileData.filename;
    end
  end
end