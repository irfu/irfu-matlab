function mms_sdc_sdp_proc( procName, varargin)
% MMS_SDC_SDP_PROC main starting point for MMS SDC processing.
%	MMS_SDC_SDP_PROC('processingType', '/pathTo/input/file1.cdf', ...
%   '/pathTo/input/file2.cdf', ...); will start the MMS processing for
%   processing types "ql", "sitl" and "scpot" which is to be performed at
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

global ENVIR MMS_CONST;

narginchk(4,4);

init_matlab_path()

HK_101_File = '';
DCV_File = '';
DCE_File = '';
HeaderInfo = [];

if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

% First argument is and should always be which mode to run.
% QuickLook, SCPOT or SITL. (Is static in the bash script that starts Matlab)
if ~ischar(procName)
    error('Matlab:MMS_SDC_SDP_PROC:Input', ...
    'MMS_SDC_SDP_PROC first argument must be a string');
end
[~,procId] = intersect( MMS_CONST.SDCProcs, lower(procName));
if isempty(procId)
  error('Matlab:MMS_SDC_SDP_PROC:Input', ...
    'MMS_SDC_SDP_PROC first argument must be one of: %s',...
    mms_constants2string('SDCProcs'));
end
procName = upper(procName);
irf.log('notice', ['Starting process: ', procName]);

%% Process inpit
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
        scNumberStr = fileIn(4);
        ENVIR = mms_sdc_sdp_init(scNumberStr);
        scId = str2double(scNumberStr);
        HeaderInfo.numberStr = scNumberStr;
        tmModeStr = fileIn(10:13);
        [~,tmMode] = intersect( MMS_CONST.TmModes, tmModeStr);
        if isempty(tmMode)
          errStr = ['Invalid file name: unrecognized tmMode( ' tmModeStr ...
            '), must be one of: ' mms_constants2string('TmModes')];
          irf.log('critical', errStr);
          error(errStr)
        end
        mms_sdc_sdp_datamanager('init',...
          struct('scId',scId,'tmMode',tmMode,'procId',procId))
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
        irf.log('notice', ['HK_101 input file: ', ...
            HK_101_File]);
        
    elseif regexpi(fileIn, '_dcv_') %DCV
        if ~isempty(DCV_File)
            err_str = ['Received multiple DC V files in input (',...
                DCV_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        DCV_File = varargin{i};
        irf.log('notice', ['DCV input file: ', DCV_File]);
        
    elseif regexpi(fileIn, '_dce_') % DCE
        if ~isempty(DCE_File)
            err_str = ['Received multiple DC E files in input (',...
                DCE_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        DCE_File = varargin{i};
        irf.log('notice', ['DCE input file: ', DCE_File]);
        
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

%% Processing for SCPOT or QL or SITL.
switch procId
  case MMS_CONST.SDCProc.scpot
    if isempty(DCV_File)
      errStr = ['missing reqired input for ' procName ': DCV_File'];
      irf.log('critical',errStr)
      error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
    end
    
    if isempty(DCE_File)
      irf.log('warning', ['MMS_SDC_SDP_PROC ' procName...
        'received no DCE file argument.']);
    else
      irf.log('notice', [procName ' proc using: ' DCE_File]);
      dce_src_fileData = mms_sdc_sdp_cdf_in_process(DCE_File,'sci','dce');
    end
    
    irf.log('notice', [procName ' proc using: ' DCV_File]);
    dcv_src_fileData = mms_sdc_sdp_cdf_in_process(DCV_File,'sci','dcv');
    update_header(dcv_src_fileData,1) % Update header with primary file dcv
    
    if(~isempty(DCE_File))
      update_header(dce_src_fileData,2); % Update header with extra file dce
    end
    
    irf.log('notice', [procName ' proc using: ' HK_101_File]);
    hk_src_fileData = mms_sdc_sdp_cdf_in_process(HK_101_File,'sci','hk_101');
    update_header(hk_src_fileData,2) % Update header with extra file hk
    
    % Write the output
    % Test the new cdf_patch:
    filename_output = mms_sdc_sdp_cdf_writing_2(HeaderInfo);
    
  case {MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2pre}
    % Check if have all the necessary input
    if isempty(DCE_File)
      errStr = ['missing reqired input for ' procName ': DCE_File'];
      irf.log('critical',errStr)
      error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
    end
    if isempty(HK_101_File)
      errStr = ['missing reqired input for ' procName ': HK_101_File'];
      irf.log('critical',errStr)
      error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
    end
    
    irf.log('notice', [procName ' proc using: ' DCE_File]);
    dce_src_fileData = mms_sdc_sdp_cdf_in_process(DCE_File,'sci','dce');
    update_header(dce_src_fileData,1) % Update header with primary file dce

    if isempty(DCV_File)
      irf.log('warning', ['MMS_SDC_SDP_PROC ' procName...
        'received no DCV file argument.']);
    else
      irf.log('notice', [procName ' proc using: ' DCV_File]);
      dcv_src_fileData = mms_sdc_sdp_cdf_in_process(DCV_File,'sci','dcv');
      update_header(dcv_src_fileData,2) % Update header with extra file dcv
    end
    
    irf.log('notice', [procName ' proc using: ' HK_101_File]);
    hk_src_fileData=mms_sdc_sdp_cdf_in_process(HK_101_File,'sci','hk_101');
    update_header(hk_src_fileData,2) % Update header with extra file hk

    % Write the output
    %filename_output = mms_sdc_sdp_cdf_writing(HeaderInfo);
    % Test the new cdf_patch:
    filename_output = mms_sdc_sdp_cdf_writing_2(HeaderInfo);
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

%% Write out filename as empty logfile so it can be easily found by SDC
% scripts.
if ~isempty(ENVIR.LOG_PATH_ROOT)
  unix(['touch', ' ', ENVIR.LOG_PATH_ROOT, filesep 'mms', ...
    HeaderInfo.numberStr, filesep, 'sdp', filesep, filename_output, ...
    '_',runTime,'.log']);
end

  function init_matlab_path()
    irfPath = [irf('path') filesep];
    % For testing and transition to new cdf patch at nasa_cdf_patch_beta
    % and skeletons (the latter is req for beta to write proper
    % multidimensional labels).
    irfDirectories = {'irf',...
      ['mission' filesep 'mms'],...
      ['mission' filesep 'cluster'],...
      ['contrib' filesep 'nasa_cdf_patch'],...
      ['contrib' filesep 'nasa_cdf_patch_beta'],...
      ['mission' filesep 'mms' filesep 'skeletons'],...
      };
    for iPath = 1:numel(irfDirectories)
      pathToAdd = [irfPath irfDirectories{iPath}];
      addpath(pathToAdd);
      irf.log('notice',['Added to path: ' pathToAdd]);
    end
  end

  function update_header(src,updateRun)
    % Update header info
    if updateRun>2 || updateRun<1, error('updateRun must be 1 or 2'), end
    if(updateRun==1) % Initialization. Store startTime and first filename as parents_1.
      HeaderInfo.startTime = src.startTime;
      HeaderInfo.parents_1 = src.filename;
      HeaderInfo.numberOfSources = 1;
    elseif(updateRun==2) % Next run.
      % Increase number of sources and new parent information.
      HeaderInfo.numberOfSources = HeaderInfo.numberOfSources + 1;
      eval(sprintf('HeaderInfo.parents_%i=''%s'';', HeaderInfo.numberOfSources, src.filename))
    end
  end

end
