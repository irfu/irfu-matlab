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
%         '/path/mms2_edp_fast_dce_20150410_v0.0.1.cdf, ...
%         '/path/mms2_edp_fast_dcv_20150410_v0.0.0.cdf', ...
%         '/path/mm2_fields_hk_101_20150410_v0.0.2.cdf');
%	Only partial scientist in the loop ("sitl") processing, with only the
%	required cdf files as input, (note for full sitl include DCV).
%       mms_sdc_sdp_proc('sitl',...
%         '/path/mms2_edp_fast_dce_20150410_v0.0.1.cdf, ...
%         '/path/mm2_fields_hk_101_20150410_v0.0.2.cdf');
%
% 	See also MMS_SDC_SDP_INIT, MMS_SDC_SDP_BITMASKING.

% Store runTime when script was called.
runTime = datestr(now,'yyyymmddHHMMSS'); % Use this time to ease for SDC to
% find output file created. An empty file is to be created with its
% filename being the same as that of the dataproduct created appended with
% _runTime.txt and placed in ENVIR.LOG_PATH_ROOT.

global ENVIR MMS_CONST;

narginchk(2,6);

init_matlab_path()

HK_101_File = ''; % HK with sunpulse, etc.
HK_105_File = ''; % HK with sweep status etc.
HK_10E_File = ''; % HK with bias guard settings etc.
ACE_File = '';
DCV_File = '';
DCE_File = '';
DEFATT_File = ''; % Defatt file used for l2pre
L2Pre_File = ''; % L2Pre file (output from L2pre process is input in L2A?)
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
    if any([isempty(pathIn), isempty(fileIn), isempty(extIn)])
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
        if (tmMode==MMS_CONST.TmMode.comm)
          % Special case, Commissioning data, identify samplerate.
          irf.log('notice',...
            'Commissioning data, trying to identify samplerate from filename.');
          if regexpi(fileIn, '_dc[ev]32_') % _dcv32_ or _dce32_
            samplerate = MMS_CONST.Samplerate.comm_32;
          elseif regexpi(fileIn, '_dc[ev]64_') % _dcv64_ or _dce64_
            samplerate = MMS_CONST.Samplerate.comm_64;
          elseif regexpi(fileIn, '_dc[ev]128_') % _dcv128_ or _dce128_
            samplerate = MMS_CONST.Samplerate.comm_128;
          else
            % Possibly try to look at "dt" from Epoch inside of file? For
            % now just default to first TmMode (slow).
            irf.log('warning',...
              ['Unknown samplerate for Commissioning data from file: ',...
              fileIn]);
            irf.log('warning', ['Defaulting samplerate to ',...
              MMS_CONST.Samplerate.(MMS_CONST.TmModes{1})]);
            samplerate = MMS_CONST.Samplerate.(MMS_CONST.TmModes{1});
          end
          mms_sdp_datamanager('init',...
            struct('scId',scId,'tmMode',tmMode,'procId',procId,...
            'samplerate',samplerate));
        else
          mms_sdp_datamanager('init',...
            struct('scId',scId,'tmMode',tmMode,'procId',procId))
        end
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

    elseif regexpi(fileIn, '_10e_') % 10E, mmsX_fields_hk_l1b_10e_20150410_v0.0.1.cdf
        if ~isempty(HK_10E_File)
            err_str = ['Received multiple HK_10E files in input (',...
                HK_10E_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        % It is the HK_10E file
        HK_10E_File = varargin{i};
        irf.log('notice', ['HK_10E input file: ', ...
            HK_10E_File]);

    elseif regexpi(fileIn, '_105_') % 105, mmsX_fields_hk_l1b_105_20150410_v0.0.1.cdf
        if ~isempty(HK_105_File)
            err_str = ['Received multiple HK_105 files in input (',...
                HK_10E_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        % It is the HK_105 file
        HK_105_File = varargin{i};
        irf.log('notice', ['HK_105 input file: ', ...
            HK_105_File]);

    elseif regexpi(fileIn, '_dcv\d{0,3}_') % _dcv_ or _dcv32_ or _dcv128_
        if ~isempty(DCV_File)
            err_str = ['Received multiple DC V files in input (',...
                DCV_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        DCV_File = varargin{i};
        irf.log('notice', ['DCV input file: ', DCV_File]);
        
    elseif regexpi(fileIn, '_dce\d{0,3}_') % _dce_ or _dce32_ or _dce128_
        if ~isempty(DCE_File)
            err_str = ['Received multiple DC E files in input (',...
                DCE_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        DCE_File = varargin{i};
        irf.log('notice', ['DCE input file: ', DCE_File]);

    elseif regexpi(fileIn, '_ace_') % _ace_
        if ~isempty(ACE_File)
            err_str = ['Received multiple AC E files in input (',...
                ACE_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        ACE_File = varargin{i};
        irf.log('notice', ['ACE input file: ', ACE_File]);

    elseif regexpi(fileIn, '_DEFATT_') % DEFATT
        if ~isempty(DEFATT_File)
            err_str = ['Received multiple DEFATT files in input (',...
                DEFATT_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        DEFATT_File = varargin{i};
        irf.log('notice', ['DEFATT input file: ', DEFATT_File]);

    elseif regexpi(fileIn, '_l2pre_') % L2Pre, used for L2A process.
        if ~isempty(L2Pre_File)
            err_str = ['Received multiple L2Pre files in input (',...
                L2Pre_File, ', ', varargin{i} ')'];
            irf.log('critical', err_str);
            error('Matlab:MMS_SDC_SDP_PROC:Input', err_str);
        end
        L2Pre_File = varargin{i};
        irf.log('notice', ['L2Pre input file: ', L2Pre_File]);
        
    else
        % Unidentified input argument
        err_str = ['MMS_SDC_SDP_PROC unrecognized input file: ',...
            varargin{i}];
        irf.log('critical', err_str);
        error('Matlab:MMS_SDC_SDP_PROC:Input', err_str)
    end
end

% All input arguments read. All files required identified correct?
% QL, SITL, SCPOT use: hk101, hk10e, hk105, dce and possibly dcv (if not incl. in dce).
% L2PRE use: defatt, hk10e, hk105, dce and possibly dcv (if not incl. in dce).
% L2A use: l2pre.
if any([(isempty(HK_101_File) && isempty(DEFATT_File)), isempty(DCE_File),...
    isempty(HK_105_File), isempty(HK_10E_File)]) && isempty(L2Pre_File)
    irf.log('warning', 'MMS_SDC_SDP_PROC missing some input.');
    for i=1:nargin-1
        irf.log('warning',...
            ['MMS_SDC_SDP_PROC received input argument: ', varargin{i}]);
    end
elseif(isempty(DCV_File) && procId~=MMS_CONST.SDCProc.l2a )
  irf.log('debug','It appears we are running with the dcv data combined into dce file.');
end


%% Processing for SCPOT or QL or SITL.
% Load and process identified files in the following order first any of the
% available files out of "DCE", "HK_10E", "HK_101", "HK_105", then lastly
% the "DCV".
% Reason: DCV in data manager calls on other subfunctions which may require
% DCE, HK_101, HK_105 and HK_10E files to already be loaded into memory.

switch procId
  case MMS_CONST.SDCProc.scpot
    if(~isempty(HK_10E_File))
      irf.log('notice', [procName ' proc using: ' HK_10E_File]);
      src_fileData = mms_sdp_load(HK_10E_File,'hk_10e');
      update_header(src_fileData); % Update header with file info.
    end
    if(~isempty(HK_105_File))
      irf.log('notice', [procName ' proc using: ' HK_105_File]);
      src_fileData = mms_sdp_load(HK_105_File,'hk_105');
      update_header(src_fileData); % Update header with file info.
    end
    if isempty(HK_101_File)
      errStr = ['missing required input for ' procName ': HK_101_File'];
      irf.log('critical', errStr);
      error('Matlab:MMS_SDC_SDP_PROC:Input', errStr);
    end
    irf.log('notice', [procName ' proc using: ' HK_101_File]);
    src_fileData = mms_sdp_load(HK_101_File,'hk_101');
    update_header(src_fileData) % Update header with file info.

    if isempty(DCE_File)
      errStr = ['missing required input for ' procName ': DCE_File'];
      irf.log('critical',errStr);
      error('Matlab:MMS_SDC_SDP_PROC:Input', errStr);
    end
    irf.log('notice', [procName ' proc using: ' DCE_File]);
    src_fileData = mms_sdp_load(DCE_File,'dce');
    update_header(src_fileData); % Update header with file info.
    
    if ~isempty(DCV_File)
      % Separate DCV file
      irf.log('notice', [procName ' proc using: ' DCV_File]);
      src_fileData = mms_sdp_load(DCV_File,'dcv');
      update_header(src_fileData) % Update header with file info.
    end

    % Write the output
    filename_output = mms_sdp_cdfwrite(HeaderInfo);
    
  case {MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2pre}
    if(~isempty(HK_10E_File))
      irf.log('notice', [procName ' proc using: ' HK_10E_File]);
      src_fileData = mms_sdp_load(HK_10E_File,'hk_10e');
      update_header(src_fileData) % Update header with file info.
    end
    if(~isempty(HK_105_File))
      irf.log('notice', [procName ' proc using: ' HK_105_File]);
      src_fileData = mms_sdp_load(HK_105_File,'hk_105');
      update_header(src_fileData); % Update header with file info.
    end
    % Phase
    if(procId == MMS_CONST.SDCProc.l2pre)
      % Defatt file => phase
      if isempty(DEFATT_File)
        errStr = ['missing required input for ' procName ': DEFATT_File'];
        irf.log('critical',errStr)
        error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
      end
      irf.log('notice', [procName ' proc using: ' DEFATT_File]);
      [dataTmp,src_fileData]=mms_load_ancillary(DEFATT_File,'defatt');
      mms_sdp_datamanager('defatt',dataTmp);
      update_header(src_fileData); % Update header with file info.
    else
      % HK101 file => phase
      if isempty(HK_101_File)
        errStr = ['missing required input for ' procName ': HK_101_File'];
        irf.log('critical',errStr)
        error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
      end
      irf.log('notice', [procName ' proc using: ' HK_101_File]);
      src_fileData=mms_sdp_load(HK_101_File,'hk_101');
      update_header(src_fileData) % Update header with file info.
    end

    if isempty(DCE_File)
      errStr = ['missing required input for ' procName ': DCE_File'];
      irf.log('critical',errStr);
      error('Matlab:MMS_SDC_SDP_PROC:Input', errStr);
    end
    irf.log('notice', [procName ' proc using: ' DCE_File]);
    src_fileData = mms_sdp_load(DCE_File,'dce');
    update_header(src_fileData) % Update header with file info.

    if ~isempty(DCV_File)
      % Separate DCV file
      irf.log('notice', [procName ' proc using: ' DCV_File]);
      src_fileData = mms_sdp_load(DCV_File,'dcv');
      update_header(src_fileData) % Update header with file info.
    end
    
    % Write the output
    filename_output = mms_sdp_cdfwrite(HeaderInfo);

  case {MMS_CONST.SDCProc.l2a}
    % L2A process with L2Pre file as input
    if isempty(L2Pre_File)
      errStr = ['missing required input for ' procName ': L2Pre_File'];
      irf.log('critical',errStr)
      error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
    end

    irf.log('notice', [procName ' proc using: ' L2Pre_File]);
    src_fileData = mms_sdp_load(L2Pre_File,'l2pre');
    update_header(src_fileData) % Update header with file info.

    % Write the output
    filename_output = mms_sdp_cdfwrite(HeaderInfo);

  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end

%% Write out filename as empty logfile so it can be easily found by SDC
% scripts.
if ~isempty(ENVIR.LOG_PATH_ROOT)
  unix(['touch', ' ', ENVIR.LOG_PATH_ROOT, filesep 'mms', ...
    HeaderInfo.numberStr, filesep, 'edp', filesep, filename_output, ...
    '_',runTime,'.log']);
end

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
      irf.log('notice',['Added to path: ' pathToAdd]);
    end
  end

  function update_header(src)
    % Update header info
    if(regexpi(src.filename,'dce')), HeaderInfo.startTime = src.startTime; end;
    % Initialization. Store startTime and first filename as parents_1.
    if(~isfield(HeaderInfo,'parents_1'))
      HeaderInfo.parents_1 = src.filename;
      HeaderInfo.numberOfSources = 1;
    else % Next run.
      % Increase number of sources and new parent information.
      HeaderInfo.numberOfSources = HeaderInfo.numberOfSources + 1;
      HeaderInfo.(sprintf('parents_%i', HeaderInfo.numberOfSources)) = ...
        src.filename;
    end
  end

end
