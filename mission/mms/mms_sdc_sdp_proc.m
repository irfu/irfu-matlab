function mms_sdc_sdp_proc( procName, varargin)
% MMS_SDC_SDP_PROC startup function for MMS processing. Rely on DEFATT
% files for phase calculations, if they are present in the DATA_PATH_ROOT
% instead of hk 101 sunpulses.
%
% 	See also MMS_SDC_SDP_INIT, MMS_SDP_DMGR.

% Store runTime when script was called. To ease for SDC, keep track of
% start time (logs are written with this time in their file name).
runTime = char(datetime("now","Format","uuuuMMdd'T'HHmmss"));

global ENVIR MMS_CONST;

% Setup irfu-matlab
init_matlab_path()

% Load global contants.
if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

HK_101_File = ''; % HK with sunpulse, etc.
HK_105_File = ''; % HK with sweep status etc.
HK_10E_File = ''; % HK with bias guard settings etc.
ASPOC_File = '';
DFG_File = ''; % B-field, L2Pre
AFG_File = ''; % B-Field, L2pre fallback (if no DFG present)
DCV_File = '';
DCE_File = '';
L2A_File = ''; % L2A file, contain offsets from fast/slow to be used by brst and for L2Pre process.
SCPOT_File = ''; % SCPOT file, used for computing DSL offsets when processing Slow L2a to Slow L2pre.
L2Pre_File = ''; % L2Pre file, containing spin residue from Fast to be used by brst for L2pre process.
DEFATT_File = ''; % Defatt file used for l2pre or reprocessing of QL.
DEFEPH_File = ''; % Defeph file used for orbit information (Slow mode offset is radius dependent).
HdrInfo = [];
procId=[]; Dmgr = [];

parse_input(procName, varargin{:});

% All input arguments read. All files required identified correct?
% QL, SCPOT req: hk101 or DEFATT, hk10e, hk105, dce and possibly dcv (if not incl. in dce).
if any([(isempty(HK_101_File) && isempty(DEFATT_File)), isempty(DCE_File),...
    isempty(HK_105_File), isempty(HK_10E_File)])
  irf.log('warning', 'MMS_SDC_SDP_PROC missing some input.');
  for i=1:nargin-1
    irf.log('warning',['Received input argument: ', varargin{i}]);
  end
elseif(isempty(DCV_File) && procId~=MMS_CONST.SDCProc.l2pre )
  irf.log('debug','It appears we are running with the dcv data combined into dce file.');
end

% When running multiple instances in parallel at SDC sometimes we cannot
% access the mat file "~./matlab_datastore_hostname" on all machines. Add a
% delay on the problematic machine and try again.
dbSuccess = false; dbCounter=0;
while(~dbSuccess && dbCounter<=5)
  try
    dbCounter = dbCounter + 1;
    mms.db_init('local_file_db', ENVIR.DATA_PATH_ROOT); % Setup mms database
    dbSuccess = true; % if previous command was successful then exit loop
  catch ME
    irf.log('warning', 'FAILED to initiate MMS db, trying again in one second.');
    irf.log('warning', ['Error was:', ME.identifier, ' with message: ', ME.message]);
    pause(floor(rand(1)*5)+1); % Sleep random value up to 6 second, then try again.
  end
end
if(~dbSuccess), irf.log('critical', 'FAILED to initiate MMS db'); end

%% Processing for SCPOT, QL or other products.
% Load and process identified files in the following order first any of the
% available files out of "DCE", "HK_10E", "HK_101", "HK_105", then lastly
% the "DCV".
% Reason: DCV in data manager calls on other subfunctions which may require
% DCE, HK_101, HK_105 and HK_10E files to already be loaded into memory.

switch procId
  case {MMS_CONST.SDCProc.scpot, MMS_CONST.SDCProc.ql, MMS_CONST.SDCProc.l2a}
    if(~isempty(HK_10E_File))
      fileSplit = strsplit(HK_10E_File,':');
      for iFile=1:size(fileSplit,2)
        irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
        src_fileData = load_file(fileSplit{iFile},'hk_10e');
        update_header(src_fileData); % Update header with file info.
      end
    end
    if(~isempty(HK_105_File))
      fileSplit = strsplit(HK_105_File,':');
      for iFile=1:size(fileSplit,2)
        irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
        src_fileData = load_file(fileSplit{iFile},'hk_105');
        update_header(src_fileData); % Update header with file info.
      end
    end
    
    %% PHASE and ASPOC information, somewhat special case.
    % Begin by loading the DCE file in order to get time interval of
    % interest, then if no DEFATT was sent go looking for it. If no DEFATT
    % is found, use HK 101. Similar for ASPOC, if no input was sent go look
    % for it. If no ASPOC is found, simply skip it.
    if isempty(DCE_File)
      errStr = ['missing required input for ' procName ': DCE_File'];
      irf.log('critical',errStr);  error(errStr);
    end
    irf.log('notice', [procName ' proc using: ' DCE_File]);
    dce_obj = dataobj(DCE_File);
    [~,tmpName, ~] = fileparts(DCE_File);
    update_header(mms_fields_file_info(tmpName)); % Update header with file info.
    if(dce_obj.data.Epoch.nrec==0)
      errStr='Empty Epoch. Possibly started processing too early..';
      irf.log('critical',errStr); error(errStr);
    end
    tint = EpochTT(sort(dce_obj.data.Epoch.data));
    % Keep only valid times (well after end of mission).
    tint(tint.tlim(irf.tint('2015-01-01T00:00:00.000000000Z/2040-12-31T23:59:59.999999999Z')));
    % Create a time interval for start and stop of dce epoch times.
    tint = irf.tint(tint.start, tint.stop);
    if(~dbSuccess)
      % Setup mms database, if it failed several times before.
      mms.db_init('local_file_db', ENVIR.DATA_PATH_ROOT);
    end
    if(isempty(DEFATT_File))
      % Go looking for DEFATT to match tint.
      list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_ancillary_defatt'],tint);
      if(isempty(list) || list(1).start >= tint.start || list(end).stop <= tint.stop)
        % If no DEFATT was found or it did not cover all of tint, use HK 101 files.
        if(~isempty(HK_101_File))
          fileSplit = strsplit(HK_101_File,':');
          for iFile=1:size(fileSplit,2)
            irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
            src_fileData = load_file(fileSplit{iFile},'hk_101');
            update_header(src_fileData) % Update header with file info.
          end
        else
          % Should not be here!
          errStr = 'No DEFATT was found and no HK 101 identified in arguments.';
          irf.log('critical',errStr); error(errStr);
        end
      else
        for ii=1:length(list)
          irf.log('notice', [procName ' proc using: ',list(ii).name]);
          [dataTmp, src_fileData] = mms_load_ancillary([list(ii).path, filesep, ...
            list(ii).name], 'defatt');
          Dmgr.set_param('defatt', dataTmp);
          update_header(src_fileData); % Update header with file info.
        end
      end
    end
    % Similar for ASPOC.
    if(isempty(ASPOC_File))
      % Go looking for ASPOC to match tint.
      list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_aspoc_srvy_l2'],tint+[-100 0]);
      if(isempty(list) || list(1).start > tint.start)
        % If no ASPOC was found or it did not cover start of tint. Simply
        % issue warning.
        irf.log('warning','No ASPOC files located or sent in input.');
      else
        for ii=1:length(list)
          irf.log('notice', [procName ' proc using: ',list(ii).name]);
          src_fileData = load_file([list(ii).path, filesep, list(ii).name],...
            'aspoc');
          update_header(src_fileData); % Update header with file info.
        end
      end
    else
      % Load input specified ASPOC
      fileSplit = strsplit(ASPOC_File,':');
      for iFile=1:size(fileSplit,2)
        irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
        src_fileData = load_file(fileSplit{iFile},'aspoc');
        update_header(src_fileData) % Update header with file info.
      end
    end
    
    %% Second type of special case, brst QL or L2A (use L2A from previously processed Fast).
    if(regexpi(DCE_File,'_brst_') )
      if( procId==MMS_CONST.SDCProc.ql || procId==MMS_CONST.SDCProc.l2a)
        if( ~isempty(L2A_File))
          irf.log('notice', [procName ' proc using: ' L2A_File]);
          src_fileData = load_file(L2A_File,'l2a');
          update_header(src_fileData); % Update header with file info.
        else
          irf.log('warning',[procName ' but no L2A file from Fast. Looking for it..']);
          list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_edp_fast_l2a_dce2d'], tint);
          if(isempty(list) || list(1).start > tint.start)
            % If no L2a dce2d was found or it did not cover start of tint.
            % Simply issue warning.
            irf.log('warning','No Fast L2a dce2d file located.');
          else
            irf.log('notice', [procName ' proc using: ',list(1).name]);
            src_fileData = load_file([list(1).path, filesep, list(1).name],...
              'l2a');
            update_header(src_fileData); % Update header with file info.
          end
        end
      end % If QL or L2A
    end % If running Brst dce
    
    %% Third type of special case, Slow mode is now dependent in radius found in DEFEPH
    if( procId==MMS_CONST.SDCProc.l2a )
      if(isempty(DEFEPH_File))
        % Go looking for DEFEPH to match tint.
        list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_ancillary_defeph'],tint);
        if(isempty(list) || list(1).start >= tint.start || list(end).stop <= tint.stop)
          % Should not be here!
          %% FIXME (go looking for Pred. Eph?)
          errStr = 'No DEFEPH was found and no DEFEPH identified in arguments when processing SLOW mode.';
          irf.log('critical',errStr);% error(errStr);
        else
          for ii=1:length(list)
            irf.log('notice', [procName ' proc using: ',list(ii).name]);
            [dataTmp, src_fileData] = mms_load_ancillary([list(ii).path, filesep, ...
              list(ii).name], 'defeph');
            Dmgr.set_param('defeph', dataTmp);
            update_header(src_fileData); % Update header with file info.
          end
        end
      else
        % Load input specified DEFEPH
        fileSplit = strsplit(DEFEPH_File,':');
        for iFile=1:size(fileSplit,2)
          irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
          [dataTmp, src_fileData] = mms_load_ancillary(fileSplit{iFile}, 'defeph');
          Dmgr.set_param('defeph', dataTmp);
          update_header(src_fileData) % Update header with file info.
        end
      end % DEFEPH special case
    end % If running L2a processing (check for DefEph).
    
    % Go on with the DCE file.
    Dmgr.set_param('dce',dce_obj);
    
    if ~isempty(DCV_File)
      % Separate DCV file (during commissioning)
      irf.log('notice', [procName ' proc using: ' DCV_File]);
      src_fileData = load_file(DCV_File,'dcv');
      update_header(src_fileData) % Update header with file info.
    end
    
  case {MMS_CONST.SDCProc.l2pre}
    % L2Pre process with L2A Fast/Slow file and DFG L2Pre as input. Or if
    % processing Brst segments then inputs are/should be DCE_File (brst),
    % L2A_File (fast), HK105_File, HK10E_File, HK101_File, DEFATT, ASPOC,
    % and the corresponding DFG L2Pre file(-s).
    
    % DFG is required for both Fast/Slow L2a->L2Pre and Brst L1b->L2Pre,
    % Mark, in e-mail dated 2020/04/18, noted that MMS2 DFG was powered off
    % for a period this week, as a result no data is available from DFG but
    % AFG is available. Use AFG if DFG is missing!
    if isempty(DFG_File)
      % IF no DFG is found, try fallback to AFG_File
      errStr = ['No DFG file found for ' procName ', trying AFG file.'];
      irf.log('warning', errStr)
      if isempty(AFG_File)
        errStr = ['missing required input for ' procName ': DFG_File & AFG_File'];
        irf.log('critical',errStr)
        error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
      else
        fileSplit = strsplit(AFG_File,':');
        for iFile=1:size(fileSplit,2)
          irf.log('notice',[procName ' proc using: ' fileSplit{iFile}]);
          src_fileData = load_file(fileSplit{iFile}, 'dfg');
          update_header(src_fileData) % Update header with file info.
        end
      end
    else
      fileSplit = strsplit(DFG_File,':');
      for iFile=1:size(fileSplit,2)
        irf.log('notice',[procName ' proc using: ' fileSplit{iFile}]);
        src_fileData = load_file(fileSplit{iFile}, 'dfg');
        update_header(src_fileData) % Update header with file info.
      end
    end % DFG
    
    if(~isempty(DCE_File))
      % L1b brst -> L2Pre
      if(~isempty(HK_10E_File))
        fileSplit = strsplit(HK_10E_File,':');
        for iFile=1:size(fileSplit,2)
          irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
          src_fileData = load_file(fileSplit{iFile},'hk_10e');
          update_header(src_fileData); % Update header with file info.
        end
      end % HK 10E
      if(~isempty(HK_105_File))
        fileSplit = strsplit(HK_105_File,':');
        for iFile=1:size(fileSplit,2)
          irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
          src_fileData = load_file(fileSplit{iFile},'hk_105');
          update_header(src_fileData); % Update header with file info.
        end
      end % HK 105
      %% PHASE and ASPOC information, somewhat special case.
      % Begin by loading the DCE file in order to get time interval of
      % interest, then if no DEFATT was sent go looking for it. If no DEFATT
      % is found, use HK 101. Similar for ASPOC, if no input was sent go look
      % for it. If no ASPOC is found, simply skip it.
      irf.log('notice', [procName ' proc using: ' DCE_File]);
      dce_obj = dataobj(DCE_File);
      [~,tmpName, ~] = fileparts(DCE_File);
      update_header(mms_fields_file_info(tmpName)); % Update header with file info.
      if(dce_obj.data.Epoch.nrec==0)
        errStr='Empty Epoch. Possibly started processing too early..';
        irf.log('critical',errStr); error(errStr);
      end
      tint = EpochTT(sort(dce_obj.data.Epoch.data));
      % Keep only valid times (well after end of mission).
      tint(tint.tlim(irf.tint('2015-01-01T00:00:00.000000000Z/2040-12-31T23:59:59.999999999Z')));
      % Create a time interval for start and stop of dce epoch times.
      tint = irf.tint(tint.start, tint.stop);
      if(~dbSuccess)
        % Setup mms database, if it failed several times before.
        mms.db_init('local_file_db', ENVIR.DATA_PATH_ROOT);
      end
      if(isempty(DEFATT_File))
        % Go looking for DEFATT to match tint.
        list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_ancillary_defatt'],tint);
        if(isempty(list) || list(1).start >= tint.start || list(end).stop <= tint.stop)
          % If no DEFATT was found or it did not cover all of tint, use HK 101 files.
          if(~isempty(HK_101_File))
            fileSplit = strsplit(HK_101_File,':');
            for iFile=1:size(fileSplit,2)
              irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
              src_fileData = load_file(fileSplit{iFile},'hk_101');
              update_header(src_fileData) % Update header with file info.
            end
          else
            % Should not be here!
            errStr = 'No DEFATT was found and no HK 101 identified in arguments.';
            irf.log('critical',errStr); error(errStr);
          end
        else
          for ii=1:length(list)
            irf.log('notice', [procName ' proc using: ',list(ii).name]);
            [dataTmp, src_fileData] = mms_load_ancillary([list(ii).path, filesep, ...
              list(ii).name], 'defatt');
            Dmgr.set_param('defatt', dataTmp);
            update_header(src_fileData); % Update header with file info.
          end
        end
      end % DEFATT special case
      % Similar for ASPOC.
      if(isempty(ASPOC_File))
        % Go looking for ASPOC to match tint.
        list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_aspoc_srvy_l2'],tint);
        if(isempty(list) || list(1).start > tint.start)
          % If no ASPOC was found or it did not cover start of tint. Simply
          % issue warning.
          irf.log('warning','No ASPOC files located or sent in input.');
        else
          for ii=1:length(list)
            irf.log('notice', [procName ' proc using: ',list(ii).name]);
            src_fileData = load_file([list(ii).path, filesep, list(ii).name],...
              'aspoc');
            update_header(src_fileData); % Update header with file info.
          end
        end
      else
        % Load input specified ASPOC
        fileSplit = strsplit(ASPOC_File,':');
        for iFile=1:size(fileSplit,2)
          irf.log('notice', [procName ' proc using: ' fileSplit{iFile}]);
          src_fileData = load_file(fileSplit{iFile},'aspoc');
          update_header(src_fileData) % Update header with file info.
        end
      end % ASPOC special case
      %% L2A Fast mode file
      if(~isempty(L2A_File))
        irf.log('notice', [procName ' proc using: ' L2A_File]);
        src_fileData = load_file(L2A_File,'l2a');
        update_header(src_fileData); % Update header with file info.
      else
        irf.log('warning',[procName ' but no L2A file from Fast. Looking for it..']);
        list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_edp_fast_l2a_dce2d'], tint);
        if(isempty(list) || list(1).start > tint.start)
          % If no L2a dce2d was found or it did not cover start of tint.
          % Simply issue warning.
          irf.log('critical','No Fast L2a dce2d file located.');
          error('NO L2A fast dce2d file in input and no file found!!');
        else
          irf.log('notice', [procName ' proc using: ',list(1).name]);
          src_fileData = load_file([list(1).path, filesep, list(1).name],...
            'l2a');
          update_header(src_fileData); % Update header with file info.
        end
      end
      %% L2Pre Fast mode file
      if(~isempty(L2Pre_File))
        irf.log('notice', [procName ' proc using: ' L2Pre_File]);
        src_fileData = load_file(L2Pre_File,'l2pre');
        update_header(src_fileData); % Update header with file info.
      else
        irf.log('warning',[procName ' but no L2Pre file from Fast. Looking for it..']);
        list = mms.db_list_files(['mms',HdrInfo.scIdStr,'_edp_fast_l2pre_dce2d'], tint);
        if(isempty(list) || list(1).start > tint.start)
          % If no L2Pre dce2d was found or it did not cover start of tint.
          % Simply issue warning.
          irf.log('critical','No Fast L2Pre dce2d file located.');
          error('NO L2Pre fast dce2d file in input and no file found!!');
        else
          for ii=1:length(list)
            irf.log('notice', [procName ' proc using: ',list(ii).name]);
            src_fileData = load_file([list(ii).path, filesep, list(ii).name],...
              'l2pre');
            update_header(src_fileData); % Update header with file info.
          end
        end
      end
      %% L1B dce file
      Dmgr.set_param('dce',dce_obj);
      Dmgr.process_l2a_to_l2pre(MMS_CONST);
    else
      % Simple L2A Fast/Slow to L2Pre file
      if isempty(L2A_File)
        errStr = ['missing required input for ' procName ': L2A_File'];
        irf.log('critical',errStr)
        error('Matlab:MMS_SDC_SDP_PROC:Input', errStr)
      end
      %% Special case, Slow mode calibration is dependent on SCpot
      if regexpi(L2A_File, '_slow_')
        if ~isempty(SCPOT_File)
          fileSplit = strsplit(SCPOT_File,':');
          for iFile=1:size(fileSplit,2)
            irf.log('notice',[procName ' proc using: ' fileSplit{iFile}]);
            src_fileData = load_file(fileSplit{iFile}, 'scpotFile');
            update_header(src_fileData) % Update header with file info.
          end
        else
          irf.log('warning',[procName ' slow mode but no SCPOT file. Looking for it..']);
          % Create a time interval for start and stop of l2a epoch times.
          dce_obj = dataobj(L2A_File);
          tint = EpochTT(dce_obj.data.(['mms', HdrInfo.scIdStr, '_edp_epoch_slow_l2a']).data);
          tint = irf.tint(tint.start, tint.stop);
          list = mms.db_list_files(['mms', HdrInfo.scIdStr, '_edp_slow_l2_scpot'], tint);
          if isempty(list)
            % If no L2 scpot was found or it did not cover start of tint.
            % Simply issue warning.
            irf.log('warning', 'No SLOW L2 scpot file located.');
          else
            for iFile=1:length(list)
              irf.log('notice', [procName, ' proc using: ', list(iFile).name]);
              src_fileData = load_file([list(iFile).path, filesep, list(iFile).name],...
                'scpotFile');
              update_header(src_fileData); % Update header with file info.
            end
          end
        end
      end % If running SLOW mode
      irf.log('notice', [procName ' proc using: ' L2A_File]);
      src_fileData = load_file(L2A_File,'l2a');
      Dmgr.process_l2a_to_l2pre(MMS_CONST);
      update_header(src_fileData) % Update header with file info.
    end % Empty L1B DCE_File
    
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
    
end

% Write the output
filename_output = mms_sdp_cdfwrite(HdrInfo, Dmgr);


%% Help functions
  function init_matlab_path()
    % Setup irfu-matlab and subdirs
    irfPath = [irf('path') filesep];
    irfDirectories = {'irf',...
      ['mission' filesep 'mms'],...
      ['mission' filesep 'cluster'],...
      ['contrib' filesep 'nasa_cdf_patch'],...
      ['contrib' filesep, 'matlab_central'],...
      };
    for iPath = 1:numel(irfDirectories)
      pathToAdd = [irfPath irfDirectories{iPath}];
      addpath(pathToAdd);
      irf.log('notice',['Added to path: ' pathToAdd]);
    end
  end

  function parse_input(procName, varargin)
    % Process input arguments
    if ~ischar(procName), error('MMS_SDC_SDP_PROC first argument must be a string.'); end
    [~,procId] = intersect( MMS_CONST.SDCProcs, lower(procName));
    if isempty(procId)
      error('MMS_SDC_SDP_PROC first argument must be one of: %s',...
        mms_constants2string('SDCProcs'));
    end
    procName = upper(procName);
    irf.log('notice', ['Starting process: ', procName]);
    
    %% Identify each input argument
    for j=1:nargin-1
      if isempty(varargin{j}), continue, end
      if(~ischar(varargin{j})), error('MMS_SDC_SDP_PROC input parameter must be string.'); end
      [pathIn, fileIn, extIn] = fileparts(varargin{j});
      if any([isempty(pathIn), isempty(fileIn), isempty(extIn)])
        errStr = sprintf('Expected files with full paths, got: %s as input argument number %i.', ...
          varargin{j}, j);
        irf.log('critical', errStr); error(errStr);
      end
      
      if j==1
        % Setup environment.
        HdrInfo.scIdStr = fileIn(4);
        ENVIR = mms_sdc_sdp_init;
        % Setup log
        mms_sdc_sdp_log_init(procName, fileIn, runTime);
      elseif(~strcmp(HdrInfo.scIdStr, fileIn(4)))
        errStr = ['MMS_SDC_SDP_PROC called using MMS S/C: ', ...
          HdrInfo.scIdStr, ' and another file from MMS S/C: ', fileIn(4),'.'];
        irf.log('critical', errStr);
        irf.log('critical', ['Argument ', varargin{j}, ' did not match ',...
          'previous s/c ',varargin{j-1},'. Aborting with error.']);
        error(errStr);
      end
      
      if regexpi(fileIn,'_dce')
        if( (procId == MMS_CONST.SDCProc.l2pre || procId == MMS_CONST.SDCProc.ql) ...
            && ~isempty(regexpi(fileIn,'_l2(a|pre)_')) && any(cell2mat(regexp(varargin(:),'_brst_'))) )
          % L2A or L2Pre file (from fast mode) for "QL Brst" or "L2A Brst" process.
        else
          % Alternative to multiple regexpi() in commissioning data.
          %           expr = ['mms(?<SCid>\d{0,1})_edp_(?<tmModeStr>(', ...
          %              strjoin(MMS_CONST.TmModes, '|'),'))', ...
          %              '_(?<dataLevel>(l1b|l2a|l2pre))', ...
          %              '_dc[ev](?<commRate>\d{0,3})', ...
          %              '_(?<dateTime>\d{8,14})', ...
          %              '_v(?<verX>\d{1,}).(?<verY>\d{1,}).(?<verZ>\d{1,})'];
          %           dceStr = regexpi(fileIn, expr, 'names');
          % This argument is the dce file, (l1b raw, l2a/pre dce2d or similar)
          % Use this file to get TMmode directly from filename, and if comm.
          % data also sample rate. And also initialize the Dmgr.
          tmModeStr = fileIn(10:13); % mmsX_edp_[TMmode]_
          [~,tmMode] = intersect( MMS_CONST.TmModes, tmModeStr);
          if isempty(tmMode)
            errStr = ['Unrecognized tmMode (',tmModeStr,'), must be one of: '...
              mms_constants2string('TmModes')];
            irf.log('critical', errStr);  error(errStr);
          end
          if (tmMode==MMS_CONST.TmMode.comm)
            % Special case, Commissioning data, identify samplerate.
            irf.log('notice',...
              'Commissioning data, trying to identify samplerate from filename.');
            if regexpi(fileIn, '_dc[ev]8_') % _dcv8_ or _dce8_
              samplerate = MMS_CONST.Samplerate.comm_8;
            elseif regexpi(fileIn, '_dc[ev]32_') % _dcv32_ or _dce32_
              samplerate = MMS_CONST.Samplerate.comm_32;
            elseif regexpi(fileIn, '_dc[ev]64_') % _dcv64_ or _dce64_
              samplerate = MMS_CONST.Samplerate.comm_64;
            elseif regexpi(fileIn, '_dc[ev]128_') % _dcv128_ or _dce128_
              samplerate = MMS_CONST.Samplerate.comm_128;
            else
              % Possibly try to look at "dt" from Epoch inside of file? For
              % now just default to first TmMode (slow).
              irf.log('warning',...
                ['Unknown samplerate for Commissioning data from file: ',fileIn]);
              if iscell(MMS_CONST.Samplerate.(MMS_CONST.TmModes{1}))
                samplerate = MMS_CONST.Samplerate.(MMS_CONST.TmModes{1}){1};
              else
                samplerate = MMS_CONST.Samplerate.(MMS_CONST.TmModes{1});
              end
              irf.log('warning', ['Defaulting samplerate to ',num2str(samplerate)]);
            end
            Dmgr = mms_sdp_dmgr(str2double(HdrInfo.scIdStr), procId, tmMode, samplerate);
          elseif(tmMode==MMS_CONST.TmMode.slow || tmMode==MMS_CONST.TmMode.brst)
            % As of 2016/08 discussion about increasing Slow mode to 32 Hz,
            % and eventually Brst mode will change from 8192 Hz to 1024 Hz.
            irf.log('warning', [tmModeStr, ' mode dce looking inside to determine sample rate.']);
            try
              % Read the first 10 records of the main "Epoch" (as int64, TT2000)
              tmpEpoch = spdfcdfread(varargin{j}, 'KeepEpochAsis', true, ...
                'Variables', 'Epoch', 'Records', 1:10);
              % Compute the sample rate of these records and compare with
              % expected values, consider match if within +/- 5%.
              tmpRate = 10^9 / median(diff(double(cell2mat(tmpEpoch))));
              for jj=1:length(MMS_CONST.Samplerate.(MMS_CONST.TmModes{tmMode}))
                if tmpRate < 1.05*MMS_CONST.Samplerate.(MMS_CONST.TmModes{tmMode}){jj} && ...
                    tmpRate > 0.95*MMS_CONST.Samplerate.(MMS_CONST.TmModes{tmMode}){jj}
                  samplerate = MMS_CONST.Samplerate.(MMS_CONST.TmModes{tmMode}){jj};
                  break
                end
              end
            catch err
              irf.log('warning', ['Failed to read sample rate. Falling back to original default for ', tmModeStr]);
            end
            if(exist('samplerate','var') && ~isempty(samplerate) )
              Dmgr = mms_sdp_dmgr(str2double(HdrInfo.scIdStr), procId, tmMode, samplerate);
            else
              Dmgr = mms_sdp_dmgr(str2double(HdrInfo.scIdStr), procId, tmMode, MMS_CONST.Samplerate.(MMS_CONST.TmModes{tmMode}){1});
            end
          else
            % Not "Comm", "Slow" or "Brst". Determine sample rate by tmMode.
            Dmgr = mms_sdp_dmgr(str2double(HdrInfo.scIdStr), procId, tmMode);
          end
        end
      end
      
      if regexpi(fileIn, '_101_') % 101, mmsX_fields_hk_l1b_101_20150410_v0.0.1.cdf
        if ~isempty(HK_101_File)
          errStr = ['Multiple HK_101 files in input (',HK_101_File,' and ',varargin{j},')'];
          irf.log('critical', errStr);  error(errStr);
        end
        HK_101_File = varargin{j};
        irf.log('notice', ['HK_101 input file: ', HK_101_File]);
      elseif regexpi(fileIn, '_10e_') % 10E, mmsX_fields_hk_l1b_10e_20150410_v0.0.1.cdf
        if ~isempty(HK_10E_File)
          errStr = ['Multiple HK_10E files in input (',HK_10E_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        HK_10E_File = varargin{j};
        irf.log('notice', ['HK_10E input file: ', HK_10E_File]);
      elseif regexpi(fileIn, '_105_') % 105, mmsX_fields_hk_l1b_105_20150410_v0.0.1.cdf
        if ~isempty(HK_105_File)
          errStr = ['Multiple HK_105 files in input (',HK_10E_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        HK_105_File = varargin{j};
        irf.log('notice', ['HK_105 input file: ', HK_105_File]);
      elseif regexpi(fileIn, '_dcv\d{0,3}_') % _dcv_ or _dcv32_ or _dcv128_
        if ~isempty(DCV_File)
          errStr = ['Multiple DC V files in input (',DCV_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        DCV_File = varargin{j};
        irf.log('notice', ['DCV input file: ', DCV_File]);
      elseif regexpi(fileIn, '_dce\d{0,3}_') % _dce_ or _dce32_ or _dce128_
        if ~isempty(DCE_File)
          errStr = ['Multiple DC E files in input (',DCE_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        DCE_File = varargin{j};
        irf.log('notice', ['DCE input file: ', DCE_File]);
      elseif regexpi(fileIn, '_l2a_') % L2A file (produced by QL Fast/slow)
        if ~isempty(L2A_File)
          errStr = ['Multiple L2A files in input (',L2A_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        L2A_File = varargin{j};
      elseif regexpi(fileIn, '_l2pre_dce2d_') % L2Pre file
        if ~isempty(L2Pre_File)
          errStr = ['Multiple L2Pre dce2d files in input (',L2Pre_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        L2Pre_File = varargin{j};
      elseif regexpi(fileIn, '_DEFATT_') % DEFATT
        if ~isempty(DEFATT_File)
          errStr = ['Multiple DEFATT files in input (',DEFATT_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        DEFATT_File = varargin{j};
        irf.log('notice', ['DEFATT input file: ', DEFATT_File]);
      elseif regexpi(fileIn, '_aspoc_') % ASPOC
        if ~isempty(ASPOC_File)
          errStr = ['Multiple ASPOC files in input (',ASPOC_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        ASPOC_File = varargin{j};
        irf.log('notice',['ASPOC input file: ',ASPOC_File]);
      elseif regexpi(fileIn, '_dfg_') % DFG - B-field
        if ~isempty(DFG_File)
          errStr = ['Multiple DFG files in input (',DFG_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        DFG_File = varargin{j};
        irf.log('notice',['DFG input file: ',DFG_File]);
      elseif regexpi(fileIn, '_afg_') % AFG - B-field (fallback for times with missing DFG
        if ~isempty(AFG_File)
          errStr = ['Multiple AFG files in input (',AFG_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        AFG_File = varargin{j};
        irf.log('notice',['AFG input file: ',AFG_File]);
      elseif regexpi(fileIn, '_DEFEPH_') % DEFEPH
        if ~isempty(DEFEPH_File)
          errStr = ['Multiple DEFEPH files in input (',DEFEPH_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        DEFEPH_File = varargin{j};
        irf.log('notice', ['DEFEPH input file: ', DEFEPH_File]);
      elseif regepxi(fileIn, '_SCPOT_') % SCPOT
        if ~isempty(SCPOT_File)
          errStr = ['Multiple SCPOT files in input (',SCPOT_File,' and ',varargin{j},')'];
          irf.log('critical', errStr); error(errStr);
        end
        SCPOT_File = varargin{j};
        irf.log('notice', 'SCPOT input file: ', SCPOT_File);
      else
        % Unidentified input argument
        errStr = ['MMS_SDC_SDP_PROC unrecognized input file: ',varargin{j}];
        irf.log('critical', errStr); error(errStr);
      end
    end % End for j=1:nargin-1
  end % End parse_input

  function [filenameData] = load_file(fullFilename, dataType)
    [~, fileName, ~] = fileparts(fullFilename);
    filenameData = mms_fields_file_info(fileName);
    Dmgr.set_param(dataType, fullFilename);
  end

  function update_header(src)
    % Update header info
    if(~isempty(regexpi(src.filename,'dce')) && ~isfield(HdrInfo,'startTime'))
      HdrInfo.startTime = src.startTime;
    end
    % Initialization. Store startTime and first filename as parents_1.
    if(~isfield(HdrInfo,'parents_1'))
      HdrInfo.parents_1 = src.filename;
      HdrInfo.numberOfSources = 1;
    else % Next run.
      % Increase number of sources and new parent information.
      HdrInfo.numberOfSources = HdrInfo.numberOfSources + 1;
      HdrInfo.(sprintf('parents_%i', HdrInfo.numberOfSources)) = ...
        src.filename;
    end
  end

end
