classdef solo_local_file_db < solo_file_db
  %SOLO_LOCAL_FILE_DB  Local file database for SOLO
  %   Class handling a database of local SOLO files
  
  properties (SetAccess = immutable)
    dbRoot
  end
  
  methods
    function obj = solo_local_file_db(rootPath)
      % Create local database for SOLO data located in rootPath.
      % Example
      %   SOLO_DB = SOLO_LOCAL_FILE_DB('/data/solo');
      if nargin == 0, rootPath = pwd; end
      if (rootPath(end)==filesep), rootPath(end)=[]; end % path only, excluding last filesep
      
      obj@solo_file_db(rootPath); obj.dbRoot = rootPath;
      if nargin == 0, return, end
      if ~ischar(rootPath)
        errStr = 'rootPath must be a directory path name';
        irf.log('critical',errStr), error(errStr)
      elseif exist(rootPath,'dir')~=7
        errStr = ['DB rootPath (',rootPath,') does not exist. Not mounted?'];
        irf.log('critical',errStr), error(errStr)
      end
    end
    
    %% LIST FILES
    function fileList = list_files(obj,filePrefix,tint)
      % fileList = list_files(obj, filePrefix, [tint]);
      % List files from Database "obj", wich match "filePrefix" and cover
      % optional time period "tint".
      % Example:
      %  SOLO_DB = solo_local_file_db('/data/solo'); % init
      %  tint1 = irf.tint('2020-09-15T00:00:00.000000000Z', '2020-09-15T23:59:59.999999999Z');
      %  fileList = SOLO.db_list_files('solo_L1R_rpw-lfr-surv-swf-b-cdag', tint1);
      narginchk(2,3)
      fileList = [];
      if nargin==3 && (~isempty(tint) && ~isa(tint,'GenericTimeArray'))
        error('Expecting TINT (GenericTimeArray)')
      elseif nargin==2, tint = [];
      end
      if length(filePrefix) < 4 || ~strcmp(filePrefix(1:4),'solo')
        errStr = 'filePrefix must begin with solo*';
        irf.log('critical',errStr), error(errStr)
      end
      C = strsplit(filePrefix,'_');
      if length(C)<3
        errStr = 'filePrefix too short';
        irf.log('critical',errStr), error(errStr)
      end
      if solo.db_index && ~isempty(obj.index)
        irf.log('notice','Using index');
        fileList = obj.index.search_files_with_dataset(filePrefix,tint);
        return
      else
        if ~isempty(tint)
          list_sci_tint();
        else
          irf.log('warning','THIS MAY TAKE SOME TIME');
          list_sci();
        end
      end
      % END LIST_FILES
      
      %% LIST_SCI_TINT
      function list_sci_tint()
        rDir = get_remotePrefix(obj, C);
        fDir = get_fileDir(obj, C);
        TStart = get_times(tint.start); TStop = get_times(tint.stop);
        dateFormat=[]; startFile=[]; stopFile=[];
        for year = TStart.year:TStop.year
          moStart = 1; moStop = 12;
          if year==TStart.year, moStart = TStart.month; end
          if year==TStop.year, moStop = TStop.month; end
          for mo = moStart:moStop
            moDir = sprintf('%s%s%s%s%d%s%02d',rDir,filesep,fDir,filesep,year,filesep,mo);
            curDir = moDir;
            if ismember(C{2}, {'L2', 'L3'})
              % L2 or L3
              dPref = sprintf('%s_%d%02d',filePrefix,year,mo);
              limited_sci_list;
            else
              % L1, L1R, HK etc have daily subfolders
              dStart = 1; dStop = 31;
              if year==TStart.year && mo==TStart.month, dStart=TStart.day; end
              if year==TStop.year && mo==TStop.month, dStop = TStop.day; end
              for day = dStart:dStop
                curDir = [moDir filesep sprintf('%02d',day)];
                dPref = sprintf('%s_%d%02d%02d',filePrefix,year,mo,day);
                limited_sci_list;
              end
            end
          end
        end

        function t = get_times(tt)
          utc = tt.toUtc();
          t.year  = str2double(utc(1:4));
          t.month = str2double(utc(6:7));
          t.day   = str2double(utc(9:10));
          t.hour  = str2double(utc(12:13));
          t.min   = str2double(utc(15:16));
          t.sec   = str2double(utc(18:end-1));
        end
        function limited_sci_list()
          listingD = dir([curDir filesep dPref '*.cdf']); % SolO have only latest file of each type
          if isempty(listingD), return, end
          if isempty(dateFormat)
            % Are we looking for files with 8 or the full 14 digits
            % in the date and time.
            fileformat = regexp(listingD(1).name, '_(?<dateFormat>\d{8,})_V','names');
            switch length(fileformat.dateFormat)
              case 8
                dateFormat = 'yyyymmdd';
              case 14
                dateFormat = 'yyyymmddHHMMSS';
              otherwise
                dateFormat = 'yyyymmdd';
            end
            % Create reconstructed file names for our interval
            startFile = [filePrefix, '_', tint.start.toUtc(dateFormat), '_V00.cdf'];
            stopFile = [filePrefix, '_', tint.stop.toUtc(dateFormat), '_V99.cdf'];
          end
          % Find index of files with names which timewise are sorted
          % between our "startFile" and "stopFile" names.
          indAfterStart = arrayfun(@(x) isequal({startFile; x.name}, sort({startFile; x.name})), listingD);
          indBeforeStop = arrayfun(@(x) isequal({x.name;  stopFile}, sort({x.name;  stopFile})), listingD);
          % Also look at files just before and after as it might be some
          % overlap.
          indLast = find(indBeforeStop, 1, 'last');
          if(indLast<length(listingD)), indBeforeStop(indLast+1) = true; end
          indFirst = find(indAfterStart, 1, 'first');
          if isempty(indFirst)
            indAfterStart(end) = true;
          elseif indFirst>1
            indAfterStart(indFirst-1) = true;
          end
          tmpIndex= find(bitand(indBeforeStop, indAfterStart));
          if isempty(tmpIndex), return, end
          listingD = listingD(tmpIndex);
          arrayfun(@(x) add2list_sci(x.name,curDir), listingD)
        end
      end
      
      %% LIST SCI
      function list_sci()
        rDir = get_remotePrefix(obj, C);
        fDir = get_fileDir(obj, C);
        fileDir = [rDir, filesep, fDir];
        if exist(fileDir,'dir')~=7, return, end
        listingY = dir(fileDir); listingY(~[listingY.isdir]) = [];
        for iDir = 1:length(listingY)
          % Loop over years
          dNameY = listingY(iDir).name;
          if length(dNameY)~=4, continue, end
          yyyy = str2double(dNameY);
          if yyyy<2020 || yyyy > 2050, continue, end
          listingM = dir([fileDir filesep dNameY]);
          listingM(~[listingM.isdir]) = [];
          for iDirMo = 1:length(listingM)
            dNameM = listingM(iDirMo).name;
            if length(dNameM)~=2, continue, end
            switch dNameM(1)
              case '0', if ~any(dNameM(2)=='123456789'), continue, end
              case '1', if ~any(dNameM(2)=='012'), continue, end
              otherwise, continue
            end
            if ismember(C{2}, {'L2', 'L3'})
              % L2 or L3
              curDir = [fileDir filesep dNameY filesep dNameM];
              listingD = dir([curDir, filesep, filePrefix, '*.cdf']);
              if isempty(listingD), continue, end
              arrayfun(@(x) add2list_sci(x.name,curDir), listingD)
            else
              % L1, L1R, HK etc have daily subfolders
              dStart = 1; dStop = 31;
              for day = dStart:dStop
                curDir = [fileDir, filesep, dNameY, filesep, dNameM, ...
                  filesep, sprintf('%02d', day)];
                if ~exist(curDir, 'dir'), continue, end
                listingD = dir([curDir, filesep, filePrefix, '*.cdf']);
                if isempty(listingD), continue, end
                arrayfun(@(x) add2list_sci(x.name,curDir), listingD)
              end
            end
          end
        end
      end % LIST_SCI
      %% ADD2LIST_SCI
      function add2list_sci(name,curDir)
        Entry = struct('name', name, 'ver', str2double(name(end-5:end-4)), ...
          'start',[], 'stop',[],...
          'path', curDir, 'dbId', obj.id);
        Entry = add_ss(Entry);
        % Check time limits of the file
        if isempty(Entry) || ~isempty(tint) && ...
            (Entry.start>tint.stop || Entry.stop<tint.start)
          return
        end
        if isempty(fileList), fileList = Entry; return, end
        fName = name(1:end-7); % Name excl. 'Vxx.cdf'
        
        hasFile = arrayfun(@(x) ~isempty(strfind(x.name,fName)), fileList);
        if ~any(hasFile), fileList = [fileList add_ss(Entry)]; return, end
        iSame = find(hasFile);
        if length(iSame)>1, error('multiple files with same name'); end
        if is_version_larger(str2double(name(end-5:end-4)), fileList(iSame).ver)
          fileList(iSame) = add_ss(Entry); % replace file
        end
        function entry = add_ss(entry)
          entryTmp = obj.cache.get_by_key(entry.name);
          if ~isempty(entryTmp)
            entry.start = entryTmp.start;
            entry.stop = entryTmp.stop;
            return
          end
          try
            info = spdfcdfinfo([entry.path filesep entry.name]);
            if ispc
              % Add a very short delay to ensure consecutive files are not
              % accessed TOO quickly as this may cause Matlab to experince a
              % hard crash on Win10 regardless of the try&catch.
              pause(0.0001);
            end
          catch
            errS = ['Cannot read: ' entry.path filesep entry.name];
            irf.log('critical',errS), error(errS)
          end
          isCdfEpochTT2000VariableArray=cellfun(@(x) strcmpi(x,'tt2000'), info.Variables(:,4));
          if ~any(isCdfEpochTT2000VariableArray)
            errS = ['no TT2000 vars in:' entry.path filesep entry.name];
            irf.log('critical',errS), error(errS)
          end
          iVar = find(isCdfEpochTT2000VariableArray,1);
          data = spdfcdfread([entry.path filesep entry.name], ...
            'Variables', info.Variables(iVar,1), 'CombineRecords', true, ...
            'KeepEpochAsIs', true, 'DataOnly', true);
          if ispc
            % Add a very short delay to ensure consecutive files are not
            % accessed TOO quickly as this may cause Matlab to experince a
            % hard crash on Win10 regardless of the try&catch.
            pause(0.0001);
          end
          if isempty(data), entry = []; return, end
          entry.start = EpochTT(data(1));
          entry.stop = EpochTT(data(end));
          % add to cache
          entryTmp.start = entry.start; entryTmp.stop = entry.stop;
          entryTmp.vars = info.Variables;
          obj.cache.add_entry(entry.name, entryTmp);
        end % ADD_SS
      end % ADD2LIST
    end % LIST_FILES
    
    %% LOAD FILES
    function res = load_file(~, fullPathFilename)
      narginchk(2,3)
      
      irf.log('notice',['loading ', fullPathFilename])
      res = dataobj(fullPathFilename);
    end % LOAD_FILES
    
    %% FILE_HAS_VAR
    function res = file_has_var(obj,fullPathFilename,varName)
      % checks if fileName includes variable name varName
      % res = true/false
      narginchk(3,3)
      res = false; if isempty(varName) || isempty(fullPathFilename), return, end
      
      entryTmp = obj.cache.get_by_key(fullPathFilename);
      if ~isempty(entryTmp)
        res = any(cellfun(@(x) strcmp(x,varName), entryTmp.vars(:,1)));
        return
      end
      if ~exist(fullPathFilename,'file')
        irf.log('warning', ['Fies does not exist: ' fullPathFilename])
        return
      end
      
      % cdf
      if solo.db_index
        res = obj.index.file_has_var(fullPathFilename,varName);
      else
        info = spdfcdfinfo(fullPathFilename);
        if ispc
          % Add a very short delay to ensure consecutive files are not
          % accessed TOO quickly as this may cause Matlab to experince a
          % hard crash on Win10 regardless of the try&catch.
          pause(0.0001);
        end
        res = any(cellfun(@(x) strcmp(x,varName), info.Variables(:,1)));
      end
    end
    
    %% Get metakernel
    function solo_metakernel = get_metakernel(obj, flown_or_predicted)
      % solo_metakernel = get_metakernel(obj, flown_or_predicted)
      % Function to get the name of a local SolO SPICE metakernel to load.
      % Due to the limitations in SPCIE it requires absolute paths to the SPICE
      % kernel files to load, this is likely different on different computers,
      % this script is meant to help with that (assuming the folder with SolO
      % data is constructed as and contains data as it does on IRFU's system).
      % That is to say: a folder subfolder "SPICE" is present in the mounted root
      % folder, and this folder is a clone of the official SPICE git repository
      % for ESA SolarOrbiter.
      %
      % Example:
      %  solo_metakernel = get_metakernel('predicted'); % Get latest predicted kernel
      %  cspice_furnsh(solo_metakernel);
      solo_metakernel = []; % Default empty result
      mkPath = [obj.dbRoot, filesep, 'SPICE', filesep, 'kernels', ...
          filesep, 'mk'];
        switch lower(flown_or_predicted)
          case 'predicted'
            dirs = dir([mkPath, filesep, '*pred-mk_v*.tm']);
          case 'flown'
            dirs = dir([mkPath, filesep, '*flown-mk_v*.tm']);
          otherwise
            error('Unexpected input');
        end
      if size(dirs, 1) > 1
        % Multiple kernels could be found if executing this script at the same
        % time as syncing new kernel files
        error('Found multiple metakernels, please check your SPICE folder.');
      elseif isempty(dirs)
        irf.log('warning', 'Did not find any SPICE metakernels.');
        return
      end
      kernelFile = [dirs.folder, filesep, dirs.name];
      srcKernel = fileread(kernelFile); % Read remote kernel file
      % Local temporary file (placed in "/tmp/" on Unix systems), in which we can
      % correct the local paths to the various SPICE kernel files.
      solo_metakernel = tempname;
      % Full local absolute path (metakernel are in subfolder "mk", so up one level)
      localSpicePath = what([dirs.folder, filesep, '..']);

      % Replace relative or remote paths with absolute path on locally mounted
      % system
      expression = 'PATH_VALUES\s{0,}=\s{0,}(\s{0,}''[a-zA-Z_0-9\.\\]*''\s{0,}';
      localKernel = regexprep(srcKernel, ...
        expression, ...
        ['PATH_VALUES = ( ''', localSpicePath.path, ''' ']);
      if ispc
        % FIXME:
        % Windows system use "\" as filesep, Unix use "/", check OS and if PC then
        % change all of the lines like:
        % KERNELS_TO_LOAD   = (   
        %           '$KERNELS/ck/solo_ANC_soc-sc-iboom-ck_20180930-21000101_V01.bc'
        % to use "\" instead.
        irf.log('critical', 'NOTE: Windows File separations in the metakernel is not yet implemented');
      end

      % Write the local metakernel to file ("cspice_furnsh" do not like it as a
      % string nor a file with relative paths)
      fileID = fopen(solo_metakernel, 'w');
      fwrite(fileID, localKernel);
      fclose(fileID);

      % Now it is ready to be loaded and used (outside of this function)
      end
  end
  
  methods (Access=private)
    function rDir = get_remotePrefix(obj, C)
      if any(contains(C, 'rpw'))
        % offical RPW data is keept in one separate sync folder at IRFU,
        % locally produced data is kept in another folder
        if exist([obj.dbRoot, filesep, 'latest'], 'dir')
          rDir = [obj.dbRoot, filesep, 'latest', filesep, 'RPW'];
        else
          rDir = [obj.dbRoot, filesep, 'remote', filesep, 'data'];
        end
      else
        % All other instruments are kept in a "soar" sync folder, sorted
        % in per instrument folder
        instr = strsplit(C{3}, '-'); % Descriptor contains instrument and dataproduct descriptor part separated by "-".
        rDir = [obj.dbRoot, filesep, 'soar', filesep, instr{1}];
      end
    end % get_remotePrefix
    
    function fileDir = get_fileDir(~, C)
      levelDir = C{2}; % "L2" (or "L1R", "L1", "L3", "HK")
      if ismember(levelDir, {'L2', 'L3'})
        switch C{3}
          case 'rpw-lfr-surv-asm-cdag'
            subDir = 'lfr_asm'; % ie combined 2nd "_" 4th
          case 'rpw-tds-surv-hist1d-cdag'
            subDir = 'hist1d';  % ie 4th
          case 'rpw-tds-surv-hist2d-cdag'
            subDir = 'hist2d';  % ie 4th
          case 'rpw-tds-surv-mamp-cdag'
            subDir = 'mamp';    % ie 4th
          case 'rpw-tds-surv-stat-cdag'
            subDir = 'stat';    % ie 4th
          case {'rpw-lfr-surv-bp1-cdag', 'rpw-lfr-surv-bp2-cdag'}
            subDir = 'lfr_bp';    % ie combined 2nd "_" 4th (excl last digit, which is unique)
          case {'rpw-lfr-surv-cwf-b-cdag', 'rpw-lfr-surv-swf-b-cdag'}
            subDir = 'lfr_wf_b';  % ie combined 2nd "_" 4th and 5th (excl first char of 4th, which is unqiue)
          case {'rpw-lfr-surv-cwf-e-cdag', 'rpw-lfr-surv-swf-e-cdag'}
            subDir = 'lfr_wf_e';  % ie combined 2nd "_" 4th and 5th (excl first char of 4th, which is unqiue)
          case {'rpw-tds-surv-rswf-b-cdag', 'rpw-tds-surf-tswf-b-cdag'}
            subDir = 'tds_wf_b';  % ie combined 2nd "_" 4th and 5th (excl first two chars of 4th, of which the first one is unqiue)
          case {'rpw-tds-surv-rswf-e-cdag', 'rpw-tds-surf-tswf-e-cdag'}
            subDir = 'tds_wf_e';  % ie combined 2nd "_" 4th and 5th (excl first two chars of 4th, of which the first one is unqiue)
          case {'rpw-hfr-surv-cdag', 'rpw-tnr-surv-cdag'}
            subDir = 'thr';  % ie combined 2nd of the two using only first and last char?
          case 'rpw-bia-scpot-10-seconds'
            subDir = 'bia-scpot-10-seconds';  % Locally produced files at IRFU (may have to change in future if official)
          case 'rpw-bia-scpot'
            subDir = 'bia-scpot'; % Locally produced files at IRFU (may have to change in future if official)
          case 'rpw-bia-efield-10-seconds'
            subDir = 'bia-efield-10-seconds';  % Locally produced files at IRFU (may have to change in future if official)
          case 'rpw-bia-efield'
            subDir = 'bia-efield'; % Locally produced files at IRFU (may have to change in future if official)
          otherwise
            % fallback to full descriptor (used for local SOAR copy at IRFU)
            subDir = C{3};
        end
        fileDir = [levelDir, filesep, subDir];
      else
        % Keep it ("HK", "L1R" etc. as these do not have separate subfolders based on descriptor)
        fileDir = levelDir;
      end
    end % get_fileDir
    
    function p = get_path_to_file(obj,fileName)
      C = strsplit(lower(fileName),'_');
      if strcmpi(fileName(end-3:end),'.cdf')
        d =  C{end-1}; p = obj.dbRoot;
        for ix=1:(length(C)-2), p = [p filesep C{ix}]; end %#ok<AGROW>
        p = [p filesep d(1:4) filesep d(5:6)];
      end
    end % get_path_to_file
  end
  
end

