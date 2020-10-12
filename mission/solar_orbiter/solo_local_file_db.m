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
      %  fileList = list_files(SOLO_DB, 'solo1_edp_comm_l1b_dce128');
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
        rDir = get_remotePrefix(obj, filePrefix);
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
            if strcmp(C{2}, 'L2')
              % L2
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
        rDir = get_remotePrefix(obj, filePrefix);
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
            if strcmp(C{2}, 'L2')
              % L2
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
        if is_version_larger(fnd.vXYZ, fileList(iSame).ver)
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
    function res = load_file(obj,fileName)
      narginchk(2,3)
      
      irf.log('notice',['loading ' fileName])
      if solo.db_index
        fileNameFullPath = fileName;
      else
        p = obj.get_path_to_file(fileName);
        fileNameFullPath = [p filesep fileName];
      end
      res = dataobj(fileNameFullPath);
    end % LOAD_FILES
    
    %% FILE_HAS_VAR
    function res = file_has_var(obj,fileName,varName)
      % checks if fileName includes variable name varName
      % res = true/false
      narginchk(3,3)
      res = false; if isempty(varName) || isempty(fileName), return, end
      
      entryTmp = obj.cache.get_by_key(fileName);
      if ~isempty(entryTmp)
        res = any(cellfun(@(x) strcmp(x,varName), entryTmp.vars(:,1)));
        return
      end
      if solo.db_index
        irf.log('notice','Using index to check if file ok');
        fullPath = fileName;
      else
        p = obj.get_path_to_file(fileName);
        fullPath = [p filesep fileName];
      end
      if ~exist(fullPath,'file')
        irf.log('warning', ['Fies does not exist: ' fullPath])
        return
      end
      
      % cdf
      if solo.db_index
        res = obj.index.file_has_var(fileName,varName);
      else
        info = spdfcdfinfo(fullPath);
        if ispc
          % Add a very short delay to ensure consecutive files are not
          % accessed TOO quickly as this may cause Matlab to experince a
          % hard crash on Win10 regardless of the try&catch.
          pause(0.0001);
        end
        res = any(cellfun(@(x) strcmp(x,varName), info.Variables(:,1)));
      end
    end
  end
  
  methods (Access=private)
    function rDir = get_remotePrefix(obj, filePrefix)
      if contains(filePrefix, 'rpw')
        % RPW data is keept in one separate sync folder at IRFU
        rDir = [obj.dbRoot, filesep, 'remote', filesep, 'data'];
      else
        % All other instruments are kept in a "soar" sync folder
        rDir = [obj.dbRoot, filesep, 'soar'];
      end
    end % get_remotePrefix
    
    function fileDir = get_fileDir(~, C)
      fileDir = C{2}; % "L2" (or "L1R", "L1", "L3", "HK")
      if isequal(fileDir, 'L2')
        switch C{3}
          case 'rpw-lfr-surv-asm-cdag'
            fileDir = [fileDir, filesep, 'lfr_asm'];
          case 'rpw-tds-surv-hist1d-cdag'
            fileDir = [fileDir, filesep, 'hist1d'];
          case 'rpw-tds-surv-hist2d-cdag'
            fileDir = [fileDir, filesep, 'hist2d'];
          case 'rpw-tds-surv-mamp-cdag'
            fileDir = [fileDir, filesep, 'mamp'];
          case 'rpw-tds-surv-stat-cdag'
            fileDir = [fileDir, filesep, 'stat'];
          case {'rpw-lfr-surv-bp1-cdag', 'rpw-lfr-surv-bp2-cdag'}
            fileDir = [fileDir, filesep, 'bp'];
          case {'rpw-lfr-surv-cwf-b-cdag', 'rpw-lfr-surv-swf-b-cdag'}
            fileDir = [fileDir, filesep, 'lfr_wf_b'];
          case {'rpw-lfr-surv-cwf-e-cdag', 'rpw-lfr-surv-swf-e-cdag'}
            fileDir = [fileDir, filesep, 'lfr_wf_e'];
          case {'rpw-tds-surv-rswf-b-cdag', 'rpw-tds-surf-tswf-b-cdag'}
            fileDir = [fileDir, filesep, 'tds_wf_b'];
          case {'rpw-tds-surv-rswf-e-cdag', 'rpw-tds-surf-tswf-e-cdag'}
            fileDir = [fileDir, filesep, 'tds_wf_e'];
          case {'rpw-hfr-surv-cdag', 'rpw-tnr-surv-cdag'}
            fileDir = [fileDir, filesep, 'thr'];
          otherwise
            % Not yet implemented
            errS = 'Not yet implemented!';
            irf.log('critical', errS);
            error(errS);
        end
      else
        % Keep it ("HK", "L1R" etc. as these do not have separate subfolders based on descriptor)
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

