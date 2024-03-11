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
      % List files from Database "obj", which match "filePrefix" and cover
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
          if year==TStop.year,  moStop  = TStop.month;  end
          for mo = moStart:moStop
            moDir = sprintf('%s%s%s%s%d%s%02d',rDir,filesep,fDir,filesep,year,filesep,mo);
            curDir = moDir;
            if ismember(C{2}, {'L2', 'L3'})
              % L2, L3 have monthly subfolders.
              dPref = sprintf('%s_%d%02d',filePrefix,year,mo);
              limited_sci_list;
            else
              % L1, L1R, HK etc have daily subfolders.
              dStart = 1; dStop = 31;
              if year==TStart.year && mo==TStart.month, dStart=TStart.day; end
              if year==TStop.year  && mo==TStop.month,  dStop = TStop.day; end
              for day = dStart:dStop
                curDir = [moDir filesep sprintf('%02d',day)];
                dPref = sprintf('%s_%d%02d%02d',filePrefix,year,mo,day);
                limited_sci_list;
              end
            end
          end
        end    % for

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
          listingD = dir([fullfile(curDir, dPref) '*.cdf']); % SolO have only latest file of each type.
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
            stopFile  = [filePrefix, '_', tint.stop.toUtc(dateFormat),  '_V99.cdf'];
          end
          % Find index of files with names which timewise are sorted
          % between our "startFile" and "stopFile" names.
          indAfterStart = arrayfun(@(x) isequal({startFile; x.name}, sort({startFile; x.name})), listingD);
          indBeforeStop = arrayfun(@(x) isequal({x.name;  stopFile}, sort({x.name;  stopFile})), listingD);
          % Also look at files just before and after as there might be some
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
        end    % limited_sci_list
      end    % list_sci_tint

      %% LIST SCI
      %
      % What does this do? Combine all available time intervals?
      % Cf list_sci_tint().
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
  end

  methods (Access=private)

    % Return ONE directory tree root depending on instrument.
    function rDir = get_remotePrefix(obj, C)
      % Descriptor contains instrument and data product descriptor part
      % separated by "-".

      temp       = strsplit(C{3}, '-');
      instrument = temp{1};
      dbRoot     = obj.dbRoot;

      if strcmp(instrument, 'rpw')
        %==============================
        % CASE: Searching for RPW data
        %==============================
        if exist(fullfile(dbRoot, 'latest'), 'dir')
          % CASE: obj.dbRoot has subdirectory "latest/".
          %       ==> Use RPW BIAS data (L2, L3) processed at IRFU.
          %
          % Ex: obj.dbRoot == /data/solo/data_irfu/
          %               ==> /data/solo/data_irfu/latest/rpw/

          rDir = fullfile(dbRoot, 'latest', 'rpw');
        else
          % CASE: obj.dbRoot DOES NOT have subdirectory "latest/".
          %       ==> Use RPW data (all subsystems) mirrored from ROC/LESIA.
          %
          % Ex: obj.dbRoot == /data/solo/
          %               ==> /data/solo/remote/data/
          rDir = fullfile(dbRoot, 'remote', 'data');
        end
      else
        %==================================
        % CASE: Searching for non-RPW data
        %==================================
        rDir = fullfile(dbRoot, 'soar', instrument);
        if exist(rDir, 'dir')
          % CASE: obj.dbRoot has subdirectory "soar".
          %       ==> Use (presumed) SOAR mirror.
          %
          % Ex: obj.dbRoot == /data/solo/
          %               ==> /data/solo/soar/<instr>/
          return
        end

        rDir = fullfile(dbRoot, instrument);
        if exist(rDir, 'dir')
          % CASE: obj.dbRoot has subdirectory named after instrument.
          %       ==> obj.dbRoot is a general folder for (multiple) non-RPW
          %           instruments.
          %
          % Ex: obj.dbRoot == /data/solo/data_manual/
          %               ==> /data/solo/data_manual/<instr>/
          return
        end
      end
    end % get_remotePrefix

    function fileDir = get_fileDir(~, C)
      levelDir = C{2}; % "L2" (or "L1R", "L1", "L3", "HK")

      % Normalization. Remove '-cdag', if present.
      descr = regexprep(C{3}, '-cdag$', '');

      if ismember(levelDir, {'L2', 'L3'})
        switch descr
          case 'rpw-lfr-surv-asm'
            subDir = 'lfr_asm'; % ie combined 2nd "_" 4th
          case 'rpw-tds-surv-hist1d'
            subDir = 'tds_hist1d';  % ie 4th
          case 'rpw-tds-surv-hist2d'
            subDir = 'tds_hist2d';  % ie 4th
          case 'rpw-tds-surv-mamp'
            subDir = 'tds_mamp';    % ie 4th
          case 'rpw-tds-surv-stat'
            subDir = 'tds_stat';    % ie 4th
          case {'rpw-lfr-surv-bp1', 'rpw-lfr-surv-bp2'}
            subDir = 'lfr_bp';    % ie combined 2nd "_" 4th (excl last digit, which is unique)
          case {'rpw-lfr-surv-cwf-b', 'rpw-lfr-surv-swf-b'}
            subDir = 'lfr_wf_b';  % ie combined 2nd "_" 4th and 5th (excl first char of 4th, which is unqiue)
          case {'rpw-lfr-surv-cwf-e', 'rpw-lfr-surv-cwf-e-1-second', 'rpw-lfr-surv-swf-e'}
            % NOTE: 'rpw-lfr-surv-cwf-e-1-second' refers to an unofficial
            % dataset (DATASET_ID) only used internally at IRF. /2021-05-21
            subDir = 'lfr_wf_e';  % ie combined 2nd "_" 4th and 5th (excl first char of 4th, which is unqiue)
          case {'rpw-tds-surv-rswf-b', 'rpw-tds-surv-tswf-b'}
            subDir = 'tds_wf_b';  % ie combined 2nd "_" 4th and 5th (excl first two chars of 4th, of which the first one is unqiue)
          case {'rpw-tds-surv-rswf-e', 'rpw-tds-surv-tswf-e'}
            subDir = 'tds_wf_e';  % ie combined 2nd "_" 4th and 5th (excl first two chars of 4th, of which the first one is unqiue)
          case {'rpw-hfr-surv', 'rpw-tnr-surv'}
            subDir = 'thr';  % ie combined 2nd of the two using only first and last char?
          case {'rpw-tnr-fp'}
            subDir = 'tnr_fp';

            % Official directory names used by ROC and that IRFU should
            % therefore also use. As per agreement with Yuri Khotyaintsev,
            % Thomas Chust, and Erik P G Johansson 2020-11-27.
            % /Erik P G Johansson 2020-12-15.
          case {'rpw-bia-density', 'rpw-bia-density-10-seconds'}
            subDir = 'lfr_density';
          case {'rpw-bia-efield',  'rpw-bia-efield-10-seconds'}
            subDir = 'lfr_efield';
          case {'rpw-bia-scpot',   'rpw-bia-scpot-10-seconds'}
            subDir = 'lfr_scpot';

          otherwise
            % Fallback to full descriptor (used for local SOAR copy at IRFU).
            subDir = descr;
        end
        fileDir = fullfile(levelDir, subDir);
      else
        % Keep it ("HK", "L1R" etc. as these do not have separate subfolders based on descriptor)
        fileDir = levelDir;
      end
    end % get_fileDir

    % UNUSED FUNCTION?!
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

