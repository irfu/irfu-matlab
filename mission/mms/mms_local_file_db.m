classdef mms_local_file_db < mms_file_db
  %MMS_LOCAL_FILE_DB  Local file database for MMS
  %   Class handling a database of local MMS files

  properties (SetAccess = immutable)
    dbRoot
  end

  methods
    function obj = mms_local_file_db(rootPath)
      % Create local database for MMS data located in rootPath.
      % Example
      %   MMS_DB = MMS_LOCAL_FILE_DB('/data/mms');
      if nargin == 0, rootPath = pwd; end
      if (rootPath(end)==filesep), rootPath(end)=[]; end % path only, excluding last filesep

      obj@mms_file_db(rootPath); obj.dbRoot = rootPath;
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
      %  MMS_DB = mms_local_file_db('/data/mms'); % init
      %  fileList = list_files(MMS_DB, 'mms1_edp_comm_l1b_dce128');
      narginchk(2,3)
      fileList = [];
      if nargin==3 && (~isempty(tint) && ~isa(tint,'GenericTimeArray'))
        error('Expecting TINT (GenericTimeArray)')
      elseif nargin==2, tint = [];
      end
      if length(filePrefix) < 3 || ~strcmp(filePrefix(1:3),'mms')
        errStr = 'filePrefix must begin with mms*';
        irf.log('critical',errStr), error(errStr)
      end
      C = strsplit(filePrefix,'_');
      if length(C)<3
        errStr = 'filePrefix too short';
        irf.log('critical',errStr), error(errStr)
      end
      if strcmp(C{2},'ancillary')
        list_ancillary();
        if isempty(fileList) || isempty(tint), return, end
        pick_ancillary();
      else
        if mms.db_index && ~isempty(obj.index)
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
      end
      % END LIST_FILES
      %% PICK ANCILLARY
      function pick_ancillary()
        fileList = fileList(arrayfun(@(x) x.start<=tint.stop,fileList) &...
          arrayfun(@(x) x.stop>=tint.start,fileList));
      end
      %% LIST ANCILLARY
      function list_ancillary()
        fileDir = [obj.dbRoot filesep C{2} filesep C{1} filesep C{3}];
        if exist(fileDir,'dir')~=7, return, end
        filePref = [upper(C{1}) '_' upper(C{3}) '_'];
        if(~isempty(tint) && ...
            ismember(upper(C{3}),{'DEFATT','DEFEPH','DEFERR','DEFQ'}) && ...
            tint.stop.ttns-tint.start.ttns<int64(86400000000000))
          % Tint is set and less than one day. Speed up by only listing and
          % reading files +/- 5 days from this interval. (Ancillary files
          % defatt, defeph, deferr cover at most 3 days)).
          t_start = irf_time(tint.start.ttns-int64(86400000000000*5), 'ttns>doy');
          t_stop = irf_time(tint.stop.ttns+int64(86400000000000*5), 'ttns>doy');
          if(t_start(1)==t_stop(1))
            % Same year
            DOY = t_start(2):t_stop(2);
            YYYY = t_start(1)*ones(size(DOY));
          else
            DOY = t_start(2):366; % Assume maximum, ie leap year
            YYYY = t_start(1)*ones(size(DOY));
            DOY = [DOY, 1:t_stop(2)];
            YYYY = [YYYY, t_stop(1)*ones(1,length(DOY)-length(YYYY))];
          end
          for ii=1:length(DOY)
            listing = dir([fileDir, filesep, filePref, ...
              sprintf('%04d%03d', YYYY(ii), DOY(ii)), '*.V*']);
            if isempty(listing), continue, end
            arrayfun(@(x) add2list(x.name), listing)
          end
        elseif (ismember(upper(C{3}), {'PREDEPH'}) && ~isempty(tint) && ...
            tint.stop.ttns - tint.start.ttns<=int64(86400000000000))
          % Tint is set and less than one day.
          t_start = irf_time(tint.start.ttns, 'ttns>doy');
          t_start_d = str2double(sprintf('%04i%03i', t_start));
          t_stop_d = str2double(sprintf('%04i%03i', irf_time(tint.stop.ttns, 'ttns>doy')));
          % Look for file starting the previous year as well as some
          % predeph cover very long time periods.
          YYYY = [t_start(1)-1, t_start(1)];
          for ii=1:length(YYYY)
            listing = dir([fileDir, filesep, filePref, ...
              sprintf('%04i', YYYY(ii)), '*.V*']);
            if isempty(listing), continue, end
            % Now look for file names which could be interesting (starting
            % at the latest on our end date and ending on or after our start date
            listing = listing(arrayfun(@(x) (str2double(x.name(14:20)) <= t_stop_d && str2double(x.name(22:28)) >= t_start_d), listing));
            if isempty(listing), continue, end
            arrayfun(@(x) add2list(x.name), listing)
          end
        else
          irf.log('warning','THIS MAY TAKE SOME TIME')
          listing = dir([fileDir filesep filePref '*.V*']);
          if isempty(listing), return, end
          arrayfun(@(x) add2list(x.name), listing)
        end

        function add2list(name)
          [~,fName,fExt] = fileparts(name);
          ver = str2double(fExt(3:4));
          Entry = struct('name',name,'ver',fExt(3:4),'start',[],'stop',[],...
            'path',fileDir,'dbId',obj.id);
          if isempty(fileList), fileList = add_ss(Entry); return, end
          hasFile = arrayfun(@(x) ~isempty(strfind(x.name,fName)),fileList);
          if ~any(hasFile), fileList = [fileList add_ss(Entry)]; return, end
          iSame = find(hasFile);
          if length(iSame) > 1, error('multiple files with same name'),end
          if ver>str2double(fileList(iSame).ver)
            fileList(iSame) = add_ss(Entry);
          end
          function e = add_ss(e)
            e.start = get_time('start');
            e.stop = get_time('stop');
            function epoch = get_time(s)
              epoch = []; sss=[];
              if isunix
                cmd = sprintf('grep -m1 -i %s_time %s/%s | awk ''{print $3}''',...
                  s,e.path,e.name);
                [sta,out] = unix(cmd); if sta>0, return, end
                if isempty(out)
                  cmd = sprintf('grep -m1 -i %stime %s/%s | awk ''{print $3}''',...
                    s,e.path,e.name);
                  [sta,out] = unix(cmd); if sta>0 || isempty(out), return, end
                end
              else
                % Windows user, findstr is slower as it "greps" all occurrances, does not stop after first match.
                cmd = sprintf('findstr /b /i %s_time "%s\\%s"', ...
                  s,e.path,e.name);
                [sta,out] = system(cmd); if sta>0, return, end
                if isempty(out)
                  cmd = sprintf('findstr /b /i %stime "%s\\%s"',...
                    s,e.path,e.name);
                  [sta,out] = system(cmd); if sta>0 || isempty(out), return, end
                end
              end
              try
                % Split up doy string YYYY-DOYThh:mm:ss.mmmuuunnn
                % works on YYYY-DOYThh:mm:ss and YYYY-DOY/hh:mm:ss or a
                % combination of these.
                doy5 = sscanf(out,'%4d-%3d%c%2d:%2d:%2f');
                sec = floor(doy5(6));
                msec = floor((doy5(6)-sec)*10^3);
                usec = floor(((doy5(6) - sec)*10^3 - msec)*10^3);
                nsec = floor((((doy5(6) - sec)*10^3 - msec)*10^3 - usec) * 10^3);
                doy8 = [doy5(1), doy5(2), doy5(4), doy5(5), sec, msec, usec, nsec]; % YYYY, DOY, hh, mm, ss, msec, usec, nsec
                sss = irf_time(doy8, 'doy8>ttns');
                epoch = EpochTT(sss);
              catch ME
                try
                  % Try for a second method. This is sligthly less accurate
                  % but should be more stable if "out" contains any
                  % erroneous characters. This is limited to DOY strings in
                  % the format of "YYYY-DOY*hh:mm:ss.mmm" where * is any
                  % char.
                  regStr = regexp(out, ...
                    '(\d{4})-(\d{3}).(\d{2}):(\d{2}):(\d{2}).(\d{1,3})', ...
                    'tokens');
                  doy8 = [str2double(regStr{1}), 0, 0];
                  sss = irf_time(doy8, 'doy8>ttns');
                  epoch = EpochTT(sss);
                catch ME2
                  errStr = ['Error reading times for ancillary file: ', ...
                    e.name, ' got: ', sss, ' from :' out];
                  irf.log('critical', errStr); rethrow(ME2);
                end
              end
            end
          end % ADD_SS
        end % ADD2LIST
      end % LOAD_ANCILLARY

      %% LIST_SCI_TINT
      function list_sci_tint()
        fDir = get_prefix();
        TStart = get_times(tint.start); TStop = get_times(tint.stop);
        dateFormat=[]; startFile=[]; stopFile=[];
        for year = TStart.year:TStop.year
          moStart = 1; moStop = 12;
          if year==TStart.year, moStart = TStart.month; end
          if year==TStop.year, moStop = TStop.month; end
          for mo = moStart:moStop
            moDir = sprintf('%s%s%d%s%02d',fDir,filesep,year,filesep,mo);
            curDir = moDir;
            dStart = 1; dStop = 31;
            if year==TStart.year && mo==TStart.month, dStart=TStart.day; end
            if year==TStop.year && mo==TStop.month, dStop = TStop.day; end
            if strcmp(C{3}, 'brst')
              for day = dStart:dStop
                % BRST files are in daily subdirs
                curDir = [moDir filesep sprintf('%02d',day)];
                dPref = sprintf('%s_%d%02d%02d',filePrefix,year,mo,day);
                limited_sci_list;
              end
            else
              %Fast, slow, srvy, comm
              dPref = sprintf('%s_%d%02d',filePrefix,year,mo);
              limited_sci_list;
            end
          end
        end

        function p = get_prefix()
          p = [obj.dbRoot, filesep, strjoin(C, filesep)];
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
          listingD = mms_find_latest_version_cdf([curDir filesep dPref '*.cdf']);
          if isempty(listingD), return, end
          if isempty(dateFormat)
            % Are we looking for files with 8 or the full 14 digits
            % in the date and time.
            fileformat = regexp(listingD(1).name, '_(?<dateFormat>\d{8,})_v','names');
            switch length(fileformat.dateFormat)
              case 8
                dateFormat = 'yyyymmdd';
              case 14
                dateFormat = 'yyyymmddHHMMSS';
              otherwise
                dateFormat = 'yyyymmddHHMMSS';
            end
            % Create reconstructed file names for our interval
            startFile = [filePrefix, '_', tint.start.toUtc(dateFormat), '_v0.0.0.cdf'];
            stopFile = [filePrefix, '_', tint.stop.toUtc(dateFormat), '_v9999999.999999.999999.cdf'];
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
        fileDir = [obj.dbRoot, filesep, strjoin(C, filesep)];
        if exist(fileDir,'dir')~=7, return, end
        listingY = dir(fileDir); listingY(~[listingY.isdir]) = [];
        for iDir = 1:length(listingY)
          % Loop over years
          dNameY = listingY(iDir).name;
          if length(dNameY)~=4, continue, end
          yyyy = str2double(dNameY);
          if yyyy<2015 || yyyy > 2050, continue, end
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
            curDir = [fileDir filesep dNameY filesep dNameM];
            listingD = mms_find_latest_version_cdf([curDir filesep filePrefix '*.cdf']);
            if isempty(listingD), continue, end
            arrayfun(@(x) add2list_sci(x.name,curDir), listingD)
          end
        end
      end % LIST_SCI
      %% ADD2LIST_SCI
      function add2list_sci(name,curDir)
        fnd = mms_fields_file_info(name);
        Entry = struct('name',name,'ver',fnd.vXYZ,'start',[],'stop',[],...
          'path',curDir,'dbId',obj.id);
        Entry = add_ss(Entry);
        % Check time limits of the file
        if isempty(Entry) || ~isempty(tint) && ...
            (Entry.start>tint.stop || Entry.stop<tint.start)
          return
        end
        if isempty(fileList), fileList = Entry; return, end
        fName = [fnd.scId '_' fnd.instrumentId '_' fnd.tmMode '_' ...
          fnd.dataLevel];
        if ~isempty(fnd.dataType), fName = [fName '_' fnd.dataType]; end
        fName = [fName '_' fnd.date '_'];

        hasFile = arrayfun(@(x) ~isempty(strfind(x.name,fName)),fileList);
        if ~any(hasFile), fileList = [fileList add_ss(Entry)]; return, end
        iSame = find(hasFile);
        if length(iSame)>1, error('multiple files with same name'); end
        if is_version_larger(fnd.vXYZ,fileList(iSame).ver)
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
      if mms.db_index
        fileNameFullPath = fileName;
      else
        p = obj.get_path_to_file(fileName);
        fileNameFullPath = [p filesep fileName];
      end
      if mms_local_file_db.is_cdf_file(fileName)
        res = dataobj(fileNameFullPath);
        return
      end

      % ancillary
      [res,~] = mms_load_ancillary(fileNameFullPath,...
        mms_local_file_db.get_anc_type(fileName));
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
      if mms.db_index
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

      if ~mms_local_file_db.is_cdf_file(fileName) % ancillary
        ANC_VARS.defatt = {'wphase','zra','zdec','zphase','lra','ldec',...
          'lphase','pra','pdec','pphase'};
        ANC_VARS.defeph = {'r','v'};
        ANS_VARS.predeph = ANC_VARS.defeph;
        ANC_VARS.defq = {'quality', 'scale'};
        ANC_VARS.predq = ANC_VARS.defq;
        if ~isempty(intersect(varName,...
            ANC_VARS.(mms_local_file_db.get_anc_type(fileName))))
          res = true;
        end
        return
      end
      % cdf
      if mms.db_index
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
    function p = get_path_to_file(obj,fileName)
      C = strsplit(lower(fileName),'_');
      if strcmpi(fileName(end-3:end),'.cdf')
        d =  C{end-1}; p = obj.dbRoot;
        for ix=1:(length(C)-2), p = [p filesep C{ix}]; end %#ok<AGROW>
        p = [p filesep d(1:4) filesep d(5:6)];
        if strcmpi(C{3},'brst'), p = [p filesep d(7:8)]; end % BRST files are in daily subdirs
      else % ancillary
        p = [obj.dbRoot filesep 'ancillary' filesep C{1} filesep C{2}];
      end
    end
  end

  methods (Static, Access=private)
    function res = is_cdf_file(fileName)
      res = false;
      if length(fileName) < 4, return, end
      res = strcmpi(fileName(end-3:end),'.cdf');
    end
    function res= get_anc_type(fileName)
      res = false;
      C = strsplit(lower(fileName),'_');
      if length(C)<3, return, end
      res = C{2};
    end
  end

end

