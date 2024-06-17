classdef mms_db < handle
  %MMS_DB Summary of this class goes here
  %   Detailed explanation goes here

  properties
    databases
    cache
    index
  end

  methods
    function obj=mms_db()
      obj.databases = [];
    end
    function obj = add_db(obj,dbInp)
      if ~isa(dbInp,'mms_file_db')
        error('expecting MMS_FILE_DB input')
      end
      if any(arrayfun(@(x) strcmpi(x.id,dbInp.id), obj.databases))
        irf.log('warning',['Database [' dbInp.id '] already added'])
        return
      end
      obj.databases = [obj.databases dbInp];
      obj.cache = mms_db_cache();
      obj.cache.enabled = false;
      obj.index.enabled = false;
    end


    function fileList = list_files(obj,filePrefix,tint,varName)
      if nargin < 4, varName = ''; end
      fileList =[];
      if nargin==2, tint =[]; end
      if isempty(obj.databases)
        irf.log('warning','No databases initialized'), return
      end
      for iDb = 1:length(obj.databases)
        db = obj.databases(iDb);
        if mms.db_index && ~isempty(obj.databases(iDb).index)
          if isempty(varName)
            flTmp = db.index.find_files('fileprefix',filePrefix,'tint',tint);
          else
            flTmp = db.index.find_files('fileprefix',filePrefix,...
              'varname',varName,'tint',tint);
          end
        else
          flTmp =  db.list_files(filePrefix,tint);
        end
        fileList = [fileList flTmp]; %#ok<AGROW>
      end
    end

    function res = get_variable(obj,filePrefix,varName,tint)
      narginchk(4,4)
      res = [];
      fileList = list_files(obj,filePrefix,tint,varName);
      if isempty(fileList), return, end

      loadedFiles = obj.load_list(fileList,varName);
      if numel(loadedFiles)==0, return, end

      flagDataobj = isa(loadedFiles{1},'dataobj');
      for iFile = 1:length(loadedFiles)
        if flagDataobj, append_sci_var(loadedFiles{iFile})
        else, append_ancillary_var(loadedFiles{iFile});
        end
      end

      function append_ancillary_var(ancData)
        if isempty(ancData), return, end
        if ~isstruct(ancData) || ~(isfield(ancData,varName) ...
            && isfield(ancData,'time'))
          error('Data does not contain %s or time',varName)
        end
        time = ancData.time; data = ancData.(varName);
        if isempty(res), res = struct('time',time,varName,data); return, end
        res.time = [res.time; time];
        res.(varName) = [res.(varName); data];
        % check for overlapping time records and remove duplicates
        [~, idxSort] = sort(res.time);
        [res.time, idxUniq] = unique(res.time(idxSort));
        irf.log('warning',...
          sprintf('Discarded %d data points',length(idxSort)-length(idxUniq)))
        res.(varName) = res.(varName)(idxSort(idxUniq),:);
      end

      function append_sci_var(sciData)
        if isempty(sciData), return, end
        if ~isa(sciData,'dataobj')
          error('Expecting DATAOBJ input')
        end
        v = get_variable(sciData,varName);

        if isempty(v)
          irf.log('waring','Empty return from get_variable()')
          return
        end
        if ~isstruct(v) || ~(isfield(v,'data') && isfield(v,'DEPEND_0'))
          error('Data does not contain DEPEND_0 or DATA')
        end

%        if v.nrec == 1 && contains(v.CATDESC, 'FPI') % this is a quick fix of the problem of having PDist burst files with nrec =1 remove when fixed by FPI people
%          irf.log('waring','PDist with nrec = 1')
%          return
%        end


        if isempty(res), res = v; return, end
        if iscell(res), res = [res {v}]; return, end
        if ~comp_struct(res,v), res = [{res}, {v}]; return, end
        if res.variance(1)=='F' % not varying in time
          if isequal(res.data,v.data), return, end
          error('Static (variance=F/) variable changing between files')
        end

        % append data
        res.data = [res.data; v.data];
        % append depend variables
        n_dep = sum(contains(fields(res),'DEPEND_'))-1;
        for idep = 0:n_dep
          DEP_str = ['DEPEND_' num2str(idep)];
          if v.(DEP_str).nrec == v.nrec % check if depend is a timeseries, if yes, then append
            res.(DEP_str).data = [res.(DEP_str).data; v.(DEP_str).data];
          end
        end

        % check for overlapping time records
        [~,idxUnique] = unique(res.DEPEND_0.data);
        idxDuplicate = setdiff(1:length(res.DEPEND_0.data), idxUnique);
        res.data(idxDuplicate, :, :, :, :, :, :, :, :, :, :, :) = [];
        for idep = 0:n_dep
          DEP_str = ['DEPEND_' num2str(idep)];
          if v.(DEP_str).nrec == v.nrec
            res.(DEP_str).data(idxDuplicate, :, :, :, :, :, :, :, :, :, :, :) = [];
          end
        end
        nDuplicate = length(idxDuplicate);
        if nDuplicate
          irf.log('warning',sprintf('Discarded %d data points',nDuplicate))
        end

        % update number of records, nrec
        res.nrec = length(res.DEPEND_0.data);
        res.DEPEND_0.nrec = res.nrec;
        for idep = 1:n_dep
          DEP_str = ['DEPEND_' num2str(idep)];
          if size(res.(DEP_str).data,1) == res.nrec
            res.(DEP_str).nrec = res.nrec;
          end
        end

        % sort data
        [res.DEPEND_0.data,idxSort] = sort(res.DEPEND_0.data);
        res.data = res.data(idxSort, :, :, :, :, :, :, :, :, :, :, :);
        for idep = 1:n_dep
          DEP_str = ['DEPEND_' num2str(idep)];
          if v.(DEP_str).nrec == v.nrec
            res.(DEP_str).data(idxSort, :, :, :, :, :, :, :, :, :, :, :);
          end
        end
        function res = comp_struct(s1,s2)
          % Compare structures
          narginchk(2,2), res = false;

          if ~isstruct(s1) ||  ~isstruct(s2), error('expecting STRUCT input'), end
          if isempty(s1) && isempty(s2), res = true; return
          elseif xor(isempty(s1),isempty(s2)), return
          end

          fields1 = fields(s1); fields2 = fields(s2);
          if ~comp_cell(fields1,fields2), return, end

          ignoreFields = {'data','nrec','Generation_date',...
            'GlobalAttributes','Logical_file_id','Data_version','Parents', ...
            'VALIDMAX'};
          for iField=1:length(fields1)
            f = fields1{iField};
            % data, nrec and the GlobalAttributes Generation_date,
            % Logical_file_id and Data_version will almost always differ
            % between files.
            if ~isempty(intersect(f,ignoreFields)), continue, end
            if isnumeric(s1.(f)) || ischar(s1.(f))
              if ~all(all(all(s1.(f)==s2.(f)))), return, end
            elseif isstruct(s1.(f)), if ~comp_struct(s1.(f),s2.(f)), return, end
            elseif iscell(s1.(f)), if ~comp_cell(s1.(f),s2.(f)), return, end
            else
              error('cannot compare : %s',f)
            end
          end
          res = true;
        end % COMP_STRUCT
        function res = comp_cell(c1,c2)
          %Compare cells
          narginchk(2,2), res = false;

          if ~iscell(c1) ||  ~iscell(c2), error('expecting CELL input'), end
          if isempty(c1) && isempty(c2), res = true; return
          elseif xor(isempty(c1),isempty(c2)), return
          end
          if ~all(size(c1)==size(c2)), return, end

          [n,m] = size(c1);
          if(m==1), c1=sort(c1); c2=sort(c2); end
          for iN = 1:n
            for iM = 1:m
              if ischar(c1{iN, iM}) && ischar(c2{iN,iM})
                if ~strcmp(c1{iN, iM},c2{iN,iM}), return , end
              elseif iscell(c1{iN, iM}) && iscell(c2{iN,iM})
                if ~comp_cell(c1{iN, iM},c2{iN,iM}), return , end
              else
                irf.log('warning','can only compare chars')
                res = true; return
              end

            end
          end
          res = true;
        end % COMP_CELL
      end % APPEND_SCI_VAR
    end % GET_VARIABLE

    function res = load_list(obj,fileList,mustHaveVar)
      narginchk(2,3), res = {};
      if isempty(fileList), return, end
      if nargin==2, mustHaveVar = ''; end

      for iFile=1:length(fileList)
        fileToLoad = fileList(iFile);
        if mms.db_index
          fileNameToLoad = fileToLoad{1};
        else
          fileNameToLoad = fileToLoad.name;
        end
        dobjLoaded = obj.cache.get_by_key(fileNameToLoad);
        if isempty(dobjLoaded)
          if mms.db_index
            db=obj.databases;
          else
            db = obj.get_db(fileToLoad.dbId);
          end
          if isempty(db) || ~db.file_has_var(fileNameToLoad,mustHaveVar)
            continue
          end
          dobjLoaded = db.load_file(fileNameToLoad);
          obj.cache.add_entry(fileNameToLoad,dobjLoaded)
        end
        res = [res {dobjLoaded}]; %#ok<AGROW>
      end
    end

    function res = get_db(obj,id)
      idx = arrayfun(@(x) strcmp(x.id,id),obj.databases);
      res = obj.databases(idx);
    end

    function res = get_ts(obj,filePrefix,varName,tint)
      narginchk(4,4)
      res = [];
      v = get_variable(obj,filePrefix,varName,tint);
      if isempty(v), return, end
      if numel(v)==1
        res = mms.variable2ts(v);
        res = res.tlim(tint);
      else
        res = cell(1,numel(v));
        for iV = 1:numel(v)
          resTmp = mms.variable2ts(v{iV});
          res{iV} = resTmp.tlim(tint);
        end
      end
    end
  end
end