classdef solo_db < handle
  %SOLO_DB Summary of this class goes here
  %   Detailed explanation goes here

  properties
    databases
    cache
    index
  end

  methods



    function obj = solo_db()
      obj.databases = [];
    end



    function obj = add_db(obj,dbInp)
      if ~isa(dbInp,'solo_file_db')
        error('expecting SOLO_FILE_DB input')
      end
      if any(arrayfun(@(x) strcmpi(x.id,dbInp.id), obj.databases))
        irf.log('warning',['Database [' dbInp.id '] already added'])
        return
      end
      obj.databases = [obj.databases dbInp];
      obj.cache = solo_db_cache();
      obj.cache.enabled = false;
      obj.index.enabled = false;
    end



    function fileList = list_files(obj,filePrefix,tint,varName)
      % Return array of data structures describing matching files (datasets).
      if nargin < 4, varName = ''; end
      fileList = [];
      if nargin==2, tint =[]; end
      if isempty(obj.databases)
        irf.log('warning','No databases initialized'), return
      end
      for iDb = 1:length(obj.databases)
        db = obj.databases(iDb);
        if solo.db_index && ~isempty(obj.databases(iDb).index)
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



    function get_metakernel(obj, flown_or_predicted)
      % Loop through all databases (if mounted local irfu-data as well).

      found = false;
      for iDb = 1:length(obj.databases)
        db = obj.databases(iDb);
        dbRoot = db.dbRoot;
        if exist([dbRoot, filesep, 'SPICE'], 'dir')
          % If root folder contains a subfolder "SPICE", then it is most
          % likely the one database we want.
          prevSpicePath = datastore('spice_paths', 'solarorbiter');
          if isempty(prevSpicePath)
            datastore('spice_paths', 'solarorbiter', [dbRoot, filesep, 'SPICE']);
          end
          found = true;
          load_mkernel('solarorbiter', flown_or_predicted);
          break;
        end
      end
      if ~found, irf.log('critical', 'Did not find any SPICE kernel paths, please try manually loading it using new function "load_mkernel(''solarorbiter'')"'); end
    end



    function res = get_variable(obj,filePrefix,varName,tint)
      narginchk(4,4)
      res = [];
      fileList = list_files(obj,filePrefix,tint,varName);
      if isempty(fileList), return, end

      loadedFiles = obj.load_list(fileList,varName);
      if numel(loadedFiles)==0, return, end

      for iFile = 1:length(loadedFiles)
        append_sci_var(loadedFiles{iFile})
      end

      % ========================================================================

      function append_sci_var(sciData)
        % Nested function
        % Append data from "sciData" to "res".

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

        if isempty(res), res = v; return, end

        % NOTE: If "res" is a cell array, then always add next time series in
        % cell array, without checking whether this or the immediately
        % preceeding time series are "compatible" (have consistent metadata).
        % Bug or feature? /EJ 2023-07-06
        if iscell(res), res = [res {v}]; return, end
        if ~comp_struct(res,v)
          % CASE: zVariable data "v" appears to be/might be different from the
          %       previous datasets.

          % Switch from res=TSeries, to cell array of multiple TSeries.
          res = [{res}, {v}];
          return
        end
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

        % ======================================================================

        function res = comp_struct(s1,s2)
          % Compare structures representing zVariable data, including zVariable
          % attributes. Return whether corresponding zVariable attributes appear
          % to be identical/the same or not.

          narginchk(2,2)

          % Default return value. Assume structs are different until (almost)
          % proven to be equal.
          res = false;

          if ~isstruct(s1) || ~isstruct(s2), error('expecting STRUCT input'), end
          if isempty(s1) && isempty(s2), res = true; return
          elseif xor(isempty(s1),isempty(s2)), return
          end

          fields1 = fields(s1); fields2 = fields(s2);
          if ~comp_cell(fields1,fields2)
            return
          end

          % "data", nrec, and the global attributes (in GlobalAttributes) named
          % Generation_date, Logical_file_id, Data_version, and Parents will
          % almost always differ between files. Therefore not comparing those
          % fields.
          ignoreFields = {'data','nrec','Generation_date',...
            'GlobalAttributes','Logical_file_id','Data_version','Parents', ...
            'SCALEMAX','SCALEMIN','VALIDMAX'};
          for iField=1:length(fields1)
            f = fields1{iField};
            if ~isempty(intersect(f,ignoreFields)), continue, end
            if isnumeric(s1.(f)) && all(isnan(s1.(f))) && all(isnan(s2.(f))), continue, end
            if isnumeric(s1.(f)) || ischar(s1.(f))
              if ~all(all(all(s1.(f)==s2.(f)))), return, end
            elseif isstruct(s1.(f)), if ~comp_struct(s1.(f),s2.(f)), return, end
            elseif iscell(s1.(f)), if ~comp_cell(s1.(f),s2.(f)), return, end
            else
              error('cannot compare : %s',f)
            end
          end
          res = true;
        end % function COMP_STRUCT

        % ======================================================================

        function res = comp_cell(c1,c2)
          % Compare cell arrays of text.

          narginchk(2,2)

          % Default return value. Assume cells are different until (almost)
          % proven to be equal.
          res = false;

          if ~iscell(c1) || ~iscell(c2), error('expecting CELL input'), end
          if isempty(c1) && isempty(c2), res = true; return
          elseif xor(isempty(c1),isempty(c2)), return
          end
          if ~all(size(c1)==size(c2))
            % CASE: Different number of zVariable attributes.
            return
          end

          [n,m] = size(c1);
          if(m==1), c1=sort(c1); c2=sort(c2); end
          for iN = 1:n
            for iM = 1:m
              if ischar(c1{iN, iM}) && ischar(c2{iN,iM})
                if ~strcmp(c1{iN, iM},c2{iN,iM}), return , end
              elseif iscell(c1{iN, iM}) && iscell(c2{iN,iM})
                % NOTE: RECURSIVE CALL
                if ~comp_cell(c1{iN, iM},c2{iN,iM})
                  return
                end
              else
                irf.log('warning','can only compare chars')
                res = true; return
              end

            end
          end
          res = true;
        end % function COMP_CELL

        % ======================================================================

      end % function APPEND_SCI_VAR

      % ========================================================================

    end % function GET_VARIABLE



    function res = load_list(obj,fileList,mustHaveVar)
      narginchk(2,3), res = {};
      if isempty(fileList), return, end
      if nargin==2, mustHaveVar = ''; end

      for iFile=1:length(fileList)
        fileToLoad = fileList(iFile);
        if solo.db_index
          fileNameToLoad = fileToLoad{1};
        else
          fileNameToLoad = fileToLoad.name;
        end
        dobjLoaded = obj.cache.get_by_key(fileNameToLoad);
        if isempty(dobjLoaded)
          if solo.db_index
            db = obj.databases;
          else
            db = obj.get_db(fileToLoad.dbId);
          end
          if isempty(db) || ~db.file_has_var([fileToLoad.path, filesep, fileNameToLoad], mustHaveVar)
            continue
          end
          dobjLoaded = db.load_file([fileToLoad.path, filesep, fileNameToLoad]);
          obj.cache.add_entry(fileNameToLoad,dobjLoaded)
        end
        res = [res {dobjLoaded}]; %#ok<AGROW>
      end
    end % function load_list



    function res = get_db(obj,id)
      idx = arrayfun(@(x) strcmp(x.id,id),obj.databases);
      res = obj.databases(idx);
    end % function get_db



    function res = get_ts(obj,filePrefix,varName,tint)
      % See solo.db_get_ts() (a wrapper around this method).

      narginchk(4,4)
      res = [];
      v = get_variable(obj,filePrefix,varName,tint);
      if isempty(v), return, end
      if numel(v)==1
        res = solo.variable2ts(v);
        res = res.tlim(tint);
      else
        res = cell(1,numel(v));
        for iV = 1:numel(v)
          resTmp = solo.variable2ts(v{iV});
          res{iV} = resTmp.tlim(tint);
        end
      end
    end % function get_ts



  end % methods

end
