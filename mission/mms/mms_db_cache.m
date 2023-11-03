classdef mms_db_cache<handle
  %MMS_DB_CACHE Caching of files

  properties
    timeout = 0          % db cache timeout in sec, off by default
    cacheSizeMax = 1024  % db cache max size in mb
  end

  properties (Dependent = true)
    enabled              % ON/OFF flag
  end

  properties (Access = private)
    enabled_ = true;
    names
    loaded
    data
  end

  methods
    function obj = mms_db_cache(szMB,tout)
      % obj = mms_db_cache([szMB,toutSEC])
      if nargin<1, return, end
      if nargin<2 && ~isempty(szMB), obj.cacheSizeMax = szMB; end
      if nargin==2 && ~isempty(tout), obj.timeout = tout; end
    end

    function res = get.enabled(obj)
      res = obj.enabled_;
    end

    function set.enabled(obj,value)
      if numel(value) ~=1 || ~islogical(value)
        error('expecting logical value (true/false)')
      end
      obj.enabled_ = value;
      if ~value % cleanup on turn off
        obj.loaded = []; obj.names = []; obj.data = [];
      end
    end

    function set.timeout(obj,value)
      if numel(value) ~=1 || ~isnumeric(value) || value<0
        error('expecting a numerical value >=0 (seconds)')
      end
      obj.timeout = value;
    end

    function set.cacheSizeMax(obj,value)
      if numel(value) ~=1 || value<=0
        error('expecting a positive numerical value (MB)')
      end
      if value>10*1024, warning('cache size > 10 GB'), end
      obj.cacheSizeMax = value;
    end

    function res = get_by_key(obj,key)
      % res = cache.get_by_key(key)
      res = [];
      if ~obj.enabled || isempty(obj.names), return, end
      idx = cellfun(@(x) strcmp(key,x),obj.names);
      if ~any(idx), return, end
      if (obj.timeout>0) && (now() > obj.loaded(idx) + obj.timeout/86400)
        obj.purge(), return
      end
      res = obj.data{idx};
    end

    function add_entry(obj,key,dataEntry)
      % cache.add_entry(key,dataEntry)
      % key - string
      % dataEntry - number or structure
      if ~obj.enabled, return, end
      if ~isempty(obj.names)
        obj.purge();
        idx = cellfun(@(x) strcmp(key,x),obj.names);
        if any(idx), return, end % file is already there
      end
      if isempty(obj.names)
        obj.loaded = now(); obj.names = {key}; obj.data = {dataEntry};
        return
      end
      obj.loaded = [obj.loaded now()];
      obj.names = [obj.names {key}];
      obj.data = [obj.data {dataEntry}];
      % check if we did not exceed cacheSizeMax
      cacheTmp = obj.data; w = whos('cacheTmp'); t0 = now(); %#ok<NASGU>
      while w.bytes > obj.cacheSizeMax*1024*1024
        if length(obj.loaded) == 1, break, end
        dt = t0 - obj.loaded; idx = dt == max(dt);
        disp(['purging ' obj.names{idx}])
        obj.loaded(idx) = []; obj.names(idx) = []; obj.data(idx) = [];
        cacheTmp = obj.data; %#ok<NASGU>
        w = whos('cacheTmp');
      end
    end

    function disp(obj)
      % Display cache memory usage and contents
      if ~obj.enabled, disp('DB caching disabled'), return, end
      if isempty(obj.names), disp('DB cache empty'), return, end
      cacheTmp = obj.data; w = whos('cacheTmp'); t0 = now(); %#ok<NASGU>
      fprintf('DB cache using %.1fMB of %dMB \n',...
        w.bytes/1024/1024,obj.cacheSizeMax)
      if obj.timeout==0, fprintf('Timeout : off sec\nEntries : \n')
      else
        fprintf('Timeout : %d sec\nEntries : \n',obj.timeout)
      end
      dt = t0 - obj.loaded; [dt,idx] = sort(dt);
      for i=1:length(idx)
        if obj.timeout==0, fprintf('''%s''\n', obj.names{idx(i)})
        else
          fprintf('''%s'' (expires in %d sec)\n', obj.names{idx(i)},...
            ceil(obj.timeout - dt(i)*86400))
        end
      end
    end

    function purge(obj)
      if ~obj.enabled || obj.timeout==0 || isempty(obj.names), return, end
      % purge old entries from cache
      t0 = now;
      idx = t0 > obj.loaded + obj.timeout/86400;
      obj.loaded(idx) = []; obj.names(idx) = []; obj.data(idx) = [];
    end
  end

end

