classdef mms_db < handle
  %MMS_DB Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    databases
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
    end
    
   function fileList = list_files(obj,filePrefix,tint)
     fileList =[];
     for iDb = 1:length(obj.databases)
       fileList = [fileList obj.databases(iDb).list_files(filePrefix,tint)]; %#ok<AGROW>
     end
   end
   
   function res = get_variable(obj,filePrefix,varName,tint)
     narginchk(4,4)
     res = [];
     
     fileList = list_files(obj,filePrefix,tint);
     if isempty(fileList), return, end
     
     loadedFiles = obj.load_list(fileList,varName);
     if numel(loadedFiles)==0, return, end
     
     flagDataobj = isa(loadedFiles{1},'dataobj');
     if flagDataobj, res = {}; end
     for iFile = 1:length(loadedFiles)
       if flagDataobj
         res = [res get_variable(loadedFiles{iFile},varName)]; %#ok<AGROW>
       else append_ancillary_var(loadedFiles{iFile});
       end
     end
     function append_ancillary_var(ancData)
       if isempty(ancData), return, end
       if ~(isfield(ancData,varName) && isfield(ancData,'time'))
         error('Data does not contain %s or time',varName)
       end
       time = ancData.time; data = ancData.(varName);
       if isempty(res), res = struct('time',time,varName,data); return, end
       res.time = [res.time; time];
       res.(varName) = [res.(varName); data];
       % check for overlapping time records
       [~,idxUnique] = unique(res.time); 
       idxDuplicate = setdiff(1:length(res.time), idxUnique);
       res.time(idxDuplicate) = []; res.(varName)(idxDuplicate) = [];
     end
   end
   
   function res = load_list(obj,fileList,mustHaveVar)
     narginchk(2,3), res = {};
     if isempty(fileList), return, end
     if nargin==2, mustHaveVar = ''; end
     
     for iFile=1:length(fileList)
       fileToLoad = fileList(iFile);
       db = obj.get_db(fileToLoad.dbId);  
       if isempty(db) || ~db.file_has_var(fileToLoad.name,mustHaveVar)
         continue
       end
       res = [res {db.load_file(fileToLoad.name)}]; %#ok<AGROW>
     end
   end
   
   function res = get_db(obj,id)
     idx = arrayfun(@(x) strcmp(x.id,id),obj.databases);
     res = obj.databases(idx);
   end
  end
  
end

