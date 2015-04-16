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
  end
  
end

