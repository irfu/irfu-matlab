classdef (Abstract) mms_file_db
  %UNTITLED2 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties (SetAccess = immutable)
    id 
  end
  
  methods
    function obj = mms_file_db(id)
      obj.id = id;
    end
    fileList = list_files(obj,filePrefix,tint)
    dataObj = load_file(obj,fileName)
  end
  
end

