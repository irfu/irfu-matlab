classdef (Abstract) mms_file_db
  %UNTITLED2 Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods
    fileList = list_files(obj,filePrefix)
    dataObj = load_file(obj,fileName)
  end
  
end

