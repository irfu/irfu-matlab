classdef (Abstract) mms_file_db
  %MMS_FILE_DB  Interface class for MMS file databses
  
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

