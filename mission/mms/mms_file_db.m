classdef (Abstract) mms_file_db
  %MMS_FILE_DB  Interface class for MMS file databses
  
  properties (SetAccess = immutable)
    id 
  end
  
  properties
    cache
		index
  end
  
  methods
    function obj = mms_file_db(id)
      obj.id = id;
      obj.cache = mms_db_cache();
      if exist([id filesep 'index_sql'],'file')
        obj.index = mms_db_sql([id filesep 'index_sql']);
      else
        obj.index = [];
      end
    end
    fileList = list_files(obj,filePrefix,tint)
    dataObj = load_file(obj,fileName)
  end
  
end

