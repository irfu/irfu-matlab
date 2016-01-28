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
      if(false)
        % DO NOT CREATE a file "index_sql" when running on SDC in the
        % data path for the entire mission (SDC eqivalent of "/data/mms" on
        % Spis/Brain at IRFU).
        obj.index = mms_db_sql([id filesep 'index_sql']); %#ok<UNRCH>
      else
        obj.index = [];
      end
    end
    fileList = list_files(obj,filePrefix,tint)
    dataObj = load_file(obj,fileName)
  end
  
end

