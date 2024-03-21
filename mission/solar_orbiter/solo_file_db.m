classdef (Abstract) solo_file_db
  %SOLO_FILE_DB  Interface class for SOLO file databses

  properties (SetAccess = immutable)
    id
  end

  properties
    cache
    index
  end

  methods
    function obj = solo_file_db(id)
      obj.id = id;
      obj.cache = solo_db_cache();
      if exist([id filesep 'index_sql'],'file') && ~verLessThan('matlab', '8.5')
        try
          obj.index = solo_db_sql([id filesep 'index_sql']);
        catch
          obj.index = [];
        end
      else
        obj.index = [];
      end
    end
  end

end

