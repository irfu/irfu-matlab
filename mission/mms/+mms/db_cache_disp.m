function db_cache_disp()
%MMS.DB_CACHE_DISP  display DB cache contents

global MMS_DB; if isempty(MMS_DB), mms.db_init(), end

MMS_DB.cache.disp()
