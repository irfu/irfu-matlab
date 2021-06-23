function db_cache_disp()
%SOLO.DB_CACHE_DISP  display DB cache contents

global SOLO_DB; if isempty(SOLO_DB), solo.db_init(), end

SOLO_DB.cache.disp()
