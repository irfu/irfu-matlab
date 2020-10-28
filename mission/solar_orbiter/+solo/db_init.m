function db_init(key, val)
%SOLO.DB_INIT  initialize SOLO file database
%
%  SOLO.DB_INIT(key,value)
%
%  Example:
%    solo.db_init('local_file_db','/data/solo')
%    solo.db_init('local_file_db','/Users/yuri/Data/solo')
%
%    solo.db_init('local_file_db_irfu','/data/solo/irfu')
%    Database with IRFU produced products
%
%    solo.db_init('db_cache_enabled',true)
%    solo.db_init('db_cache_size_max',4096) % set cache to 4GB
%    solo.db_init('db_index_enabled',true) % enable use of database index


global SOLO_DB;
if isempty(SOLO_DB), SOLO_DB = solo_db; end

narginchk(0,2)

if nargin==0
  localFileDbRoot = datastore('solo_db','local_file_db_root');
  if isempty(localFileDbRoot)
    irf.log('warining','local_file_db_root empty - run solo.db_init()')
  else
    localFileDb = solo_local_file_db(localFileDbRoot);
    SOLO_DB.add_db(localFileDb);
    irf.log('notice','Initialized local_file_db')
  end
  localFileDbRoot = datastore('solo_db','local_file_db_irfu_root');
  if ~isempty(localFileDbRoot)
    localFileDb = solo_local_file_db(localFileDbRoot);
    SOLO_DB.add_db(localFileDb);
    irf.log('notice','Initialized local_file_irfu_db')
  end
  cacheEnabled = datastore('solo_db','db_cache_enabled');
  if ~isempty(cacheEnabled)
    SOLO_DB.cache.enabled = cacheEnabled;
    if cacheEnabled, s = 'ON'; else, s = 'OFF'; end
    irf.log('notice',['db_cache: ' s])
  end
  cacheTimeout = datastore('solo_db','db_cache_timeout');
  if ~isempty(cacheTimeout)
    SOLO_DB.cache.timeout = cache.timeout;
    irf.log('notice',sprintf('db_cache_timeout: %d sec',cacheTimeout))
  end
  cacheSizeMax = datastore('solo_db','db_cache_size_max');
  if ~isempty(cacheSizeMax)
    SOLO_DB.cache.cacheSizeMax = cacheSizeMax;
    irf.log('notice',sprintf('db_cache_size_max: %d MB',cacheSizeMax))
  end
  return
end

if nargin ~=2, error('need a pair KEY,VALUE'), end

switch key
  case 'local_file_db'
    localFileDb = solo_local_file_db(val);
    SOLO_DB.add_db(localFileDb);
    datastore('solo_db','local_file_db_root',val);
  case 'local_file_db_irfu'
    localFileDb = solo_local_file_db(val);
    SOLO_DB.add_db(localFileDb);
    datastore('solo_db','local_file_db_irfu_root',val);
  case 'db_index_enabled'
    SOLO_DB.index.enabled = val;
    datastore('solo_db',key,val)
  case 'db_cache_enabled'
    SOLO_DB.cache.enabled = val;
    datastore('solo_db',key,val)
  case 'db_cache_timeout'
    SOLO_DB.cache.timeout = val;
    datastore('solo_db',key,val)
  case 'db_cache_size_max'
    SOLO_DB.cache.cacheSizeMax = val;
    datastore('solo_db',key,val)
  otherwise
    error('Unrecognized KEY')
end
