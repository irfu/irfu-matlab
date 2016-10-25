function db_init(key, val)
%MMS.DB_INIT  initialize MMS file database
%
%  MMS.DB_INIT(key,value)
%  
%  Example:
%    mms.db_init('local_file_db','/data/mms')
%    mms.db_init('local_file_db','/Users/yuri/Data/mms')
%
%    mms.db_init('local_file_db_irfu','/data/mms/irfu')
%    Database with IRFU produced products
%
%    mms.db_init('db_cache_enabled',true)
%    mms.db_init('db_cache_size_max',4096) % set cache to 4GB
%    mms.db_init('db_index_enabled',true) % enable use of database index


global MMS_DB;
if isempty(MMS_DB), MMS_DB = mms_db; end

narginchk(0,2)

if nargin==0
  localFileDbRoot = datastore('mms_db','local_file_db_root');
  if isempty(localFileDbRoot)
    irf.log('warining','local_file_db_root empty - run mms.db_init()')
  else
    localFileDb = mms_local_file_db(localFileDbRoot);
    MMS_DB.add_db(localFileDb);
    irf.log('notice','Initialized local_file_db')
  end
  localFileDbRoot = datastore('mms_db','local_file_db_irfu_root');
  if ~isempty(localFileDbRoot)
    localFileDb = mms_local_file_db(localFileDbRoot);
    MMS_DB.add_db(localFileDb);
    irf.log('notice','Initialized local_file_irfu_db')
  end
  cacheEnabled = datastore('mms_db','db_cache_enabled');
  if ~isempty(cacheEnabled)
    MMS_DB.cache.enabled = cacheEnabled;
    if cacheEnabled, s = 'ON'; else, s = 'OFF'; end
    irf.log('notice',['db_cache: ' s])
  end
  cacheTimeout = datastore('mms_db','db_cache_timeout');
  if ~isempty(cacheTimeout)
    MMS_DB.cache.timeout = cache.timeout;
    irf.log('notice',sprintf('db_cache_timeout: %d sec',cacheTimeout))
  end
  cacheSizeMax = datastore('mms_db','db_cache_size_max');
  if ~isempty(cacheSizeMax)
    MMS_DB.cache.cacheSizeMax = cacheSizeMax;
    irf.log('notice',sprintf('db_cache_size_max: %d MB',cacheSizeMax))
  end
  return
end

if nargin ~=2, error('need a pair KEY,VALUE'), end
  
switch key
  case 'local_file_db'
    localFileDb = mms_local_file_db(val);
    MMS_DB.add_db(localFileDb);
    datastore('mms_db','local_file_db_root',val);
  case 'local_file_db_irfu'
    localFileDb = mms_local_file_db(val);
    MMS_DB.add_db(localFileDb);
    datastore('mms_db','local_file_db_irfu_root',val);
  case 'db_index_enabled'
    MMS_DB.index.enabled = val;
    datastore('mms_db',key,val)
  case 'db_cache_enabled'
    MMS_DB.cache.enabled = val;
    datastore('mms_db',key,val)
  case 'db_cache_timeout'
    MMS_DB.cache.timeout = val;
    datastore('mms_db',key,val)
  case 'db_cache_size_max'
    MMS_DB.cache.cacheSizeMax = val;
    datastore('mms_db',key,val)
  otherwise
    error('Unrecoglized KEY')
end
