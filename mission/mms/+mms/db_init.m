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
  otherwise
    error('Unrecoglized KEY')
end
