function res = db_get_ts(filePrefix,varName,tint)
%MMS.DB_GET_VARIABLE  get variable TSeries from a file DB
%
% res = mms.db_get_ts(filePrefix,varName,tint)

narginchk(3,3)

global MMS_DB; if isempty(MMS_DB), mms.db_init(), end

res = MMS_DB.get_ts(filePrefix,varName,tint);

if isempty(res), res = TSeries([]); end
