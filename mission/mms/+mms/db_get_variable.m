function res = db_get_variable(filePrefix,varName,tint)
%MMS.DB_GET_VARIABLE  get variable from a file DB
%
% res = mms.db_get_variable(filePrefix,varName,tint)

narginchk(3,3)

global MMS_DB; if isempty(MMS_DB), mms.db_init(), end

res = MMS_DB.get_variable(filePrefix,varName,tint);