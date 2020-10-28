function res = db_get_variable(filePrefix,varName,tint)
%SOLO.DB_GET_VARIABLE  get variable from a file DB
%
% res = solo.db_get_variable(filePrefix,varName,tint)

narginchk(3,3)

global SOLO_DB; if isempty(SOLO_DB), solo.db_init(), end

res = SOLO_DB.get_variable(filePrefix,varName,tint);