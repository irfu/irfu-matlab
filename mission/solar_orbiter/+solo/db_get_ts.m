function res = db_get_ts(filePrefix,varName,tint)
%SOLO.DB_GET_VARIABLE  get variable TSeries from a file DB
%
% res = solo.db_get_ts(filePrefix,varName,tint)
%
% ARGUMENTS
% ---------
% filePrefix
%     Beginning of file name for the datasets that shall be searched, up until
%     (but excluding) the third underscore.
%     Note: Should include "-cdag" when present.
%     Ex: solo_L1_rpw-hfr-surv-cdag_20230101_V06.cdf
%         ==> solo_L1_rpw-hfr-surv-cdag
% varName
%     zVariable name
% tint
%     Time interval for which data shall be retrieved.
%
% RETURN VALUE
% ------------
% (1) TSeries (nominally), or
% (2) Empty array (double), if there is no data, or
% (2) 1D cell array of multiple TSeries, if the zVariables with the same name
%     in the underlying datasets (data files) appear to be *not* equivalent e.g.
%     due to having different zVariable attributes.
%     NOTE: The underlying algorithm for finding the underlying datasets (data
%     files) may include datasets for the day before the begin timestamp, and/or
%     possibly the day after the end timestamp. If the zVariable in such dataset
%     is deemed not equivalent to the other datasets, then the function will
%     still return a cell array of TSeries, despite the change of e.g. zVariable
%     metadata occurring outside the specified time range.

narginchk(3,3)

global SOLO_DB; if isempty(SOLO_DB), solo.db_init(), end

res = SOLO_DB.get_ts(filePrefix,varName,tint);
end
