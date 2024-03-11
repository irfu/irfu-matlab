function fileList = db_list_files(filePrefix,tint)
%fileList = db_list_files(filePrefix,tint)
%
% Obtain list of matching SolO datasets (files).
%
% NOTE: This is the same function which underlies retreiving zVariable data via
% solo.db_get_ts().
%
%
% ARGUMENTS
% =========
% filePrefix
%     Beginning of file name for the datasets that shall be searched, up until
%     (but excluding) the third underscore.
%     NOTE: Must include "-cdag" when present. Is otherwise equivalent to
%     "dataset ID".
%     Ex: solo_L1_rpw-hfr-surv-cdag_20230101_V06.cdf
%         ==> solo_L1_rpw-hfr-surv-cdag
% tint
%     2x1 GenericTimeArray. Represents time interval.
%
%
% RETURN VALUE
% ============
% fileList
%     Either
%     (1) 1D array of structs. One struct per file, including path.
%     (2) Empty array (double), if no matching files were found.
%
%
% NOTE: If no dataset is found (e.g. due to adding/omitting "-cdag" to argument
% "filePrefix"), then no exception is raised.

narginchk(1,2)
if nargin<2, tint=[]; end
global SOLO_DB; if isempty(SOLO_DB), solo.db_init(), end

fileList = SOLO_DB.list_files(filePrefix,tint);
end
