function fileList = db_list_files(filePrefix,tint)

narginchk(1,2)
if nargin<2, tint=[]; end
global MMS_DB; if isempty(MMS_DB), mms.db_init(), end

fileList = MMS_DB.list_files(filePrefix,tint);

