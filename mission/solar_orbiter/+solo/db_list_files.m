function fileList = db_list_files(filePrefix,tint)

narginchk(1,2)
if nargin<2, tint=[]; end
global SOLO_DB; if isempty(SOLO_DB), solo.db_init(), end

fileList = SOLO_DB.list_files(filePrefix,tint);
end
