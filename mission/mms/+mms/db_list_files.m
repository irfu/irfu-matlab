function fileList = db_list_files(filePrefix,tint)

narginchk(1,2)

global MMS_DB;
if isempty(MMS_DB)
  errS = 'MMS_DB not initialized, run mms.db_init()';
  irf.log('critical', errS), error(errS)
end

fileList = MMS_DB.list_files(filePrefix,tint);

