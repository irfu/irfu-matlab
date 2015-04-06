function filenameData = mms_fields_file_info(fileName)

pos = strfind(fileName,'_');
if length(pos)<5
  irf.log('warning',...
    'mms_sdp_load sci filename has too few parts seperated by underscore.');
end

filenameData = [];
filenameData.scId = fileName(1:pos(1)-1);
filenameData.instrumentId = fileName(pos(1)+1:pos(2)-1);
filenameData.tmMode = fileName(pos(2)+1:pos(3)-1);
filenameData.dataLevel = fileName(pos(3)+1:pos(4)-1);
filenameData.dataType = '';
if length(pos)==6, filenameData.dataType = fileName(pos(4)+1:pos(5)-1); end
filenameData.startTime = fileName(pos(end-1)+1:pos(end)-1);
filenameData.date = filenameData.startTime;
filenameData.vXYZ = fileName(pos(end)+2:end-4);
% Also store the filename, used as Parents reference for output files.
filenameData.filename = fileName;

lST = length(filenameData.startTime);
if lST == 14, % = yyyymmddhhmmss, fine
elseif lST<14 && lST>7
  % If length has been shortened, padd it so we can find day of year and
  % proper startTime.
  filenameData.startTime((lST+1):14) = '0';
else
  % Something went wrong in reading startTime from filename.
  errStr = ['unexpected format of startTime from ', fullFilename];
  irf.log('critical', errStr);
  error('MATLAB:MMS_FIELDS_FILE_INFO:INPUTFILE', errStr);
end