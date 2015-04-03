function filenameData = mms_sdp_load( fullFilename, dataType )
% MMS_SDP_LOAD reads a SDP data files and store them in MMS_SDP_DATAMANAGER
%
%	MMS_SDP_LOAD( fullFilename, , datatype)
%   read the CDF file fullFilename and send the data along to DATAMANAGER 
%   to be stored as dataType. Returns information about the file.
%
%	Example:
%		filenameData = mms_sdp_load(...
%     '/full/path/2015/04/10/mms2_sdp_fast_l1b_20150410_v0.0.0.cdf', 'dce');
%
% 	See also MMS_SDC_SDP_PROC, MMS_SDP_DATAMANAGER, MMS_LOAD_ANCILLARY

narginchk(2,2);

% Convert input file name to parameters, use this to go down through the
% data folder structure. (if not "HK" then "science" folder, mms1 folder,
% instrument folder, mode folder, datalevel folder, start time (as year and
% day of year) folder.
irf.log('notice',...
  ['MMS_SDP_LOAD received input filename: ', ...
  fullFilename]);
if(~exist(fullFilename,'file'))
  irf.log('critical', ...
    ['mms_sdp_load CDF file not found: ',fullFilename]);
  error('MATLAB:SDCcode','184');
end
[~, filename, ~] = fileparts(fullFilename);
pos = strfind(filename,'_');
if length(pos)<5
  irf.log('warning',...
    'mms_sdp_load sci filename has too few parts seperated by underscore.');
end

filenameData = [];
filenameData.scId = filename(1:pos(1)-1);
filenameData.instrumentId = filename(pos(1)+1:pos(2)-1);
filenameData.tmMode = filename(pos(2)+1:pos(3)-1);
filenameData.dataLevel = filename(pos(3)+1:pos(4)-1);
% filenameData.optional, ignored.
filenameData.startTime = filename(pos(end-1)+1:pos(end)-1);
filenameData.vXYZ = filename(pos(end)+2:end-4);
% Also store the filename, used as Parents reference for output files.
filenameData.filename = filename;

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
  error('MATLAB:MMS_SDP_LOAD:INPUTFILE', errStr);
end

tmpDataObj = dataobj(fullFilename,'KeepTT2000');
% Store it using the DataManager as dataType
mms_sdp_datamanager(dataType,tmpDataObj);

end
