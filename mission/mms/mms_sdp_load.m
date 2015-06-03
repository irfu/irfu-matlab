function filenameData = mms_sdp_load( fullFilename, dataType )
% MMS_SDP_LOAD reads a SDP data files and store them in MMS_SDP_DATAMANAGER
%
%	MMS_SDP_LOAD( fullFilename, datatype)
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

[~, fileName, ~] = fileparts(fullFilename);
filenameData = mms_fields_file_info(fileName);

tmpDataObj = dataobj(fullFilename);
% Store it using the DataManager as dataType
mms_sdp_datamanager(dataType,tmpDataObj);

end
