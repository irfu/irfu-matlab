function filenameData = mms_sdc_sdp_cdf_in_process( fullFilename, sci_or_ancillary, dataType )
% MMS_SDC_SDP_CDF_IN_PROCESS reads a MMS CDF files or ancillary data files
% and store them in memory using the MMS_SDC_SDP_DATAMANAGER.
%	MMS_SDC_SDP_CDF_IN_PROCESS( fullFilename, sci_or_ancillary, datatype) 
%   read the CDF file or (Ancillary ASCII?) found as fullFilename and send
%   the data along to the DATAMANAGER to be stored as datatype. And return
%   some decoded information about the file in filenameData.
%
%	Example:
%		filenameData = mms_sdc_sdp_cdf_in_process(...
%     '/full/path/2015/04/10/mms2_sdp_fast_l1b_20150410_v0.0.0.cdf',...
%     'sci', 'dce');
%
% 	See also MMS_SDC_SDP_PROC, MMS_SDC_SDP_DATAMANAGER.

narginchk(2,3);

if(strcmp(sci_or_ancillary,'sci'))
    % Convert input file name to parameters, use this to go down through the
    % data folder structure. (if not "HK" then "science" folder, mms1 folder,
    % instrument folder, mode folder, datalevel folder, start time (as year and
    % day of year) folder.
    % 
    
    irf.log('debug',...
        ['mms_sdc_sdp_cdf_in_process received input filename: ', ...
        fullFilename]);
    
    filenameData = [];
    
    [~, filename, ~] = fileparts(fullFilename); 
    
    pos = strfind(filename,'_');
    if length(pos)<5
        irf.log('warning',...
            ['mms_sdc_sdp_cdf_process sci filename has too few parts seperated by underscore.']);
    end

    
    filenameData.scId = filename(1:pos(1)-1);
    filenameData.instrumentId = filename(pos(1)+1:pos(2)-1);
    filenameData.dataMode = filename(pos(2)+1:pos(3)-1);
    filenameData.dataLevel = filename(pos(3)+1:pos(4)-1);
    % filenameData.optional, ignored.
    filenameData.startTime = filename(pos(end-1)+1:pos(end)-1);
    filenameData.vXYZ = filename(pos(end)+2:end-4);
    % Also store the filename, used as Parents reference for output files.
    filenameData.filename = filename;

    switch(length(filenameData.startTime))
        % If length has been shortened, padd it so we can find day of year and
        % proper startTime.
        case 8 % = yyyymmdd
            filenameData.startTime = [filenameData.startTime,'000000'];
        case 9 % = yyyymmddh ?
            filenameData.startTime = [filenameData.startTime,'00000'];
        case 10 % = yyyymmddhh
            filenameData.startTime = [filenameData.startTime,'0000'];
        case 11 % = yyyymmddhhm
            filenameData.startTime = [filenameData.startTime,'000'];
        case 12 % = yyyymmddhhmm
            filenameData.startTime = [filenameData.startTime,'00'];
        case 13 % = yyyymmddhhmms
            filenameData.startTime = [filenameData.startTime,'0'];
        case 14 % = yyyymmddhhmmss, fine
        otherwise
            % Something went wrong in reading startTime from filename.
            irf.log('warning',...
                'mms_sdc_sdp_cdf_process something went wrong in reading startTime from filename input.');
    end
    
    if(~exist(fullFilename,'file'))
        irf.log('critical', ...
            ['mms_sdc_sdp_cdf_process CDF file not found: ',fullFilename]);
        error('MATLAB:SDCcode','184');
    else
        % FIXME: tint is ignored with 4 arguments, but KeepTT2000 for MMS.
        tmpDataObj = dataobj(fullFilename,'tint',0,true);
        % Store it using the DataManager as dataType
        mms_sdc_sdp_datamanager(tmpDataObj, dataType);
    end

elseif(strcmp(sci_or_ancillary,'ancillary'))
    % Ancillary data, named differently and stored differently. At first
    % glance it appears it is stored as simply ASCII files in folders by
    % DOY.
    % FIXME: Implement reading of these files and store the data in
    % DATAMANAGER.

else
    % Processing cdf files req. either SCIENCE or ANCILLARY data files.
    err_str = 'Input must be either "sci" or "ancillary" information';
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_CDF_IN_PROCESSS', err_str);

end

