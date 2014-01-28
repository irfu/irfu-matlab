function [outObj, filenameData] = mms_cdf_in_process( filename, sci_or_ancillary )
% Process input arguments to find data file (in CDF format) and read this
% file as dataobj, return the dataobj containing all data and info from the
% CDF file.
% Latest modification: 2014/01/16

global ENVIR;

if(strcmp(sci_or_ancillary,'sci'))
    % Convert input file name to parameters, use this to go down through the
    % data folder structure. (if not "HK" then "science" folder, mms1 folder,
    % instrument folder, mode folder, datalevel folder, start time (as year and
    % day of year) folder.
    % 
    
    irf.log('debug',['mms_cdf_in_process recived input filename: ',filename]);
    
    outObj = [];
    filenameData = [];
    
    pos = strfind(filename,'_');
    if length(pos)<5
        irf.log('warning',['mms_cdf_process sci filename has too few parts seperated by underscore.']);
    end

    
    filenameData.scId = filename(1:pos(1)-1);
    filenameData.instrumentId = filename(pos(1)+1:pos(2)-1);
    filenameData.dataMode = filename(pos(2)+1:pos(3)-1);
    filenameData.dataLevel = filename(pos(3)+1:pos(4)-1);
    % filenameData.optional, ignored.
    filenameData.startTime = filename(pos(end-1)+1:pos(end)-1);
    filenameData.vXYZ = filename(pos(end)+2:end-4);

    switch(length(filenameData.startTime))
        % If length has been shortened, padd it so we can find day of year and
        % proper startTime.
        case 8 % = yyyymmdd
            filenameData.startTime = strcat(filenameData.startTime,'000000');
        case 9 % = yyyymmddh ?
            filenameData.startTime = strcat(filenameData.startTime,'00000');
        case 10 % = yyyymmddhh
            filenameData.startTime = strcat(filenameData.startTime,'0000');
        case 11 % = yyyymmddhhm
            filenameData.startTime = strcat(filenameData.startTime,'000');
        case 12 % = yyyymmddhhmm
            filenameData.startTime = strcat(filenameData.startTime,'00');
        case 13 % = yyyymmddhhmms
            filenameData.startTime = strcat(filenameData.startTime,'0');
        case 14 % = yyyymmddhhmmss, fine
        otherwise
            % Something went wrong in reading startTime from filename.
            irf.log('warning','mms_cdf_process something went wrong in reading startTime from filename input.');
    end

    dirDOY = strcat(filenameData.startTime(1:4),'/',filenameData.startTime(5:6),'/',filenameData.startTime(7:8),'/');

    dirToInput = strcat(ENVIR.DATA_PATH_ROOT,'/','science','/',filenameData.scId,'/',filenameData.instrumentId,'/',filenameData.dataMode,'/',filenameData.dataLevel,'/',dirDOY);
    
    % Add debug message to know where we are trying to look for file.
    irf.log('debug',['mms_cdf_process CDF file decoded as located here: ',strcat(dirToInput,filename)]);

    if(~exist(strcat(dirToInput,filename),'file'))
        irf.log('critical',['mms_cdf_process CDF file not found: ',strcat(dirToInput,filename)]);
        error('MATLAB:MMS:mms_cdf_in_process',['inputfile not found:', strcat(dirToInput,filename)]);
    else
        % FIXME: tint is ignored with 4 arguments, but KeepTT2000 for MMS.
        outObj = dataobj(strcat(dirToInput,filename),'tint',0,true);
    end

elseif(strcmp(sci_or_ancillary,'ancillary'))
    % HouseKeeping and Ancillary data, named differently and stored
    % differently.
    

else
    % Processing cdf files req. either SCIENCE or ANCILLARY data files.
    irf.log('warning','Input to mms_cdf_process must contain "sci" or "ancillary" information.');

end

