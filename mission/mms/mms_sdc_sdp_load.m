function filenameData = mms_sdc_sdp_load( fullFilename, sci_or_ancillary, dataType )
% MMS_SDC_SDP_LOAD reads a MMS CDF files or ancillary data files
% and store them in memory using the MMS_SDC_SDP_DATAMANAGER.
%	MMS_SDC_SDP_LOAD( fullFilename, sci_or_ancillary, datatype) 
%   read the CDF file or (Ancillary ASCII?) found as fullFilename and send
%   the data along to the DATAMANAGER to be stored as datatype. And return
%   some decoded information about the file in filenameData.
%
%	Example:
%		filenameData = mms_sdc_sdp_load(...
%     '/full/path/2015/04/10/mms2_sdp_fast_l1b_20150410_v0.0.0.cdf',...
%     'sci', 'dce');
%
% 	See also MMS_SDC_SDP_PROC, MMS_SDC_SDP_DATAMANAGER.

narginchk(2,3);

switch lower(sci_or_ancillary)
  case 'sci'
    % Convert input file name to parameters, use this to go down through the
    % data folder structure. (if not "HK" then "science" folder, mms1 folder,
    % instrument folder, mode folder, datalevel folder, start time (as year and
    % day of year) folder.
    irf.log('notice',...
      ['MMS_SDC_SDP_LOAD received input filename: ', ...
      fullFilename]);
    if(~exist(fullFilename,'file'))
      irf.log('critical', ...
        ['mms_sdc_sdp_cdf_process CDF file not found: ',fullFilename]);
      error('MATLAB:SDCcode','184');
    end
    [~, filename, ~] = fileparts(fullFilename); 
    pos = strfind(filename,'_');
    if length(pos)<5
        irf.log('warning',...
            'mms_sdc_sdp_cdf_process sci filename has too few parts seperated by underscore.');
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
      error('MATLAB:MMS_SDC_SDP_LOAD:INPUTFILE', errStr);
    end
    
    tmpDataObj = dataobj(fullFilename,'KeepTT2000');
    % Store it using the DataManager as dataType
    mms_sdc_sdp_datamanager(dataType,tmpDataObj);

  case 'ancillary'
    % Ancillary data, named differently and stored differently. At first
    % glance it appears it is stored as simply ASCII files.
    % MMS SDC Developer Guide list these as stored in folders structure:
    % $DATA_PATH_ROOT/ancillary/mmsX/defatt/YYYY/

    % FIXME: Call upon DataManager to store the values in global memory (DATAC).

    irf.log('notice', ['Received input filename: ', fullFilename, ' for ancillary data.']);

    if(~exist(fullFilename,'file'))
      errStr = ['File not found. ', fullFilename];
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDC_SDP_LOAD:INPUTFILE', errStr);
    end
    % DEFATT File start with header, number of lines with header is not
    % constant nor do all header lines beging with "COMMENT", but the
    % last header line does. Also last line in file contain only one
    % column, "DATA_STOP". In order to import the data use textscan with
    % format specified (doesn't match "DATA_STOP") and skip all headers.
    %
    % Get number of last line with COMMENT using unix commands grep, tail and cut.
    if ismac
      [~, numHeaders] = unix(['grep -onr COMMENT ',fullFilename,...
        ' | tail -n1 | cut -d'':'' -f2']);
    else
      [~, numHeaders] = unix(['grep -onr COMMENT ',fullFilename,...
        ' | tail -n1 | cut -d'':'' -f1']);
    end
    numHeaders = str2double(numHeaders);
    
    fileID = fopen(fullFilename, 'r');
    % DEFATT file:
    % Column 1 time in format YYYY-DOYTHH:mm:ss.SSS (where DOY is day of year and SSS is milliseconds)
    % Column 10 Z-Phase (in degrees).
    formatSpec='%f-%f%s %*f %*f %*f %*f %*f %*f %*f %*f %f %*[^\n]';
    tmpData = textscan( fileID, formatSpec,...
      'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', numHeaders );
    fclose(fileID);
    
    % Convert time to format YYYY-MM-DDTHH:mm:ss.mmmuuunnn (where mmm = ms, uu = us, nnn = ns)
    timeFullStr=[irf_time([tmpData{1,1}, tmpData{1,2}],'doy2yyyy-mm-dd'), ...
      cell2mat(tmpData{1,3}), repmat('000000',size(tmpData{1,3},1),1)];
    % Convert time to TT2000 and store as struct (for now)
    DEFATT.time = spdfparsett2000(timeFullStr);
    DEFATT.zphase = tmpData{1,4};
    
    % Return filename (to be stored in CDF GATTRIB Parents)
    filenameData = [];
    [~, filenameData.filename, ~] = fileparts(fullFilename);
    
    % Store it using the DataManager as dataType
    mms_sdc_sdp_datamanager(dataType,DEFATT);
  otherwise
    % Processing cdf files req. either SCIENCE or ANCILLARY data files.
    err_str = 'Input must be either "sci" or "ancillary" information';
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDC_SDP_LOADS', err_str);
end
end
