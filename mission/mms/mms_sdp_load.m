function filenameData = mms_sdp_load( fullFilename, sci_or_ancillary, dataType )
% MMS_SDP_LOAD reads a MMS CDF files or ancillary data files
% and store them in memory using the MMS_SDP_DATAMANAGER.
%	MMS_SDP_LOAD( fullFilename, sci_or_ancillary, datatype) 
%   read the CDF file or (Ancillary ASCII?) found as fullFilename and send
%   the data along to the DATAMANAGER to be stored as dataType. And return
%   some decoded information about the file in filenameData.
%
%	Example:
%		filenameData = mms_sdp_load(...
%     '/full/path/2015/04/10/mms2_sdp_fast_l1b_20150410_v0.0.0.cdf',...
%     'sci', 'dce');
%		filenameData = mms_sdp_load(...
%     '/full/path/2015/MMS2_DEFEPH_2015074_2015075.V00',...
%     'ancillary', 'defeph');
%
% 	See also MMS_SDC_SDP_PROC, MMS_SDP_DATAMANAGER.

narginchk(3,3);

switch lower(sci_or_ancillary)
  case 'sci'
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

  case 'ancillary'
    % Ancillary data, named differently and stored differently. At first
    % glance it appears it is stored as simply ASCII files.
    % MMS SDC Developer Guide list these as stored in folders structure:
    % $DATA_PATH_ROOT/ancillary/mmsX/defatt/YYYY/

    irf.log('notice', ['Received input filename: ', fullFilename, ' for ancillary data.']);

    if(~exist(fullFilename,'file'))
      errStr = ['File not found. ', fullFilename];
      irf.log('critical', errStr);
      error('MATLAB:MMS_SDP_LOAD:INPUTFILE', errStr);
    end

    if(strcmp(dataType,'defatt'))
      % DEFATT file:
      % The DEFATT files start with a header, the number of lines with
      % header is not constant nor do all the header lines begin with
      % "COMMENT", but the last header line does. Also the last line in
      % DEFATT contain only one column with string "DATA_STOP". This last
      % line will not match the specified format and will not be an issue,
      % if the number of headers are specified to textscan.
      % Last header line is identified by the existence of "COMMENT".
      headerGrep = 'COMMENT';

      % Column 1 time in format yyyy-doyTHH:MM:SS.mmm
      % (where doy is day of year and mmm is milliseconds)
      % Column 10 Z-Phase (in degrees).
      formatSpec='%f-%f%s %*f %*f %*f %*f %*f %*f %*f %*f %f %*[^\n]';
    elseif(strcmp(dataType,'defeph'))
      % DEFEPH file:
      % The DEFEPH files start with a header, the number of lines with
      % header in not constant and this file does not contain things like
      % "COMMENT" for header lines. However the last line of header will
      % contain the units of each column, for DEFEPH this means we can
      % identify the last header line by looking for a string which only
      % should appear as a unit.
      % Last header is identified by the existence of "Km/Sec"
      headerGrep = 'Km/Sec';

      % Column 1 time in format yyyy-doy/HH:MM:SS.mmm
      % (where doy is day of year and mmm is milliseconds),
      % Column 3, 4, 5 is position in X,Y,Z (in some ref.frame, TBC which)
      formatSpec='%f-%f%s %f %f %f %*[^\n]';
    else
      errStr = ['Unknown dataType: ',dataType,' for ancillary data. ', ...
        'Valid values are "defatt" and "defeph".'];
      irf.log('critical', errStr); error(errStr);
    end

    % Get number of last header line using unix commands grep, tail and cut.
    if ismac
      [~, numHeaders] = unix(['grep -onr ',headerGrep,' ',fullFilename,...
        ' | tail -n1 | cut -d'':'' -f2']);
    else
      [~, numHeaders] = unix(['grep -onr ',headerGrep,' ',fullFilename,...
        ' | tail -n1 | cut -d: -f1']);
    end
    numHeaders = str2double(numHeaders);
    
    fileID = fopen(fullFilename, 'r');
    tmpData = textscan( fileID, formatSpec,...
      'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', numHeaders );
    fclose(fileID);
    
    % Return filename (to be stored in CDF GATTRIB Parents)
    filenameData = [];
    [~, filename, fileext] = fileparts(fullFilename);
    filenameData.filename = [filename, fileext]; % DEFATT/DEFEPH has extension
    %'V00', 'V01' etc being the version number, store this as well.

    if(strcmp(dataType,'defatt'))
      % Convert time to format TT2000 using the irf_time function by first
      % converting [yyyy, doy] to a 'yyyy-mm-dd' string, then adding the
      % remaining 'THH:MM:SS.mmm' which was read into a string previously.
      dataIN.time = irf_time([irf_time([tmpData{1,1}, tmpData{1,2}],'doy>utc_yyyy-mm-dd'), ...
        cell2mat(tmpData{1,3})],'utc>ttns');
      % And simply include the Z-phase.
      dataIN.zphase = tmpData{1,4};

    elseif(strcmp(dataType,'defeph'))
      % Convert time to fromat TT2000 using the irf_time function by first
      % converting [yyyy, doy] to a 'yyyy-mm-dd' string, then add the
      % remaining 'HH:MM:SS.mmm' string (excluding the "/") which was read
      % into a string previously.
      tmpStr = cell2mat(tmpData{1,3});
      dataIN.time = irf_time([irf_time([tmpData{1,1}, tmpData{1,2}],'doy>utc_yyyy-mm-dd'), ...
        tmpStr(:,2:end)],'utc>ttns');
      % And simply include the position in X, Y, Z (which ref.syst. used is
      % still TBC)
      dataIN.Pos_X = tmpData{1,4};
      dataIN.Pos_Y = tmpData{1,5};
      dataIN.Pos_Z = tmpData{1,6};
    end

    % Store it using the DataManager as dataType
    mms_sdp_datamanager(dataType,dataIN);

  otherwise
    % Processing cdf files req. either SCIENCE or ANCILLARY data files.
    err_str = 'Input must be either "sci" or "ancillary" information';
    irf.log('critical', err_str);
    error('MATLAB:MMS_SDP_LOADS', err_str);
end
end
