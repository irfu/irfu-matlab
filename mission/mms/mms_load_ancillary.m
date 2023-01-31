function [dataIN,filenameData] = mms_load_ancillary(fullFilename,dataType)
%MMS_LOAD_ANCILLARY  load DEFATT/DEFEPH file
%
% [dataIN,filenameData] = mms_load_ancillary(fullFilename,dataType)
%
% dataType is one of : 'defatt','defeph','defq', 'predq'

if(~exist(fullFilename,'file'))
  errStr = ['File not found. ', fullFilename];
  irf.log('critical', errStr);
  error('MATLAB:MMS_LOAD_ANCILLARY:INPUTFILE', errStr);
end

% Return filename (to be stored in CDF GATTRIB Parents)
filenameData = [];
[~, filename, fileext] = fileparts(fullFilename);
filenameData.filename = [filename, fileext]; % DEFATT/DEFEPH/DEFQ has
% extension 'V00', 'V01' etc being the version number, store this as well.

switch lower(dataType)
  case 'defatt'
    % DEFATT file:
    % The DEFATT files start with a header, the number of lines with
    % header is not constant nor do all the header lines begin with
    % "COMMENT", but the last header line does. Also the last line in
    % DEFATT contain only one column with string "DATA_STOP". This last
    % line will not match the specified format and will not be an issue,
    % if the number of headers are specified to textscan.
    % Last header line is identified by the existence of "COMMENT".
    headerGrep = 'COMMENT';
    nHeaders = 49;
    % Column 1 time in format yyyy-doyTHH:MM:SS.mmm
    % (where doy is day of year and mmm is milliseconds)
    % Column 10 Z-Phase (in degrees).
    formatSpec='%f-%f%s %*f %*f %*f %*f %*f %*f %*f %*f %f %f %f %f %f %f %f %f %f %f %*[^\n]';
  case 'defeph'
    % DEFEPH file:
    % The DEFEPH files start with a header, the number of lines with
    % header in not constant and this file does not contain things like
    % "COMMENT" for header lines. However the last line of header will
    % contain the units of each column, for DEFEPH this means we can
    % identify the last header line by looking for a string which only
    % should appear as a unit.
    % Last header is identified by the existence of "Km/Sec"
    headerGrep = 'Kg';
    nHeaders = 14;
    % Column 1 time in format yyyy-doy/HH:MM:SS.mmm
    % (where doy is day of year and mmm is milliseconds),
    % Column 3, 4, 5 is position in X,Y,Z (in some ref.frame, TBC which)
    formatSpec='%f-%f%s %f %f %f %f %f %f %f %*[^\n]';
  case {'defq', 'predq'}
    % DEFQ file:
    % The DEFQ files start with a header, the number of lines with header
    % is not constant and this files does not contain things like "COMMENT"
    % for header lines. However the last line contains "". We can identify
    % number of header lines by looking for this string.
    headerGrep = 'Scale';
    nHeaders = 11;
    % Column 1 time in format yyyy-doy/HH:MM:SS.mmm
    % (where doy is day of year and mmm is miliseconds.
    % Column 3 is Quality factor, column 4 is scale.
    formatSpec='%f-%f%s %*f %f %f %*[^\n]';
    
  otherwise
    errStr = ['Unknown dataType: ',dataType,' for ancillary data. ', ...
      'Valid values are "defatt" and "defeph".'];
    irf.log('critical', errStr); error(errStr);
end

% Get number of last header line using unix commands grep, tail and cut.
if ispc
else
  if ismac, cutArgs = '-d'':'' -f2'; else, cutArgs = '-d: -f1'; end
  nTries = 1;
  [status, nHeaders] = unix(['grep -onr "' headerGrep '" "' fullFilename,...
    '" | tail -n1 | cut ' cutArgs]);
  nHeaders = str2double(nHeaders);
  while (status || ~isa(nHeaders, 'double') || nHeaders<=0) && nTries<=5
    % For some reason it failed, log it and try again.
    irf.log('warning', ['Failed with "grep" on ', fullFilename,'. Trying again.']);
    nTries = nTries+1;
    [status, nHeaders] = unix(['grep -onr "' headerGrep '" "' fullFilename,...
      '" | tail -n1 | cut ' cutArgs]);
    irf.log('warning', ['This time around "grep" got nHeaders:', nHeaders,'.']);
    nHeaders = str2double(nHeaders);
  end
end

fileID = fopen(fullFilename, 'r');
tmpData = textscan( fileID, formatSpec,...
  'delimiter', ' ', 'MultipleDelimsAsOne', 1, 'HeaderLines', nHeaders );
fclose(fileID);

switch lower(dataType)
  case 'defatt'
    % Convert time to format TT2000 using the irf_time function by first
    % converting [yyyy, doy] to a 'yyyy-mm-dd' string, then adding the
    % remaining 'THH:MM:SS.mmm' which was read into a string previously.
    dataIN.time = irf_time([irf_time([tmpData{1,1}, tmpData{1,2}],'doy>utc_yyyy-mm-dd'), ...
      cell2mat(tmpData{1,3})],'utc>ttns');
    % And simply include the Z-phase.
    dataIN.wphase = tmpData{1,4};
    dataIN.zra =    tmpData{1,5};
    dataIN.zdec =   tmpData{1,6};
    dataIN.zphase = tmpData{1,7};
    dataIN.lra =    tmpData{1,8};
    dataIN.ldec =   tmpData{1,9};
    dataIN.lphase = tmpData{1,10};
    dataIN.pra =    tmpData{1,11};
    dataIN.pdec =   tmpData{1,12};
    dataIN.pphase = tmpData{1,13};
    % As per e-mail discussion of 2015/04/07, duplicated timestamps can
    % occur in Defatt (per design). If any are found, use the last data
    % point and disregard the first duplicate.
    idxBad = diff(dataIN.time)==0; % Identify first duplicate
    fs = fields(dataIN);
    for idxFs=1:length(fs), dataIN.(fs{idxFs})(idxBad) = []; end
    
  case 'defeph'
    % Convert time to fromat TT2000 using the irf_time function by first
    % converting [yyyy, doy] to a 'yyyy-mm-dd' string, then add the
    % remaining 'HH:MM:SS.mmm' string (excluding the "/") which was read
    % into a string previously.
    tmpStr = cell2mat(tmpData{1,3});
    dataIN.time = irf_time([irf_time([tmpData{1,1}, tmpData{1,2}],'doy>utc_yyyy-mm-dd'), ...
      tmpStr(:,2:end)],'utc>ttns');
    % And simply include the position in X, Y, Z (which ref.syst. used is
    % still TBC)
    dataIN.r = [tmpData{1,5} tmpData{1,6} tmpData{1,7}];
    dataIN.v = [tmpData{1,8} tmpData{1,9} tmpData{1,10}];
  case {'defq', 'predq'}
    % Convert time to TT2000 using irf_time
    tmpStr = cell2mat(tmpData{1,3});
    dataIN.time = irf_time([irf_time([tmpData{1,1}, tmpData{1,2}],'doy>utc_yyyy-mm-dd'),...
      tmpStr(:,2:end)], 'utc>ttns');
    dataIN.quality = tmpData{1,4};
    dataIN.scale   = tmpData{1,5};
    
  otherwise
end
