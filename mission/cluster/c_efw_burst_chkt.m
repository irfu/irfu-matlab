function startSAT=c_efw_burst_chkt(database,filename)
%C_EFW_BURST_CHKT Get the correct start time for the efw internal burst
%
% STARTTIME=C_EFW_BURST_CHKT(DATABASE,FILENAME) reads the EFW burst
% header from the file FILENAME, and tries to transform the
% starting time from the EFW clock to the satellite clock, using
% linear interpolation.
%
% STARTTIME is the start time for the burst in epoch according to
% the spacecraft clock.
%
% There are several problems that can occur here, so don't be
% scared if the calculated time seems strange.
%

% By Anders Tjulin, last update 17/5-2004

% Extract the times from the burst data file

times=c_efw_burst_geth(filename);
playbackEFW=times(1);
playbackSAT=times(2);
startEFW=times(3);
startSAT=[];

% Choose the earliest possible starting time, preferably a short
% time before what ISDAT tells you

shorttime=150; % Seconds before ISDATs startingtime

% Choose time interval to be longer than the short time used before
duration=5*shorttime;

% Get EFW time as function of satellite time
spacecraft = str2double(filename(end));
while shorttime >= 0
  startdate=(playbackSAT-(playbackEFW-startEFW) - shorttime );
  
  [sctime,temp]=caa_is_get(database,startdate, ...
    duration,spacecraft,'efw','DSC');
  if ~isempty(temp), break, end
  shorttime = shorttime - 5;
end
if isempty(temp) || size(sctime,1)<3
  irf_log('proc','Cannot get EFW time for IB');
  return;
end
efwtime=(temp(81,:)+temp(82,:)*256+temp(83,:)*65536+ temp(84,:)* ...
  16777216+temp(85,:)*4294967296)/1000;

% Find the time-index before which the burst was collected

tempindex=find((efwtime-startEFW)>0);

if isempty(tempindex)
  irf_log('proc','Cannot get EFW time for IB. Will extrapolate.')
  tempindex = length(efwtime);
else
  tempindex = tempindex(1);
end
if tempindex==1
  idx=2;
  % adjust if time gap in efwtime or sctime
  if (efwtime(idx) - efwtime(idx-1))>60 || (sctime(idx) - sctime(idx-1))>60
    tempindex=tempindex+1;
    idx=idx+1;
  end
else
  % adjust if time gap in efwtime or sctime
  if (efwtime(tempindex) - efwtime(tempindex-1))>60 || (sctime(tempindex) - sctime(tempindex-1))>60
    if tempindex==length(efwtime)
      tempindex=tempindex-1;
    else
      tempindex=tempindex+1;
    end
  end
  idx = tempindex;
end

% Interpolate to get the exact start time
startSAT = sctime(tempindex) + (startEFW-efwtime(tempindex))* ...
  ( sctime(idx) - sctime(idx-1) )/( efwtime(idx) - efwtime(idx-1) );