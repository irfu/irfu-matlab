function startSAT=checkbursttime(database,filename)

%CHECKBURSTTIME Get the correct start time for the efw internal burst
%
% STARTTIME=CHECKBURSTLIST(DATABASE,FILENAME) reads the EFW burst
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
%
% By Anders Tjulin, last update 17/5-2004

  % Extract the times from the burst data file

  times=getburstheader(filename);
  playbackEFW=times(1);
  playbackSAT=times(2);
  startEFW=times(3);
  
  % Choose the earliest possible starting time, preferably a short
  % time before what ISDAT tells you
  
  shorttime=150; % Seconds before ISDATs startingtime
  startdate=(playbackSAT-(playbackEFW-startEFW)-150);

  % Choose time interval to be longer than the short time used before
  
  duration=4*shorttime;

  % Get EFW time as function of satellite time
  
  spacecraft=filename(end);
  [sctime,temp]=isGetDataLite(database,fromepoch(startdate), ...
			    duration,'Cluster',spacecraft,'efw','DSC');
  efwtime=(temp(81,:)+temp(82,:)*256+temp(83,:)*65536+ temp(84,:)* ...
	 16777216+temp(85,:)*4294967296)/1000;

  % Find the time-index before which the burst was collected
  
  tempindex=find((efwtime-startEFW)>0);

  % Interpolate to get the exact start time
  
  startSAT=sctime(tempindex(1))+(startEFW-efwtime(tempindex(1)))* ...
	 (sctime(tempindex(1))-sctime(tempindex(1)-1))/ ...
	 (efwtime(tempindex(1))-efwtime(tempindex(1)-1));
