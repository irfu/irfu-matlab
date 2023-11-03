function [time_of_events, dist_t] = search_events(start_time, dt, dist_to_MP, p_solarwind)
%
%Input:
% start_time -when to start the calculation [yyyy mm dd hh mm ss]
% dt -the duration of the calculation in seconds.
% dist_to_MP -the max distance to the MP for beeing called an event, in Re.
% p_solarwind -the preasure of the solarwind in nPa.
%
%Output:
% time_of_events -the timestamps when the satellites where closer
% than "dist_to_MP" to the MP and when the satellites leaves the
% close area. (The times are in epochs)
% dist_t - a table with distances to MP for all the samples when
%        when the satellites are within distance to MP.
%
%Descrition of the function:
% This functions finds the event when the satelliets are within a
% certain distance from the magnetopause, MP, and then time when
% it leaves the MP.
%
%Using:
% get pos_from_ISDAT
% distance_to_MP()
% gse2gsm()
%
%Work method:
%% This functions downloads, all the positions for ISDAT, within the time
% frame for "start_time" and "dt" seconds forward. The distance to the
% MP is calculated and compared with the given limit "dist_to_MP". When the
% position is in GSE so it has to be changed to GSM to fit the MP-modell
% If "dist" falls within "dist_to_MP", then the timestamp of this event
% is stored in "time_of_events". If no events are found, -1 is returned.
% The function gse2gsm requires that the time is in a special time format
% This makes that only one day at the time can be loaded.
%
%Error:
% If no events are found, -1 is returned
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

inside = 0;
dist_pos = 0;
dist_t = -1;

%an error check should be here to avoid reading outside the
%datafield when calling get position.

% get the position of the refernce satellite
pos = get_pos_from_ISDAT(3, start_time, dt);
nr_of_positions = length(pos);

% This section converting the date to fit the function gse2gsm
% Convertes yyyy -> yy, keeping the last to digits
y = num2str(start_time(1));
y = strread(y(3:4));
% Convertes month and day -> number of day in the year
d = datenum(start_time(1:3)) - datenum(start_time(1),1 ,1) +1;

if pos(1,1) == -1 %errorcheck. if -1 the there is no data
  time_of_events = -1

else
  b=0; %starting with zero events
  for i = 1:nr_of_positions
    %the seconds of the day
    s = pos(i,1) - toepoch([start_time(1:3) 0 0 0]);
    xgsm = pos(i,2);
    [ygsm, zgsm]= irf_gse2gsm(y,d,s,pos(i,3), pos(i,4),0);
    xgsm_re = xgsm/6378;
    ygsm_re = ygsm/6378;
    zgsm_re = zgsm/6378;

    dist = distance_to_MP(p_solarwind, -1, xgsm_re, ygsm_re, zgsm_re);

    %logging the distance
    if dist < dist_to_MP
      dist_pos = dist_pos + 1;
      dist_t(dist_pos,1) = pos(i,1);
      dist_t(dist_pos,2) = dist;
    end

    % gets the timetag when the satellite enters the MP
    if dist < dist_to_MP && inside == 0
      b=b+1;
      inside = 1;
      time_of_events(b,1) = pos(i,1);
    end

    % gets the timetag when the satellite leaves the MP
    if dist > dist_to_MP && inside == 1
      inside = 0;
      time_of_events(b,2) = pos(i,1);
    end

  end

  %when the duration time, dt, has ended and the satellites are still within
  %the event area.
  if inside == 1
    time_of_events(b,2) = pos(i,1);
  end

  %of no events where found then -1 is returned
  if b == 0
    time_of_events = -1;
  end

end




