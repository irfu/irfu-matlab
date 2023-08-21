function [passing_MP,dist_t] = c_ri_auto_event_search(start_time_m,end_time_m,dist2MP, p_solarwind, output_path)
%
% passing_MP = c_ri_auto_event_search(start_time,end_time,dist2MP, p_solarwind, output_path)
%
%Input:
% start_time, end_time - the time range in which the satellites position
%                      - compared to the MP is intressting. [YYYY MM DD HH MM SS.SSS]
% p_solarwind -the preassure from the solarwind
% dist2MP -the distance to the magnetopause
% output_path -where to save the outputfile. Intervall for MP-crossing and distance to MP during
%              the crossings
% input can be:(start_time,end_time) OR (start_time,end_time,dist2MP, p_solarwind)
%               OR (start_time,end_time,dist2MP, p_solarwind, output_path)
%
%Output:
% passing_MP -the time when the satellite enters the MP and when it leaves the MP
%            -the time is in epoch
% save to file: passing_MP and dist_t. (the distance to the MP during a crossing)
%
%Descrition of the function:
% Finds the times when the satellites enters and leaves the magnetopause. The magnetopause
% is defined by a modell using the solarwind preassure as the input.
%
%Using:
% Mat_DbOpen
% Mat_DbClose
% isGetContentLite
% toepoch
% search_events
% add_A2M
% save /share/robert/MP/MP_(from date)_to_(to_date) !!!saving to file
%
%Work method:
% The function loads a matrix containing all the data entries for the satellite positions.
% Then the function finds the best matching start and end times in this matrix. This is
% done to avoid trying to read outside the matrix. Then there is a position calculation
% for all the valid times.
%
%Error:
% [-1,-1] is retured of there are no events found.
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
if nargin == 2
  dist2MP = 3;
  p_solarwind = 2;
  output_path = [pwd '/'];
end

if nargin == 4
  output_path = [pwd '/'];
end

[nr_intervalls,c] =size(start_time_m);
[nr_of_end_intervalls,c] = size(end_time_m);

if nr_intervalls ~= nr_of_end_intervalls
  disp(['error in input. Nr of start points: ' int2str(nr_intervalls) ' nr of endpoints: ' int2str(nr_of_end_intervalls)]);
end

for start_time_intervalls=1:nr_intervalls
  start_time = start_time_m(start_time_intervalls, :);
  end_time = end_time_m(start_time_intervalls, :);


  db = Mat_DbOpen('disco:10');

  [pos_time, dur_time] = isGetContentLite(db,'Cluster','3','ephemeris','position',' ', ' ', ' ');

  Mat_DbClose(db);

  %p_and_f_n = sprintf('%spos_time',output_path)
  %save(p_and_f_n ,'pos_time', 'dur_time');

  passing_MP = 0;
  dist_t = 0;
  pos =0;
  start_time_e = toepoch(start_time);
  end_time_e = toepoch(end_time);
  [k_max,col] = size(pos_time);
  i_start = -1;
  i_end = -1;

  % if any of the values are out of range
  if toepoch(pos_time(1,:)) > start_time_e
    i_start = 1;
  end
  if toepoch(pos_time(k_max,:)) < end_time_e
    i_end = k_max;
  end

  %if start and endtime are in timerange of postime
  if i_start == -1 || i_end == -1

    for k = 1:k_max-1

      t_temp = toepoch(pos_time(k,:));
      n_t_temp = toepoch(pos_time(k+1,:));
      if ((t_temp <= start_time_e) && (t_temp+dur_time(k)) > start_time_e) || ( t_temp > start_time_e && i_start == -1)
        i_start = k;
      end

      t_temp = toepoch(pos_time(k,:));
      n_t_temp = toepoch(pos_time(k+1,:));
      if (t_temp+dur_time(k) > end_time_e && i_end == -1) || ( n_t_temp > end_time_e && i_end == -1)
        i_end = k;
      end

    end

  end

  % searches for events for the specified time
  for i = i_start:i_end

    if start_time_e > toepoch(pos_time(i,:))
      s_time = start_time;
    else
      s_time = pos_time(i,:);
    end

    if end_time_e < toepoch(pos_time(i,:))+dur_time(i)
      dur_t = end_time_e - start_time_e;
    else
      dur_t = dur_time(i);
    end

    disp(['processing: ' R_datestring(pos_time(i,:)) ' with length ' num2str(dur_time(i)/3600) ' hr'])
    tpm = -1;
    d_temp = 0;
    [tpm, d_temp] = c_ss_search_events(s_time,dur_t,dist2MP,p_solarwind);

    if tpm ~= -1
      passing_MP = add_A2M(passing_MP,tpm);
      dist_t = add_A2M(dist_t, d_temp);
    end

  end

  %if passing_MP ~= 0
  s_t = deblank(R_datestring(start_time));
  e_t = deblank(R_datestring(end_time));
  file_name = sprintf('MP_%s_to_%s',s_t,e_t);
  path_and_name = sprintf('%s%s',output_path, file_name);
  save(path_and_name, 'passing_MP', 'dist_t', 'dist2MP','p_solarwind')
  %end

end
