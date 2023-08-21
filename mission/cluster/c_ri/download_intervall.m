function [B,download_logg]  = download_intervall(Btdt, start_time,start_row, end_time, end_row, cl_nr)
%
%Input:
% start_time -in epochs
% end_time -in epochs
% Btdt - matrix with [time, dt] which are the start of an event
%        and the duration of that event.
% start_row -which row that corresponds to the start_time
% end_row -which row that corresponds to the end_time
% cl_nr -which clustersatellite
% block_time -the length of the timeblock to be downloaded... not in use
%
%Output:
%
%Descrition of the function:
%
%Using:
% get_B_from_ISDAT
% add_A2M
%
%Work method:
%
%Error:
%
%Description of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
B = 0;
B_temp = 0;
download_logg = 0;

for k = start_row:end_row
  %disp(['download intervall, row ' int2str(k)])

  time = Btdt(k,1);
  if start_time > time && k == start_row
    time = start_time;
  end

  time_e = Btdt(k,1) + Btdt(k,2);
  if end_time < time_e && k == end_row
    time_e = end_time;
  end


  dt = floor(time_e - time);
  time_U = fromepoch(ceil(time));

  if dt < 0
    disp(['start_time : ' R_datestring(fromepoch(time)) ' end_time: ' R_datestring(fromepoch(time_e)) 'start row ' int2str(start_row) ' end row: ' int2str(end_row)])
    input('continue')
  end

  %disp(['download intervall ' int2str(time_U) ' dt ' int2str(dt)])
  clear B_temp
  B_temp = 0;
  B_temp = get_B_from_ISDAT(cl_nr, time_U, dt);
  B = add_A2M(B,B_temp);

  d_times = [time dt -1];
  if B_temp(1,1) ~= -1
    d_times(1,3) = 1;
  end
  download_logg = add_A2M(download_logg, d_times);

end
