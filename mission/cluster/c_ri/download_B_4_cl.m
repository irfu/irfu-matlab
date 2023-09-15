function [B, download_logg] = download_B_4_cl(passing_MP, cl_nr, B_t_dt)
%
%Input:
% passing_MP - the information
%         	about the crossing of the MP. [entering | leaving], where
%        	the times are in epoch.
% cl_nr -nr of which satellite the download is for
% B_t_dt -the time and duration of the downloads
%
%Output:
% B -a matrix with all the B-data for the MP-crossing.
%
%Descrition of the function:
%
%Using:
% find_row
% download_intervall
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
[nr_MP_cross,col] = size(passing_MP);
%this loop is for every MP-crossing
download_logg = 0;

for i = 1:nr_MP_cross
  start_time = passing_MP(i,1);
  end_time = passing_MP(i,2);

  disp( ' ')
  disp(['download_B_4_cl, MP-crossing ' R_datestring(fromepoch(start_time)) ' to '  R_datestring(fromepoch(end_time))])

  [start_row, end_row] = find_row(start_time, end_time, B_t_dt,1);
  if end_row < start_row || start_time > end_time
    %something to avoid
    disp(['start_time : ' R_datestring(fromepoch(start_time)) ' end_time: ' R_datestring(fromepoch(end_time)) ' start row ' int2str(start_row) ' end row: ' int2str(end_row)])
    input('continue')
  end

  clear B_temp
  B_temp = 0;

  if start_row ~= -1 && end_row  ~= -1
    [B_temp,t_download_logg] = download_intervall(B_t_dt, start_time, start_row, end_time, end_row, cl_nr);
    B = add_A2M(B,B_temp);
    download_logg = add_A2M(download_logg, t_download_logg);
  else
    t_download_logg = [start_time end_time -1];
    download_logg = add_A2M(download_logg, t_download_logg);
  end

end % end of rows in "passing_MP"
