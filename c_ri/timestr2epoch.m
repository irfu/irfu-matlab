function time = timestr2epoch(strline) 
%
%time = timestr2epoch(strline) 
%
%Input:
% strline -ex 2002-03-25T04:11:00.871Z
%
%Output:
% time -time in epoch
%
%Descrition of the function:
% convertes the time string in to epoch time. Removes the time Z and the .xxx
% in the end of the time string.
%
%Using:
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
str_length = max(size(strline));
str_length = str_length-5;
t_line = strline(1:str_length);

yyyy_mm_dd = t_line(1:10);
mm_dd_yyyy(4:5) = yyyy_mm_dd(9:10);
mm_dd_yyyy(1:2) = yyyy_mm_dd(6:7);
mm_dd_yyyy(7:10) = yyyy_mm_dd(1:4);
mm_dd_yyyy(3) = '-';
mm_dd_yyyy(6) = '-';

yyyy_mm_dd = datevec(mm_dd_yyyy);

hr_mm_ss = t_line (12:length(t_line));
hr_mm_ss = datevec(hr_mm_ss);

time_UT = [yyyy_mm_dd(1:3) hr_mm_ss(4:6)];
time = toepoch(time_UT);
