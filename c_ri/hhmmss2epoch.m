function epoch_time = hhmmss2epoch(t_start, t)
%
% epoch_time = hhmmss2epoch(t_start, t
%
%Input:
% t_start -the epoch time for: 2002-03-01
% t_line -the time in format: 2002-03-01T07:37:28.096Z
%
%Output:
% epoch_time -the time in epoch
%
%Descrition of the function:
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
[r, c] = size(t);

t_L = c-5;
tic
hr_mm_ss = t(:,12:t_L);
t_str  = toc;
disp(['string cut: ' int2str(t_str)])

tic
hr = hr_mm_ss(:,1:2);
hr(1,:)
hr = str2num(hr);
hr(1,1)
mm = hr_mm_ss(:,4:5);
mm(1,:)
mm = str2num(mm);
mm(1,1) 
ss = hr_mm_ss(:,7:8);
ss(1,:)
ss = str2num(ss);
ss(1,1)
t_hemma = toc
disp(['home conversion: ' int2str(t_hemma)])

%tic
%hr_mm_ss = datevec(hr_mm_ss);
%t_datevec  = toc;
%disp(['datevec: ' int2str(t_datevec)])

tic
sec_hr_mm_ss = 3600*hr_mm_ss(:,4) + 60*hr_mm_ss(:,5) + hr_mm_ss(:,6);
epoch_time = sec_hr_mm_ss + t_start;
t_calc  = toc;
disp(['string calc: ' int2str(t_calc)])
