function within = c_ri_timestr_within_tint(time_str,time_interval)
%
%within = c_ri_timestr_within_intervall(time_str,s_t,e_t)
%time_interval = c_ri_timestr_within_intervall(time_str)
%
%Input:
% time_str -timestring ex: "..x20020302x030301_x20030302x030306.mat"
% time_interval - interval in isdat_epoch [start_time end_time]
%
%Output:
% within -1 if there is a union of the two timeintervalls
%        -0 of there is no union
% time_interval - if only one input assume it is time string and return time interval in isdat epoch [t_start _end]
%
%Descrition of the function:'
% Finds if the to intervalls intersects each other
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
ts=time_str;

if nargin ==2
  st_e = time_interval(1);
  et_e = time_interval(2);
end

f_t=toepoch([str2num(ts(end-35+[0:3])) str2num(ts(end-31+[0:1])) ...
             str2num(ts(end-29+[0:1])) str2num(ts(end-26+[0:1])) ...
             str2num(ts(end-24+[0:1])) str2num(ts(end-22+[0:1]))]) ;

t_t=toepoch([str2num(ts(end-18+[0:3])) str2num(ts(end-14+[0:1])) ...
             str2num(ts(end-12+[0:1])) str2num(ts(end-9+[0:1])) ...
             str2num(ts(end-7+[0:1])) str2num(ts(end-5+[0:1]))]) ;

switch nargin
case 1 % return time interval
 within =[f_t t_t];
case 2 % return flag if time string is within specified time interval
  if st_e <= t_t && f_t <= et_e
   within = 1;
  else
   within = 0;
  end
end

