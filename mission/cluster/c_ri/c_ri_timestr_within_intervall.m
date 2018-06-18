function within = c_ri_timestr_within_intervall(time_str,s_t,e_t)
%
%within = c_ri_timestr_within_intervall(time_str,s_t,e_t)
%
%Input:
% time_str -timestring ex: "Ba_20020302_F030301_T030306"
% s_t -beginning of intervall ex: [2002 03 02 0 0 0]
% e_t -end of intervall ex:[2002 03 02 0 0 0]
%
%Output:
% within -1 if there is a union of the two timeintervalls
%        -0 of there is no union
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
st_e = toepoch(s_t);
et_e = toepoch(e_t);

mmddyy = sprintf('%s/%s/%s',time_str(6:7), time_str(8:9), time_str(4:5));
v_ymd = datevec(mmddyy);
fhms = sprintf('%s:%s:%s', time_str(12:13), time_str(14:15), time_str(16:17));
v_fhms = datevec(fhms);
thms = sprintf('%s:%s:%s', time_str(20:21), time_str(22:23), time_str(24:25));
v_thms = datevec(thms);

f_t = toepoch([v_ymd(1:3) v_fhms(4:6)]);
t_t = toepoch([v_ymd(1:3) v_thms(4:6)]);

if st_e <= t_t && f_t <= et_e
within = 1;
else
within = 0;
end
