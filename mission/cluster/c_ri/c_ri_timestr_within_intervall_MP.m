function within = c_ri_timestr_within_intervall_MP(time_str,s_t,e_t)
%
%within = c_ri_timestr_within_intervall_MP(time_str,s_t,e_t)
%
%Input:
% time_str -timestring ex: "MP_20020302_03:03:01_to_20030302_03:03:06.mat"
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

mmddyy = sprintf('%s/%s/%s',time_str(8:9), time_str(10:11), time_str(6:7));
v_ymd = datevec(mmddyy);
fhms = sprintf('%s:%s:%s', time_str(13:14), time_str(16:17), time_str(19:20));
v_fhms = datevec(fhms);

tmmddyy = sprintf('%s/%s/%s',time_str(29:30), time_str(31:32), time_str(27:28));
v_tymd = datevec(tmmddyy);
thms = sprintf('%s:%s:%s', time_str(34:35), time_str(37:38), time_str(40:41));
v_thms = datevec(thms);

f_t = toepoch([v_ymd(1:3) v_fhms(4:6)]);
t_t = toepoch([v_tymd(1:3) v_thms(4:6)]);

if st_e <= t_t && f_t <= et_e
within = 1;
else
within = 0;
end
