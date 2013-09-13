function view&filter_angles(path,filename) 
%
%Input:
%
%Output:.
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

disp('these files are possible to read')
!ls /share/robert/angle/A*.*
fn = input('Write the file to load, ' , 's');

p_and_fn = sprintf('/share/robert/angle/%s',fn);
load(p_and_fn);

[a_max,c] = size(angles);

from = fromepoch(angles(1,1));
to = fromepoch(angles(a_max,1));

disp(['There are data in the range: ' R_datestring(from) ' to 'R_datestring(to)])

start_time = input('from [yyyy mm dd hh mm ss.mss] (0 -from beginning): ');
end_time = input('to [yyyy mm dd hh mm ss.mss] (0 -from beginning): ');

%sets the start time in epochs, start_time goes from UT to epoch
if start_time == 0| toepoch(start_time) < toepoch(angles(1,1))
start_time = angles(1,1)
else
start_time = toepoch(start_time);
end

%sets the end time in epochs, end_time goes from UT to epoch
if end_time == 0 | toepoch(end_time) > toepoch(angles(a_max,1))
end_time = angles(a_max,1)
else
end_time = toepoch(start_time);
end

s_row = time2row(start_time, angles,1)
e_row = time2row(start_time, angles,s_row)
