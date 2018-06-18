function view_and_filter_angles() 
%
%view_and_filter_angles()
%
%Input:
%
%Output:.
%
%Descrition of the function:
% Loads a file with variable "angles" and "ampl". Filters the angle with min_ampl to get
% the max angle. Then plots the max angles.
% 
%
%Using:
% fromepoch
% toepoch
% time2row
% find_max_angles
% plot_max_angles
% class_angle_as_event 
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
ls /share/robert/angle/*.*
fn = input('Write the file to load, ' , 's');

p_and_fn = sprintf('/share/robert/angle/%s',fn);
load(p_and_fn);

[a_max,c] = size(angles);

from = fromepoch(angles(1,1));
to = fromepoch(angles(a_max,1));

disp(['There are data in the range: ' R_datestring(from) ' to ' R_datestring(to)])

start_time = input('from [yyyy mm dd hh mm ss.mss] (0 -from beginning): ');
end_time = input('to [yyyy mm dd hh mm ss.mss] (0 -from beginning): ');

%sets the start time in epochs, start_time goes from UT to epoch
if start_time == 0 || toepoch(start_time) < toepoch(angles(1,1))
start_time = angles(1,1);
else
start_time = toepoch(start_time);
end

%sets the end time in epochs, end_time goes from UT to epoch
if end_time == 0 || toepoch(end_time) > toepoch(angles(a_max,1))
end_time = angles(a_max,1);
else
end_time = toepoch(start_time);
end

s_row = time2row(start_time, angles,1);
e_row = time2row(end_time, angles,s_row);

min_ampl = input('Give the smallest amplitude of the B-vector to be classified as non-zero, ');

t_angles = angles(s_row:e_row,:);
t_ampl = ampl(s_row:e_row,:);
max_angles = find_max_angles(t_angles,t_ampl, min_ampl);

fp = plot_max_angles(max_angles);

min_angle = input('Give the smallest angle to class as an event: (0-180) ');
time_of_events = class_angle_as_event(t_angles,t_ampl, min_angle, min_ampl);

figure(fp)
hold on
plot(time_of_events(:,1), time_of_events(:,2),'rx')
hold off

save_or_not = input('Do you want to save the results?(y/n) ','s')
if save_or_not == 'y'
[a ,b] = size(fn);
tf = sprintf('/share/robert/events/E%s',fn(2:b));
disp(['save to: ' tf])
save(tf,'max_angles','time_of_events');
end





