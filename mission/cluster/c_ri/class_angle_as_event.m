function [time_of_events,angles_out,ampl_out] = class_angle_as_event(angles,ampl, min_angle, min_ampl,mode) 
%function [time_of_events,angles,ampl] = class_angle_as_event(angles,ampl, min_angle, min_ampl,mode) 
%
%Input: 
% angles -[time | ang1 -> ang6] and ampl[ampl1 -> ampl 4]
% ampl -[ampl1 -> ample 4] 
% min_angle -the minimum angle to be called a event
% min_ampl -the minimum ampl to be called a event
%
%Output:
% time_of_events - [time in epoch | angle | amplitude | amplitude | mode]
%
%Descrition of the function:
% search angle and ampl for events. Saves the time, angle and ampl of the event
%
%Using:
% ind2nr
% 
%Work method:
%
%Error:
% if no events are found 0 is returned
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
%turning a zero into a zero-vector
time_of_events = [];
angles_out = [];
ampl_out = [];

[r_ang, c_ang] = size(angles);
if r_ang == 1 && c_ang == 1
  angles = [-1 0 0 0 0 0 0 0];
  disp('no data in the file')
end

[b_angle, b_angles_pos] = max(angles(:,2:7));
[biggest_angle,biggest_angle_pos] = max(b_angle);
ar_pos = b_angles_pos(biggest_angle_pos);
tba = angles(ar_pos,1);
tt = isnan(tba);
if tt == 1 || tba == -1
  tba = 0;
end

[a1_p,a2_p] = ind2nr(biggest_angle_pos);
a1 = ampl(ar_pos,a1_p);
a1t = isnan(a1);
if a1t == 1
  a1 = 0;
end
a2 = ampl(ar_pos,a2_p);
a2t = isnan(a2);
if a2t == 1
  a2 = 0;
end

ainfo = sprintf('max angle : %f, ampl1: %f ampl2: %f at time %s', biggest_angle,a1,a2,datestring(fromepoch(tba))); 
disp(ainfo)
[nr_angles, col] = size(angles);
%starting with zero events
b = 0;

[ang_i, ang_j] = find(angles(:,2:7) >= min_angle);
[ampl_i, ampl_j] = find(ampl(:,2:4) >= min_ampl);
ii = intersect(ang_i,ampl_i);

angles_out = angles(ii,:);
ampl_out = ampl(ii,:);
time_of_events = angles(ii,1);




