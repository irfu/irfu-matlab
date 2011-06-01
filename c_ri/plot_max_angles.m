function fp = plot_max_angles(max_angles) 
%
%Input:
% max_angles -[time|angle] the time in epochs
%
%Output:
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

fp = figure
plot(max_angles(:,1),max_angles(:,2))
irf_timeaxis;
title('The maximum of the rotation of the magnetic-field')
xlabel('time in UT')
ylabel('the rotation of the magnetic-field')

