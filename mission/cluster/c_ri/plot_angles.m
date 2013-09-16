function u = plot_angles(angle) 
%
%Input:
% One matrix containing the angle between the satelites magnetic vectors and the
% time. 
% |time|1 to 2| 1 to 3|1 to 4|2 to 3|2 to 4|3 to 4|
%
%Output:
% One figure with 6 plots
%
%Descrition of the function:
%
%Using: irf_timeaxis
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
%calculating the first and last point
first = angle(1,1);
last = angle(max(size(angle)),1);
axis_size= [first last 0 180];

figure

zoom on
subplot(3,2,1)
plot(angle(:,1),angle(:,2))
title('The angle between the magnetic vectors at cluster 1 and cluster 2');
ylabel('angle');
axis(axis_size);
irf_timeaxis;

subplot(3,2,2)
plot(angle(:,1),angle(:,3))
title('The angle between the magnetic vectors at cluster 1 and cluster 3');
ylabel('angle');
axis(axis_size);
irf_timeaxis;

subplot(3,2,3)
plot(angle(:,1),angle(:,4))
title('The angle between the magnetic vectors at cluster 1 and cluster 4');
ylabel('angle');
axis(axis_size);
irf_timeaxis;

subplot(3,2,4)
plot(angle(:,1),angle(:,5))
title('The angle between the magnetic vectors at cluster 2 and cluster 3');
ylabel('angle');
axis(axis_size);
irf_timeaxis;

subplot(3,2,5)
plot(angle(:,1),angle(:,6))
title('The angle between the magnetic vectors at cluster 2 and cluster 4');
ylabel('angle');
axis(axis_size);
irf_timeaxis;

subplot(3,2,6)
plot(angle(:,1),angle(:,7))
title('The angle between the magnetic vectors at cluster 3 and cluster 4');
ylabel('angle');
axis(axis_size);
irf_timeaxis;

hold off
