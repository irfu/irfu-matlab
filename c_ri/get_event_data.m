function get_event_data(fp) 
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
g=50;

fp = figure;
h1 = uicontrol('Position', [10,10,100,100]);
set(h1,'CallBack','get(h3,''value''),');
set(h1,'Style','slider');
set(h1,'Max',g);
set(h1,'Min',1);
set(h1,'SliderStep',[2 3])

get(h1)
