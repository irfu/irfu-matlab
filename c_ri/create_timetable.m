function [timetable_b, timetable_n] = create_timetable(from,to,cl_nr) 
%
%Input:
% from -time in epoch
% to -time in epoch
% from and to must be within the same day
% cl_nr -which cluster satellite
%
%Output:
% timetable - [from | to | mode] , where mode (b/n), burst or normal
% for which times there are data
%
%Descrition of the function:
% This function creates a timetable for which times the are 
% data possible to download. And it is in cronological order.
%
%Using:
% create_file
% get_timetable
% 
%Work method:
%
%Error:
% 
%Discription of variables:
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------
file_p_n_b = create_file(from,'b',1);
file_p_n_n = create_file(from,'n',1);

timetable_b = get_timetable(from,to,file_p_n_b);
timetable_n = get_timetable(from,to,file_p_n_n);

%unix_com = sprintf('rm %s %s', file_p_n_b ,file_p_n_n);
%unix(unix_com);
