function [nr_1, nr_2] = ind2nr(ind) 
%
%Input:
% ind -index of the position of and angle
%
%Output:
% nr_1, nr_2 -positions in the B_ampl.
%
%Descrition of the function:
% finding the positions in the B_ampl who correspond to a certain
% position in time_and_angles
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

if ind == 1
nr_1 = 1;
nr_2 = 2;
end

if ind == 2
nr_1 = 1;
nr_2 = 3;
end

if ind == 3
nr_1 = 1;
nr_2 = 4;
end

if ind == 4
nr_1 = 2;
nr_2 = 3;
end

if ind == 5
nr_1 = 2;
nr_2 = 4;
end

if ind == 6
nr_1 = 3;
nr_2 = 4;
end

