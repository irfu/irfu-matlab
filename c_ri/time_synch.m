function [B_vectors, n_positions] = time_synch(t, dt, B1, n1, B2, n2, B3, n3, B4, n4) 
%
%Input:
% t, is the time at the previous (n-1) sample
% dt, is the time between samples, where t(n-1)+dt = t(n)
% n is the current position in the vector
% b1, b2, b3 ,b4 are the vectors containg the time i row 1 and B three component for
% sampling of magnetic field.
% n1, n2, n3, n4 are the last postion (n-1) in each of the time matrix 
%
%Output:
% A matrix containing vector synchronized in time. One row is one vector.
% If one row is zero the this matrix har no valid value.
% N-positions is where in the matrix this vector i positioned.
%
%Descrition of the function:
% If one of the cluster satellites
% doesnt have delivered a valid value, the value in the vector will
% be -1. Otherwise the value will be >0  
%
%Using: find_vector_to_time
% 
%Work method: calls find_vector_to_time four times
%
%Error: Zero-vector i return if one 
% 
%Discription of variables:
% t, is the time at the previous (n-1) sample
% dt, is the time between samples, where t(n-1)+dt = t(n)
% n, is the current position in the vector
% B1, B2, B3, B4 are the vectors containing the time and the 
% sampling of magnetic field.
%

%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

new_time = t + dt;

[B_vectors(1,:) ,n_positions(1)] = find_vector_to_time(B1, n1, new_time, dt);

[B_vectors(2,:) ,n_positions(2)] = find_vector_to_time(B2, n2, new_time, dt);

[B_vectors(3,:) ,n_positions(3)] = find_vector_to_time(B3, n3, new_time, dt);

[B_vectors(4,:) ,n_positions(4)] = find_vector_to_time(B4, n4, new_time, dt);



