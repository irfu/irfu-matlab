function angle_table = four_vector_angles(vector_table) 
%
%Input:
% A table containing four, non-zero, vectors. They must all have the same 
% length. Each row is a vector.  
%
%Output:
% A table containing all the six angles in this order:
% |v_1 to v_2 |v_1 to v_3 |v_1 to v_4 |v_2 to v_3 |v_2 to v_4 |v_3 to v_4 |
% The output is in absolut angle in the range 0-180 degrees
% 
%Descrition of the function:
% Calculates all the six angles between all the four vectors.
%
%Using: vector_angles
% 
%Work method:
% Six functioncall to vector_angles.
%
%Error:
% If one of the vector i a zero-vector then all the angles compared with
% that vector will be zero.
% 
%Discription of variables:
% vector_table is an array containg four vectors. Each row contains a vector.
% i is a loop index
% angle_table is the output table
%
%Written by Robert Isaksson in the summer of -03

%--------------------- the beginning --------------------------

angle_table(1) = vector_angles(vector_table(1,:), vector_table(2,:));
angle_table(2) = vector_angles(vector_table(1,:), vector_table(3,:));
angle_table(3) = vector_angles(vector_table(1,:), vector_table(4,:));
angle_table(4) = vector_angles(vector_table(2,:), vector_table(3,:));
angle_table(5) = vector_angles(vector_table(2,:), vector_table(4,:));
angle_table(6) = vector_angles(vector_table(3,:), vector_table(4,:));
