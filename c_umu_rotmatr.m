function [R] = c_umu_rotmatr(alpha, beta, gamma)

% function [R] = c_umu_rotmatr(alpha, beta, gamma)
% 
%   (alpha, beta, gamma) are angles in radians.
%
%  Computes the matrix for rotation through
%  1. the angle alpha about the z axis to give (x', y', z'),
%  2. the angle beta  about the y' axis to give (x'', y'', z''),
%  3. and the angle gamma about the z'' axis to give (x_, y_, z_)
%  which are the coordinates in the new system.

R = ...
[ cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)  sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma)  -sin(beta)*cos(gamma); 
 -cos(alpha)*cos(beta)*sin(gamma)-sin(alpha)*cos(gamma) -sin(alpha)*cos(beta)*sin(gamma)+cos(alpha)*cos(gamma)   sin(beta)*sin(gamma);
  cos(alpha)*sin(beta)                                   sin(alpha)*sin(beta)                                    cos(beta);           ];
 



