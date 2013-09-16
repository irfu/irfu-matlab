function [volTensor,R_Center,dR1,dR2,dR3,dR4,L,E,P]=c_4_r(r1,r2,r3,r4,flag)
%C_4_R  Calculate volumetric tensor R
%
% volTensor=c_4_r(r1,r2,r3,r4) returns the volumetric tensor.
%	r1..r4 are row vectors with the satellite positions,
%	if more than 3 columns assumes that first is time and 
%	for x,y,z components uses columns 2,3,4
%
% volTensor=c_4_r(r1,r2,r3,r4,-1) returns the inverse of volumetric tensor.
%
% [volTensor,R_Center,dR1,dR2,dR3,dR4,L,E,P]=c_4_r(..)
%	Returns also the coordinates of barycenter and the satellite relative 
%	position with respect to that. In addition the are returned
%	L - characgeristic size
%	E - elongation
%	P - planarity
%
%  The volumetric tensor is calculated with respect to the center of r1..r4.
%  The units of volumetric tensor is the same as r.
%
%  Reference: ISSI book
%  R=1/N Sum[sc=1->N] (r_sc -r_Center) (r_sc-r_Center)^T
%  r_Center = 1/N Sum[sc=1->N] r_sc
%

% TODO: vectorize

% $Id$


if nargin<4
	disp('ERROR! 4 s/c positions should be given. c_4_r()');
	disp('See usage:')
	help c_4_r
	return
end

returnInverseVolumetricTensor = false; % default return volumetric tensor

if nargin>4 && flag == -1,
	returnInverseVolumetricTensor = true;
end

indXYZ_r1=min(size(r1,2),4)+[-2 -1 0];
indXYZ_r2=min(size(r2,2),4)+[-2 -1 0];
indXYZ_r3=min(size(r3,2),4)+[-2 -1 0];
indXYZ_r4=min(size(r4,2),4)+[-2 -1 0];

R1=r1(1,indXYZ_r1);
R2=r2(1,indXYZ_r2);
R3=r3(1,indXYZ_r3);
R4=r4(1,indXYZ_r4);

R_Center=(R1+R2+R3+R4)/4;

dR1=R1-R_Center;
dR2=R2-R_Center;
dR3=R3-R_Center;
dR4=R4-R_Center;

volTensor=(dR1'*dR1+dR2'*dR2+dR3'*dR3+dR4'*dR4)/4; % 1e6 to get into SI units km->m

if nargout >6 % calculate also configuration parameters
	eigenValues = eig(volTensor);                   % eigenvalues of volumetric tensor
	eigenValuesSorted = sort(eigenValues);
	
	a = sqrt(eigenValuesSorted(3));
	b = sqrt(eigenValuesSorted(2));
	c = sqrt(eigenValuesSorted(1));
	
	L = 2*a;        %characteristic size
	E = 1 - b/a;    %elongation
	P = 1 - c/b;    %planarity
end

if returnInverseVolumetricTensor
	volTensor=inv(volTensor);
end


