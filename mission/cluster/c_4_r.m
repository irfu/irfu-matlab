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

if nargin<4
  disp('ERROR! 4 s/c positions should be given. c_4_r()');
  disp('See usage:')
  help c_4_r
  return
end

returnInverseVolumetricTensor = false; % default return volumetric tensor

if nargin>4 && flag == -1, returnInverseVolumetricTensor = true; end

indXYZ_r1 = min(size(r1,2),4)+[-2 -1 0];
indXYZ_r2 = min(size(r2,2),4)+[-2 -1 0];
indXYZ_r3 = min(size(r3,2),4)+[-2 -1 0];
indXYZ_r4 = min(size(r4,2),4)+[-2 -1 0];

R1 = r1(:, indXYZ_r1);
R2 = r2(:, indXYZ_r2);
R3 = r3(:, indXYZ_r3);
R4 = r4(:, indXYZ_r4);

R_Center = (R1 + R2 + R3 + R4)/4;

dR1 = R1-R_Center;
dR2 = R2-R_Center;
dR3 = R3-R_Center;
dR4 = R4-R_Center;

n = size(R1,1);
j = repmat(1:3, [3, 1]);
S_sum = repmat(dR1', 3, 1) .* dR1(:,j).' + ...
 repmat(dR2', 3, 1) .* dR2(:,j).' + ...
 repmat(dR3', 3, 1) .* dR3(:,j).' + ...
 repmat(dR4', 3, 1) .* dR4(:,j).';
volTensor = reshape(S_sum/4, [3 3 n]); % 1e6 to get into SI units km->m

% TODO vectorize (if possible)
if returnInverseVolumetricTensor || nargout>6
  computedVolTensor = volTensor;
  L = zeros(n,1); E = zeros(n,1); P = zeros(n,1);
  for i=1:n
    if(returnInverseVolumetricTensor)
      volTensor(:,:,i) = inv(computedVolTensor(:,:,i));
    end
    if nargout>6
      eigenValues = eig(computedVolTensor(:,:,i)); % eigenvalues of volumetric tensor
      eigenSqrt = sqrt(sort(eigenValues));
      L(i) = 2*eigenSqrt(3);                % characteristic size
      E(i) = 1 - eigenSqrt(2)/eigenSqrt(3); % elongation
      P(i) = 1 - eigenSqrt(1)/eigenSqrt(2); % planarity
    end
  end
end

end
