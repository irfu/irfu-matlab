function gradPhi = c_4_gradphi(r1, r2, r3, r4, u1, u2, u3, u4)
%C_4_GRADPHI  Calculate spatial gradient of scalar
%
%  gradPhi = c_4_gradphi(r1, r2, r3, r4, u1, u2, u3, u4)  calculate spatial
%  gradient of scalar. The gradient is in the same units as u.
%
%  r1..r4 are position vectors in TSeries or nx2 array (1st column time)
%  u1..u4 are scalars in TSeries or nx2 array (1st column time)
%
%  Reference: ISSI book SR-001 Eq.12.16
%  http://www.issibern.ch/PDF-Files/analysis_methods_1_1a.pdf
%  k_l=S_l R_kl^-1
%  where S_l=(Sum a!=b dx_ab dRab)/N^2   N=4 number of satellites

if nargin<8
  disp('Too few parameters. See usage:');
  help c_4_gradPhi;
  return
end
% check input TSeries or array
if all([isa(r1,'TSeries'), isa(r2,'TSeries'), isa(r3,'TSeries'), ...
    isa(r4,'TSeries'), isa(u1,'TSeries'), isa(u2,'TSeries'), ...
    isa(u3,'TSeries'), isa(u4,'TSeries')])
  output_as_TSeries = true;
  tref = u1.time;
  name = u1.name;
  units = u1.units;
  % everything is resampled to u1
  c_eval('r?=resample(r?,u1);')
  c_eval('u?=resample(u?,u1);')
  c_eval('r?=r?.data(:,1:3);')
  c_eval('u?=u?.data(:,1);')
elseif all([isnumeric(r1), isnumeric(r2), isnumeric(r3), isnumeric(r4), ...
    isnumeric(u1), isnumeric(u2), isnumeric(u3), isnumeric(u4)])
  output_as_TSeries = false;
  % everything is resampled to u1
  c_eval('r?=irf_resamp(r?,u1);')
  c_eval('u?=irf_resamp(u?,u1);')
  c_eval('r?=r?(:,2:4);')
  c_eval('u?=u?(:,2);')
else
  error('Input type not consistent.')
end

c_eval('r?=double(r?); u?=double(u?);');
R_Center = (r1 + r2 + r3 + r4)/4;
dR1 = r1 - R_Center; dR2 = r2 - R_Center;
dR3 = r3 - R_Center; dR4 = r4 - R_Center;

S = (u2-u1) .* (dR2-dR1) + ...
  (u3-u1) .* (dR3-dR1) + ...
  (u3-u2) .* (dR3-dR2) + ...
  (u4-u1) .* (dR4-dR1) + ...
  (u4-u2) .* (dR4-dR2) + ...
  (u4-u3) .* (dR4-dR3);
S = S/16; % with "/16 = /N^2 for N=4"

Rinv = c_4_r(dR1, dR2, dR3, dR4, -1); % inverse of volumetric tensor

n = size(r1,1);
gradPhi = zeros(n, 3);

for i=1:n % TO BE VECTORIZED AT SOME POINT.
  gradPhi(i,:) = S(i,:)*Rinv(:,:,i);
end

% Output
if output_as_TSeries
  gradPhi = irf.ts_vec_xyz(tref, gradPhi);
  gradPhi.name = ['Grad of ' name];
  gradPhi.units = units;
end

end