function [varargout]=c_4_gradvec(r1,r2,r3,r4,u1,u2,u3,u4)
% C_4_GRADVEC calculate spatial gradient tensor of vector
%
%  [gradvec]=c_4_gradvec(r1,r2,r3,r4,u1,u2,u3,u4)  calculate spatial gradient
%  tensor of vector. The gradient is in the same units as u.
%
%  [div,curl]=c_4_gradvec(r1,r2,r3,r4,u1,u2,u3,u4)  calculate divergence and curl of vector
%
%  r1..r4 are position vectors in TSeries or nx4 array (1st column time)
%  u1..u4 are vectors in TSeries or nx4 array (1st column time)
%
%  Reference: ISSI book  Eq.12.18
%  k_l=S_l R_kl^-1
%  where S_l=(Sum a!=b dx_ab dR_ab)/N^2   N=4 number of satellites
%

if nargin<8; disp('Too few parameters. See usage:'); help c_4_gradVec; return; end

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
  c_eval('u?=u?.data(:,1:3);')
elseif all([isnumeric(r1), isnumeric(r2), isnumeric(r3), isnumeric(r4), ...
    isnumeric(u1), isnumeric(u2), isnumeric(u3), isnumeric(u4)])
  output_as_TSeries = false;
  % everything is resampled to u1
  c_eval('r?=irf_resamp(r?,u1);')
  c_eval('u?=irf_resamp(u?,u1);')
  c_eval('r?=r?(:,2:4);')
  c_eval('u?=u?(:,2:4);')
else
  error('Input type not consistent.')
end

n = size(r1,1);
gradVec = zeros(n,3,3);
for i=1:n % to be vectorized at some point
  r1i=r1(i,1:3); r2i=r2(i,1:3); r3i=r3(i,1:3); r4i=r4(i,1:3);
  u1i=u1(i,1:3); u2i=u2(i,1:3); u3i=u3(i,1:3); u4i=u4(i,1:3); %#ok<NASGU>
  R_Center = (r1i+r2i+r3i+r4i)/4;
  dR1 = r1i-R_Center; dR2 = r2i-R_Center;
  dR3 = r3i-R_Center; dR4 = r4i-R_Center;
  S = 0;
  for a=1:4
    for b=1:a
      S = S + eval(irf_ssub('(u?i-u!i)''*(dR?-dR!)',a,b));
    end
  end
  S = S/16;
  Rinv = c_4_r(dR1, dR2, dR3, dR4, -1);  % inverse of volumetric tensor
  gradVec(i,:,:) = S*Rinv;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0||nargout==1
  if output_as_TSeries
    gradVec = TSeries(tref, gradVec, 'TensorOrder', 2, ...
      'repres', {'x','y','z'}, 'repres',{'x','y','z'});
    gradVec.name = ['Grad of ' name];
    gradVec.units = units;
    varargout(1) = {gradVec};
  else
    varargout(1) = {gradVec};
  end
elseif nargout == 2
  div = trace(gradVec);
  curl = [gradVec(2,3)-gradVec(3,2) -gradVec(1,3)+gradVec(3,1) gradVec(1,2)-gradVec(2,1)];
  if output_as_TSeries
    div = TSeries(tref, div, 'TensorOrder', 0);
    div.name = ['div of ' name];
    div.units = units;
    varargout(1) = {div};
    curl = TSeries(tref, curl, 'TensorOrder', 1, 'repres', {'x','y','z'});
    curl.name = ['Curl of ' name];
    curl.units = units;
    varargout(2) = {curl};
  else
    varargout(1) = {div};
    varargout(2) = {curl};
  end
end

end