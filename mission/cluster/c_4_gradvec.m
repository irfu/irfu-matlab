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
%  Reference: ISSI book SR-001 Eq.12.18
%  http://www.issibern.ch/PDF-Files/analysis_methods_1_1a.pdf
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

c_eval('r?=double(r?); u?=double(u?);');
R_Center = (r1 + r2 + r3 + r4)/4;
dR1 = r1-R_Center; dR2 = r2-R_Center; dR3 = r3-R_Center; dR4 = r4-R_Center;
% Differences in dR (aka "dR_ab, for all n(n-1)/2 independent terms with a ~= b")
dR_21 = dR2 - dR1;
dR_31 = dR3 - dR1; dR_32 = dR3 - dR2;
dR_41 = dR4 - dR1; dR_42 = dR4 - dR2; dR_43 = dR4 - dR3;
j = repmat(1:3, [3, 1]);
S_sum = repmat((u2-u1)', 3, 1) .* dR_21(:,j).' + ...
 repmat((u3-u1)', 3, 1) .* dR_31(:,j).' + ...
 repmat((u3-u2)', 3, 1) .* dR_32(:,j).' + ...
 repmat((u4-u1)', 3, 1) .* dR_41(:,j).' + ...
 repmat((u4-u2)', 3, 1) .* dR_42(:,j).' + ...
 repmat((u4-u3)', 3, 1) .* dR_43(:,j).';
S = reshape(S_sum/16, [3 3 n]); % with "/16 = /N^2 for N=4"

Rinv = c_4_r(dR1, dR2, dR3, dR4, -1);

for i=1:n % TO BE VECTORIZED AT SOME POINT.
  gradVec(i,:,:) = S(:,:,i)*Rinv(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%
if nargout==0||nargout==1
  if output_as_TSeries
    gradVec = TSeries(tref, gradVec, 'TensorOrder', 2, ...
      'repres', {'x','y','z'}, 'repres',{'x','y','z'});
    gradVec.name = ['Grad of ' name];
    gradVec.units = units;
  end
  varargout(1) = {gradVec};
elseif nargout == 2
  % Sum of the diagonal elements for each 3D "slice"
  div = sum(gradVec(bsxfun(@plus,(1:n)',[0, 4*n, 8*n])), 2);
  curl = [gradVec(:,2,3)-gradVec(:,3,2), -gradVec(:,1,3)+gradVec(:,3,1), gradVec(:,1,2)-gradVec(:,2,1)];
  if output_as_TSeries
    div = TSeries(tref, div, 'TensorOrder', 0);
    div.name = ['div of ' name];
    div.units = units;
    curl = TSeries(tref, curl, 'TensorOrder', 1, 'repres', {'x','y','z'});
    curl.name = ['Curl of ' name];
    curl.units = units;
  end
  varargout(1) = {div};
  varargout(2) = {curl};
end

end
