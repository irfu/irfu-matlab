function [z]=irf_vec_x_scal(x,y,p)
%IRF_VEC_X_SCAL   Calculate vec*scal^p
%
% function [z]=irf_vec_x_scal(x,y,p)
% estimates z=x*y^p;
% x - vector time series
% y - scalar time series
% p - power
%

if nargin==2,p=1;end
if isscalar(y), y=[x(:,1) ones(size(x(:,1)))*y]; end

if size(x,1) ~= size(y,1)
  irf_log('fcal','interpolating y to x')
  y = irf_resamp(y,x);
end

z=x;
for j=2:size(x,2)
  z(:,j)=x(:,j).*y(:,2).^p;
end
