function [z]=irf_vec_x_scal(x,y,p)
%IRF_VEC_X_SCAL   Calculate vec*scal^p
%
% function [z]=irf_vec_x_scal(x,y,p)
% estimates z=x*y^p;
% x - vector time series
% y - scalar time series
% p - power
%
% $Id$

global AV_DEBUG;if isempty(AV_DEBUG), debug=0; else debug=AV_DEBUG;end

if nargin==2,p=1;end

if size(x,1) ~= size(y,1),
 if debug ==1, disp('interpolating y to x, irf_vec_x_scal()');end
 y=av_interp(y,x);
end

z=x;
for j=2:size(x,2);
 z(:,j)=x(:,j).*y(:,2).^p;
end
