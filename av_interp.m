function [z]=av_interp(x,y,method)
% AV_INTERP Interpolate time series
%    [z]=av_interp(x,y,method) interpolates x to the time line of y
%    x,y,z - row vectors where first column is time
%    z has y-time and x-column number
%    method - method of interpolation 'spline', 'linear' etc. (default 'linear')
%
% See also: INTERP1
%

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_resamp')

if nargin < 3, method='linear';end
if nargin<2;  disp('Not enough arguments. See usage:');help av_interp;      return;end

if size(x,1)==size(y,1), if x(:,1)==y(:,1),  z=x; return;  end, end  % no interpolation necessary as time series agree

if size(x,1) > 1,
  z=[y(:,1) interp1(x(:,1),x(:,2:end),y(:,1),method,'extrap')];
else
  z=[y(:,1) (y(:,1)*0+1)*x(:,2:end)];
end
