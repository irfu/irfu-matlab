function out = irf_resamp(x,y,method)
%IRF_RESAMP   Resample X to the time line of Y
%
% if sampling of X is more than two times higher than Y, we average X,
% otherwise we interpolate X.
%
% out = irf_resamp(x,y,[method])
% method - method of interpolation 'spline', 'linear' etc. (default 'linear')
%          if method is given the interpolate independant of sampling
%
% See also INTERP1
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)
%

error(nargchk(2,3,nargin))
if nargin==3,
  flag_do='interpolation'; % if method is given do interpolation
else
  flag_do='check'; % if no method check if interpolate or average
end

% return in case no time axis
if prod(size(y))==0,
    out=[];return;
end

% construct output time axis
if size(y,2)==1, t = y(:); % y is only time
else, t = y(:,1); t = t(:);   % first column of y is time
end
ndata = length(t);

if size(x,1) == 1,            % if X has only one point, this is a trivial
                              % case and we return directly
  out = [t (t*0+1)*x(:,2:end)];
  return
end

if strcmp(flag_do,'check'), % check if interpolation or average
  if ndata>1, % if more than one output time check sampling frequencies 
              % to decide interpolation/average
    sfy = ndata/(t(end) - t(1)); % guess samplings frequency y
    if length(x(:,1))/(x(end,1) - x(1,1)) > 2*sfy
      flag_do='average';
    else
      flag_do='interpolation';
    end
  else
    flag_do='interpolation';  % if one output time then do interpolation
  end
end

if strcmp(flag_do,'average'),
  out = zeros(ndata,size(x,2));
  out(:,1) = t;
  dt2 = .5/sfy; % half interval
  for j=1:ndata
    ii = find(x(:,1) <  t(j) + dt2 & x(:,1) >  t(j) - dt2);
    if isempty(ii), out(j,2:end) = NaN;
    else, out(j,2:end) = mean(x(ii,2:end));
    end
  end
elseif strcmp(flag_do,'interpolation'),
  if nargin < 3, method = 'linear'; end

  % If time series agree, no interpolation is necessary.
  if size(x,1)==size(y,1), if x(:,1)==y(:,1), out = x; return, end, end

  out = [y(:,1) interp1(x(:,1),x(:,2:end),y(:,1),method,'extrap')];
end
