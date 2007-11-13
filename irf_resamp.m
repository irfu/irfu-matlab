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

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,3,nargin))
if nargin==3,
	flag_do='interpolation'; % if method is given do interpolation
else
	flag_do='check'; % if no method check if interpolate or average
end

% return in case inputs are empty
if numel(x)==0 || numel(y)==0
    out=[];return;
end

% construct output time axis
if size(y,2)==1, t = y(:); % y is only time
else t = y(:,1); t = t(:);   % first column of y is time
end
ndata = length(t);

if size(x,1) == 1,            % if X has only one point, this is a trivial
                              % case and we return directly
  out = [t (t*0+1)*x(:,2:end)];
  return
end

if strcmp(flag_do,'check'), % Check if interpolation or average
	if ndata>1 
		% If more than one output time check sampling frequencies
		% to decide interpolation/average
		
		% Guess samplings frequency y
		sfy1 = 1/(t(2) - t(1));
		if ndata==2, sfy = fsy1;
		else
			not_found = 1; cur = 3; MAXTRY = 10;
			while (not_found && cur<=ndata && cur-3<MAXTRY)
				sfy = 1/(t(cur) - t(cur-1));
				if abs(sfy-sfy1)<sfy*0.001
					not_found = 0;
					sfy = (sfy+sfy1)/2;
					break
				end
				cur = cur + 1;
			end
			if not_found
				sfy = fsy1;
				irf_log('proc',	sprintf(...
					'Cannot guess sampling frequency. Tried %d times',MAXTRY));
			end
		end
		clear sfy1
		if length(x(:,1))/(x(end,1) - x(1,1)) > 2*sfy
			flag_do='average';
		else
			flag_do='interpolation';
		end
	else
		flag_do='interpolation';  % If one output time then do interpolation
	end
end

if strcmp(flag_do,'average')
	dt2 = .5/sfy; % Half interval
	if exist('irf_average_mx','file')~=3
		irf_log('fcal','cannot find mex file, defaulting to Matlab code.')

		out = zeros(ndata,size(x,2));
		out(:,1) = t;
		for j=1:ndata
			ii = find(x(:,1) <  t(j) + dt2 & x(:,1) >  t(j) - dt2);
			if isempty(ii), out(j,2:end) = NaN;
			else out(j,2:end) = mean(x(ii,2:end));
			end
		end
	else
		out = irf_average_mx(x,t,dt2);
	end
elseif strcmp(flag_do,'interpolation'),
  if nargin < 3, method = 'linear'; end

  % If time series agree, no interpolation is necessary.
  if size(x,1)==size(y,1), if x(:,1)==y(:,1), out = x; return, end, end

  out = [y(:,1) interp1(x(:,1),x(:,2:end),y(:,1),method,'extrap')];
end
