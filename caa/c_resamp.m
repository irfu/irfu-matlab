function out = c_resamp(x,y,varargin)
%C_RESAMP resample X to the time line of Y
% if sampling of X is more than two times higher than Y, we average X, 
% otherwise we interpolate X.
%
% out = c_resamp(x,y,varargin)
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)
%

error(nargchk(2,2,nargin))

% guess the sampling frequency
if min(size(y))==1, t = y(:); 
else, t = y(:,1); t = t(:);
end

ndata = length(t);
sfy = ndata/(t(end) - t(1));
if length(x(:,1))/(x(end,1) - x(1,1)) > 2*sfy
	% we average
	out = zeros(ndata,size(x,2));
	out(:,1) = t;
	dt2 = .5/sfy; % half interval
	for j=1:ndata
		ii = find(x(:,1) <  t(j) + dt2 & x(:,1) >  t(j) - dt2);
		if isempty(ii), out(j,2:end) = NaN;
		else, out(j,2:end) = mean(x(ii,2:end));
		end
	end
else
	% we interpolate
	out = av_interp(x,y,'linear');
end
