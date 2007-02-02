function [st,dt] = irf_stdt(st_in,dt_in)
%CAA_STDT  parse input start and stop times
%
% [st,dt] = irf_stdt(start_t,dt)
% [st,dt] = irf_stdt(start_t,stop_t)
% [st,dt] = irf_stdt(iso_start_t,iso_stop_t)
%
%    Parse input parameters to produce interval start and length
%
% See also ISO2EPOCH
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev

error(nargchk(2,2,nargin))
error(nargoutchk(2,2,nargout))

if isnumeric(st_in), st = st_in;
elseif ischar(st_in), st = iso2epoch(st_in); 
else error('ST must be eather ISDAT epoch or ISO string')
end

if isnumeric(dt_in)
	if dt_in > 788918400, dt = dt_in -st; % 788918400 == 1995-01-01
	else dt = dt_in;
	end
elseif ischar(dt_in), dt = iso2epoch(dt_in) -st; 
else error('DT/ET must be eather numeric or ISO string')
end

if dt<=0, error('negative of zero DT'), end
