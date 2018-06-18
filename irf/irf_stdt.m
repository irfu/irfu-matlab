function [st,dt] = irf_stdt(st_in,dt_in)
%IRF_STDT  parse input start and stop times
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

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)
nargoutchk(2,2)

if isnumeric(st_in), st = st_in;
elseif ischar(st_in), st = iso2epoch(st_in); 
else, error('ST must be eather ISDAT epoch or ISO string')
end

if isnumeric(dt_in)
	if dt_in > 788918400, dt = dt_in -st; % 788918400 == 1995-01-01
	else, dt = dt_in;
	end
elseif ischar(dt_in), dt = iso2epoch(dt_in) -st; 
else, error('DT/ET must be eather numeric or ISO string')
end

if dt<=0, error('negative or zero DT'), end
