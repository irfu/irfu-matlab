function [st,dt] = irf_stdt(st_in,dt_in)
%IRF_STDT  parse input start and stop times
%
% [st,dt] = irf_stdt(start_t,dt)
% [st,dt] = irf_stdt(start_t,stop_t)
% [st,dt] = irf_stdt(iso_start_t,iso_stop_t)
%
%           If first argument is numeric, then it is used as is.
%           If second argument is numeric and >788918400,
%           then it is interpreted as an absolut time
%           If second argument is numeric and <=788918400,
%           then it is interpreted as a length of time.
%           If an argument is a string, then it is interpreted as an absolute
%           time through iso2epoch().
%
%    Parse input parameters to produce interval start and length
%
% See also ISO2EPOCH
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,2)
nargoutchk(2,2)

if     isnumeric(st_in), st = st_in;
elseif ischar(st_in),    st = iso2epoch(st_in);
else, error('ST must be either ISDAT epoch or ISO string')
end

if isnumeric(dt_in)
  if dt_in > 788918400, dt = dt_in -st; % 788918400 == 1995-01-01
  else, dt = dt_in;
  end
elseif ischar(dt_in), dt = iso2epoch(dt_in) -st;
else, error('DT/ET must be either numeric or ISO string')
end

if dt<=0, error('negative or zero DT'), end
