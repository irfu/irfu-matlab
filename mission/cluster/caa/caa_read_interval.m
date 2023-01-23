function [st,dt] = caa_read_interval(sp)
%CAA_READ_INTERVAL  read caa interval information from .interval file
%
% [iso_t,dt] = caa_read_interval([sp])
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<1, sp=pwd; end

old_pwd = pwd;
cd(sp);
if exist('./.interval','file')
  [st_s,dt_s] = textread('./.interval','%s %s',-1);
  st_s = st_s{1}; dt_s = dt_s{1};
  if nargout==0
    disp([st_s ' : ' dt_s ' sec'])
  elseif nargout==1
    st = st_s;
  elseif nargout==2
    dt = str2double(dt_s);
    st = st_s;
  else
    error('wrong number of input arguments')
  end
else
  if nargout==0
    irf_log('load','cannot find .interval')
  elseif nargout==1
    st = '';
  elseif nargout==2
    dt = [];
    st = '';
  else
    error('wrong number of input arguments')
  end
end
cd(old_pwd)
