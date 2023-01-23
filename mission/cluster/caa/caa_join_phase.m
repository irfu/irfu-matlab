function caa_join_phase(cl_id,gap)
%CAA_JOIN_PHASE  refetch missing phase
%
% caa_join_phase(cl_id,iso_t)
% caa_join_phase(cl_id,gap_len)
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

[iso_t,dt] = caa_read_interval;
st = iso2epoch(iso_t);

irf_log('proc',['interval: ' epoch2iso(st,1) ' -- ' epoch2iso(st+dt,1) ])

[AOK,A_old] = c_load('A?',cl_id);
[AtwoOK,Atwo_old] = c_load('Atwo?',cl_id);

if isnumeric(gap)
  if AtwoOK, old_end = Atwo_old(end,1);
  else, old_end = st;
  end
  st_refetch = old_end + gap;
elseif ischar(gap), st_refetch = iso2epoch(gap);
else, error('invalid input')
end

if st_refetch < st || st_refetch > st+dt
  error([epoch2iso(st_refetch,1) ' is outside the interval'])
end

irf_log('proc',['refetching phase from: ' epoch2iso(st_refetch,1) ])
getData(ClusterDB,st_refetch,st+dt - st_refetch,cl_id,'a')



[AOK,A] = c_load('A?',cl_id);
[AtwoOK,Atwo] = c_load('Atwo?',cl_id);

if ~AtwoOK, irf_log('load','refetch failed'), return, end

if ~isempty(A_old), A =[A_old; A]; end
if ~isempty(Atwo_old), Atwo =[Atwo_old; Atwo]; end


irf_log('save','saving updated phase')
c_eval('A?=A; Atwo?=Atwo; figure, irf_plot(Atwo?), save mA.mat A? Atwo? -append',cl_id)