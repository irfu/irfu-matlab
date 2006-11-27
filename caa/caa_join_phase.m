function caa_join_phase(cl_id,gap)
%CAA_JOIN_PHASE  refetch missing phase
%
% caa_join_phase(cl_id,iso_t)
% caa_join_phase(cl_id,gap_len)
%
% $Id$

% Copyright 2006 Yuri Khotyaintsev


[iso_t,dt] = caa_read_interval;
st = iso2epoch(iso_t);

irf_log('proc',['interval: ' epoch2iso(st,1) ' -- ' epoch2iso(st+dt,1) ])

A_old = c_load('A?',cl_id,'var');
Atwo_old = c_load('Atwo?',cl_id,'var');

if isnumeric(gap)
    if isempty(Atwo_old), old_end = st;
    else old_end = Atwo_old(end,1);
    end
    st_refetch = old_end + gap;
elseif ischar(gap), st_refetch = iso2epoch(gap);
else error('invalid input')
end

if st_refetch < st || st_refetch > st+dt
    error([epoch2iso(st_refetch,1) ' is outside the interval'])
end

irf_log('proc',['refetching phase from: ' epoch2iso(st_refetch,1) ])
getData(ClusterDB,st_refetch,st+dt - st_refetch,cl_id,'a')



A = c_load('A?',cl_id,'var');
Atwo = c_load('Atwo?',cl_id,'var');

if isempty(Atwo), irf_log('load','refetch failed'), return, end

if ~isempty(A_old), A =[A_old; A]; end
if ~isempty(Atwo_old), Atwo =[Atwo_old; Atwo]; end


irf_log('save','saving updated phase')
c_eval('A?=A; Atwo?=Atwo; figure, irf_plot(Atwo?), save mA.mat A? Atwo? -append',cl_id)