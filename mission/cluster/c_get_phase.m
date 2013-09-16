function [t,data] = c_get_phase(db,start_time,dt,cl_id,phasetype)
%C_GET_PHASE fetches phase, skipping over a potential gap
%              (as caused by a maneuver, for example)
%
% [t,data] = c_get_phase(db,start_time,dt,cl_id,phasetype)
%
% Primarily intended as a subroutine for getData. Skips at most one gap.
%  db: ISDAT database name
%  start_time: epoch start time
%  dt: length of data interval
%  cl_id: Cluster satellite number
%  phasetype: which phase variable you want ('phase' or 'phase_2')

[t,data] = caa_is_get(db, start_time, dt,cl_id, 'ephemeris', phasetype); 

if isempty(t)
	low=start_time+60;
else
	if t(end)+4-start_time > dt, return, end
	low=t(end)+60;
end
up=start_time+dt;
if up<=(low+60), return, end
irf_log('dsrc',[phasetype ' for C' num2str(cl_id) ' does not cover the full interval.'])

while up>(low+60)
	test_t=floor((up+low)/2);
	[t1,data1] = caa_is_get(db, test_t, dt-(test_t-start_time),cl_id, 'ephemeris', phasetype); 
	if isempty(t1)
		low=test_t;
	else
		up=test_t;
	end
end
if isempty(t1)
	[t1,data1] = caa_is_get(db, up, dt-(up-start_time),cl_id, 'ephemeris', phasetype); 
end

if isempty(t1), return, end
if isempty(t)
	t=t1;
	data=data1;
	return;
end

irf_log('dsrc',['Skipped over ' num2str(t1(1)-t(end)) ' second phase gap in ' ...
	phasetype ' for C' num2str(cl_id)])
t=[t' t1']';
data=[data' data1']';