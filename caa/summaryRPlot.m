function summaryRPlot(cl_id,sp,st,dt,flag)
%summaryRPlot make a simple summary plot with P, E and Es DSI
% This function is useful, for example, when we do not have B PP.
% All data needs to be loaded to SP. 
%
% summaryRPlot(cl_id,sp,st,dt,flag)
% 
% Input:
% cl_id - SC#
% sp -storage path (data directory)
% st - start time, ISDAT epoch
% dt - interval length in sec
% flag - 0 (only Es, default), 1 (E and Es) [optional]
%
% $Revision$  $Date$
%
% See also ClusterProc/summaryPlot

% Copyright 2004 Yuri Khotyaintsev

if nargin<5, flag=0; end

eval(irf_ssub('load mP P?;P=P?;',cl_id))
if flag
	eval(irf_ssub('load mEDSI diE?p1234 diEs?p34 Ddsi?;E=diE?p1234;Es=diEs?p34;Dx=Ddsi?;',cl_id))
	E(:,2) = E(:,2) - Dx;
	Es(:,2) = Es(:,2) - Dx;
	n_plots = 3;
else
	eval(irf_ssub('load mEDSI diEs?p34 Ddsi?;Es=diEs?p34;Dx=Ddsi?;',cl_id))
	Es(:,2) = Es(:,2) - Dx;
	n_plots = 2;
end

subplot(n_plots,1,1)
irf_plot(P);
ylabel('SC pot [-V]')
title(['EFW, Cluster ' num2str(cl_id,'%1d')]),
set(gca,'XLim',st+[0 dt])
xlabel('')

if flag
	subplot(n_plots,1,2)
	irf_plot(E);
	ylabel('E DSI [mV/m]')
	xlabel('')
	set(gca,'XLim',st+[0 dt])
	legend('X','Y','Z')
end

subplot(n_plots,1,n_plots)
irf_plot(Es);
ylabel('E DSI [mV/m]')
set(gca,'XLim',st+[0 dt])
legend('X','Y','Z')

irf_pl_add_info
