function caa_comp_efw_edi_corr(cl_id)
%CAA_COMP_EFW_EDI_CORR  compare EFW, EDI with CORROTATION
%
% CAA_COMP_EFW_EDI_CORR(CL_ID)
% 
% Compare E-fileds measured by  EFW and EDI with 
% CORROTATION E-filed.
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

getData(ClusterProc,cl_id,'edbs')
getData(ClusterProc,cl_id,'iedbs')

diEs = c_load('diEs?',cl_id,'var');
diEDI = c_load('diEDI?',cl_id,'var');
diBr = c_load('diBr?',cl_id,'var');
P = c_load('P?',cl_id,'var');
R = c_load('R?',cl_id,'var');
SAX = c_load('SAX?',cl_id,'var');
diV = c_load('diV?',cl_id,'var');

diR = c_gse2dsi(R,SAX);
diRr = irf_resamp(diR,diBr);

% Earth rotation vector in DSI
geiOM = diBr;
geiOM(:,2:3) = 0; 
geiOM(:,4) = 1;
om = irf_gse2gei(geiOM,-1);
om(:,2:4) = om(:,2:4)*2*pi/86400;
diOMr = c_gse2dsi(om,SAX);
%diOMr = irf_resamp(diOM,diBr);

% Corrotation E-filed
vCorr = irf_cross(diRr,diOMr);
vCorr(:,2:4) = -vCorr(:,2:4);
idiECorr = irf_e_vxb(vCorr,diBr);

% SC motion induced E-filed in the inertial frame
diEi = irf_tappl(irf_cross(diBr,irf_resamp(diV,diBr)),'*1e-3*(-1)');

% Corrotation E-filed in the SC frame
diECorr = idiECorr;
diECorr(:,2:4) = diECorr(:,2:4) + diEi(:,2:4);

%figure
clf
h = irf_plot({P,P,P,P});
for ax=1:3
	axes(h(ax)); cla
	irf_plot({diEDI(:,[1 (ax+1)]),diEs(:,[1 (ax+1)]),diECorr(:,[1 (ax+1)])},...
		'linestyle',{'.','-','-'},'comp');
end
legend(h(1),'EDI','EFW','corrotation')
ylabel(h(1),'Ex [mV/m]')
ylabel(h(2),'Ey [mV/m]')
ylabel(h(3),'Ez [mV/m]')
ylabel(h(4),'-ScPot [V]')
title(h(1),['Cluster ' num2str(cl_id), 'DSI SC-frame'])
irf_zoom([diBr(1,1) diBr(end,1)],'x',h)
orient tall
