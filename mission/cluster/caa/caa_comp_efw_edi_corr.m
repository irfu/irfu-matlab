function caa_comp_efw_edi_corr(cl_id)
%CAA_COMP_EFW_EDI_CORR  compare EFW, EDI with CORROTATION
%
% CAA_COMP_EFW_EDI_CORR(CL_ID)
%
% Compare E-fileds measured by  EFW and EDI with
% CO-ROTATION E-filed.
%
% See also: C_EFW_CORROT
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


% Control parameters
DE_LIM = 1;       % Limit on deviation from corrotation in mV/m
SCPOT_LIM = -1.5; % Limit for spacecraft potential
TAV = 60;         % Step in sec
WIN = 10;         % Averaging window WIN*TAV

getData(ClusterProc,cl_id,'edbs')
getData(ClusterProc,cl_id,'iedbs')

diEs = c_load('diEs?',cl_id,'var');
if diEs(1,1) == -157e8, error('No E-field'), end
diEs(isnan(diEs(:,2)),:) = [];
diEDI = c_load('diEDI?',cl_id,'var');
if diEDI(1,1) == -157e8
  diEDI = [];
  disp('No EDI data')
end

diBr = c_load('diBr?',cl_id,'var');
if diBr(1,1) == -157e8
  diBr = c_load('diBrs?',cl_id,'var');
  if diBr(1,1) == -157e8, error('No B-field'), end
end
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
h = 1:4;
for ax=1:3
  h(ax) = irf_subplot(4,1,-ax);
  if isempty(diEDI)
    irf_plot({diEs(:,[1 (ax+1)]),diECorr(:,[1 (ax+1)])},...
      'linestyle',{'-','-'},'comp');
  else
    irf_plot({diEDI(:,[1 (ax+1)]),diEs(:,[1 (ax+1)]),diECorr(:,[1 (ax+1)])},...
      'linestyle',{'.','-','-'},'comp');
  end
  sfit_probe = caa_sfit_probe(cl_id);
  [ok, lowake] = c_load(sprintf('LOWAKE?p%d',sfit_probe),cl_id);
  if ok && ~isempty(lowake)
    diEs_wake = caa_rm_blankt(diEs,lowake,1);
    hold on
    irf_plot(diEs_wake(:,[1 (ax+1)]),'g*')
  end
end
if isempty(diEDI), legend(h(1),'EFW','corrotation')
else, legend(h(1),'EDI','EFW','corrotation')
end
ylabel(h(1),'Ex [mV/m]')
ylabel(h(2),'Ey [mV/m]')
ylabel(h(3),'Ez [mV/m]')
h(4) = irf_subplot(4,1,-4);
irf_plot(P)
ylabel(h(4),'-ScPot [V]')
title(h(1),['Cluster ' num2str(cl_id), 'DSI SC-frame'])
irf_zoom([diBr(1,1) diBr(end,1)],'x',h)
orient tall

ndata = ceil((diEs(end,1) - diEs(1,1))/TAV);
t = diEs(1,1) + (1:ndata)*TAV - TAV/2; t = t';

diEr = irf_resamp(diEs,t,'fsample',.1/TAV);
diECr = irf_resamp(diECorr,t,'fsample',.1/TAV);
Pr = irf_resamp(P(~isnan(P(:,2)),:),t,'fsample',1/TAV);

dE = sqrt( (diEr(:,2) - diECr(:,2)).^2 + (diEr(:,3) - diECr(:,3)).^2 );
idx = find( Pr(:,2) > SCPOT_LIM & dE > DE_LIM );
clear diEr diECr Pr

if ~isempty(idx)
  irf_log('proc',sprintf('max diff from corrotation is %.2f mV/m',max(dE)))
  for j=idx'
    diEs( diEs(:,1)>=t(j)-WIN*TAV/2 & diEs(:,1)<=t(j)+WIN*TAV/2,...
      2:end ) = NaN;
  end
  for ax=1:2
    axes(h(ax)); hold on
    irf_plot(diEs(:,[1 (ax+1)]),'g.');
    hold off
  end
  set(gca,'XTickLabel',[])
  xlabel('')
else
  irf_log('proc','data is good')
end

