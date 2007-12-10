function c_efw_lobewake(cl_id)
%C_EFW_LOBEWAKE  detect lobe wakes
%
% c_efw_lobewake(cl_id)
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

pp = caa_sfit_probe(cl_id);
[ok,diE] = c_load(sprintf('diEs%dp%d',cl_id,pp));
if ~ok, error('cannot load E'), end
[Ddsi,Damp] = c_efw_dsi_off(diE(1,1),cl_id);
diE = caa_corof_dsi(diE,Ddsi,Damp);
 
[ok,diB] = c_load('diBPP?',cl_id);
if ~ok, error('cannot load B'), end
[ok,Ps] = c_load('Ps?',cl_id);
[ok,diEDI] = c_load('diEDI?',cl_id);

nplots = 6;
h = 1:nplots;
h(1) = irf_subplot(nplots,1,-1);
irf_plot(Ps);
ylabel('Sc pot [-V]')
st_iso = epoch2iso(diE(end,1));
title(['Cluster ' num2str(cl_id) ', ' st_iso(1:10)])
h(2) = irf_subplot(nplots,1,-2);
irf_plot(diE(:,[1 2]));
set(gca,'YLim',[-9 9])
ylabel('Ex [mV/m]')
h(3) = irf_subplot(nplots,1,-3);
irf_plot(diE(:,[1 3]));
set(gca,'YLim',[-9 9])
ylabel('Ey [mV/m]')
h(4) = irf_subplot(nplots,1,-4);

TAV = 60; % 1 minute window

ndata = ceil((diE(end,1) - diE(1,1))/TAV);
t = diE(1,1) + (1:ndata)*TAV - TAV/2; t = t';

diEr = irf_resamp(diE,t,'fsample',1/TAV);
if ok
	diEDI = diEDI(~isnan(diEDI(:,2)),:);
	diEDIr = irf_resamp(diEDI,t,'fsample',1/TAV,'thresh',1.3);
	axes(h(2)), hold on
	irf_plot(diEDIr(:,[1 2]),'k.');
	axes(h(3)), hold on
	irf_plot(diEDIr(:,[1 3]),'k.');
	diEDIr = diEDIr(:,1:3);
	diEDIr(:,2:3) = diEr(:,2:3) - diEDIr(:,2:3);
	axes(h(4))
	irf_plot(diEDIr)
	ylabel('E diff [mV/m]')
end

diBr = irf_resamp(diB,t,'fsample',1/TAV);
bele = (180/pi)*asin(...
	diBr(:,4)./sqrt(diBr(:,2).^2 + diBr(:,3).^2 + diBr(:,4).^2) );
h(5) = irf_subplot(nplots,1,-5);
irf_plot([diEr(:,1) bele]);
ylabel('\theta [deg]')

ANG_LIM = 15;
EPAR_LIM = 1; % 1 mV/m
EZ_LIM = 2;   % 2 mV/m
EPAR_EPERP_RATIO_LIM = .2;
EZ_EPERP_RATIO_LIM = 2;

% Look for strong apparent Epar if B field is within 15 deg of spin plane:
ind = find( abs(bele) < ANG_LIM );
Epar = ( diEr(:,2).*diBr(:,2) + diEr(:,3).*diBr(:,3) )...
	./sqrt( diBr(:,2).^2 + diBr(:,3).^2 );
Eperp = abs( diEr(:,2).*diBr(:,3) - diEr(:,3).*diBr(:,2) )...
	./sqrt( diBr(:,2).^2 + diBr(:,3).^2 );
wind = ind(abs(Epar(ind)) > EPAR_LIM &...
	abs(Epar(ind)./Eperp(ind)) > EPAR_EPERP_RATIO_LIM);

if ~isempty(wind)
	axes(h(2)), hold on
	irf_plot(diEr(wind,[1 2]),'g.');
	axes(h(3)), hold on
	irf_plot(diEr(wind,[1 3]),'g.');
end

% Look for strong apparent Ez if B field is above 15 deg of spin plane:
ind = find( abs(bele) >= ANG_LIM );
Ez = -(diEr(:,2).*diBr(:,2)+diEr(:,3).*diBr(:,3))./diBr(:,4);
ind = ind(abs(Ez(ind)) > EZ_LIM &...
	abs(Ez(ind))./Eperp(ind) > EZ_EPERP_RATIO_LIM);

if ~isempty(ind)
	axes(h(2)), hold on
	irf_plot(diEr(ind,[1 2]),'r.');
	axes(h(3)), hold on
	irf_plot(diEr(ind,[1 3]),'r.');
end

wind = sort([wind; ind]);

for j=wind'
	%diE( diE(:,1)>=t(j)-TAV/2 & diE(:,1)<=t(j)+TAV/2, 2:end) = NaN;
	diE( diE(:,1)>=t(j)-TAV & diE(:,1)<=t(j)+TAV, 2:end) = NaN;
end

h(6) = irf_subplot(nplots,1,-6);
irf_plot(diE(:,1:3))
ylabel('E cor [mV/m]')

irf_zoom([t(1)-TAV/2 t(end)+TAV/2],'x',h)

[ok,r] = c_load('R?',cl_id);
if ok, add_position(h(6),r), end
