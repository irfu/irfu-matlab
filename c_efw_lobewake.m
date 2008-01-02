function res = c_efw_lobewake(cl_id,diE,diB,Ps)
%C_EFW_LOBEWAKE  detect lobe wakes
%
% wake = c_efw_lobewake(cl_id,[diEs,diBrs,Ps])
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% Control parameters
TAV = 60; % window in sec
ANG_LIM = 15;
EPAR_LIM = 1;          % limit on E|| in mV/m
EPERP_LIM = 1.5;       % limit on full Eperp (E*B=0) in mV/m
EDEV_LIM = .8;         % limit on s-dev of E
EZ_LIM = 1.0;          % limit on Ez in mV/m
EPAR_EPERP_RATIO_LIM = .2;
EPERP_RATIO_LIM = 1.8; % Ratio between full Eperp (E*B=0) and 
                       % E along the projecttion of B to spin plane
SCPOT_LIM = -5; % Limit for spacecraft potential

DOPLOT = 0;

% Load data
probe_p = caa_sfit_probe(cl_id);
if nargin==1
	[ok,diE] = c_load(sprintf('diEs%dp%d',cl_id,probe_p));
	if ~ok, error('cannot load E'), end

	[ok,diB] = c_load('diBrs?',cl_id);
	if ~ok, error('cannot load B'), end
	[ok,Ps] = c_load('Ps?',cl_id);
	if ~ok, error('cannot load P'), end
	[ok,diEDI] = c_load('diEDI?',cl_id);
	if ~ok, disp('cannot load EDI'), end
end

ndata = ceil((diE(end,1) - diE(1,1))/TAV);
t = diE(1,1) + (1:ndata)*TAV - TAV/2; t = t';

diEr = irf_resamp(diE,t,'fsample',1/TAV);
Psr = irf_resamp(Ps(~isnan(Ps(:,2)),:),t,'fsample',1/TAV);
diBr = irf_resamp(diB,t,'fsample',1/TAV);
bele = (180/pi)*asin(...
	diBr(:,4)./sqrt(diBr(:,2).^2 + diBr(:,3).^2 + diBr(:,4).^2) );
if DOPLOT
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
	h(5) = irf_subplot(nplots,1,-5);
	irf_plot([diEr(:,1) bele]);
	ylabel('\theta [deg]')
end

% Look for strong apparent Epar if B field is within 15 deg of spin plane:
ind = find( abs(bele) < ANG_LIM & Psr(:,2) < SCPOT_LIM);
Epar = ( diEr(:,2).*diBr(:,2) + diEr(:,3).*diBr(:,3) )...
	./sqrt( diBr(:,2).^2 + diBr(:,3).^2 );
Eperp = abs( diEr(:,2).*diBr(:,3) - diEr(:,3).*diBr(:,2) )...
	./sqrt( diBr(:,2).^2 + diBr(:,3).^2 );
wind = ind(abs(Epar(ind)) > EPAR_LIM & abs(Eperp(ind)) > EPERP_LIM &...
	abs(Epar(ind)./Eperp(ind)) > EPAR_EPERP_RATIO_LIM);

if DOPLOT && ~isempty(wind)
	axes(h(2)), hold on
	irf_plot(diEr(wind,[1 2]),'g.');
	axes(h(3)), hold on
	irf_plot(diEr(wind,[1 3]),'g.');
end

% Look for strong apparent Ez if B field is above 15 deg of spin plane:
ind = find( abs(bele) >= ANG_LIM & Psr(:,2) < SCPOT_LIM);
if ~isempty(ind)
	diEr(:,4) = -(diEr(:,2).*diBr(:,2)+diEr(:,3).*diBr(:,3))./diBr(:,4);
	diEr = irf_abs(diEr); % diEr(:,6) now contains abs(E)
	ind = ind(abs(diEr(ind,4)) > EZ_LIM & diEr(ind,6) > EPERP_LIM &...
		diEr(ind,6)./Eperp(ind) > EPERP_RATIO_LIM & diEr(ind,5) < EDEV_LIM);
end

if DOPLOT && ~isempty(ind)
	axes(h(2)), hold on
	irf_plot(diEr(ind,[1 2]),'r.');
	axes(h(3)), hold on
	irf_plot(diEr(ind,[1 3]),'r.');
end

wind = sort([wind; ind]);
clear ind

if DOPLOT
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
elseif ~isempty(wind)
	DT2 = TAV;
	res = t(wind(1)) + [-DT2 DT2];
	wind(1) = [];
	if ~isempty(wind)
		for j=wind'
			if res(end,2)>=t(j)-DT2*2, res(end,2) = t(j) + DT2; % throw away one good point in between bad (DT2*2)
			else res = [res; t(j)-DT2 t(j)+DT2];
			end
		end
	end
else
	irf_log('proc','no lobe wakes')
	res = [];
end
