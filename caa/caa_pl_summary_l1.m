function caa_pl_summary_l1(iso_t,dt,sdir)
%CAA_PL_SUMMARY_L1 CAA summary plot for L1 & L2 P data & EF
%
% caa_pl_summary_l1(iso_t,dt,sdir)
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

if ~exist(sdir,'dir'), error(['directory ' sdir ' does not exist']), end

old_pwd = pwd;
st = iso2epoch(iso_t);

% Save the screen size
sc_s = get(0,'ScreenSize');
if sc_s(3)==1600 & sc_s(4)==1200, scrn_size = 2;
else, scrn_size = 1;
end
	
figure(75)
if scrn_size==1 ,set(gcf,'position',[91  40 909 640])
else, set(gcf,'position',[691   159   909   916])
end
clf

for pl=1:6
	h(pl) = irf_subplot(6,1,-pl);
end

c_eval('p?=[];spec?=[];')
for cli=1:4
	cdir = [sdir '/C' num2str(cli)];
	p = []; spec = [];
	
	if exist(cdir, 'dir') & exist([cdir '/mINTER.mat'], 'file')
		% Load data
		c_eval(['load ' cdir '/mINTER.mat INTERVALS?'],cli)
		
		if exist(irf_ssub('INTERVALS?',cli),'var')
			c_eval('inter=INTERVALS?;',cli)

			for jj=1:size(inter,1)
				cd([cdir '/' irf_fname(inter(jj,1))])
				% Load P
				p_tmp = c_load('P?',cli,'var');
				if ~isempty(p_tmp), p = [p; p_tmp]; end
				% Load spec
				
				cd(old_pwd)
			end
			if ~isempty(p), c_eval('p?=p;',cli), end
		else
			irf_log('load', irf_ssub(['No INTERVALS? in ' cdir '/mINTER.mat'],cli))
		end
	end
end

% Plot everything
axes(h(1))
title(['EFW L1 & Co summary for ' iso_t])
axes(h(5))
c_pl_tx('p?')
ylabel('P L2 [-V]')

irf_zoom(st +[0 dt],'x',h)
