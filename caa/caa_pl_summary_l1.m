function caa_pl_summary_l1(iso_t,dt,sdir,options)
%CAA_PL_SUMMARY_L1 CAA summary plot for L1 & L2 P data & EF
%
% caa_pl_summary_l1(iso_t,dt,sdir,[options])
%   options:
%           saveps  - save PS and PDF
%           savepng - save PNG
%           save    - save PNG, PS and PDF
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

if ~exist(sdir,'dir'), error(['directory ' sdir ' does not exist']), end

savePS = 0;
savePNG = 0;
if nargin>3
	if strcmp(options,'save')
		savePS = 1;
		savePNG = 1;
	elseif strcmp(options,'saveps')
		savePS = 1;
	elseif strcmp(options,'savepng')
		savePNG = 1;
	else
		irf_log('fcal','unknown option')
	end
end

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
irf_zoom([0 5],'y',h(6))
hold(h(6),'on')
set(h(6),'YTick',1:4)
ylabel(h(6),'SC')
krgb = 'krgb';
r = [];
ri = [];

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
				% Load R
				if isempty(r) | ri==cli
					r_tmp = c_load('R?',cli,'var');
					if ~isempty(r_tmp), r = [r; r_tmp]; end
					if isempty(ri), ri = cli; end
				end
				% Load P
				p_tmp = c_load('P?',cli,'var');
				if ~isempty(p_tmp)
					% Remove Sweeps from plot, not from data
					[ok,sweep] = c_load('SWEEP?',cli);
					if ok
						if ~isempty(sweep)
							irf_log('proc','blanking sweeps')
							p_tmp = caa_rm_blankt(p_tmp,sweep);
							clear sweep
						end
					else
						irf_log('load',...
							irf_ssub(['No SWEEP?. Use getData(CP,cl_id,''sweep'')'],cl_id))
					end
					p = [p; p_tmp]; 
				end
				% Load spectrum
				spec = c_load('diESPEC?p1234',cli,'var');
				if ~isempty(spec)
					axes(h(cli))
					if jj>1, hold on, end
					caa_spectrogram(h(cli),spec)
					if jj==size(inter,1)
						hold off
						caxis([-6 1])
					end
				end
				% Load intervals & TM mode
				[st_s,dt1] = caa_read_interval;
				t1 = iso2epoch(st_s);
				st_s = st_s([12:19]);
				tm = c_load('mTMode?',cli,'var');
				axes(h(6))
				ud = get(gcf,'userdata');
				if isfield(ud,'t_start_epoch'), 
					t_start_epoch = ud.t_start_epoch;
				else
					% Set start_epoch if time is in isdat epoch, warn about changing t_start_epoch
					t_start_epoch = t1;
					ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
					irf_log('proc',['user_data.t_start_epoch is set to ' epoch2iso(t_start_epoch,1)]);
				end
				pp = plot(t1-t_start_epoch + [0 dt1],[cli cli],krgb(cli));
				set(pp,'Marker','+');
				if ~isempty(tm), if tm(1), set(pp,'LineWidth',3); end, end
				text(t1-t_start_epoch+60,cli+0.2,st_s)
				cd(old_pwd)
			end
			if ~isempty(p), c_eval('p?=p;',cli), end
		else
			irf_log('load', irf_ssub(['No INTERVALS? in ' cdir '/mINTER.mat'],cli))
		end
	end
end

for cli=1:4
	axes(h(cli))
	ylabel(sprintf('Ex C%d freq [Hz]',cli))
	set(gca,'YTick',[.25 .5 1 10],'YScale','log')
	grid
end
hold(h(6),'off')
grid(h(6),'on')

% Plot the rest
axes(h(1))
ds = irf_fname(st);
tit = ['EFW EF and P L2 (' ds(1:4) '-' ds(5:6) '-' ds(7:8) ' ' ds(10:11) ':' ds(12:13) ')'];
title(tit)
axes(h(5))
c_pl_tx('p?')
ylabel('P L2 [-V]')
a = get(gca,'YLim');
if a(1)<-70
	a(1)=-70;
	set(gca,'YLim',a);
end

irf_zoom(st +[0 dt],'x',h)
irf_zoom([0 12.5],'y',h(1:4))

if ~isempty(r)
	r = irf_abs(r);
	add_timeaxis(h(6),'usefig',[r(:,1) r(:,2:end)/6300],...
		{'X [Re]','Y [Re]','Z [Re]','R [Re]'})
	axes(h(1)), title([tit ', Position C' num2str(ri)])
end

orient tall
fn = sprintf('EFW_SPLOT_L1__%s',irf_fname(st));
if savePS
	irf_log('save',['saving ' fn '.[ps,pdf]'])
	print( gcf, '-dpsc2', fn) 
	[s,w] = unix(['/usr/local/bin/ps2pdf12 ' fn '.ps']);
	if s~=0, irf_log('save','problem with ps2pdf'), end
end
if savePNG
	irf_log('save',['saving ' fn '.png'])
	print( gcf, '-depsc2', fn) 
	[s,w] = unix(['/usr/local/bin/eps2png ' fn '.eps; rm -f ' fn '.eps']);
	if s~=0, irf_log('save','problem with eps2png'), end
end
