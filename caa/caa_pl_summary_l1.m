function caa_pl_summary_l1(iso_t,dt,sdir,varargin)
%CAA_PL_SUMMARY_L1 CAA summary plot for L1 & L2 P data & EF
%
% caa_pl_summary_l1(iso_t,dt,sdir,[options])
%   options:
%           saveps    - save PS and PDF
%           savepng   - save JPG
%           savepng   - save PNG
%           save      - save PNG, PS and PDF
%           nosave
%           fullscale - use full scale (up to 180 Hz) on spectrograms
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

if ~exist(sdir,'dir'), error(['directory ' sdir ' does not exist']), end

savePS = 0;
savePNG = 0;
saveJPG = 0;
fullscale = 0;

if nargin > 3, have_options = 1; args = varargin;
else, have_options = 0;
end
while have_options
	l = 1;
	switch(args{1})
	case 'nosave'
		savePS = 0;
		savePNG = 0;
	case 'save'
		savePS = 1;
		savePNG = 1;
	case 'saveps'
		savePS = 1;
	case 'savepng'
		savePNG = 1;
	case 'savejpg'
		saveJPG = 1;
	case 'fullscale'
		fullscale = 1;
	otherwise
		irf_log('fcal,',['Option ''' args{1} '''not recognized'])
	end
	if length(args) > l, args = args(l+1:end);
	else break
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
set(h(6),'YTick',1:4,'YTickLabel',4:-1:1)
ylabel(h(6),'proc intrerv/SC')
krgb = 'krgb';
r = [];
ri = [];
fmax = 12.5;
cli_pos = [4 3 2 1];

c_eval('p?=[];spec?=[];')
for cli=1:4
	cdir = [sdir '/C' num2str(cli)];
	p = []; spec = [];
	
	if exist(cdir, 'dir')
		d = dir([cdir '/2*_*']);
		if isempty(d), continue, end
		
		for jj=1:length(d)
			curdir = [cdir '/' d(jj).name];
			%{
			if ~(exist([curdir '/.interval'],'file') & ...
				(exist([curdir '/mP.mat'],'file') | ...
				exist([curdir '/mEDSIf.mat'],'file'))), continue, end
			%}
			
			if ~exist([curdir '/.interval'],'file'), continue, end
			
			cd(curdir)
			% Load R
			if isempty(r) | ri==cli
				r_tmp = c_load('R?',cli,'var');
				if ~isempty(r_tmp), r = [r; r_tmp]; end
				if isempty(ri), ri = cli; end
			end
			% Load P
			p_tmp = c_load('P?',cli,'var');
			if ~isempty(p_tmp), p = [p; p_tmp]; end
			% Load spectrum
			spec = c_load('diESPEC?p1234',cli,'var');
			if ~isempty(spec)
				axes(h(cli))
				if jj>1, hold on, end
				caa_spectrogram(h(cli),spec)
				if spec.f(end)>fmax, fmax = spec.f(end); end
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
			pp = plot(t1-t_start_epoch + [0 dt1],[cli_pos(cli) cli_pos(cli)],krgb(cli));
			set(pp,'Marker','+');
			if ~isempty(tm), if tm(1), set(pp,'LineWidth',3); end, end
			text(t1-t_start_epoch+60,cli_pos(cli)+0.2,st_s)
			cd(old_pwd)
		end
		if ~isempty(p), c_eval('p?=p;',cli), end
	end
end

ytick =  [.25 .5 1 10];
if fullscale & fmax>100, ytick = [ytick 100];, end
for cli=1:4
	axes(h(cli))
	ylabel(sprintf('Ex C%d freq [Hz]',cli))
	set(gca,'YTick',ytick,'YScale','log')
	grid
	caxis([-4 1])
	hold off
end
hold(h(6),'off')
grid(h(6),'on')

% Plot the rest
axes(h(1))
ds = irf_fname(st);
tit = ['EFW E and P 5Hz (' ds(1:4) '-' ds(5:6) '-' ds(7:8) ' ' ds(10:11) ':'...
	ds(12:13) ', produced ' date ')'];
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

if fullscale, irf_zoom([0 fmax],'y',h(1:4))
else, irf_zoom([0 12.5],'y',h(1:4))
end

if ~isempty(r)
	r = irf_abs(r);
	add_timeaxis(h(6),'usefig',[r(:,1) r(:,2:end)/6371.2],...
		{'X [Re]','Y [Re]','Z [Re]','R [Re]'})
	axes(h(1)), title([tit ', GSE Position C' num2str(ri)])
end

orient tall
if fullscale,fn = sprintf('EFW_SPLOT_L1FULL__%s',irf_fname(st));
else,fn = sprintf('EFW_SPLOT_L1__%s',irf_fname(st));
end
if savePS
	irf_log('save',['saving ' fn '.[ps,pdf]'])
	print( gcf, '-dpsc2', fn) 
	[s,w] = unix(['/usr/local/bin/ps2pdf12 ' fn '.ps']);
	if s~=0, irf_log('save','problem with ps2pdf'), end
end
if savePNG
	irf_log('save',['saving ' fn '.png'])
	print( gcf, '-depsc2', fn) 
	[s,w] = unix(['/usr/local/bin/eps2png -res 150 ' fn '.eps; rm -f ' fn '.eps']);
	if s~=0, irf_log('save','problem with eps2png'), end
end
if saveJPG
	irf_log('save',['saving ' fn '.jpg'])
	print( gcf, '-depsc2', fn) 
	[s,w] = unix(['/usr/local/bin/eps2png -jpg -res 150 ' fn '.eps; rm -f ' fn '.eps']);
	if s~=0, irf_log('save','problem with eps2png'), end
end
