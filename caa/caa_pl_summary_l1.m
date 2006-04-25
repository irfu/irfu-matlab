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
% In iso_t='-1' and dt=-1, they will be determined automatically
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

if ~exist(sdir,'dir'), error(['directory ' sdir ' does not exist']), end

savePS = 0;
savePNG = 0;
saveJPG = 0;
fullscale = 0;

int_s = realmax;
int_e = -1;

if nargin > 3, have_options = 1; args = varargin;
else have_options = 0;
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

% Save the screen size
sc_s = get(0,'ScreenSize');
if sc_s(3)==1600 && sc_s(4)==1200, scrn_size = 2;
else scrn_size = 1;
end

% Load data
r = [];
ri = [];
fmax = 12.5;
c_eval('p?=[];spec?={};es?=[];rspec?=[];in?={};')
for cli=1:4
	cdir = [sdir '/C' num2str(cli)];
	p = []; spec = {}; es = []; rspec = []; in = [];
	
	if exist(cdir, 'dir')
		d = dir([cdir '/2*_*']);
		if isempty(d), continue, end
		
		for jj=1:length(d)
			curdir = [cdir '/' d(jj).name];
			if ~exist([curdir '/.interval'],'file'), continue, end
			cd(curdir)
			
			% Load R
			if isempty(r) || ri==cli
				r_tmp = c_load('R?',cli,'var');
				if ~isempty(r_tmp) && r_tmp(1,1)~=-157e8, r = [r; r_tmp]; end
				if isempty(ri), ri = cli; end
			end
			clear r_tmp
			
			% Load P
			p_tmp = c_load('P?',cli,'var');
			if ~isempty(p_tmp) && p_tmp(1,1)~=-157e8, p = [p; p_tmp]; end
			clear p_tmp
			
			% Load spectrum
			spec_tmp = c_load('diESPEC?p1234',cli,'var');
			if ~isempty(spec_tmp) && isstruct(spec_tmp)
				spec = [spec; {spec_tmp}];
				if spec_tmp.f(end)>fmax, fmax = spec_tmp.f(end); end
			end
			clear spec_tmp
			
			% Load Es
			pp = caa_sfit_probe(cli);
			es_tmp = c_load(['diEs?p' num2str(pp)],cli,'var');
			if ~isempty(es_tmp) && es_tmp(1,1)~=-157e8
				dsiof = c_ctl(cli,'dsiof');
				if isempty(dsiof), dsiof = [1+0i 1]; end
				[ok,Ddsi] = c_load('Ddsi?',cli); if ~ok, Ddsi = dsiof(1); end
				[ok,Damp] = c_load('Damp?',cli); if ~ok, Damp = dsiof(2); end
				clear dsiof
				es_tmp = caa_corof_dsi(es_tmp,Ddsi,Damp);
				es = [es; es_tmp];
			end
			clear es_tmp
			
			% Load RSPEC
			rspec_tmp = c_load(['RSPEC?p' num2str(pp)],cli,'var');
			if ~isempty(rspec_tmp) && rspec_tmp(1,1)~=-157e8 
				rs = rspec_tmp;
				rs(:,2) = sqrt(rspec_tmp(:,2).^2+rspec_tmp(:,3).^2);
				rs(:,3) = sqrt(rspec_tmp(:,4).^2+rspec_tmp(:,5).^2);
				rs(:,4) = sqrt(rspec_tmp(:,6).^2+rspec_tmp(:,7).^2);
				rs(:,5) = sqrt(rspec_tmp(:,8).^2+rspec_tmp(:,9).^2);
				rs(:,6) = sqrt(rspec_tmp(:,10).^2+rspec_tmp(:,11).^2);
				rs(:,7:end) = [];
				rspec = [rspec; rs];
				clear rs
			end
			clear rspec_tmp
			
			% Load intervals & TM mode
			[st_s,dt1] = caa_read_interval;
			t1 = iso2epoch(st_s);
			if t1<int_s, int_s = t1; end
			if t1+dt1>int_e, int_e = t1+dt1; end
			in_tmp.interv = [t1 dt1];
			in_tmp.st_s = st_s(12:16);
			tm = c_load('mTMode?',cli,'var');
			if ~isempty(tm) && tm(1,1)~=-157e8
				if tm(1), in_tmp.tm = 1; else in_tmp.tm = 0; end 
			end
			in = [in; {in_tmp}];
			clear in_tmp
			
			cd(old_pwd)
		end
		if ~isempty(p), c_eval('p?=p;',cli), end, clear p
		if ~isempty(es), c_eval('es?=es;',cli), end, clear es
		if ~isempty(rspec), c_eval('rspec?=rspec;',cli), end, clear rspec
		if ~isempty(spec), c_eval('spec?=spec;',cli), end, clear spec
		if ~isempty(in), c_eval('in?=in;',cli), end, clear in
	end
end

if strcmp(iso_t,'-1') && dt==-1
	st = int_s;
	dt = int_e - int_s;
else st = iso2epoch(iso_t);
end
ds = irf_fname(st);
tit = ['EFW E and P 5Hz (' ds(1:4) '-' ds(5:6) '-' ds(7:8) ' ' ds(10:11) ':'...
	ds(12:13) ', produced ' date ')'];
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrum figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(75)
if scrn_size==1 ,set(gcf,'position',[91  40 909 640])
else set(gcf,'position',[7   159   790   916])
end
clf

h = 1:6;
for pl=1:6, h(pl) = irf_subplot(6,1,-pl); end

ytick =  [.25 .5 1 10];
if fullscale && fmax>100, ytick = [ytick 100]; end
for cli=1:4
	axes(h(cli))
	hold on
	c_eval('spec=spec?;',cli)
	if ~isempty(spec), for k=1:length(spec)
		caa_spectrogram(h(cli),spec{k})
	end, end
	ylabel(sprintf('Ex C%d freq [Hz]',cli))
	set(gca,'YTick',ytick,'YScale','log')
	grid
	caxis([-4 1])
	hold off
	if fullscale, set(h(cli),'YLim',[0 fmax])
	else set(h(cli),'YLim',[0 12.5])
	end
	if cli==1
		if isempty(r), title(h(1),tit)
		else title(h(1),[tit ', GSE Position C' num2str(ri)])
		end
	end
	if dt>0, irf_zoom(st +[0 dt],'x',h(cli)), end
end

% Plot P
axes(h(5))
c_pl_tx('p?')
ylabel('P L2 [-V]')
a = get(gca,'YLim');
if a(1)<-70, a(1)=-70; set(gca,'YLim',a); end
irf_zoom(st +[0 dt],'x',h(5))

if dt>0 
	plot_intervals(h(6),{in1,in2,in3,in4},st)
	irf_zoom(st +[0 dt],'x',h(6))
	if ~isempty(r), add_position(h(6),r), end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E-field figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(76)
if scrn_size==1 ,set(gcf,'position',[91  40 909 640])
else set(gcf,'position',[807   159   790   916])
end
clf

he = 1:8;
for pl=1:8,	he(pl) = irf_subplot(8,1,-pl); end

% Plot E
axes(he(1)), c_pl_tx('es?',2), ylabel('Ex [mV/m]'), axis tight
if isempty(r), title(he(1),tit)
else title(he(1),[tit ', GSE Position C' num2str(ri)])
end

axes(he(2)), c_pl_tx('es?',3), ylabel('Ey [mV/m]'), axis tight

% Plot RSPEC
c_eval('axes(he(2+?)),if exist(''rspec?'',''var'') && ~isempty(rspec?),irf_plot(rspec?),axis tight,end, ylabel(''Rspec C?''), grid on')

% Plot P
axes(he(7))
c_pl_tx('p?')
ylabel('P L2 [-V]')
a = get(gca,'YLim');
if a(1)<-70, a(1)=-70; set(gca,'YLim',a); end

if dt>0 
	plot_intervals(he(8),{in1,in2,in3,in4},st)
	irf_zoom(st +[0 dt],'x',he)
	if ~isempty(r), add_position(he(8),r), end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orient(75,'tall'), orient(76,'tall')
if fullscale,fn = sprintf('EFW_SPLOT_L1FULL__%s',irf_fname(st));
else fn = sprintf('EFW_SPLOT_L1ESPEC__%s',irf_fname(st));
end
fne = sprintf('EFW_SPLOT_L1ERSPEC__%s',irf_fname(st));
fone = sprintf('EFW_SPLOT_L1__%s',irf_fname(st));

if savePS
	irf_log('save',['saving ' fn '.[ps,pdf]'])
	irf_log('save',['saving ' fne '.[ps,pdf]'])
	print( 75, '-dpsc2', fn), print( 76, '-dpsc2', fne)
	[s,w] = unix(['/usr/local/bin/ps2pdf12 ' fn '.ps']);
	if s~=0, irf_log('save',['problem with ps2pdf: ' w]), end
	[se,w] = unix(['/usr/local/bin/ps2pdf12 ' fne '.ps']);
	if se~=0, irf_log('save',['problem with ps2pdf: ' w]), end
	if s==0 && se==0
		irf_log('save',['joining to ' fone '.pdf'])
		[s,w] = unix(['LD_LIBRARY_PATH="" /usr/local/bin/pdfjoin ' fn '.pdf ' fne '.pdf --outfile ' fone '.pdf']);
		if s~=0, irf_log('save','problem with pdfjoin'), end
	end
end
if savePNG
	irf_log('save',['saving ' fn '.png'])
	irf_log('save',['saving ' fne '.png'])
	print( 75, '-depsc2', fn), print( 76, '-depsc2', fn)
	[s,w] = unix(['/usr/local/bin/eps2png -res 150 ' fn '.eps; rm -f ' fn '.eps']);
	if s~=0, irf_log('save','problem with eps2png'), end
	[s,w] = unix(['/usr/local/bin/eps2png -res 150 ' fne '.eps; rm -f ' fne '.eps']);
	if s~=0, irf_log('save','problem with eps2png'), end
end
if saveJPG
	irf_log('save',['saving ' fn '.jpg'])
	irf_log('save',['saving ' fne '.jpg'])
	print( 75, '-depsc2', fn), print( 76, '-depsc2', fn)
	[s,w] = unix(['/usr/local/bin/eps2png -jpg -res 150 ' fn '.eps; rm -f ' fn '.eps']);
	if s~=0, irf_log('save','problem with eps2png'), end
	[s,w] = unix(['/usr/local/bin/eps2png -jpg -res 150 ' fne '.eps; rm -f ' fne '.eps']);
	if s~=0, irf_log('save','problem with eps2png'), end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Help functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function add_position(h,r)
% Add position labels
r = irf_abs(r);
add_timeaxis(h,'usefig',[r(:,1) r(:,2:end)/6371.2],...
	{'X [Re]','Y [Re]','Z [Re]','R [Re]'})

function plot_intervals(h,ints,st)
% Plot intervals

krgb = 'krgb';
cli_pos = [4 3 2 1];

axes(h)
ud = get(gcf,'userdata');
if isfield(ud,'t_start_epoch'), 
	t_start_epoch = ud.t_start_epoch;
else
	% Set start_epoch if time is in isdat epoch, 
	% warn about changing t_start_epoch
	t_start_epoch = st;
	ud.t_start_epoch = t_start_epoch; set(gcf,'userdata',ud);
	irf_log('proc',['user_data.t_start_epoch is set to '...
		epoch2iso(t_start_epoch,1)]);
end

hold(h,'on')
for cli=1:4
	in = ints{cli};
	if ~isempty(in), for k=1:length(in)
		in_tmp = in{k};
		pp = plot(in_tmp.interv(1)-t_start_epoch + [0 in_tmp.interv(2)],...
			[cli_pos(cli) cli_pos(cli)],krgb(cli));
		set(pp,'Marker','+');
		if in_tmp.tm, set(pp,'LineWidth',3); end
		text(in_tmp.interv(1)-t_start_epoch+60,cli_pos(cli)+0.2,in_tmp.st_s)
	end, end
end
hold(h,'off')
set(h,'YLim',[0 5],'YTick',1:4,'YTickLabel',4:-1:1)
ylabel(h,'proc intrerv/SC')
grid(h,'on')
