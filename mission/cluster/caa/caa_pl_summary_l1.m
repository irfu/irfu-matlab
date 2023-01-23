function caa_pl_summary_l1(iso_t,dt,sdir,varargin)
%CAA_PL_SUMMARY_L1  CAA summary plot for L1 & L2 P data & EF
%
% caa_pl_summary_l1(iso_t,dt,sdir,[options])
%   options:
%           savepdf   - save PDF
%           saveps    - same as 'savepdf' (deprecated)
%           savepng   - save JPG
%           savepng   - save PNG
%           save      - save PNG, PS and PDF
%           nosave
%           fullscale - use full scale (up to 180 Hz) on spectrograms
%           nospec    - do not plot spectrum
%           usextra   - use and plot DdsiX offsets
%
% In iso_t='-1' and dt=-1, they will be determined automatically
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin==0, sdir = pwd; iso_t = -1; dt = -1; end

if ~exist(sdir,'dir'), error(['directory ' sdir ' does not exist']), end

savePDF = 0;
savePNG = 0;
saveJPG = 0;
fullscale = 0;
plotspec = 1;
plotbitmask = 1;
usextra = 0;
data_level = 2;
QUALITY = 3;

int_s = realmax;
int_e = -1;

if nargin > 3, have_options = 1; args = varargin;
else, have_options = 0;
end
while have_options
  l = 1;
  switch(args{1})
    case 'nosave'
      savePDF = 0;
      savePNG = 0;
    case 'save'
      savePDF = 1;
      savePNG = 1;
    case 'saveps'
      savePDF = 1;
    case 'savepdf'
      savePDF = 1;
    case 'savepng'
      savePNG = 1;
    case 'savejpg'
      saveJPG = 1;
    case 'fullscale'
      fullscale = 1;
    case 'nospec'
      plotspec = 0;
    case 'usextra'
      usextra = 1;
    otherwise
      irf_log('fcal,',['Option ''' args{1} '''not recognized'])
  end
  if length(args) > l, args = args(l+1:end);
  else, break
  end
end

old_pwd = pwd;

% Save the screen size
sc_s = get(0,'ScreenSize');
if sc_s(3)==1600 && sc_s(4)==1200, scrn_size = 2;
elseif sc_s(3)==1920 && sc_s(4)==1200, scrn_size = 2;
elseif sc_s(3)==2560 && sc_s(4)==1440, scrn_size = 2;
else, scrn_size = 1;
end

% Load data
r = [];
ri = [];
fmax = 12.5;
c_eval('p?=[];ps?=[];spec?={};es?=[];rspec?=[];in?={};wamp?=[];pswake?=[];lowake?=[];edi?=[];')

for cli=1:4
  cdir = [sdir '/C' num2str(cli)];
  in = {};
  
  if ~exist(cdir, 'dir'), continue, end
  d = dir([cdir '/2*_*']);
  if isempty(d), continue, end
  
  for jj=1:length(d)
    curdir = [cdir '/' d(jj).name];
    if ~exist([curdir '/.interval'],'file'), continue, end
    cd(curdir)
    
    % Load intervals & TM mode
    [st_s,dt1] = caa_read_interval;
    t1 = iso2epoch(st_s);
    if t1<int_s, int_s = t1; end
    if t1+dt1>int_e, int_e = t1+dt1; end
    in_tmp.interv = [t1 dt1];
    in_tmp.st_s = st_s(12:16);
    [tt,ttt,in_tmp.probeID] = caa_sfit_probe(cli); %#ok<ASGLU>
    tm = c_load('mTMode?',cli,'var');
    if ~isempty(tm) && tm(1,1)~=-157e8
      if tm(1), in_tmp.tm = 1; else, in_tmp.tm = 0; end
    else, in_tmp.tm = -1;
    end
    in = [in; {in_tmp}]; %#ok<AGROW>
    clear in_tmp
    
    cd(old_pwd)
  end
  if ~isempty(in), c_eval('in?=in;',cli), end, clear in
end

if ( strcmp(iso_t,'-1') || (isnumeric(iso_t) && iso_t==-1) ) && dt==-1
  st = int_s;
  dt = int_e - int_s;
else, st = iso2epoch(iso_t);
end

dEx = cell(4,1);
for cli=1:4
  cdir = [sdir '/C' num2str(cli)];
  p = []; ps = [];spec = {}; es = []; rspec = []; wamp = [];
  pswake = []; lowake = []; edi = [];
  
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
        if ~isempty(r_tmp) && r_tmp(1,1)~=-157e8, r = [r; r_tmp]; end %#ok<AGROW>
        if isempty(ri), ri = cli; end
      end
      clear r_tmp
      
      % Load EDI
      edi_tmp = c_load('diEDI?',cli,'var');
      if ~isempty(edi_tmp) && edi_tmp(1,1)~=-157e8
        edi = [edi; edi_tmp]; %#ok<AGROW>
      end
      clear edi_tmp
      
      % Load P
      p_tmp = c_load('P?',cli,'var');
      if ~isempty(p_tmp) && p_tmp(1,1)~=-157e8, p = [p; p_tmp]; end %#ok<AGROW>
      clear p_tmp
      
      % Load Ps
      p_tmp = c_load('Ps?',cli,'var');
      [ok,probe_info,msg] = c_load('Ps?_info',cli);
      if ~isempty(p_tmp) && p_tmp(1,1)~=-157e8 && ok && ~isempty(probe_info)
        % Extend data array to accept bitmask and quality flag (2 columns at the end)
        p_tmp = [p_tmp zeros(size(p_tmp, 1), 2)];
        p_tmp(:, end) = QUALITY;    % Default quality column to best quality, i.e. good data/no problems.
        p_quality_column = size(p_tmp, 2);
        p_bitmask_column = p_quality_column - 1;
        
        p_tmp = caa_identify_problems(p_tmp, 3, num2str(probe_info.probe), cli, p_bitmask_column, p_quality_column, 1);
        ps = [ps; p_tmp];
      else
        irf_log('load',msg)
      end
      clear p_tmp
      
      % Load SW WAKE amplitude
      for pp=[12 32 34]
        wamp_tmp = c_load(['WAKE?p' num2str(pp)],cli,'var');
        if ~isempty(wamp_tmp) && wamp_tmp(1,1)~=-157e8
          if isempty(wamp) || wamp_tmp(1,1) > wamp(end,1)
            wamp = [wamp; wamp_tmp(:,[1 3])]; %#ok<AGROW>
          else
            % Add a NaN to avoid lines going backward
            wamp = [wamp; [wamp(end,1) NaN]; wamp_tmp(:,[1 3])]; %#ok<AGROW>
          end
        end
        clear wamp_tmp
      end
      
      % Load PS/LO WAKEs
      for pp=[12 32 34 42]
        pswake_tmp = c_load(['PSWAKE?p' num2str(pp)],cli,'var');
        if ~isempty(pswake_tmp) && pswake_tmp(1,1)~=-157e8
          pswake = [pswake; pswake_tmp]; %#ok<AGROW>
        end
        clear pswake_tmp
        lowake_tmp = c_load(['LOWAKE?p' num2str(pp)],cli,'var');
        if ~isempty(lowake_tmp) && lowake_tmp(1,1)~=-157e8
          lowake = [lowake; lowake_tmp]; %#ok<AGROW>
        end
        clear lowake_tmp
      end
      
      % Load spectrum
      spec_tmp = c_load('diESPEC?p1234',cli,'var');
      if ~isstruct(spec_tmp) && (isempty(spec_tmp) || spec_tmp==-157e8)
        spec_tmp = c_load('diELXSPEC?p1234',cli,'var');
      end
      if ~isempty(spec_tmp) && isstruct(spec_tmp)
        spec = [spec; {spec_tmp}]; %#ok<AGROW>
        if spec_tmp.f(end)>fmax, fmax = spec_tmp.f(end); end
      end
      clear spec_tmp
      
      % Load Es
      spinFits = caa_sfit_load(cli);
      
      if ~isempty(spinFits)
        
        if spinFits.flagLX
          probe_numeric = spinFits.probePair;
        else
          E_info = c_load('diESPEC?p1234_info', cli, 'var');    % Load info; need list of probe pairs!
          if isempty(E_info) || ~isfield(E_info, 'probe')
            irf_log('load','Could not load probe pair info!')
            probe_numeric = spinFits.probePair;
          else
            probe_numeric=str2double(E_info.probe);
          end
        end
        % Remove saturation due to too high bias current
        if probe_numeric<50, probepair_list=probe_numeric;
        else, probepair_list=[12 32 34];end
        for probepair=probepair_list
          [ok,hbias,msg] = c_load(irf_ssub('HBIASSA?p!',cli,probepair));
          if ok
            % Special trick for C2 after April 2011
            if ~isempty(hbias) && cli==2
              hbias(hbias(:,1)==iso2epoch('2011-04-30T06:00:00.00Z'),:)=[];
            end
            if ~isempty(hbias)
              irf_log('proc','blanking HB saturation')
              spinFits.diEs = caa_rm_blankt(spinFits.diEs,hbias);
            end
          else, irf_log('load',msg)
          end
          clear ok hbias msg
        end
        
        % Remove saturation
        if probe_numeric<50, probepair_list=[mod(probe_numeric,10),fix(probe_numeric/10)];
        else, probepair_list=[1 2 3 4];end
        for probe=probepair_list
          [ok,hbias,msg] = c_load(irf_ssub('PROBESA?p!',cli,probe));
          if ok
            if ~isempty(hbias)
              irf_log('proc','blanking probe saturation')
              spinFits.diEs = caa_rm_blankt(spinFits.diEs,hbias);
            end
          else, irf_log('load',msg)
          end
          clear ok hbias msg
        end
        
        % Remove whisper pulses
        [ok,whip,msg] = c_load('WHIP?',cli);
        if ok
          if ~isempty(whip)
            irf_log('proc','blanking Whisper pulses')
            spinFits.diEs = caa_rm_blankt(spinFits.diEs,whip);
          end
        else, irf_log('load',msg)
        end
        clear ok whip msg
        
        % Remove ns_ops intervals
        ns_ops = c_ctl('get', cli, 'ns_ops');
        if isempty(ns_ops)
          c_ctl('load_ns_ops', [c_ctl('get', 5, 'data_path') '/caa-control'])
          ns_ops = c_ctl('get', cli, 'ns_ops');
        end
        if ~isempty(ns_ops)
          ns_ops_intervals = [caa_get_ns_ops_int(spinFits.diEs(1,1), spinFits.diEs(end,1)-spinFits.diEs(1,1), ns_ops, 'bad_data')' ...
            caa_get_ns_ops_int(spinFits.diEs(1,1), spinFits.diEs(end,1)-spinFits.diEs(1,1), ns_ops, 'bad_tm')'...
            caa_get_ns_ops_int(spinFits.diEs(1,1), spinFits.diEs(end,1)-spinFits.diEs(1,1), ns_ops, 'high_bias')']';
          if ~isempty(ns_ops_intervals)
            ns_ops_intervals(:,1)=ns_ops_intervals(:,1)-4;
            ns_ops_intervals(:,2)=ns_ops_intervals(:,2)+4;
            irf_log('proc', 'blanking NS_OPS')
            spinFits.diEs = caa_rm_blankt(spinFits.diEs,ns_ops_intervals);
          end
          clear ns_ops ns_ops_intervals
        end
        
        %%% Fill gaps in data???? %%%
        
        
        % Extend data array to accept bitmask and quality flag (2 columns at the end)
        spinFits.diEs = [spinFits.diEs zeros(size(spinFits.diEs, 1), 2)];
        spinFits.diEs(:, end) = QUALITY;    % Default quality column to best quality, i.e. good data/no problems.
        e_quality_column = size(spinFits.diEs, 2);
        e_bitmask_column = e_quality_column - 1;
        
        % Identify and flag problem areas in data with bitmask and quality factor:
        if probe_numeric == 120 || probe_numeric == 320 || probe_numeric == 340 || probe_numeric == 420
          probe_numeric = probe_numeric/10;
        end
        spinFits.diEs = caa_identify_problems(spinFits.diEs, data_level, sprintf('%d',probe_numeric), cli, e_bitmask_column, e_quality_column);
        
        % Delta offsets
        Del_caa = c_efw_delta_off(spinFits.diEs(1,1),cli);
        if ~isempty(Del_caa)
          [ok,Delauto] = c_load('D?p12p34',cli);
          if ~ok || isempty(Delauto)
            irf_log('load',irf_ssub('Cannot load/empty D?p12p34',cli))
          else
            spinFits.diEs = caa_corof_delta(spinFits.diEs,spinFits.probePair,Delauto,'undo');
            spinFits.diEs = caa_corof_delta(spinFits.diEs,spinFits.probePair,Del_caa,'apply');
          end
        end
        
        
        % DSI offsets
        dsiof = c_ctl(cli,'dsiof');
        if isempty(dsiof)
          [ok,Ps,msg] = c_load('Ps?',cli,'var');
          if ~ok, irf_log('load',msg), end
          if caa_is_sh_interval
            [dsiof_def, dam_def] = c_efw_dsi_off(spinFits.diEs(1,1),cli,[]);
          else
            [dsiof_def, dam_def] = c_efw_dsi_off(spinFits.diEs(1,1),cli,Ps);
          end
          clear ok Ps msg
          
          if usextra % Xtra offset
            [ok1,Ddsi] = c_load('DdsiX?',cli);
            if ~ok1
              [ok1,Ddsi] = c_load('Ddsi?',cli);
              if ~ok1, Ddsi = dsiof_def; end
            else
              iso_t = caa_read_interval;
              dEx(cli)={[dEx{cli}, {[iso2epoch(iso_t) real(Ddsi)]}]};
            end
          else
            [ok1,Ddsi] = c_load('Ddsi?',cli);
            if ~ok1, Ddsi = dsiof_def; end
          end
          [ok2,Damp] = c_load('Damp?',cli); if ~ok2, Damp = dam_def; end
          
          if ok1 || ok2, irf_log('calb',...
              ['Saved DSI offsets on C' num2str(cli)])
            %else irf_log('calb','Using default DSI offsets')
          end
          clear dsiof_def dam_def
        else
          Ddsi = dsiof(1); Damp = dsiof(2);
          irf_log('calb',['User DSI offsets on C' num2str(cl_id)])
        end
        clear dsiof
        
        spinFits.diEs = caa_corof_dsi(spinFits.diEs,Ddsi,Damp); clear Ddsi Damp
        es = [es; spinFits.diEs]; %#ok<AGROW>
        
        % Load RSPEC
        if spinFits.probePair == 120 || spinFits.probePair == 320 || spinFits.probePair == 340 || spinFits.probePair == 420
          spinFits.probePair = spinFits.probePair/10;
        end
        rspec_tmp = c_load(['RSPEC?p' num2str(spinFits.probePair)],cli,'var');
        if ~isempty(rspec_tmp) && rspec_tmp(1,1)~=-157e8
          rs = rspec_tmp;
          rs(:,2) = sqrt(rspec_tmp(:,2).^2+rspec_tmp(:,3).^2);
          rs(:,3) = sqrt(rspec_tmp(:,4).^2+rspec_tmp(:,5).^2);
          rs(:,4) = sqrt(rspec_tmp(:,6).^2+rspec_tmp(:,7).^2);
          rs(:,5) = sqrt(rspec_tmp(:,8).^2+rspec_tmp(:,9).^2);
          rs(:,6) = sqrt(rspec_tmp(:,10).^2+rspec_tmp(:,11).^2);
          rs(:,7:end) = [];
          rspec = [rspec; rs]; %#ok<AGROW>
          clear rs
        end
        clear rspec_tmp
      end
      
      cd(old_pwd)
    end
    if ~isempty(edi), c_eval('edi?=edi;',cli), end, clear edi
    if ~isempty(p), c_eval('p?=p;',cli), end, clear p
    if ~isempty(ps), c_eval('ps?=ps;',cli), end, clear ps
    if ~isempty(wamp), c_eval('wamp?=wamp;',cli), end, clear wamp
    if ~isempty(pswake)
      c_eval('pswake?=pswake;',cli)
      %			if ~isempty(es), es = caa_rm_blankt(es,pswake); end
    end
    clear pswake
    if ~isempty(lowake)
      c_eval('lowake?=lowake;',cli)
      %			if ~isempty(es), es = caa_rm_blankt(es,lowake); end
    end
    clear lowake
    
    if ~isempty(es), c_eval('es?=es;',cli), end, clear es
    if ~isempty(rspec), c_eval('rspec?=rspec;',cli), end, clear rspec
    if ~isempty(spec), c_eval('spec?=spec;',cli), end, clear spec
  end
end
ds = irf_fname(st);
tit = ['EFW E and P 5Hz (' ds(1:4) '-' ds(5:6) '-' ds(7:8) ' ' ds(10:11) ':'...
  ds(12:13) ', produced ' date ')'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrum figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotspec
  
  figure(75)
  if scrn_size==1 ,set(gcf,'position',[91  40 909 640])
  else, set(gcf,'position',[7   159   790   916])
  end
  clf
  
  h = 1:7;
  for pl=1:7, h(pl) = irf_subplot(7,1,-pl); end
  
  figure_start_epoch(st);
  
  ytick =  [.25 .5 1 10];
  if fullscale && fmax>100, ytick = [ytick 100]; end
  for cli=1:4
    hca = h(cli);
    hold(hca,'on')
    c_eval('spec=spec?;',cli)
    if ~isempty(spec)
      for k=1:length(spec), irf_spectrogram(hca,spec{k}), end
    end
    ylabel(hca,sprintf('Ex C%d freq [Hz]',cli))
    set(hca,'YTick',ytick,'YScale','log')
    grid(hca,'on')
    caxis(hca,[-4 1])
    hold(hca,'off')
    if fullscale, set(hca,'YLim',[0 fmax])
    else, set(hca,'YLim',[0 12.5])
    end
    if cli==1
      if isempty(r), title(h(1),tit)
      else, title(h(1),[tit ', GSE Position C' num2str(ri)])
      end
    end
    set(hca,'XTickLabel',[])
  end
  
  % Plot quality
  plot_quality(h(5), {{es1,ps1}, {es2,ps2}, {es3,ps3}, {es4,ps4}}, st)
  
  % Plot P
  hca = h(6);
  c_pl_tx(hca,'ps?',2)
  ylabel(hca,'P L3 [-V]'), xlabel(hca,'')
  a = get(hca,'YLim');
  if a(1)<-70, a(1)=-70; set(hca,'YLim',a); end
  
  if dt>0
    plot_intervals(h(7),{in1,in2,in3,in4},st)
    irf_zoom(h,'x',st +[0 dt])
    xlabel(h(end),'')
    if ~isempty(r), add_position(h(7),r), end
  end
  
  orient(75,'tall')
  
end % if plotspec

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% E-field figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resample EDI
if dt>0
  TAV = 180;
  ndata = ceil(dt/TAV);
  t = st + (1:ndata)*TAV - TAV/2; t = t'; %#ok<NASGU>   '
  c_eval('if ~isempty(edi?), edi?_tmp=irf_resamp(edi?,t,''fsample'',1/TAV,''thresh'',1.3); if size(edi?,1)>1,edi?=irf_tlim(edi?_tmp,[edi?(1,1) edi?(end,1)]); else edi?=edi?_tmp; end, clear edi?_tmp, end')
end
% Limit EDI
c_eval('if ~isempty(edi?) && ~isempty(es?) && any(~isnan(es?(:,2))), for c=2:3, edi?( edi?(:,c)>max(es?(~isnan(es?(:,c)),c)) & edi?(:,c)<min(es?(~isnan(es?(:,c)),c)), c) = NaN; end, end')

figure(76)
if scrn_size==1 ,set(gcf,'position',[91  40 909 640])
else, set(gcf,'position',[807   159   790   916])
end
clf

he = 1:8;
for pl=1:8,	he(pl) = irf_subplot(8,1,-pl); end

figure_start_epoch(st);

% Plot E
hca = he(1); c_pl_tx(hca,'edi?',2,'.'), hold(hca,'on')
c_pl_tx(hca,'es?',2), hold(hca,'off'), ylabel(hca,'Ex [mV/m]'), axis(hca,'tight')
if isempty(r), title(he(1),tit)
else, title(he(1),[tit ', GSE Position C' num2str(ri)])
end

hca = he(2); c_pl_tx(hca,'edi?',3,'.'), hold(hca,'on')
c_pl_tx(hca,'es?',3), hold(hca,'off'), ylabel(hca,'Ey [mV/m]'), axis(hca,'tight')

% Plot RSPEC
for cli=1:4
  c_eval('hca=he(2+?);,if ~isempty(rspec?),irf_plot(hca,rspec?), if ~isempty(lowake?),hold(hca,''on''),irf_plot(hca,caa_rm_blankt(rspec?(:,1:2),lowake?,1),''rO''),end,if ~isempty(pswake?),hold(hca,''on''),irf_plot(hca,caa_rm_blankt(rspec?(:,1:2),pswake?,1),''gd''),end,axis(hca,''tight''),end, ylabel(hca,''Rspec C?''), grid(hca,''on''), hold(hca,''off'')',cli);
  c_eval('in = in?;',cli);
  if isempty(in), continue, end
  for k=1:length(in)
    in_tmp = in{k};
    yrange = get(hca,'YLim');
    text(in_tmp.interv(1)-figure_start_epoch(st)+60,yrange(2)*0.8,['p' in_tmp.probeID],'parent',hca)
  end
  
end

t_start_epoch = figure_start_epoch(st);
for cli=1:4
  if ~isempty(dEx{cli})
    hca = he(2+cli);
    yy=get(hca,'YLim');
    yy=yy(1)+0.7*(yy(2)-yy(1));
    for in=1:length(dEx{cli})
      text(dEx{cli}{in}(1) - t_start_epoch, yy, ...
        sprintf('%.2f',dEx{cli}{in}(2)),'color','g','parent',hca)
    end
  end
end

% Plot P
hca = he(7);
if ~isempty(wamp1) || ~isempty(wamp2)|| ~isempty(wamp3)|| ~isempty(wamp4)
  c_pl_tx(hca,'wamp?')
  ylabel(hca,'Wake [mV/m]')
else
  c_pl_tx(hca,'p?')
  ylabel(hca,'P L2 [-V]')
  a = get(hca,'YLim');
  if a(1)<-70, a(1)=-70; set(hca,'YLim',a); end
end

if dt>0
  plot_intervals(he(8),{in1,in2,in3,in4},st)
  irf_zoom(he,'x',st +[0 dt])
  xlabel(he(end),'')
  if ~isempty(r), add_position(he(8),r), end
end

orient(76,'tall')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quality bitmask figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotbitmask
  
  figure(77)
  if scrn_size==1 ,set(gcf,'position',[91  40 909 640])
  else, set(gcf,'position',[7   159   790   916])
  end
  clf
  
  h = 1:6;
  for pl=1:6, h(pl) = irf_subplot(6,1,-pl); end
  
  t_start_epoch = figure_start_epoch(st);
  
  % Plot bits in bitmask
  ytick_label =  1:16;
  ytick_offset = fix( (ytick_label-1) / 4 );
  ytick = ytick_label + ytick_offset;
  
  for cli=1:4
    hca = h(cli);
    hold(hca,'on')
    c_eval('es=es?;',cli)
    if ~isempty(es)
      for k = 1:16
        plot(hca,[es(1,1) es(end,1)] - t_start_epoch, [k k] + fix((k-1)/4), 'k')
        index = find( bitget(es(:,e_bitmask_column), k) );
        if ~isempty(index)
          ind_d = find( diff(index) > 1 );
          if ~isempty(ind_d)
            ind_start = 1;
            for ii = 1:length(ind_d)
              ind_stop = ind_d(ii);
              plot(hca,[es(index(ind_start),1) es(index(ind_stop),1)] - t_start_epoch, ...
                [ytick(k) ytick(k)], 'r', 'LineWidth', 2);
              ind_start = ind_d(ii)+1;
            end
            plot(hca,[es(index(ind_start),1) es(index(end),1)] - t_start_epoch, ...
              [ytick(k) ytick(k)], 'r', 'LineWidth', 2);
          else
            plot(hca,[es(index(1),1) es(index(end),1)] - t_start_epoch, ...
              [ytick(k) ytick(k)], 'r', 'LineWidth', 2);
          end
        end
      end
    end
    if cli==1
      if isempty(r), title(h(1),tit)
      else, title(h(1),[tit ', GSE Position C' num2str(ri)])
      end
    end
    ylabel(hca,sprintf('Bitmask C%d [bit no.]',cli))
    set(hca,'YTick',ytick, 'YTickLabel', ytick_label, 'FontSize', 6)
    grid(hca,'on')
    hold(hca,'off')
    
  end
  
  % Plot quality
  plot_quality(h(5), {{es1,ps1}, {es2,ps2}, {es3,ps3}, {es4,ps4}}, st)
  
  if dt>0
    plot_intervals(h(6),{in1,in2,in3,in4},st)
    irf_zoom(h,'x',st +[0 dt])
    xlabel(h(end),'')
    if ~isempty(r), add_position(h(6),r), end
  end
  
  orient(77, 'tall')
end % if plotbitmask

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fullscale,fn = sprintf('EFW_SPLOT_L1FULL__%s',irf_fname(st));
else, fn = sprintf('EFW_SPLOT_L1ESPEC__%s',irf_fname(st));
end
fne = sprintf('EFW_SPLOT_L1ERSPEC__%s',irf_fname(st));
fnq = sprintf('EFW_SPLOT_L1QUAL__%s',irf_fname(st));
fone = sprintf('EFW_SPLOT_L1__%s',irf_fname(st));

if savePDF
  irf_log('save',['saving ' fn '.pdf'])
  irf_log('save',['saving ' fne '.pdf'])
  irf_log('save',['saving ' fnq '.pdf'])
  print( 75, '-dpdf', fn), print( 76, '-dpdf', fne), print( 77, '-dpdf', fnq)
  if exist('/usr/local/bin/pdfjoin','file')
    irf_log('save',['joining to ' fone '.pdf'])
    s = unix(['LD_LIBRARY_PATH="" /usr/local/bin/pdfjoin ' fn '.pdf ' fne '.pdf ' fnq '.pdf --outfile ' fone '.pdf']);
    if s~=0, irf_log('save','problem with pdfjoin'), end
  else
    irf_log('proc',...
      'cannot join PDFs: /usr/local/bin/pdfjoin does not exist')
  end
end
if savePNG
  if exist('/usr/local/bin/eps2png','file')
    irf_log('save',['saving ' fn '.png'])
    irf_log('save',['saving ' fne '.png'])
    irf_log('save',['saving ' fnq '.png'])
    print( 75, '-depsc2', fn), print( 76, '-depsc2', fn), print( 77, '-depsc2', fnq)
    s = unix(['/usr/local/bin/eps2png -res 150 ' fn '.eps; rm -f ' fn '.eps']);
    if s~=0, irf_log('save','problem with eps2png'), end
    s = unix(['/usr/local/bin/eps2png -res 150 ' fne '.eps; rm -f ' fne '.eps']);
    if s~=0, irf_log('save','problem with eps2png'), end
    s = unix(['/usr/local/bin/eps2png -res 150 ' fnq '.eps; rm -f ' fnq '.eps']);
    if s~=0, irf_log('save','problem with eps2png'), end
  else
    irf_log('proc',...
      'cannot save JPG: /usr/local/bin/eps2png does not exist')
  end
end
if saveJPG
  if exist('/usr/local/bin/eps2png','file')
    irf_log('save',['saving ' fn '.jpg'])
    irf_log('save',['saving ' fne '.jpg'])
    irf_log('save',['saving ' fnq '.jpg'])
    print( 75, '-depsc2', fn), print( 76, '-depsc2', fn), print( 77, '-depsc2', fnq)
    s = unix(['/usr/local/bin/eps2png -jpg -res 150 ' fn '.eps; rm -f ' fn '.eps']);
    if s~=0, irf_log('save','problem with eps2png'), end
    s = unix(['/usr/local/bin/eps2png -jpg -res 150 ' fne '.eps; rm -f ' fne '.eps']);
    if s~=0, irf_log('save','problem with eps2png'), end
    s = unix(['/usr/local/bin/eps2png -jpg -res 150 ' fnq '.eps; rm -f ' fnq '.eps']);
    if s~=0, irf_log('save','problem with eps2png'), end
  else
    irf_log('proc',...
      'cannot save JPG: /usr/local/bin/eps2png does not exist')
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Help functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_intervals(h,ints,st)
% Plot intervals

krgb = 'krgb';
cli_pos = [4 3 2 1];

t_start_epoch = figure_start_epoch(st);

hold(h,'on')
for cli=1:4
  in = ints{cli};
  if isempty(in), continue, end
  for k=1:length(in)
    in_tmp = in{k};
    pp = plot(h,in_tmp.interv(1)-t_start_epoch + [0 in_tmp.interv(2)],...
      [cli_pos(cli) cli_pos(cli)],krgb(cli));
    set(pp,'Marker','+');
    if in_tmp.tm==1, set(pp,'LineWidth',3)
    elseif in_tmp.tm==-1, set(pp,'LineStyle','--')
    end
    text(in_tmp.interv(1)-t_start_epoch+60,cli_pos(cli)+0.2,in_tmp.st_s,...
      'parent',h)
  end
end
hold(h,'off')
set(h,'YLim',[0 5],'YTick',1:4,'YTickLabel',4:-1:1)
ylabel(h,'proc intrerv/SC')
grid(h,'on')
end


function plot_quality(h, dataset, st)
% Plot quality factor

linecolor = 'mrgb';
cli_pos = [4 3 2 1];

t_start_epoch = figure_start_epoch(st);
t_end_epoch=t_start_epoch+10800;

hold(h,'on')
for cli=1:4
  % Start by plotting nsops intervals.
  ns_ops = c_ctl('get',cli,'ns_ops');
  if isempty(ns_ops)
    c_ctl('load_ns_ops',[cdb.dp '/caa-control'])
    ns_ops = c_ctl('get',cl_id,'ns_ops');
  end
  if isempty(ns_ops),error('Nonstandard operations table not found!'),end
  ii = find( ns_ops(:,1)<=t_end_epoch & ns_ops(:,1)+ns_ops(:,2)>t_start_epoch);
  for j=1:length(ii)
    plot(h,[ns_ops(ii(j),1) ns_ops(ii(j),1)+ns_ops(ii(j),2)] - t_start_epoch, ...
      [cli_pos(cli) cli_pos(cli)]-.2, 'm','LineWidth', 3.5);
  end
  % Then plot quality lines.
  data_c = dataset{cli};
  if isempty(data_c), continue, end
  for iVar=1:length(data_c)
    data = data_c{iVar};
    if isempty(data), continue, end
    quality = data(:, end);
    indexes = [];
    start_ind = 1;
    finished = 0;
    
    while ~finished
      if isnan(quality(start_ind))
        next_ind = find(~isnan(quality(start_ind:end)), 1);
      else
        next_ind = find(quality(start_ind:end) ~= quality(start_ind), 1);
      end
      if isempty(next_ind)
        next_ind = length(quality(start_ind:end)) + 1;
        finished = 1;
      end
      
      end_ind = start_ind + next_ind - 2;
      indexes = [indexes; [start_ind end_ind] quality(start_ind)]; %#ok<AGROW>
      
      if ~isnan(quality(start_ind))
        plot(h,[data(start_ind, 1)-2 data(end_ind, 1)+2] - t_start_epoch, ...
          [cli_pos(cli) cli_pos(cli)]+0.2*(iVar-1), linecolor(quality(start_ind)+1), ...
          'LineWidth', 3-quality(start_ind)+0.5);
      end
      
      start_ind = start_ind + next_ind - 1;
    end
  end
end
hold(h,'off')
set(h,'YLim',[0 5],'YTick',1:4,'YTickLabel',4:-1:1)
ylabel(h,'Quality/SC')
grid(h,'on')
end


function t_start_epoch = figure_start_epoch(st)
ud = get(gcf,'userdata');
if isfield(ud,'t_start_epoch')
  t_start_epoch = ud.t_start_epoch;
else
  t_start_epoch = st;
  ud.t_start_epoch = t_start_epoch;
  set(gcf,'userdata',ud);
  irf_log('proc',['user_data.t_start_epoch is set to '...
    epoch2iso(t_start_epoch,1)]);
end
end