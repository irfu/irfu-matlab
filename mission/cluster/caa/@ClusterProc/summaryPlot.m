function out=summaryPlot(cp,cl_id,varargin)
% summaryPlot make EFW summary plot
%
% h = summaryPlot(cp,cl_id,[options])
%
% Input:
%   cp - ClusterProc object
%   cl_id - SC#
%   Options: go in pair 'option', value
%    'cs'        - coordinate system : 'dsi' [default] of 'gse'
%    'st', 'dt'  - start time (ISDAT epoch) and interval length (sec)
%    'fullb'     - use full resolution B FGM
%    'withwhip' - plot time intervals with Whisper pulses
%    'ib'        - Internal Burst format
%    'wo_r'      - Do not add position labels
%    'vars'      - parameters for title
%
% Output:
%   h - axes handles // can be omitted
%
% Example:
%   summaryPlot(ClusterProc('/home/yuri/caa-data/20020304'),1,'cs','gse')
%   summaryPlot(ClusterProc('.'),2,'st',toepoch([2004 1 4 12 47 0]),'dt',60,'fullb')
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(2,10)

if nargin>2, have_options = 1; args = varargin;
else, have_options = 0;
end

% Default values
cs = 'dsi';
st = 0;
dt = 0;
have_tint = 0;
use_fullb = 'rs';
flag_rmwhip = 1;
flag_ib = 0;
flag_r = 1;
flag_ibs = 0;

while have_options
  l = 2;
  if length(args)>=1
    switch(args{1})
      case 'cs'
        if ischar(args{2})
          cs = args{2};
          if ~strcmp(cs,'dsi') && ~strcmp(cs,'gse')
            irf_log('fcal','unknown CS. defaulting to DSI')
            cs = 'dsi';
          end
        else, irf_log('fcal','wrongArgType : CS must be a string')
        end
      case 'st'
        if isnumeric(args{2}), st = args{2};
        else, irf_log('fcal','wrongArgType : ST must be numeric')
        end
      case 'dt'
        if isnumeric(args{2}), dt = args{2};
        else, irf_log('fcal','wrongArgType : DT must be numeric')
        end
      case 'fullb'
        use_fullb = '';	l = 1;
      case {'leavewhip','withwhip'}
        flag_rmwhip = 0; l = 1;
      case 'ib'
        flag_ib = 1; l = 1;
      case 'wo_r'
        flag_r = 0; l = 1;
      case 'vars'
        vars=args{2};
        l=size(vars,1);
        vars_names=[];
        if l>0
          flag_ibs = 1;
          vars_names = strtrim(vars(1,:));
          for i=2:l
            vars_names = [vars_names ', ' strtrim(vars(i,:))];
          end
        end
      otherwise
        irf_log('fcal',['Option ''' args{1} '''not recognized'])
    end
    if length(args) > l, args = args(l+1:end);
    else, break
    end
  else
    error('caa:wrongArgType','use summaryPlot(..,''option'',''value'')')
  end
end

if flag_ib && ~strcmp(cs,'dsi')
  irf_log('func','IB can only be used in DSI. Defaulting to DSI')
  cs = 'dsi';
end

if st && dt, have_tint = 1; end

% Define variables we want to plot
if flag_ib
  q_list = {'P?',['diB' use_fullb '?'],'diE?p1234','dibE?p1234','diBSC4kHz?','diVExBs?'};
  l_list = {'SC pot [-V]','B DSI [nT]','E DSI [mV/m]','E DSI [mV/m]','B DSI [nT]','V=ExB DSI [km/s]'};
else
  if strcmp(cs,'dsi')
    q_list = {'P?',['diB' use_fullb '?'],'diE?p1234','diEs?','diVExBs?'};
    l_list = {'SC pot [-V]','B DSI [nT]','E DSI [mV/m]','E DSI [mV/m]','V=ExB DSI [km/s]'};
  else
    q_list = {'P?',['B' use_fullb '?'],'E?','Es?','VExBs?'};
    l_list = {'SC pot [-V]','B GSE [nT]','E GSE [mV/m]','E GSE [mV/m]','V=ExB GSE [km/s]'};
  end
end

old_pwd = pwd;
cd(cp.sp);

n_plots = 0;
data = {};
labels = {};

if flag_rmwhip, c_load('WHIP?',cl_id), end

% Load data
for k=1:length(q_list)
  %q_list{k}
  
  if c_load(q_list{k},cl_id)
    if have_tint
      c_eval([q_list{k} '=irf_tlim(' q_list{k} ',st+[0 dt]);'],cl_id);
    end
    
    if flag_rmwhip && (k==1 || k==3)
      if exist(irf_ssub('WHIP?',cl_id),'var')
        irf_log('proc','not using times with Whisper pulses');
        c_eval([q_list{k} '=caa_rm_blankt(' q_list{k}  ',WHIP? );'],cl_id);
      end
    end
    
    n_plots = n_plots + 1;
    
    switch k
      case 1 % P
        d_t = [];
        c_eval(['d_t=' q_list{k} ';'],cl_id)
        labels{n_plots} = l_list{k};
        if flag_ib
          [ok,P_ib] = c_load('bP?',cl_id);
          if ok
            data{n_plots} = {d_t, P_ib};
          else
            data{n_plots} = d_t;
          end
          clear ok P_ib;
        else
          data{n_plots} = d_t;
        end
        
      case 2 % B-field
        c_eval(['data{n_plots}=irf_abs(' q_list{k} '(:,1:4));'],cl_id)
        labels{n_plots} = l_list{k};
        if size(data{n_plots},2)<4  % data satity check
          continue;
        end
        elev_ang = [data{n_plots}(:,1) atan2(data{n_plots}(:,4),...
          sqrt(data{n_plots}(:,2).^2+data{n_plots}(:,3).^2))*180/pi];
        
      case 3 % E-field
        if strcmp(cs,'dsi')
          % correct DSI offsets
          dsiof = c_ctl(cl_id,'dsiof');
          if isempty(dsiof)
            if ~have_tint
              c_eval(['st=' q_list{k} '(:,1);'],cl_id)
              st = st(~isnan(st));
              if isempty(st), st = 0;
              else, st = st(1);
              end
            end
            
            [ok,Ps,msg] = c_load('Ps?',cl_id);
            if ~ok, irf_log('load',msg), end
            [dsiof_def, dam_def] = c_efw_dsi_off(st,cl_id,Ps);
            
            [ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end %#ok<NASGU>
            [ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end %#ok<NASGU>
            
            if ok1 || ok2, irf_log('calb','Using saved DSI offsets')
            else, irf_log('calb','Using default DSI offsets')
            end
            clear dsiof_def dam_def;
          else
            Ddsi = dsiof(1); Damp = dsiof(2); %#ok<NASGU>
            irf_log('calb','Using user specified DSI offsets')
          end
          clear dsiof
          c_eval([q_list{k} '=caa_corof_dsi(' q_list{k} ',Ddsi,Damp);'],cl_id);
        end
        
        c_eval(['data{n_plots}=' q_list{k} '(:,1:3);'],cl_id);
        labels{n_plots} = l_list{k};
        
      case 4 % Es-field /IB
        c_eval(['d_t=' q_list{k} ';'],cl_id);
        if flag_ib
          E_hx = data{n_plots-1};
          E_ib = [];
          c_eval(['E_ib=' q_list{k} ';'],cl_id);
          % correct DSI offsets
          dsiof = c_ctl(cl_id,'dsiof');
          if isempty(dsiof)
            if ~have_tint
              st = E_ib(:,1);
              st = st(~isnan(st));
              if isempty(st), st = 0;
              else, st = st(1);
              end
            end
            
            [dsiof_def, dam_def] = c_efw_dsi_off(st,cl_id);
            
            [ok1,Ddsi] = c_load('Ddsi?',cl_id); if ~ok1, Ddsi = dsiof_def; end
            [ok2,Damp] = c_load('Damp?',cl_id); if ~ok2, Damp = dam_def; end
            
            if ok1 || ok2, irf_log('calb','Using saved DSI offsets');
            else, irf_log('calb','Using default DSI offsets');
            end
            clear dsiof_def dam_def;
          else
            Ddsi = dsiof(1); Damp = dsiof(2);
            irf_log('calb','Using user specified DSI offsets')
          end
          clear dsiof
          E_ib=caa_corof_dsi(E_ib,Ddsi,Damp);
          if size(E_hx,2)<3   % fix for cell E_hx C3 201210->
            n_plots=n_plots+1;
            data{n_plots-1} = {E_hx{1,1}, E_ib(:,[1 2])};
            data{n_plots} = {E_hx{1,2}, E_ib(:,[1 3])};
          else
            data{n_plots-1} = {E_hx(:,[1 2]), E_ib(:,[1 2])};
            data{n_plots}   = {E_hx(:,[1 3]), E_ib(:,[1 3])};
          end
          labels{n_plots-1} = 'Ex DSI [mV/m]';
          labels{n_plots}   = 'Ey DSI [mV/m]';
          if exist('elev_ang','var')
            n_plots = n_plots + 1;
            data{n_plots} = elev_ang;
            labels{n_plots} = '\theta (B,spin) [deg]';
          end
          clear E_hx E_ib elev_ang
        else
          if exist('elev_ang','var')
            data{n_plots} = elev_ang;
          else
            data{n_plots} = d_t(:,[1 5]);
          end
          labels{n_plots} = '\theta (B,spin) [deg]';
          n_plots = n_plots + 1;
          data{n_plots} = d_t(:,1:4);
          labels{n_plots} = l_list{k};
        end
        clear d_t Ddsi Damp
        
      otherwise
        d_t = [];
        c_eval(['d_t=' q_list{k} ';'],cl_id);
        labels{n_plots} = l_list{k};
        if min(size(d_t))> 4
          data{n_plots} = d_t(:,1:4);
        else
          data{n_plots} = d_t; %#ok<*AGROW>
        end
        clear d_t
    end
    
  end
end

if n_plots==0, return, end % Nothing to plot

% Define time limits
if have_tint
  t_st = st;
  t_end = st + dt;
else
  t_st = 1e32;
  t_end = 0;
  for k=1:n_plots
    if iscell(data{k})
      for col=1:length(data{k})
        t_st = min(t_st,data{k}{col}(1,1));
        t_end = max(t_end,data{k}{col}(end,1));
      end
    else
      t_st = min(t_st,data{k}(1,1));
      t_end = max(t_end,data{k}(end,1));
    end
  end
end

% Plotting
clf;
orient tall;
h = irf_figure(747,n_plots,'reset');

for k=1:n_plots
  if isempty(data{k})
    if k==1
      st_s = epoch2iso(t_st);
      if flag_ibs == 1
        title(h(k),['EFW, Cluster ' num2str(cl_id,'%1d') ', ' vars_names ' - '  st_s(1:10)]);
      else
        title(h(k),['EFW, Cluster ' num2str(cl_id,'%1d') ', ' st_s(1:10)]);
      end
      clear st_s;
    end
    continue
  end
  if iscell(data{k}), irf_plot(h(k),data{k},'comp');
  else, irf_plot(h(k),data{k})
  end
  %axis(h(k),'tight');
  irf_zoom(h(k),'x',[t_st t_end]);
  YL = get(h(k),'YLim'); DYL = (YL(2) -YL(1))*0.05;
  set(h(k),'YLim',YL +[-DYL DYL]);
  clear YL DYL;
  ylabel(h(k),labels{k});
  if k==1
    st_s = epoch2iso(t_st);
    if flag_ibs == 1
      title(h(k),['EFW, Cluster ' num2str(cl_id,'%1d') ', ' vars_names ' - '  st_s(1:10)]);
    else
      title(h(k),['EFW, Cluster ' num2str(cl_id,'%1d') ', ' st_s(1:10)]);
    end
    clear st_s;
  end
  if k<n_plots, xlabel(h(k),''),set(h(k),'XTickLabel',[])
  else
    irf_timeaxis;
    % This magic is needed for correct location of panels on printouts
    ttt = get(h(k),'XTickLabel');
    ttt(end) = {' '};
    set(h(k),'XTickLabel',ttt)
  end
  if k==n_plots
    YLIM = 990;
    yl = get(h(k),'YLim');
    if yl(1)<-YLIM, yl(1) = -YLIM; end
    if yl(2)>YLIM, yl(2) = YLIM; end
    set(h(k),'YLim',yl);
    clear yl YLIM;
  end
end

irf_pl_add_info

lyy = 0;
for k=n_plots:-1:1
  ncol = size(data{k},2) -1;
  if iscell(data{k}), ncol = ncol +1; end
  if ncol>1
    switch ncol
      case 2
        if flag_ib
          legend(h(k),'nm','IB','Location','NorthEastOutside');
        else
          legend(h(k),'X','Y','Location','NorthEastOutside');
        end
      case 3
        legend(h(k),'X','Y','Z','Location','NorthEastOutside');
      case 4
        legend(h(k),'X','Y','Z','tot','Location','NorthEastOutside');
      otherwise
        error('too many columns');
    end
    if lyy==0, pos = get(h(k),'Position'); lyy = pos(3); clear pos, end
  end
end
clear ncol;

if lyy
  for k=n_plots:-1:1
    pos = get(h(k),'Position');
    set(h(k),'Position', [pos(1) pos(2) lyy pos(4)])
  end
end

if flag_r
  [ok,r] = c_load('R?',cl_id);
  if ok, add_position(h(n_plots),r), end
end

cd(old_pwd);

if nargout>0, out=h; end
