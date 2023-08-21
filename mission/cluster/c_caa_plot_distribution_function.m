function [ax, cb] = c_caa_plot_distribution_function(varargin)
% C_CAA_PLOT_DISTRIBUTION_FUNCTION  Plots particle pitch angle distribution with
%                                   data prepared by c_caa_distribution_data.m.
%   Plots cross-section or polar pitch angle (angle with respect to the
%   magnetic field) distributions for various PEACE, CIS and RAPID products.
%
%   [ax ax_cb] = C_CAA_PLOT_DSITRIBUTION_FUNTION(ax,'tint',tint,plot_type,...
%           data_structure,'pitchangle',pitchangles,data_structure)
%   Needed input:
%       data_structure - obtained from c_caa_distribution_data.m
%            can plot two products simultaneously for polar plots
%   Additional input:
%       ax - axis handle
%       ax_cb - colorbar handle if plot_type is 'polar'
%       'tint' - time interval or single time in epoch
%           if [t1 t2] - averages over time for each energy and pitch angle
%           if t1 - takes the closest energy sweep
%       't_display' - 'given': displays given time(s) (default)
%                     'tags': displays tagged time for relevant energy sweep(s)
%       'plot_type' - 'polar' or 'cross-section'
%           'polar' supports up to two data_structure, 'cross-section' supports one
%       'pitchangle' - pitch angles to plot if 'cross-section' is chosen,
%           [0 90 180] is default
%       'emin' - will only plot values above given energy, in eV
%       'emin_scale' - the scale will start at this value, in eV, must be
%                      above emin, if emin_scale>emin then emin_scale=emin
%
%   Examples:
%       data_structure = c_caa_distribution_data('C3_CP_PEA_3DXPH_PSD');
%       h=C_CAA_PLOT_DISTRIBUTION_FUNCTION('tint',tint,'polar',data_structure);
%
% See also c_caa_distribution_data.m

% Check for axes
[ax,args,nargs] = axescheck(varargin{:});
original_args=args;
original_nargs=nargs;

% Default values
emin=[];
emin_scale=[];
plot_type='polar';
pitch_angles=[0 90 180];
t_display='given';

% Read input
n_toplot=0;
while ~isempty(args)
  if isstruct(args{1})
    n_toplot=n_toplot+1;
    to_plot{n_toplot}=args{1};
    args=args(2:end);
  elseif isstr(args{1})
    switch lower(args{1})
      case 'tint'
        tint=args{2};
        args=args(3:end);
      case 'polar'
        plot_type='polar';
        args=args(2:end);
      case 'cross-section'
        plot_type='cross-section';
        args=args(2:end);
      case 'emin'
        emin=args{2};
        args=args(3:end);
      case 'emin_scale'
        emin_scale=args{2};
        args=args(3:end);
      case 'pitchangle'
        pitch_angles=args{2};
        args=args(3:end);
      case 't_display'
        t_display=args{2};
        args=args(3:end);
      otherwise % additional input
        eval([args{1},'=args{2};'])
        args=args(2:end);
    end
  else
    args=args(2:end);
  end
end

% Return if not enough input is given.
% Reduce products if too much input is given.
if ~exist('tint','var')
  disp('No time interval was given!')
  disp('Taking the largest possible (common) time interval of the product(s)!')
  tmp_tstart = [];
  tmp_tstop = [];
  for oo=1:n_toplot
    tmp_tstart = [tmp_tstart to_plot{1}.t(1)];
    tmp_tstop = [tmp_tstop  to_plot{1}.t(end)];
  end
  tint = [max(tmp_tstart) min(tmp_tstop)];
  clear tmp_tstart tmp_tstop
  %return;
end
switch n_toplot
  case 0
    disp('No particle product was given!')
    return;
  case 1
  case 2
    if strcmp(plot_type,'cross-section')
      disp('Cross-section distributions are only plotted for single products.')
      disp('Plotting first product given.')
      to_plot=toplot(1);
    end
  otherwise
    switch plot_type
      case 'polar'
        disp('Polar distributions are plotted for at most two products.')
        disp('Plotting two first products given.')
        to_plot=to_plot(1:2);
        n_toplot=2;
      case 'cross-section'
        disp('Cross-section distributions are only plotted for single products.')
        disp('Plotting first product given.')
        to_plot=to_plot(1);
        n_toplot=1;
    end
end

% take out time interval
for k=1:n_toplot
  if length(tint)==2 % start and stop interval
    [~,ind_t{k}]=irf_tlim(to_plot{k}.t,tint);
  elseif length(tint)==1 % only one time, take closest energy sweep
    [tt,ind_t{k}]=irf_tlim(to_plot{k}.t,[tint-1 tint+1]);
    ind_min{k}=find(abs((tt-tint))==min(abs((tt-tint))));
    ind_t{k}=ind_t{k}(ind_min{k});
  end
end

% Scale must start below lowest value
if emin_scale>emin
  emin_scale=emin;
  % could also do emin=emin_scale;, but that possible cuts away data
end

% If no axes is given, initialize figure.
if isempty(ax)
  ax=irf_plot(1);
end
axes(ax);

% Plot
switch plot_type
  case 'polar'
    if n_toplot==1 % Mirror plot if only one product is given
      to_plot{2}=to_plot{1};
      ind_t{2}=ind_t{1};
    end
    for k=1:2
      % Put energy in log eV
      rlog{k} = log10(double(to_plot{k}.en_pol))';
      % Pitch angles, turn so that pitch angle 0 is on top
      theta{k} = double(to_plot{k}.f_pol)+90;
    end

    % Take away all values below emin
    if isempty(emin) % then set to lowest energy bin
      emin=min([rlog{1};rlog{1}]);
    else
      emin=log10(emin);
    end
    if emin_scale==0
      emin_scale=1; % set to 1eV, because dont want log(0)
    end
    if isempty(emin_scale) %|| log10(emin)>min([rlog{1};rlog{1}]); % Take out r0
      r0log = emin;%min(min([rlog{1};rlog{2}]));
    else
      r0log = log10(emin_scale);
    end

    for k=1:2 % find all values > emin
      ind_r{k}=find(log10(double(to_plot{k}.en_cs))>=emin);
    end

    for k=1:2 % Create surf grids
      r{k} = tocolumn(rlog{k}(ind_r{k})-r0log);
      X{k} = r{k}*cosd(theta{k});
      Y{k} = r{k}*sind(theta{k});
      C{k} = squeeze(log10(irf.nanmean(to_plot{k}.p(ind_t{k},:,ind_r{k}),1)))';
    end

    % Plot data
    surf(ax,X{1},Y{1},X{1}*0,C{1}); hold(ax,'on');
    surf(ax,-flipdim(X{2},2),Y{2},X{2}*0,C{2});
    view(ax,2);
    axis(ax,'equal','tight');
    shading(ax,'flat');
    grid(ax,'off');
    if isa(ax,'handle'), cb = colorbar(ax); % HG2
    else, cb = colorbar('peer',ax);
    end
    ylabel(cb,to_plot{1}.p_label)

    if 1 % Energy ticks
      xticks=log10([1e-1 1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7]*1e-3)-r0log;
      xticks=xticks(xticks>0);
      xticks=xticks(xticks<max([max(r{1}) max(r{2})]));
      xticklabels=cell(size(xticks));
      for k=1:length(xticklabels)
        xticklabels{k}=num2str(10.^(xticks(k)+r0log)*1e-3);
      end
      xticks=[-flipdim(xticks,2) 0 xticks];
      xticklabels=[flipdim(xticklabels,2) num2str(1e-3*10^r0log,'%.3f') xticklabels];
      yticks=xticks;
      yticklabels=xticklabels;
      set(ax,'xtick',xticks,'xticklabel',xticklabels,'TickDir','in',...
        'XMinorTick','off','ytick',yticks,'yticklabel',yticklabels)
      xlabel(ax,'Energy  [keV]'); ylabel(ax,'Energy  [keV]')
    end
    if 1 % Pitch angle labels
      rmax=max([max(r{1}) max(r{2})]);
      text(0-0.2,rmax-0.5,0,'0^o')
      text(0-0.2,-rmax+0.5,0,'180^o')
      text(-0.2-rmax+0.5,0,0,'90^o')
      text(-0.2+rmax-0.5,0,0,'90^o')
    end

  case 'cross-section'
    k=1; % only one product supported
    % Pick out angles
    n_pa=length(pitch_angles);

    if isempty(emin)
      emin=min(to_plot{k}.en_cs);
    end
    ind_r{k}=find(to_plot{k}.en_cs>emin);

    for p=1:n_pa
      diff_pa=abs(to_plot{k}.f_cs-pitch_angles(p));
      ind_pa{k,p}=find(diff_pa==min(diff_pa));
      pa_toplot{k,p}=squeeze(mean(irf.nanmean(to_plot{k}.p(ind_t{k},ind_pa{p},ind_r{k}),1),2));
      pa_legends{k,p}=num2str(pitch_angles(p),'%.0f');
      disp(['Plotting average of bins: ',num2str(to_plot{k}.f_cs(ind_pa{k,p}))])
    end
    if ~isempty(to_plot{k}.p_bg)
      PAbg=squeeze(irf.nanmean(irf.nanmean(to_plot{k}.p_bg(ind_t{k},:,ind_r{k}),1),2)); % One common zero-count level for all levels
    else
      PAbg=NaN(size(to_plot{k}.en_cs));
    end
    pa_legends{end+1}='Bg';

    % Plotting data, making string to adapt to varying # of pitch angles
    plt_str='loglog(ax';
    for p=1:n_pa; eval(['plt_str=[plt_str,'',to_plot{1}.en_cs(ind_r{k}),pa_toplot{',num2str(p),'}''];']); end
    plt_str=[plt_str,',to_plot{1}.en_cs(ind_r{k}),PAbg,''--'');'];
    eval(plt_str);

    if ~isempty(emin_scale)
      xlim(1)=emin_scale;
    else
      xlim(1)=to_plot{1}.en_cs(1)*0.8;
    end
    xlim(2)=to_plot{1}.en_cs(end)*1.2;

    set(ax,'xlim',xlim);%[to_plot{1}.en_cs(1)*0.8 to_plot{1}.en_cs(end)*1.2])
    % irf_legend(ax,{'0','90','180','-- Zero count'},[0.94 0.94])
    % irf_legend gets too big for box if box is small
    legend(ax,pa_legends,'edgecolor','w')
    ylabel(ax,to_plot{k}.p_label)
    xlabel(ax,'Energy  [eV]')
    grid(ax,'off');
    cb=[];
end

% Title
% Time
switch length(tint)
  case 1 % Only one time
    switch [t_display,num2str(n_toplot)]
      case {'given1','given2'}
        t1str=datestr(epoch2date(tint),'dd-mmm-yyyy  HH:MM:SS.FFF');
        titleStr{1}=[t1str,' UT'];
      case 'tags1'
        t1str=datestr(epoch2date(to_plot{1}.t(ind_t{1}(1))),'dd-mmm-yyyy  HH:MM:SS.FFF');
        titleStr{1}=[t1str,' UT'];
      case 'tags2'
        t1str=datestr(epoch2date(to_plot{1}.t(ind_t{1})),'dd-mmm-yyyy  HH:MM:SS.FFF');
        t2str=datestr(epoch2date(to_plot{2}.t(ind_t{2})),'dd-mmm-yyyy  HH:MM:SS.FFF');
        titleStr{1}=['Left: ',t1str,' UT'];
        titleStr{2}=['Right: ',t2str,' UT'];
    end
  case 2 % A time interval given
    switch [t_display,num2str(n_toplot)]
      case {'given1','given2'}
        t1str=datestr(epoch2date(tint(1)),'dd-mmm-yyyy  HH:MM:SS.FFF');
        t2str=datestr(epoch2date(tint(2)),'HH:MM:SS.FFF');
        titleStr{1}=[t1str,'-',t2str,' UT'];
      case 'tags1'
        t1str=datestr(epoch2date(to_plot{1}.t(ind_t{1}(1))),'dd-mmm-yyyy  HH:MM:SS.FFF');
        t2str=datestr(epoch2date(to_plot{1}.t(ind_t{1}(2))),'HH:MM:SS.FFF');
        titleStr{1}=[t1str,'-',t2str,' UT'];
      case 'tags2'
        for k=1:2
          t1str{k}=datestr(epoch2date(to_plot{k}.t(ind_t{k}(1))),'dd-mmm-yyyy  HH:MM:SS.FFF');
          t2str{k}=datestr(epoch2date(to_plot{k}.t(ind_t{k}(2))),'HH:MM:SS.FFF');
        end
        titleStr{1}=['Left: ',t1str{1},'-',t2str{1},' UT'];
        titleStr{2}=['Right: ',t1str{2},'-',t2str{2},' UT'];
    end
end

% Product
switch n_toplot % Next line is the product plotted
  case 1 % Only one product
    titleStr{end+1}=[to_plot{1}.product,to_plot{1}.detector];
  case 2 % Two products
    titleStr{end+1}=['Left: ',to_plot{1}.product,to_plot{1}.detector];
    titleStr{end+1}=['Right: ', to_plot{2}.product,to_plot{1}.detector];
end
% Change underscores to spaces
for k=1:length(titleStr), titleStr{k}(strfind(titleStr{k},'_'))=' '; end
title(ax,titleStr);
end