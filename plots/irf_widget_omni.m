function hout=irf_widget_omni(varargin)
% IRF_WIDGET_OMNI plot solar wind and indices
%
% IRF_WIDGET_OMNI
% IRF_WIDGET_OMNI('tint',tint) plot for give time interval, tint=[tmin tmax]; tmin and tmax in epoch
%
% uses OMNI data base to obtain the data
% if interval larger than 1day uses hourly database
Units=irf_units;
%% Check input
if nargin==0 % initialize
  action='initialize';
  have_options=0;
elseif nargin==1 % should be action
  if ischar(varargin{1})
    action=varargin{1};
    have_options=0;
  else
    return;
  end
else % there are additional options
  have_options=1;
  args=varargin;
end

%% Check input options
while have_options
  l = 1;
  switch(lower(args{1}))
    case {'t','tstart'}
      if nargs>1 && isnumeric(args{2})
        tstart = args{2};
        data=get(gcf,'userdata');
        data.t=tstart;
        set(gcf,'userdata',data);
        l = 2;
      else, irf.log('warning','wrongArgType : tint must be numeric')
      end
  end
  args = args(l+1:end);
  if isempty(args), break, end
end


%% Check actions
switch lower(action)
  case 'initialize' % read in all data and open figure
    initialize_figure(6); % default 5 subplots
    data=get(gcf,'userdata');
    if ~isfield(data,'t')
      time=irf_time([2010 12 31 01 01 01]);dt=24*3600;
      data.t=time;
      data.dt=dt;
      set(gcf,'userdata',data);
    end
    irf_widget_omni('read_data');
    irf_widget_omni('plot');

  case 'read_data'
    data=get(gcf,'userdata');
    tint=[data.t data.t+data.dt];
    omni2=irf_get_data(tint,'dst,f10.7','omni2');
    if diff(tint)< 48*3600 % interval larger than 48 h use 1h resolution
      disp(['Reading OMNI_MIN 1min data :' irf_time(tint,'tint>utc')]);
      ff=irf_get_data(tint,'b,bx,bygsm,bzgsm,T,n,v,P,beta,pc,ae,al,au','omni_min');
    else
      disp(['Reading OMNI2 1h data :' irf_time(tint,'tint>utc')]);
      ff=irf_get_data(tint,'b,bx,bygsm,bzgsm,T,n,v,P,beta,pc,ae,al,au','omni2');
    end
    data.ff=ff;
    data.omni2=omni2;
    set(gcf,'userdata',data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% action plot %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'plot'
    data=get(gcf,'userdata');

    %%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
    h=data.subplot_handles;
    ff=data.ff;
    tint=[data.t data.t+data.dt];
    for j=1:numel(h)
      ud=get(h(j),'userdata');
      if isstruct(ud), ud=rmfield(ud,'zoom_x');end % remove zoom_x if exists
      set(h(j),'userdata',ud);
    end
    %%%%%%%%%%%%
    % B field
    hca=h(1);
    irf_plot(hca,ff(:,[1 3 4 5 2]));
    ylabel(hca,'B [nT] GSM');
    irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.05]);
    title(hca,'OMNI solar wind parameters');

    %%%%%%%%%%%%
    % Velocity
    hca=h(2);
    irf_plot(hca,ff(:,[1 8]));
    ylabel(hca,'V [km/s]');

    %%%%%%%%%%%%
    % Density
    hca=h(3);
    irf_plot(hca,ff(:,[1 7]));
    ylabel(hca,'density [cc]');

    %%%%%%%%%%%%
    % Proton temperature
    hca=h(4);
    Tp=ff(:,[1 6]);
    Tp(:,2)=Tp(:,2)*Units.kB/Units.e;
    irf_plot(hca,Tp);
    ylabel(hca,'T proton [eV]');

    %%%%%%%%%%%%
    % Pressure
    hca=h(5);
    irf_plot(hca,ff(:,[1 9]));
    ylabel(hca,'pressure [nPa]');

    %%%%%%%%%%%%
    % A indices
    hca=h(6);
    irf_plot(hca,ff(:,[1 12 13 14]));
    ylabel(hca,'AE,AL,AU [nT]');
    irf_legend(hca,{'AE','AL','AU'},[0.02 0.05]);

    irf_zoom(h,'x',tint);
    irf_zoom(h,'y');
    irf_timeaxis(h);

  case 'new_start_time'
    data=get(gcf,'userdata');
    xx=inputdlg('Enter new start time. [yyyy mm dd hh mm ss]','**',1,{mat2str(irf_time(data.t,'vector'),4)});
    if ~isempty(xx)
      variable_str=xx{1};
      data.t=irf_time(eval(variable_str));
      set(gcf,'userdata',data);
      irf_widget_omni('read_data');
      irf_widget_omni('plot');
    end
  case 'new_time_interval'
    data=get(gcf,'userdata');
    xx=inputdlg('Enter time interval in hours.','**',1,{mat2str(data.dt/3600,3)});
    if ~isempty(xx)
      variable_str=xx{1};
      data.dt=eval(variable_str)*3600;
      set(gcf,'userdata',data);
      irf_widget_omni('read_data');
      irf_widget_omni('plot');
    end
end
if nargout
  hout=h;
else
  clear hout;
end

function initialize_figure(number_of_subplots)
figure
set(gcf,'color','white'); % white background for figures (default is grey)
set(gcf,'renderer','zbuffer'); % opengl has problems on Mac (no log scale in spectrograms)
set(gcf,'PaperUnits','centimeters');
set(gcf,'defaultlinelinewidth',1.0);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0])

xSize = 10;
ySize = 5+5*sqrt(number_of_subplots);
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
sz=get(0,'screensize');
xx=min(min(600,sz(3))/xSize,min(900,sz(4))/ySize); % figure at least 600 wide or 900 height but not outside screen
set(gcf,'Position',[10 10 xSize*xx ySize*xx])
clear xSize sLeft ySize yTop

for j=1:number_of_subplots
  c(j)=irf_subplot(number_of_subplots,1,-j);
  cla(c(j));
  set(c(j),'tag','');
end
user_data = get(gcf,'userdata');
user_data.subplot_handles = c;
user_data.current_panel=0;
set(gcf,'userdata',user_data);
figure(gcf); % bring figure to front

% generate menus
if isempty(findobj(gcf,'type','uimenu','label','&Options'))
  hcoordfigmenu=uimenu('Label','&Options');
  uimenu(hcoordfigmenu,'Label','Start time','Callback','irf_widget_omni(''new_start_time'')')
  uimenu(hcoordfigmenu,'Label','Time interval','Callback','irf_widget_omni(''new_time_interval'')')
  user_data = get(gcf,'userdata');
  user_data.coordfigmenu=1;
  set(gcf,'userdata',user_data);
end

