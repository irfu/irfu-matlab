function status = c_4_v_interacitve(x1,x2,x3,x4,column)
% C_4_V_INTERACTIVE interactively estimate the discontinuity velocity
%  status = c_4_v_interactive(x1,x2,x3,x4,column)
%

persistent ud
tlim = [];

if       (nargin==1 & isstr(x1)), action=x1;disp(['action=' action]);
elseif   (nargin == 5)                   , action='initialize';
else      help c_4_v_interactive
end

switch action,
  case 'initialize'
  dgh=figure;clf;av_figmenu;
  h(1)=subplot(3,1,1);
  c_pl_tx(x1,x2,x3,x4,column);zoom on;
  av_pl_info(['c\_4\_v\_interactive() ' datestr(now)]); % add information to the plot

  h(2)=subplot(3,1,2);
  c_pl_tx(x1,x2,x3,x4,column);


  hh=h(1,1);  % use the first subplot to estimate available time interval
  xl=get(hh,'XLim');
  hc=get(hh,'Children');
  xd=get(hc(end),'XData');
  avail=[min([xl xd]) max([xl xd])];
  presel=xl;

  dt = 0.02*diff(avail);
  xlim = [avail(1)-dt avail(2)+dt];
  ttics = timeaxis(xlim);

  ud.tlim = avail;
  ud.h=h;
  ud.data={x1 x2 x3 x4 column};

  xp=0.05;yp=0.25;
  uch1 = uicontrol('style', 'text', 'string', '[dt1 dt2 dt3 dt4] =','position', [xp yp 0.15 0.03],'units','normalized');
  ud.dt = uicontrol('style', 'edit', ...
        'string', '[0 0 0 0]', ...
      'callback', 'c_4_v_interactive(''dt'')', ...
      'backgroundcolor','white','position', [xp+0.15 yp 0.29 0.05],'units','normalized');

  xp=0.05;yp=0.2;
  uch1 = uicontrol('style', 'text', 'string', '[vx vy vz] GSE km/s =','position', [xp yp 0.15 0.03],'units','normalized');
  ud.v = uicontrol('style', 'edit', ...
        'string', '1*[0 0 0]', ...
      'callback', 'c_4_v_interactive(''v'')', ...
      'backgroundcolor','white','position', [xp+0.15 yp 0.29 0.05],'units','normalized');

  xp=0.05;yp=0.15;
  uch1 = uicontrol('style', 'text', 'string', 'Low pass filter f/Fs = ','position', [xp yp 0.15 0.03],'units','normalized');
  ud.filter = uicontrol('style', 'edit', ...
        'string', '1', ...
      'callback', 'c_4_v_interactive(''dt'')', ...
      'backgroundcolor','white','position', [xp+0.15 yp 0.1 0.05],'units','normalized');

  uimenu('label','Auto &YLim','accelerator','y','callback','c_4_v_interactive(''autoY'')');
  uimenu('label','&Distance','accelerator','d','callback','c_4_v_interactive(''distance'')');

  %set(dgh, 'userdata', dgud);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case 'dt'
    hh=ud.h(1);  % use the first subplot to estimate available time interval
    xl=get(hh,'XLim');yl=get(hh,'YLim');
    hc=get(hh,'Children');
    dt=eval(['[' get(ud.dt,'string') ']']);
    t=xl(1)+dt;
    v=c_v(t);
    tstr=['[' num2str(dt,'%7.2f') '] s'];
    vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),2) ']'];
    set(ud.v,'string',vstr);
    axes(ud.h(2));
    if eval(get(ud.filter,'string'))<1,
      x1=ud.data{1};Fs=1/(x1(2,1)-x1(1,1));flim=Fs*eval(get(ud.filter,'string'));
      for ic=1:4,eval(av_ssub('x?=av_t_lim(ud.data{?},xl+[-20/Fs 20/Fs]);x?=av_tsfilt(x?,0,flim,Fs,5);',ic)),end
      c_pl_tx(x1,x2,x3,x4,ud.data{5},1,dt);
    else,
      c_pl_tx(ud.data{1},ud.data{2},ud.data{3},ud.data{4},ud.data{5},1,dt);
    end
      axis([xl yl]);add_timeaxis;
      text(.5,-.6,['dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s GSE'],'units','normalized');
  case 'v'
    hh=ud.h(1);  % use the first subplot to estimate available time interval
    xl=get(hh,'XLim');yl=get(hh,'YLim');
    hc=get(hh,'Children');
    v=eval(['[' get(ud.v,'string') ']']),
    t=xl(1)+dt;
    dt=c_v([xl(1) v]);
    tstr=['[' num2str(dt,'%7.2f') ']'];
    vstr=[num2str(norm(v),3) ' * [' num2str(v./norm(v),'%7.2f') ']'];
    set(ud.dt,'string',tstr);
    axes(ud.h(2));
    if eval(get(ud.filter,'string'))<1,
      x1=ud.data{1};Fs=1/(x1(2,1)-x1(1,1));flim=Fs*eval(get(ud.filter,'string'));
      for ic=1:4,eval(av_ssub('x?=av_t_lim(ud.data{?},xl+[-20/Fs 20/Fs]);x?=av_tsfilt(x?,0,flim,Fs,5);',ic)),end
      c_pl_tx(x1,x2,x3,x4,ud.data{5},1,dt);
    else,
      c_pl_tx(ud.data{1},ud.data{2},ud.data{3},ud.data{4},ud.data{5},1,dt);
    end
      axis([xl yl]);add_timeaxis;
      text(.5,-.6,['dt = ' tstr '\newline V_{discontinuity}=' vstr ' km/s GSE'],'units','normalized');
  case 'distance'
    axes(ud.h(2));
      hh=ud.h(1);  % use the first subplot to estimate available time interval
      xl=get(hh,'XLim');
      v=eval(['[' get(ud.v,'string') ']']);
      tcenter = mean(xl);distance=norm(v)*diff(xl)/2;logd=log10(distance);
      if logd>round(logd), dx=10^(round(logd))/2;
      else dx=10^(round(logd))/5;
      end
      xticks=[-30:30]*dx/norm(v)/5+tcenter;
      xticklabels=cell(size(-30:30));
      for j=-30:30, xticklabels{j+31}=' ';end
      for j=-6:6, xticklabels{j*5+31}=num2str(j*dx);end
      set(ud.h(2),'xtick',xticks,'xticklabel',xticklabels);
      xlabel('km');
  case 'autoY'
    for h=ud.h(1:2),
      set(h,'YLimMode','auto');
    end
  case 'ylabel'
    h_select=get(ud.ylabpanel,'Value')-1;
    i_select=get(ud.ylab,'Value');
    strings=get(ud.ylab,'String');
    label=strings{i_select};
    if strcmp(label,'enter'); label=input('Input Y label>','s');end
    axes(ud.h(h_select));
    ylabel(label,'verticalalignment','bottom');
  case 'toggle'
    if ud.from
      set(gcbf, 'pointer', 'right');
      set(ud.helph, 'string', 'Click on axis selects ''To'' time');
      ext = get(ud.helph, 'extent');
      set(ud.helph, 'position', [395-ext(3) 55 ext(3:4)]);
      ud.from = 0;
    else
      set(gcbf, 'pointer', 'left');
      set(ud.helph, 'string', 'Click on axis selects ''From'' time');
      ext = get(ud.helph, 'extent');
      set(ud.helph, 'position', [5 55 ext(3:4)]);
      ud.from = 1;
    end
    set(gcbf, 'userdata', ud);
end

