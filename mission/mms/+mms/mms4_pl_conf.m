function hout=mms4_pl_conf(time,R_init)
%MMS.MMS4_PL_CONF   Plot the configuration of MMS in XYZ coordinates
%
%   h = MMS.MMS4_PL_CONF([time]);

mms_marker={{'ks','markersize',12},{'rd','markersize',12},...
  {'go','markersize',12,'color',[0 0.6 0]},{'bv','markersize',12}};
mms_marker_small={{'ks','markersize',8},{'rd','markersize',8},...
  {'go','markersize',8,'color',[0 0.6 0]},{'bv','markersize',8}};
mms_marker_shaded={{'ks','color',[0.3 0.3 0.3]},...
  {'rd','color',[1 0.3 0.3]},{'go','color',[.3 1 .3]},{'bv','color',[.3 .3 1]}};
R1=[];R2=[];R3=[];R4=[]; %#ok<NASGU>
R.C1=[];R.C2=[];R.C3=[];R.C4=[];R.R=[]; % positions of each s/c and mass centrum
tr=[];r=[]; %#ok<NASGU>
XRe=cell(1,4);rr=cell(1,4);
action = 'initialize'; % default
plot_type='default';

if       (nargin==1 && ischar(time))
  action=time;
  irf.log('debug',['action=' action]);
end
if nargin==0 % default time
  time=[2015 09 01 01 01 01];
  t=irf_time(time);
end
if (nargin==1 && ischar(time))
  action=time;
  irf.log('debug',['action=' action]);
end

if nargin ==2
  R=R_init;
end

sc_list=1:4;
if exist('coord_label','var') % define coord label if not defined so far
  if isempty(coord_label)
    coord_label='GSE';
  end
else % in case coord_sys not specified
  coord_label='GSE';
end

switch lower(action)
  case 'initialize' % read in all data and open figure
    if isa(time,'GenericTimeArray')
      t = time.start.epochUnix;
    elseif isnumeric(time) && isscalar(time) % time given in epoch
      t=time;
    elseif isnumeric(time) && length(time)==6 % time given as vector
      t=irf_time(time);
    else
      irf.log('critical','Wrong input format of time.');
      return;
    end
    % Open new figure
    figNumber=figure( ...
      'Name','MMS s/c configuration in XYZ', ...
      'Tag','mmsconfXYZ');
    set(figNumber,'color','white'); % white background for figures (default is grey)
    set(figNumber,'defaultLineLineWidth', 1.5);
    set(figNumber,'defaultAxesFontSize', 12);
    set(figNumber,'defaultTextFontSize', 12);
    set(figNumber,'defaultAxesFontUnits', 'pixels');
    menus;
    data.t=t;
    data.figNumber=figNumber;
    data.coord_label=coord_label;
    data.plot_type=plot_type;
    data.sc_list=sc_list;
    c_eval('data.r.C?=[];data.R.C?=[];');
    if nargin ==2
      t0 = EpochUnix(data.t);
      timeLine = R.time.epochUnix;
      R.C1 = [timeLine R.gseR1];
      R.C2 = [timeLine R.gseR2];
      R.C3 = [timeLine R.gseR3];
      R.C4 = [timeLine R.gseR4];
      data.R =R;
      set(gcf,'userdata',data);
    else
      set(gcf,'userdata',data);
      mms.mms4_pl_conf('read_position');
    end
    mms.mms4_pl_conf(plot_type);
  case 'read_position'
    data=get(gcf,'userdata');
    R=data.R;
    %		tint = [data.t-120 data.t+120]; % interval to read position
    if ~is_R_ok     % try reading from disk mat files
      t0 = EpochUnix(data.t);
      RR = mms.get_data('R_gse',irf.tint(EpochTT(t0+(-120)),240));
      timeLine = RR.time.epochUnix;
      R.C1 = [timeLine RR.gseR1];
      R.C2 = [timeLine RR.gseR2];
      R.C3 = [timeLine RR.gseR3];
      R.C4 = [timeLine RR.gseR4];
    end
    data.R=R;
    set(gcf,'userdata',data);
    return;
  case 'gse'
    data=get(gcf,'userdata');
    data.coord_label='GSE';
    data.r=data.R;
    set(gcf,'userdata',data);
    if strcmp(data.plot_type,'lmn') % need to redraw lmn text
      mms.mms4_pl_conf('lmn');
    else
      mms.mms4_pl_conf('plot');
    end
  case 'gsm'
    data=get(gcf,'userdata');
    data.coord_label='GSM';
    c_eval('data.r.C?=irf_gse2gsm(data.R.C?);',data.sc_list);
    set(gcf,'userdata',data);
    if strcmp(data.plot_type,'lmn') % need to redraw lmn text
      mms.mms4_pl_conf('lmn');
    else
      mms.mms4_pl_conf('plot');
    end
  case 'default'
    data=get(gcf,'userdata');
    ss=get(0,'screensize');
    sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
    set(gcf,'Position',[10 10 600/sfactor 1000/sfactor]);
    delete(findall(gcf,'Type','axes'))
    data.h=[];h=gobjects(1,8);
    xsize=.35;ysize=.195;dx=.13;dy=.05;
    for ix=1:2
      for iy=1:4
        h(iy*2-2+ix)=axes('position',[dx*ix+(ix-1)*xsize dy*(5-iy)+(4-iy)*ysize xsize ysize]); %#ok<LAXES>
      end
    end
    axis(h(8),'off');
    axis(h(1:3),[-19.99 19.99 -19.99 19.99]);
    axis(h(4),[-19.99 19.99 0 19.99]);
    for ii=1:4, hold(h(ii),'on');daspect(h(ii),[1 1 1]);end
    data.h=h;
    data.showClusterDescription = true; % show cluster description
    data.plot_type='default';
    set(gcf,'userdata',data);
    mms.mms4_pl_conf(data.coord_label);
  case 'compact'
    data=get(gcf,'userdata');
    clf;menus;
    set(gcf,'userdata',data);
    ss=get(0,'screensize');
    sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
    set(gcf,'Position',[10 ss(4)-80-700/sfactor 700/sfactor 700/sfactor]);
    initialize_figure;
    h=gobjects(0);
    h(1)=axes('position',[0.1  0.56 0.3 0.36]); % [x y dx dy]
    h(2)=axes('position',[0.59 0.56 0.3 0.36]); % [x y dx dy]
    h(3)=axes('position',[0.1  0.06 0.3 0.36]); % [x y dx dy]
    h(4)=axes('position',[0.59 0.06 0.3 0.36]); % [x y dx dy]
    h(21) = axes('Position',get(h(1),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
    h(22) = axes('Position',get(h(2),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
    h(23) = axes('Position',get(h(3),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
    axis(h(4),'off');hold(h(4),'on');
    data.h=h;
    data.showClusterDescription = true; % show cluster description
    data.plot_type='compact';
    set(gcf,'userdata',data);
    mms.mms4_pl_conf(data.coord_label);
  case 'config3d'
    data=get(gcf,'userdata');
    clf;menus;
    set(gcf,'userdata',data);
    ss=get(0,'screensize');
    sfactor=max([1 700/(ss(3)-80) 1000/(ss(4)-80)]);
    set(gcf,'Position',[10 ss(4)-80-500/sfactor 500/sfactor 500/sfactor]);
    initialize_figure;
    h=gobjects(0);
    h(1)=axes('position',[0.15  0.16 0.7 0.7]); % [x y dx dy]
    h(2)=axes('position',[0.5 0.8 0.5 0.2]);    % for legends
    data.h=h;
    data.showClusterDescription = true; % show cluster description
    data.plot_type='config3d';
    data.showClusterDescription = false;
    set(gcf,'userdata',data);
    mms.mms4_pl_conf(data.coord_label);
  case 'lmn'
    mms.mms4_pl_conf('compact');
    data=get(gcf,'userdata');
    delete(data.h(21)); % remove secondary axis
    delete(data.h(22));
    delete(data.h(23));
    data.h(5:end)=[];
    data.plot_type='lmn';
    hca=data.h(4);
    cla(hca);hold(hca,'on');
    % L vector
    callbackStr='mms.mms4_pl_conf(''plot'')';
    data.LMN_text_hndl=uicontrol('string',['LMN vectors in ' data.coord_label '. One of L/M/N can be zero.'],'style','text','units','normalized','Position',[0.5 0.25 .35 .05]);
    uicontrol('string','L','style','text','units','normalized','Position',[0.5 0.2 .05 .04])
    if isfield(data,'Lstr'), Lstr=data.Lstr; else, Lstr='[1 0 0]';end
    data.L_hndl=uicontrol('Style','edit','Units','normalized', ...
      'Position',[0.55 0.2 .3 .05],'String',Lstr,'Callback',callbackStr);
    uicontrol('string','M','style','text','units','normalized','Position',[0.5 0.15 .05 .04])
    if isfield(data,'Lstr'), Mstr=data.Mstr; else, Mstr='[0 1 0]';end
    data.M_hndl=uicontrol('Style','edit','Units','normalized', ...
      'Position',[0.55 0.15 .3 .05],'String',Mstr,'Callback',callbackStr);
    uicontrol('string','N','style','text','units','normalized','Position',[0.5 0.1 .05 .04])
    if isfield(data,'Nstr'), Nstr=data.Nstr; else, Nstr='0';end
    data.N_hndl=uicontrol('Style','edit','Units','normalized', ...
      'Position',[0.55 0.1 .3 .05],'String',Nstr,'Callback',callbackStr);
    set(gcf,'userdata',data);
    mms.mms4_pl_conf('plot');
  case 'supercompact'
    data=get(gcf,'userdata');
    ss=get(0,'screensize');
    sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
    set(gcf,'Position',[10 ss(4)-80-350/sfactor 750/sfactor 350/sfactor]);
    initialize_figure;
    h=gobjects(0);
    h(1)=axes('position',[0.09 0.13 0.32 0.74]); % [x y dx dy]
    h(2)=axes('position',[0.59 0.13 0.32 0.74]); % [x y dx dy]
    h(3)=axes('position',[0 0 1 1]);
    h(21)=axes('Position',get(h(1),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
    h(22)=axes('Position',get(h(2),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
    axis(h(3),'off');
    data.h=h;
    data.showClusterDescription = false;
    data.plot_type='supercompact';
    set(gcf,'userdata',data);
    mms.mms4_pl_conf(data.coord_label);
  case 'supercompact2'
    data=get(gcf,'userdata');
    ss=get(0,'screensize');
    sfactor=max([1 600/(ss(3)-80) 1000/(ss(4)-80)]);
    set(gcf,'Position',[10 ss(4)-80-650/sfactor 350/sfactor 650/sfactor]);
    initialize_figure;
    h=gobjects(0);
    h(2)=axes('position',[0.18 0.06 0.63 0.36]); % [x y dx dy]
    h(1)=axes('position',[0.18 0.57 0.63 0.36]); % [x y dx dy]
    h(3)=axes('position',[0 0 1 1]);axis off;
    h(21)=axes('Position',get(h(1),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
    h(22)=axes('Position',get(h(2),'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k');
    axis(h(3),'off');
    data.h=h;
    data.showClusterDescription = false;
    data.plot_type='supercompact2';
    set(gcf,'userdata',data);
    mms.mms4_pl_conf(data.coord_label);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% action plot %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'plot'
    data=get(gcf,'userdata');
    R=data.R;
    flag_using_omni_data=[]; % default that there is no plot involving OMNI data
    % estimate mass center
    drref=0;      % scale for plots
    R.R=[0 0 0 0];% default mass centrum in origo
    if is_R_ok
      for ic = data.sc_list
        rr{ic}=irf_resamp(data.r.(['C' num2str(ic)]),data.t);
        rr{ic}=rr{ic}(1:4); % remove magnitude in 5th col if present
        R.R=R.R+rr{ic}/length(data.sc_list);
      end
      % estimate relative position wrt mass center
      x=cell(1,4); % relative position
      for ic = data.sc_list
        x{ic}=rr{ic}-R.R;
        x{ic}(1)=data.t;
        x{ic}=irf_abs(x{ic});
        drref=max([drref x{ic}(5)]);
        XRe{ic}=irf_tappl(rr{ic},'/6372');
      end
    else
      %
      disp('');
    end
    if drref==0, drref=1; end % in case 1 satellite or satellites in the same location:)
    set(gcf,'userdata',data);
    %%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%
    h=data.h;
    coord_label=data.coord_label;
    sc_list=data.sc_list;
    switch data.plot_type
      case 'default'
        plotAxes = 'XZ';
        plot_position(h(1));
        [flag_using_omni_data,omni]=add_magnetopause(h(1));
        add_bowshock(h(1));
        add_Earth(h(1));

        plotAxes = 'YZ';
        plot_position(h(2));
        add_Earth(h(2),'day');

        plotAxes = 'XY';
        plot_position(h(3));
        add_magnetopause(h(3));
        add_bowshock(h(3));
        add_Earth(h(3));

        cla(h(4));
        for ic=sc_list
          plot(h(4),XRe{ic}(2),sqrt(XRe{ic}(3)^2+XRe{ic}(4)^2),mms_marker_small{ic}{:});
          hold(h(4),'on');
        end
        xlabel(h(4),['X [R_E] ' coord_label]);
        ylabel(h(4),['sqrt (Y^2+Z^2) [R_E] ' coord_label]);
        grid(h(4),'on');
        set(h(4),'xdir','reverse')
        add_magnetopause(h(4));
        add_bowshock(h(4));
        add_Earth(h(4));

        plotAxes = 'XZ';
        plot_relative_position(h(5));

        plotAxes = 'YZ';
        plot_relative_position(h(6));

        plotAxes = 'XY';
        plot_relative_position(h(7));
      case 'compact'

        cla(h(1));
        plotAxes = 'XZ';
        plot_relative_position(h(1));
        fix_RE_axis(h(1),h(21));

        cla(h(2));
        plotAxes = 'YZ';
        plot_relative_position(h(2));
        fix_RE_axis(h(2),h(22));

        cla(h(3));
        plotAxes = 'XY';
        plot_relative_position(h(3));
        fix_RE_axis(h(3),h(23));
        %
      case 'config3d'

        hold(h(1),'off');
        for ic=sc_list
          % put Cluster markers
          plot3(h(1),x{ic}(2),x{ic}(3),x{ic}(4),mms_marker{ic}{:});
          hold(h(1),'on');
          % lines from sc to projection planes
          lineProperties = {'parent',h(1),'linestyle',':','linewidth',0.6};
          line([x{ic}(2) -drref],[x{ic}(3) x{ic}(3)],[x{ic}(4) x{ic}(4)],lineProperties{:});
          line([x{ic}(2) x{ic}(2)],[x{ic}(3) -drref],[x{ic}(4) x{ic}(4)],lineProperties{:});
          line([x{ic}(2) x{ic}(2)],[x{ic}(3) x{ic}(3)],[x{ic}(4) -drref],lineProperties{:});
          % put Cluster markers on projections
          plot3(h(1),-drref,x{ic}(3),x{ic}(4),mms_marker_shaded{ic}{:});
          plot3(h(1),x{ic}(2),-drref,x{ic}(4),mms_marker_shaded{ic}{:});
          plot3(h(1),x{ic}(2),x{ic}(3),-drref,mms_marker_shaded{ic}{:});
        end
        axis(h(1),[-drref drref -drref drref -drref drref ]);
        for ii=1:4
          for jj=ii+1:4
            if any(find(sc_list==ii)) && any(find(sc_list==jj))
              line([x{ii}(2) x{jj}(2)],...
                [x{ii}(3) x{jj}(3)],...
                [x{ii}(4) x{jj}(4)],...
                'parent',h(1),'linewidth',2,'linestyle','-',...
                'color',[0.6 0.6 0.6]);
            end
          end
        end
        % draw origo axis
        line([-drref drref],[-drref -drref],[-drref -drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
        line([-drref -drref],[-drref drref],[-drref -drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
        line([-drref -drref],[-drref -drref],[-drref drref],'parent',h(1),'linestyle','-','color','k','linewidth',0.6);
        text(0.1,1,0,irf_time(data.t,'utc_yyyy-mm-dd HH:MM:SS.mmm'),'parent',h(1),'units','normalized','horizontalalignment','center','fontsize',9);
        xlabel(h(1),['{\Delta}X [km] ' coord_label]);
        ylabel(h(1),['{\Delta}Y [km] ' coord_label]);
        zlabel(h(1),['{\Delta}Z [km] ' coord_label]);
        set(h(1),'xdir','reverse');
        set(h(1),'ydir','reverse');
        grid(h(1),'on');
        axis(h(1),[-drref drref -drref drref]);
        hold(h(1),'off');

        hca=h(2);
        hold(hca,'on');
        axis(hca,[0 1 0 1]);
        yy=.5;dxx=.25;xs=.01;
        plot(hca,xs,yy,'ks',xs+1*dxx,yy,'rd',xs+2*dxx,yy,'go',xs+3*dxx,yy,'bv','LineWidth',1.5);
        text(xs+0.03,yy,'MMS1','parent',hca);
        text(xs+1*dxx+0.03,yy,'MMS2','parent',hca);
        text(xs+2*dxx+0.03,yy,'MMS3','parent',hca);
        text(xs+3*dxx+0.03,yy,'MMS4','parent',hca);
        axis(hca,'off');

      case 'lmn'
        cla(h(1));
        x=get_in_lmn(x);
        for ic=1:numel(sc_list)
          plot(h(1),x{ic}(2),x{ic}(4),mms_marker{sc_list(ic)}{:});
          hold(h(1),'on');
        end
        xlabel(h(1),'L [km]');
        ylabel(h(1),'N [km] ');
        set(h(1),'xdir','reverse');
        grid(h(1),'on');
        axis(h(1),[-drref drref -drref drref]);

        cla(h(2));
        for ic=1:numel(sc_list)
          plot(h(2),x{ic}(3),x{ic}(4),mms_marker{sc_list(ic)}{:});
          hold(h(2),'on');
        end
        xlabel(h(2),'M [km] ');
        ylabel(h(2),'N [km] ');
        grid(h(2),'on');
        axis(h(2),[-drref drref -drref drref]);

        cla(h(3));
        for ic=1:numel(sc_list)
          plot(h(3),x{ic}(2),x{ic}(3),mms_marker{sc_list(ic)}{:});
          hold(h(3),'on');
        end
        xlabel(h(3),'L [km] ');
        ylabel(h(3),'M [km] ');
        grid(h(3),'on');
        axis(h(3),[-drref drref -drref drref]);
        set(h(3),'xdir','reverse')
        set(h(3),'ydir','reverse')

      case 'supercompact'

        plotAxes = 'XZ';
        plot_relative_position(h(1));
        fix_RE_axis(h(1),h(21));
        irf_legend(h(1),irf_time(data.t,'epoch>utc_yyyy-mm-dd HH:MM:SS'),[0.02 0.98],'fontsize',9);

        plotAxes = 'YZ';
        plot_relative_position(h(2));
        text_Cluster_markers(h(2));
        fix_RE_axis(h(2),h(22));
        %
      case 'supercompact2'

        plotAxes = 'XZ';
        plot_relative_position(h(1));
        fix_RE_axis(h(1),h(21));
        irf_legend(h(1),irf_time(data.t,'epoch>utc_yyyy-mm-dd HH:MM:SS'),[0.02 0.98],'fontsize',9);

        plotAxes = 'XY';
        plot_relative_position(h(2));
        text_Cluster_markers(h(2));
        fix_RE_axis(h(2),h(22));
    end
    if data.showClusterDescription
      if strcmpi(data.plot_type,'compact') || ...
          strcmpi(data.plot_type,'lmn')  % decide in which axes write labels
        hca=h(4);
      else
        hca=h(8);
      end
      cla(hca);
      axis(hca,[0 1 0 1]);
      hold(hca,'on');
      yy=1;
      plot(hca,0,yy,'ks',.27,yy,'rd',.54,yy,'go',.81,yy,'bv','LineWidth',1.5);
      text(0.03,yy,'MMS1','parent',hca);
      text(.30,yy,'MMS2','parent',hca);
      text(.57,yy,'MMS3','parent',hca);
      text(.84,yy,'MMS4','parent',hca);
      axis(hca,'off');
      ht=irf_legend(hca,['mms.mms4_pl_conf() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))],[0,0],'fontsize',8);
      set(ht,'interpreter','none');
      htime=irf_legend(hca,['MMS configuration\newline ' irf_time(data.t,'utc_yyyy-mm-dd HH:MM:SS.mmm')],[0,.95]);
      set(htime,'fontsize',12);
      if ~isempty(flag_using_omni_data)
        if ~strcmpi(data.plot_type,'lmn')
          if flag_using_omni_data==1 % succeeded downloading OMNI
            irf_legend(hca,['IMF from OMNI 1h database:\newline P=' num2str(omni.Dp,'%6.1f') '[nPa],\newline Bx=' num2str(omni.Bx,'%6.1f') ',By=' num2str(omni.By,'%6.1f') ',Bz=' num2str(omni.Bz,'%6.1f') '[nT] GSM' ],[0,0.7]);
          elseif flag_using_omni_data==0 % did not succeeded downloading OMNI
            irf_legend(hca,['IMF using assumed model:\newline P=' num2str(omni.Dp,'%6.1f') '[nPa],\newline Bz=' num2str(omni.Bz,'%6.1f') '[nT] GSM' ],[0,0.7]);
          end
        end
      end
    end
  case 'new_time'
    data=get(gcf,'userdata');
    xx=inputdlg('Enter new time. [yyyy mm dd hh mm ss]','**',1,{mat2str(irf_time(data.t,'vector'))});
    if ~isempty(xx)
      variable_str=xx{1};
      data.t=irf_time(eval(variable_str));
      set(gcf,'userdata',data);
      mms.mms4_pl_conf('read_position');
      mms.mms4_pl_conf(data.coord_label);
    end
  case 'new_sc_list'
    xx=inputdlg('Enter new sc_list. ex. [1 3 4]','**',1,{mat2str(sc_list)});
    if ~isempty(xx)
      variable_str=xx{1};
      sc_list=eval(variable_str);
      data=get(gcf,'userdata');
      data.sc_list=sc_list;
      set(gcf,'userdata',data);
      mms.mms4_pl_conf('read_position');
      mms.mms4_pl_conf(data.coord_label);
    end
end
if nargout
  data=get(gcf,'userdata');
  hout = data.h;
else
  clear hout;
end

  function plot_position(h)
    ax1=h;
    cla(ax1);
    colX = plotAxes(1)-'W'+1;
    colY = plotAxes(2)-'W'+1;
    for iSc=data.sc_list
      if is_R_ok(iSc)
        plot(ax1,XRe{iSc}(colX),XRe{iSc}(colY),mms_marker_small{iSc}{:},'LineWidth',1.5);
        hold(ax1,'on');
      end
    end
    xlabel(ax1,[plotAxes(1) ' [RE] ' coord_label]);
    ylabel(ax1,[plotAxes(2) ' [RE] ' coord_label]);
    if colX == 2 , set(ax1,'xdir','reverse'); end
    if colY == 3 , set(ax1,'ydir','reverse'); end
    grid(ax1,'on');
  end
  function plot_relative_position(h)
    ax1=h;
    cla(ax1);
    colX = plotAxes(1)-'W'+1;
    colY = plotAxes(2)-'W'+1;
    hold(ax1,'off');
    for iSc=data.sc_list
      if is_R_ok(iSc)
        plot(ax1,x{iSc}(colX),x{iSc}(colY),mms_marker{iSc}{:},'LineWidth',1.5);
        hold(ax1,'on');
      end
    end
    xlabel(ax1,['{\Delta}' plotAxes(1) ' [km] ' coord_label]);
    ylabel(ax1,['{\Delta}' plotAxes(2) ' [km] ' coord_label]);
    if colX == 2 , set(ax1,'xdir','reverse'); end
    if colY == 3 , set(ax1,'ydir','reverse'); end
    grid(ax1,'on');
    axis(ax1,[-drref drref -drref drref]);
  end
  function fix_RE_axis(axis1,axis2)
    colX = plotAxes(1)-'W'+1;
    colY = plotAxes(2)-'W'+1;
    if drref>10000, REform='%6.1f';
    elseif drref<100, REform='%6.3f';
    else, REform='%6.2f';
    end
    xlim_ax1=get(axis1,'XLim');ylim_ax1=get(axis1,'YLim');
    xtick_ax1=get(axis1,'XTick');ytick_ax1=get(axis1,'YTick');
    xlabel(axis2,[plotAxes(1) ' [R_E] ' coord_label]);
    ylabel(axis2,[plotAxes(2) ' [R_E] ' coord_label]);
    xtlax2=num2str((xtick_ax1'+R.R(colX))/6372,REform);
    ytlax2=num2str((ytick_ax1'+R.R(colY))/6372,REform);
    set(axis2,'xdir',get(axis1,'xdir'));
    set(axis2,'ydir',get(axis1,'ydir'));
    set(axis2,'xlim',xlim_ax1,'xticklabel',xtlax2);
    set(axis2,'ylim',ylim_ax1,'yticklabel',ytlax2);
  end
  function text_Cluster_markers(ax1)
    colX = plotAxes(1)-'W'+1;
    colY = plotAxes(2)-'W'+1;
    signX = sign(~strcmp(get(ax1,'xdir'),'reverse')-0.5);
    signY = sign(~strcmp(get(ax1,'ydir'),'reverse')-0.5);
    for iSc=data.sc_list
      if is_R_ok(iSc)
        coordX = x{iSc}(colX);
        coordY = x{iSc}(colY);
        if coordX*signX > 0
          horAl = 'right';
        else
          horAl = 'left';
        end
        if coordY*signY > 0
          verAl = 'top';
        else
          verAl = 'bottom';
        end
        text(coordX,coordY,[' M' num2str(iSc) ' '],'parent',ax1,...
          'HorizontalAlignment',horAl,'VerticalAlignment',verAl);
      end
    end
  end
  function initialize_figure
    delete(findall(gcf,'Type','axes'))
    set(gcf,'PaperUnits','centimeters')
    pos=get(gcf,'position');
    factor=min(16./pos(3:4));
    xSize = pos(3)*factor; ySize = pos(4)*factor; % should not exceed 16cm in x or y direction
    xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
    set(gcf,'defaultLineLineWidth', 1.5);
    set(gcf,'defaultAxesFontUnits', 'pixels');
    set(gcf,'defaultAxesFontSize', 12);
  end
  function menus
    % generate menus
    if isempty(findobj(gcf,'type','uimenu','label','&Options'))
      hcoordfigmenu=uimenu('Label','&Options');
      uimenu(hcoordfigmenu,'Label','New time','Callback','mms.mms4_pl_conf(''new_time'')')
      uimenu(hcoordfigmenu,'Label','&GSE','Callback','mms.mms4_pl_conf(''GSE'')','Accelerator','G')
      uimenu(hcoordfigmenu,'Label','GS&M','Callback','mms.mms4_pl_conf(''GSM'')','Accelerator','M')
      uimenu(hcoordfigmenu,'Label','normal','Callback','mms.mms4_pl_conf(''default'')')
      uimenu(hcoordfigmenu,'Label','config3D','Callback','mms.mms4_pl_conf(''config3D'')')
      uimenu(hcoordfigmenu,'Label','compact','Callback','mms.mms4_pl_conf(''compact'')')
      uimenu(hcoordfigmenu,'Label','supercompact','Callback','mms.mms4_pl_conf(''supercompact'')')
      uimenu(hcoordfigmenu,'Label','supercompact2','Callback','mms.mms4_pl_conf(''supercompact2'')')
      uimenu(hcoordfigmenu,'Label','LMN','Callback','mms.mms4_pl_conf(''lmn'')')
      uimenu(hcoordfigmenu,'Label','New sc_list','Callback','mms.mms4_pl_conf(''new_sc_list'')')
      user_data = get(gcf,'userdata');
      user_data.coordfigmenu=1;
      set(gcf,'userdata',user_data);
    end
  end
  function answer=is_R_ok(sc)
    % check if position data are ok for spacecraft number 'sc'
    % if input argument not given check if ok for all spacecraft that needs
    % to be plotted.
    if nargin == 0
      scList = data.sc_list;
    else
      scList = sc;
    end
    for iSc=scList
      strSc = ['C' num2str(iSc)];
      if numel(R.(strSc)) < 8 % less than 2 time points
        answer=false;
        return;
      else
        tintR=[R.(strSc)(1,1) R.(strSc)(end,1)];
        if (tintR(1)>data.t) || (tintR(2)<data.t)
          answer=false;
          return;
        end
      end
    end
    answer=true;
  end
  function [flag_omni,omni]=add_magnetopause(h)
    % flag_omni=1 - using OMNI, flag_omni=0 - using default values
    tMP=getfield(get(gcf,'userdata'),'t');
    [xMP,yMP,omni]=irf_magnetosphere('mp_shue1998',tMP);
    if isempty(xMP)
      flag_omni=0;
      [xMP,yMP,omni]=irf_magnetosphere('mp_shue1998');
    else
      flag_omni=1;
    end
    xMP=[fliplr(xMP) xMP];
    yMP=[fliplr(yMP) -yMP];
    line(xMP,yMP,'parent',h,'linewidth',0.5,'linestyle','-','color','k');
  end
  function [flag_omni,omni]=add_bowshock(h)
    % flag_omni=1 - using OMNI, flag_omni=0 - using default values
    t=getfield(get(gcf,'userdata'),'t');
    [xBS,yBS,omni]=irf_magnetosphere('bs',t);
    if isempty(xBS)
      flag_omni=0;
      [xBS,yBS,omni]=irf_magnetosphere('bs');
    else
      flag_omni=1;
    end
    xBS=[fliplr(xBS) xBS];
    yBS=[fliplr(yBS) -yBS];
    line(xBS,yBS,'parent',h,'linewidth',0.5,'linestyle','-','color','k');
  end
  function add_Earth(h,flag)
    if nargin == 1
      flag='terminator';
    end
    switch flag
      case 'terminator'
        theta=0:pi/20:pi;
        xEarth=sin(theta);yEarth=cos(theta);
        patch(-xEarth,yEarth,'k','edgecolor','none','parent',h)
        patch(xEarth,yEarth,'w','edgecolor','k','parent',h)
      case 'day'
        theta=0:pi/20:2*pi;
        xEarth=sin(theta);yEarth=cos(theta);
        patch(xEarth,yEarth,'w','edgecolor','k','parent',h)
      otherwise
        error('unknown flag');
    end
  end
  function y=get_in_lmn(x)
    % get x in LMN reference frame
    data=get(gcf,'userdata');
    y=x; % intialize
    Lstr=get(data.L_hndl, 'string');
    Mstr=get(data.M_hndl, 'string');
    Nstr=get(data.N_hndl, 'string');
    L=eval(Lstr);
    M=eval(Mstr);
    N=eval(Nstr);
    if numel(L)~=3, L=0; end
    if numel(M)~=3, M=0; end
    if numel(N)~=3, N=0; end
    if numel(L)+numel(M)+numel(N) <= 5 % not enough information
      irf.log('critical','LMN not correctly defined');
      return
    end
    data.Lstr=Lstr;
    data.Mstr=Mstr;
    data.Nstr=Nstr;
    for iSc=1:numel(data.sc_list)
      y{iSc}=irf_newxyz(x{iSc},L,M,N);
    end
    if L==0, data.Lstr=['[' num2str(irf_norm(cross(M,N)),'%6.2f') ']']; set(data.L_hndl,'string',data.Lstr);end
    if M==0, data.Mstr=['[' num2str(irf_norm(cross(N,L)),'%6.2f') ']']; set(data.M_hndl,'string',data.Mstr);end
    if N==0, data.Nstr=['[' num2str(irf_norm(cross(L,M)),'%6.2f') ']']; set(data.N_hndl,'string',data.Nstr);end
    set(gcf,'userdata',data);
  end
end
