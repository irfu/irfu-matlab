function irf_figmenu(action,x1)
%IRF_FIGMENU   Add to the current figure a menu with useful commands
%

if nargin < 1, action = 'initialize'; end
if nargin < 2, x1 = ''; end % x1 is used when printing figures

switch lower(action)
  case 'initialize'
    % execute irf_figmenu if there is no such menu
    if isempty(findobj(gcf,'type','uimenu','label','&irf'))
      hfigmenu=uimenu('Label','&irf');
      uimenu(hfigmenu,'Label','&Update time axis','Callback','irf_timeaxis(gca,''date'')','Accelerator','t')
      uimenu(hfigmenu,'Label','Fit &Y axis','Callback','set(gca,''YLimMode'',''auto'')','Accelerator','y')
      uimenu(hfigmenu,'Label','&irf_tm','Callback','irf_figmenu(''irf_tm'')','Accelerator','i')
      %			not anymore supported in MATLAB 2015
      %			uimenu(hfigmenu,'Label','Pointer &Crosshair','Callback','set(gcbf,''pointer'',''fullcrosshair'')') %
      %			not needed as long as full cross hair is not supported
      %			uimenu(hfigmenu,'Label','&Pointer Arrow','Callback','set(gcbf,''pointer'',''arrow'')')
      uimenu(hfigmenu,'Label','&Align axis','Callback','irf_plot_axis_align','Accelerator','a')
      hmc  = uimenu(hfigmenu,'Label','&Cluster');
      hmcp = uimenu(hmc,'Label','satellite position','Callback','irf_figmenu(''c_position'')');
      hmct = uimenu(hmc,'Label','4 s/c timing','Callback','irf_figmenu(''c_timing'')');
      hmcm = uimenu(hmc,'Label','MVA','Callback','irf_figmenu(''mva'')');

      hmenu_zoom=uimenu(hfigmenu,'Label','Zoom all panels OFF','Callback','irf_figmenu(''zoomall'')');
      hmenu_zoom=uimenu(hfigmenu,'Label','Zoom &back','Callback','irf_figmenu(''zoomback'')');
      hmenu_zoom=uimenu(hfigmenu,'Label','Zoom &whole interval','Callback','irf_figmenu(''zoom_whole_interval'')');
      hmenu_print=uimenu(hfigmenu,'Label','Export figure as...');
      set(zoom(gcf),'ActionPostCallback', @adaptiveDateTicks);
      set(zoom(gcf),'ActionPreCallback', @saveZoomStack);
      user_data = get(gcf,'userdata');
      user_data.irf_figmenu=1;
      user_data.hmenu_zoom=hmenu_zoom;

      fileTypes = {'eps','png'};
      hmenu_print2 = gobjects(size(fileTypes));
      for i = 1:length(fileTypes)
        hmenu_print2(i) = uimenu(hmenu_print,'Label',fileTypes{i});
      end
      set(hmenu_print2(1),'Callback','irf_figmenu(''print_fig'',''eps'')')
      set(hmenu_print2(2),'Callback','irf_figmenu(''print_fig'',''png'')')

      set(gcf,'userdata',user_data);
    end
    % reset zoomStack
    ud = get(gcf,'userdata');
    ud.zoomStack = [];
    set(gcf,'userdata',ud);
  case 'irf_tm'
    user_data = get(gcf,'userdata');
    if ~isfield(user_data,'suplot_handles')
      user_data.subplot_handles=irf_plot_get_subplot_handles;
      set(gcf,'userdata',user_data);
    end
    irf_tm(user_data.subplot_handles);
  case 'zoomall'
    hmenu_zoom = getfield(get(gcf,'userdata'),'hmenu_zoom');
    lab=get(hmenu_zoom,'label');
    if strcmpi(lab,'Zoom all panels OFF') % switch off auto zoom
      set(zoom(gcf),'ActionPostCallback', '');
      set(hmenu_zoom,'label','Zoom all panels ON');
    else % switch on auto zoom
      set(hmenu_zoom,'label','Zoom all panels OFF');
      set(zoom(gcf),'ActionPostCallback', @adaptiveDateTicks);
    end
  case 'zoomback'
    ud = get(gcf,'userdata');
    if ~isfield(ud,'zoomStack') || isempty(ud.zoomStack)
      return
    end
    xlim = ud.zoomStack(end,:);
    figure_zoom_in(gcf,xlim)
    ud.zoomStack(end,:)=[];
    set(gcf,'userdata',ud);
  case 'zoom_whole_interval'
    ud = get(gcf,'userdata');
    if ~isfield(ud,'zoomStack') || isempty(ud.zoomStack)
      return
    end
    xlim = ud.zoomStack(1,:);
    figure_zoom_in(gcf,xlim)
    ud.zoomStack=[];
    set(gcf,'userdata',ud);
  case 'c_position'
    xlim=get(gca,'xlim');
    tStart= getfield(get(gcf,'userdata'),'t_start_epoch');
    c_pl_sc_conf_xyz(tStart+xlim(1));
  case 'c_timing'
    xlim=get(gca,'xlim');
    tStart= getfield(get(gcf,'userdata'),'t_start_epoch');
    tint = tStart + xlim;
    c_eval('v?=[];h?=findobj(gca,''Tag'',''C?'');');
    c_eval('if h?, v?     =get(h?,''XData'')+tStart;v?=v?(:); end;');
    c_eval('if h?, v?(:,2)=get(h?,''YData'');                 end;');
    h=c_4_v_gui(v1,v2,v3,v4,2);
    irf_zoom(h,'x',tint);
  case 'mva'
    xlim=get(gca,'xlim');
    tStart= getfield(get(gcf,'userdata'),'t_start_epoch');
    v=[];h=findobj(gca,'type','line');
    for jh = 1:numel(h)
      t = get(h(jh),'xdata') + tStart;
      y = get(h(jh),'ydata');
      if jh == 1
        v = [t(:) y(:)];
      else
        if all(t(:) == v(:,1))
          v = [v y(:)];
        else
          irf_log('dsrc','Lines do not have the same time sampling!');
        end
      end
    end
    irf_minvar_gui(v);
  case 'print_fig'
    fileType = x1;
    fileName = inputdlg('File name:',['Export as ',fileType,'...'],1,{'fig'});
    fileName = fileName{1};
    irf_print_fig(fileName,fileType)
end

  function adaptiveDateTicks(figureHandle,eventObjectHandle)
    % Resetting x axis to automatic tick mark generation
    xlim=get(eventObjectHandle.Axes,'xlim');
    figure_zoom_in(figureHandle,xlim,eventObjectHandle.Axes);
  end

  function saveZoomStack(figureHandle,eventObjectHandle)
    % save zoom levels for zooming out
    ud=get(figureHandle,'userdata');
    xlim=get(eventObjectHandle.Axes,'xlim');
    if isfield(ud,'zoomStack')
      ud.zoomStack(end+1,1:2) = xlim;
    else
      ud.zoomStack(1,1:2) = xlim;
    end
    set(figureHandle,'userdata',ud);
  end

  function figure_zoom_in(figureHandle,xlim,varargin)
    % figure_zoom_in(figureHandle,xlim,[currentAxisHandle])
    if nargin == 3
      currentAxisHandle = varargin{1};
    else
      currentAxisHandle = NaN;
    end
    ud=get(figureHandle,'userdata');
    hsubplots=irf_plot_get_subplot_handles(figureHandle);
    if isfield(ud,'t_start_epoch')
      irf_zoom(hsubplots,'x',xlim+double(ud.t_start_epoch));
    else
      set(hsubplots,'xlim',xlim);
    end
    set(figureHandle,'userdata',ud);
    for ih=1:numel(hsubplots)
      h=hsubplots(ih);
      if h~=currentAxisHandle
        % if not spectrogram do also smart y adjustment
        if ~isempty(findobj(h,'tag','irf_pl_mark')) || ...
            ~any(~isempty([findobj(h,'Type','surface') ...
            findobj(h,'Type','patch')]))
          irf_zoom(h,'y');
        end
      end
    end
    irf_timeaxis(hsubplots);
  end
end