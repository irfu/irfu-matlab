function irf_figmenu(action)
%IRF_FIGMENU   Add to the current figure a menu with some useful commands
%
% $Id$

if nargin < 1, action = 'initialize'; end

if strcmp(action,'initialize'),
    % execute irf_figmenu if there is no such menu
    if isempty(findobj(gcf,'type','uimenu','label','&irf'))
        hfigmenu=uimenu('Label','&irf');
        uimenu(hfigmenu,'Label','&Update time axis','Callback','add_timeaxis(gca,''date'')','Accelerator','t')
        uimenu(hfigmenu,'Label','Fit &Y axis','Callback','set(gca,''YLimMode'',''auto'')','Accelerator','y')
        uimenu(hfigmenu,'Label','&irf_tm','Callback','irf_figmenu(''irf_tm'')','Accelerator','i')
        uimenu(hfigmenu,'Label','Pointer &Crosshair','Callback','set(gcbf,''pointer'',''fullcrosshair'')')
        uimenu(hfigmenu,'Label','&Pointer Arrow','Callback','set(gcbf,''pointer'',''arrow'')')
        uimenu(hfigmenu,'Label','&Align axis','Callback','irf_plot_axis_align','Accelerator','a')
        hmenu_zoom=uimenu(hfigmenu,'Label','Zoom all panels OFF','Callback','irf_figmenu(''zoomall'')');
        set(zoom(gcf),'ActionPostCallback', @adaptiveDateTicks);
        user_data = get(gcf,'userdata');
        user_data.irf_figmenu=1;
        user_data.hmenu_zoom=hmenu_zoom;
        set(gcf,'userdata',user_data);
    end
elseif strcmpi(action,'irf_tm'),
    user_data = get(gcf,'userdata');
    user_data.subplot_handles=irf_plot_get_subplot_handles;
    set(gcf,'userdata',user_data);
    irf_tm(user_data.subplot_handles);
elseif strcmpi(action,'zoomall'),
    hmenu_zoom = getfield(get(gcf,'userdata'),'hmenu_zoom');
    lab=get(hmenu_zoom,'label');
    if strcmpi(lab,'Zoom all panels OFF'), % switch off auto zoom
        set(zoom(gcf),'ActionPostCallback', '');
        set(hmenu_zoom,'label','Zoom all panels ON');
    else % switch on auto zoom
        set(hmenu_zoom,'label','Zoom all panels OFF');
        set(zoom(gcf),'ActionPostCallback', @adaptiveDateTicks);
    end
end

function adaptiveDateTicks(figureHandle,eventObjectHandle)
% Resetting x axis to automatic tick mark generation
hsubplots=irf_plot_get_subplot_handles(figureHandle);
xlim=get(eventObjectHandle.Axes,'xlim');
set(hsubplots,'xlim',xlim);
add_timeaxis(hsubplots);