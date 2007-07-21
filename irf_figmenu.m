function s=irf_figmenu(action)
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
		user_data = get(gcf,'userdata');
		user_data.irf_figmenu=1;
		set(gcf,'userdata',user_data);
	end
elseif strcmp(action,'irf_tm'),
    user_data = get(gcf,'userdata');
    user_data.subplot_handles=irf_plot_get_subplot_handles;
    set(gcf,'userdata',user_data);
    irf_tm(user_data.subplot_handles);
end

