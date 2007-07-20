function s=irf_figmenu(action)
%IRF_FIGMENU   Add to the current figure a menu with some useful commands
%
% $Id$

if nargin < 1, action = 'initialize'; end

if strcmp(action,'initialize'),
	% execute irf_figmenu if there is no such menu
	need_init = 1;
	ch = get(gcf,'children');
	if ~isempty(ch)
		for j=1:length(ch)
			if strcmp(get(ch(j),'Type'),'uimenu')
				if strcmp(get(ch(j),'Label'),'&irf'), need_init = 0; break, end
			end
		end
	end

	if need_init
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
  h=findobj(gcf,'type','axes','-not','tag','Colorbar');
  hmax=1;
  for ih=1:length(h),
    ax=get(h(ih),'position');
    axy(ih)=ax(2);axx(ih)=ax(1);
  end
  ind_ax=find(axx<0.2); % find 
  hax=h(ind_ax);
  [xx,ind]=sort(axy(ind_ax));
  ind=fliplr(ind);
  hh=hax(ind);
		user_data = get(gcf,'userdata');
		user_data.subplot_handles=hh;
		set(gcf,'userdata',user_data);
  irf_tm(hh);
end

