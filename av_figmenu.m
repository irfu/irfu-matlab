function s=av_figmenu(action)
% add to the current figure a menu with some useful commands
warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'irf_figmenu')

if nargin < 1, action = 'initialize'; end

if strcmp(action,'initialize'),
  hfigmenu=uimenu('Label','&av');
  uimenu(hfigmenu,'Label','&Update time axis','Callback','add_timeaxis(gca,''date'')','Accelerator','t')
  uimenu(hfigmenu,'Label','Fit &Y axis','Callback','set(gca,''YLimMode'',''auto'')','Accelerator','y')
  uimenu(hfigmenu,'Label','&irf_tm','Callback','av_figmenu(''irf_tm'')','Accelerator','i')
  uimenu(hfigmenu,'Label','Pointer &Crosshair','Callback','set(gcbf,''pointer'',''fullcrosshair'')')
  uimenu(hfigmenu,'Label','&Pointer Arrow','Callback','set(gcbf,''pointer'',''arrow'')')
  user_data=get(gcf,'userdata'); user_data.av_figmenu=1; set(gcf,'userdata',user_data);

elseif strcmp(action,'irf_tm'),
  h=findobj(gcf,'type','axes');
  hmax=1;
  for ih=1:length(h),
    ax=get(h(ih),'position');
    axy(ih)=ax(2);
  end
  [xx,ind]=sort(axy);
  ind=fliplr(ind);
  hh=h(ind);
  irf_tm(hh);
end

