function irf_fromto(fromto)
% UI_FROMTO - callback for the buttons in the time selection window
%   fromto - action to do
%

ud = get(gcbf, 'userdata'); % get userdata of time manager window

ud_fig=get(ud.figure,'userdata'); % get userdata of figures windows

% check if t_start_epoch is necessary
if isfield(ud_fig,'t_start_epoch')
  t_start_epoch=ud.t_start_epoch;
else
  t_start_epoch=0;
end

% check if autoY field exists, if not put it to 1
if ~isfield(ud_fig,'autoY')
  ud_fig.autoY=1;set(ud.figure,'userdata',ud_fig);
end


if isfield(ud_fig,'subplot_handles') % check the handles to subplots
  SUBPLOT_HANDLES=ud_fig.subplot_handles;
else
  SUBPLOT_HANDLES=get(ud.figure,'Children');
  % Do not assume that subplots are the only children of figure
  ii = [];
  for j=1:length(SUBPLOT_HANDLES)
    if ~isgraphics(SUBPLOT_HANDLES( j ),'axes'), ii = [ii j]; end
  end
  if ~isempty(ii), SUBPLOT_HANDLES(ii) = []; clear ii, end
end

switch fromto
  case 'cancel'
    uiresume;
  case 'ok'
    uiresume;
  case 'ax'
    tlim = get(ud.lnh, 'xdata');
    p = get(gca, 'currentpoint');
    if ud.from
      tlim(1) = max(ud.tlim(1), p(1));
      if tlim(1)>tlim(2), tlim(2)=tlim(1);end % start time should be smaller or equal than end time
      update_fromto(ud,tlim);
      set(gcbf, 'pointer', 'right');
      set(ud.helph, 'string', 'Click on axis selects ''To'' time');
      ext = get(ud.helph, 'extent');
      set(ud.helph, 'position', [395-ext(3) 55 ext(3:4)]);
      ud.from = 0;
    else
      tlim(2) = min(ud.tlim(2), p(1));
      if tlim(2)<tlim(1), tlim(1)=tlim(2);end % start time should be smaller or equal than end time
      update_fromto(ud,tlim);
      set(gcbf, 'pointer', 'left');
      set(ud.helph, 'string', 'Click on axis selects ''From'' time');
      ext = get(ud.helph, 'extent');
      set(ud.helph, 'position', [5 55 ext(3:4)]);
      ud.from = 1;
    end
  case 'from'
    [tlim, step]= get_fromto(ud);
    if tlim(1) > tlim(2)
      tlim(2)=tlim(1)+step;
    else
      tlim(2)=tlim(1)+step;
    end
    update_fromto(ud,tlim);
  case 'to'
    [tlim, step]= get_fromto(ud);
    if tlim(2)< tlim(1) % end time smaller than start time
      tlim(1)=tlim(2);
      step=0;
    else
      step=diff(tlim);
    end
    update_fromto(ud,tlim);
  case 'step'
    [tlim, step]=get_fromto(ud);
    tlim(2) = tlim(1)+step;
    update_fromto(ud,tlim);
  case 'update'
    [tlim, step]=get_fromto(ud);
    irf_zoom(SUBPLOT_HANDLES,'x',tlim);
    irf_timeaxis(SUBPLOT_HANDLES);
  case 'prev'
    [tlim, step]=get_fromto(ud);
    tlim = tlim-step;
    update_fromto(ud,tlim)
    irf_fromto('update')
  case 'next'
    [tlim, step]=get_fromto(ud);
    tlim = tlim+step;
    update_fromto(ud,tlim)
    irf_fromto('update')
  case 'all_interval'
    hc=get(SUBPLOT_HANDLES(1,1),'Children'); % use the first subplot to estimate available time interval
    %xd=get(hc(1),'XData');minx=min(xd);maxx=max(xd);clear xd;
    %if length(hc)>1, for ii=2:length(hc), xd=get(hc(1),'XData');minx=min([minx xd]);maxx=max([maxx xd]);clear xd;end;end
    minx=[];maxx=[];
    for ii=1:length(hc)
      if isgraphics(hc( ii ),'axes') || isgraphics(hc( ii ),'line')
        minx=[minx min(get(hc(ii),'xdata'))];
        maxx=[maxx max(get(hc(ii),'xdata'))];
      end
    end
    minx=min(minx);maxx=max(maxx);
    avail=[minx maxx]+t_start_epoch;
    tlim=avail;step=diff(tlim);
    update_fromto(ud,tlim);
    irf_fromto('update')
  case 'zoom_to_1'
    xl=get(SUBPLOT_HANDLES(1,1),'XLim'); % use the first subplot to estimate available time interval
    tlim=[xl(1) xl(2)]+t_start_epoch;step=diff(tlim);
    update_fromto(ud,tlim);
  case 'autoY'
    for h=SUBPLOT_HANDLES
      set(h,'YLimMode','auto');
    end
  case 'ylabel'
    h_select=get(ud.ylabpanel,'Value')-1;
    i_select=get(ud.ylab,'Value');
    strings=get(ud.ylab,'String');
    label_and_legend=tokenize(strings{i_select},';');
    label=label_and_legend{1};
    if strcmp(label,'enter'); label=input('Input Y label>','s');end
    axes(ud.h(h_select));
    ylabel(label,'verticalalignment','bottom');
    if size(label_and_legend,2)>1, legend(label_and_legend(2:end));end
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
end
set(gcbf, 'userdata', ud);

if ud.autoY==1
  for h=SUBPLOT_HANDLES
    set(h,'YLimMode','auto');
  end
end

end

function update_fromto(ud,tlim)
% update input strings
t_from_str = epoch2iso(tlim(1));
t_from_str = strrep(t_from_str(1:end-3), 'T', '_');
t_to_str = epoch2iso(tlim(2));
t_to_str = strrep(t_to_str(1:end-3), 'T', '_');
set(ud.fromh, 'string', t_from_str);
set(ud.toh, 'string', t_to_str);
set(ud.step, 'string', regexp(epoch2iso(diff(tlim)),'\d+:\d+:\d+\.\d{4}','match'));
% update the red line showing the data interval
set(ud.lnh, 'xdata', tlim);
end

function [tlim, step]=get_fromto(ud)
xx=get(ud.fromh, 'string');tstr=[strrep(xx,'_','T') 'Z'];
tlim(1) = iso2epoch(tstr);
xx=get(ud.toh, 'string');tstr=[strrep(xx,'_','T') 'Z'];
tlim(2) = iso2epoch(tstr);
xx=get(ud.step, 'string');
step = 86400*(datenum(xx)-datenum('00:00:00'));
%step = (datenum(get(ud.step, 'string'))-datenum('00:00:00'))*86400;
end
