function ui_fromto(fromto)
% UI_FROMTO - callback for the buttons in the time selection window
%   fromto - action to do

global SUBPLOT_HANDLES

ud = get(gcbf, 'userdata');

switch fromto
  case 'cancel'
    ud.cancel = 1; set(gcbf, 'userdata', ud), uiresume;
  case 'ok'
    uiresume;
  case 'ax'
    tlim = get(ud.lnh, 'xdata');
    p = get(gca, 'currentpoint');
    if ud.from
      tlim(1) = max(ud.tlim(1), p(1));
      set(ud.lnh, 'xdata', tlim);
      set(gcbf, 'pointer', 'right');
      strfrom=datestr(datenum(fromepoch(tlim(1))), 0);
      set(ud.fromh, 'string', strfrom);
      set(ud.helph, 'string', 'Click on axis selects ''To'' time');
      ext = get(ud.helph, 'extent');
      set(ud.helph, 'position', [395-ext(3) 55 ext(3:4)]);
      ud.from = 0;
    else
      tlim(2) = min(ud.tlim(2), p(1));
      set(ud.lnh, 'xdata', tlim);
      set(gcbf, 'pointer', 'left');
      strto=datestr(datenum(fromepoch(tlim(2))), 0);
      set(ud.toh, 'string', strto);
      set(ud.helph, 'string', 'Click on axis selects ''From'' time');
      ext = get(ud.helph, 'extent');
      set(ud.helph, 'position', [5 55 ext(3:4)]);
      ud.from = 1;
    end
    set(gcbf, 'userdata', ud);
  case 'from'
    tlim(1) = toepoch(datevec(strrep(get(ud.fromh, 'string'),'_',' ')));
    step = datenum(get(ud.step, 'string'))*86400;
    tlim(2)=tlim(1)+step;
    set(ud.lnh, 'xdata', tlim);
    set(ud.toh, 'string', strrep(datestr(datenum(fromepoch(tlim(2))), 0),' ','_'));
  case 'to'
    tlim = get(ud.lnh, 'xdata');
    tlim(2) = toepoch(datevec(strrep(get(ud.toh, 'string'),'_',' ')));
    set(ud.lnh, 'xdata', tlim);
    set(ud.step, 'string', datestr(diff(tlim)/86400 ,13));
  case 'step'
    step = datenum(get(ud.step, 'string'))*86400;
    tlim = get(ud.lnh, 'xdata');
    tlim(2) = tlim(1)+step;
    set(ud.lnh, 'xdata', tlim);
    set(ud.toh, 'string', strrep(datestr(datenum(fromepoch(tlim(2))), 0),' ','_'));
  case 'update'
    step = datenum(get(ud.step, 'string'))*86400;
    tlim = get(ud.lnh, 'xdata');
    av_zoom(tlim,'x',SUBPLOT_HANDLES);
    add_timeaxis(SUBPLOT_HANDLES);
  case 'prev'
    step = datenum(get(ud.step, 'string'))*86400;
    tlim = get(ud.lnh, 'xdata');
    tlim = tlim-step;
    set(ud.lnh, 'xdata', tlim);
    set(ud.fromh, 'string', strrep(datestr(datenum(fromepoch(tlim(1))), 0),' ','_'));
    set(ud.toh, 'string', strrep(datestr(datenum(fromepoch(tlim(2))), 0),' ','_'));
    av_zoom(tlim,'x',SUBPLOT_HANDLES);
    add_timeaxis(SUBPLOT_HANDLES);
  case 'next'
    step = datenum(get(ud.step, 'string'))*86400;
    tlim = get(ud.lnh, 'xdata');
    tlim = tlim+step;
    set(ud.lnh, 'xdata', tlim);
    set(ud.fromh, 'string', strrep(datestr(datenum(fromepoch(tlim(1))), 0),' ','_'));
    set(ud.toh, 'string', strrep(datestr(datenum(fromepoch(tlim(2))), 0),' ','_'));
    av_zoom(tlim,'x',SUBPLOT_HANDLES);
    add_timeaxis(SUBPLOT_HANDLES);
  case 'all_interval'
    hh=SUBPLOT_HANDLES(1,1);  % use the first subplot to estimate available time interval
    xl=get(hh,'XLim');
    hc=get(hh,'Children');
    xd=get(hc(1),'XData');
    avail=[min([xl xd]) max([xl xd])];
    tlim=avail;step=diff(tlim);
    set(ud.lnh, 'xdata', tlim);
    set(ud.fromh, 'string', strrep(datestr(datenum(fromepoch(tlim(1))), 0),' ','_'));
    set(ud.toh, 'string', strrep(datestr(datenum(fromepoch(tlim(2))), 0),' ','_'));
    set(ud.step, 'string', datestr(diff(tlim)/86400 ,13));
    av_zoom(tlim,'x',SUBPLOT_HANDLES);
    add_timeaxis(SUBPLOT_HANDLES);
  case 'zoom_to_1'
    hh=SUBPLOT_HANDLES(1,1);  % use the first subplot to estimate available time interval
    xl=get(hh,'XLim');
    tlim=[xl(1) xl(2)];step=diff(tlim);
    set(ud.lnh, 'xdata', tlim);
    set(ud.fromh, 'string', strrep(datestr(datenum(fromepoch(tlim(1))), 0),' ','_'));
    set(ud.toh, 'string', strrep(datestr(datenum(fromepoch(tlim(2))), 0),' ','_'));
    set(ud.step, 'string', datestr(diff(tlim)/86400 ,13));
    av_zoom(tlim,'x',SUBPLOT_HANDLES);
    add_timeaxis(SUBPLOT_HANDLES);
  case 'autoY'
    for h=SUBPLOT_HANDLES,
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
