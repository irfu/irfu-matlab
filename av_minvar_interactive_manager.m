%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
function status = av_minvar_interactive_manager(arg)
%
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
global ud
persistent tlim
if gcbf,
  ud = get(gcbf, 'userdata');
else
  par=get(gca,'parent');ud = get(par, 'userdata');
end

switch arg
  case 'cancel'
    ud.cancel = 1; set(gcbf, 'userdata', ud), uiresume;
  case 'ok'
    uiresume;
  case 'ax'
    tlim = get(ud.mvar_intervals, 'xdata'); tlim(3:4)=[];
    p = get(gca, 'currentpoint');
    tlim_interval=get(gca,'xlim');
    if ud.from
      tlim(1) = max(tlim_interval(1), p(1));
      tlim(2) = max(p(1),tlim(2));
      set(ud.fromtext,'backgroundcolor','w');
      set(ud.totext,'backgroundcolor','r');
%      set(gcbf, 'pointer', 'right');
      ud.from = 0;
    else
      tlim(2) = min(tlim_interval(2), p(1));
      tlim(1) = min(tlim(1), p(1));
      set(ud.totext,'backgroundcolor','w');
      set(ud.fromtext,'backgroundcolor','r');
%      set(gcbf, 'pointer', 'left');
      ud.from = 1;
    end
    strfrom=datestr(datenum(fromepoch(tlim(1))), 0);
    set(ud.fromh, 'string', strfrom);
    strto=datestr(datenum(fromepoch(tlim(2))), 0);
    set(ud.toh, 'string', strto);
    set(ud.mvar_intervals,'xdata',[tlim(1) tlim(2) tlim(2) tlim(1)]);
    set(gcbf, 'userdata', ud);
  case 'from'
    tlim(1) = toepoch(datevec(strrep(get(ud.fromh, 'string'),'_',' ')));
    tlim(2) = toepoch(datevec(strrep(get(ud.toh, 'string'),'_',' ')));
    set(ud.mvar_intervals,'xdata',[tlim(1) tlim(2) tlim(2) tlim(1)]);
  case 'to'
    tlim(1) = toepoch(datevec(strrep(get(ud.fromh, 'string'),'_',' ')));
    tlim(2) = toepoch(datevec(strrep(get(ud.toh, 'string'),'_',' ')));
    set(ud.mvar_intervals,'xdata',[tlim(1) tlim(2) tlim(2) tlim(1)]);
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

% update the minvar plot
X=av_t_lim(ud.X,tlim);
clear ud.Xminvar;
[ud.Xminvar, l, v]=av_minvar(X);
ud.l=l;ud.v=v;ud.v1=v(1,:);ud.v2=v(2,:);ud.v3=v(3,:);
axes(ud.h(2));
av_tplot([ud.Xminvar X(:,5)]);
set(ud.h(2),    'buttondownfcn', 'av_minvar_interactive_manager(''ax'')');
axis tight;add_timeaxis(ud.h(2),'date');
legend('max','interm','min','abs');
axes(ud.h(3));
plot(ud.Xminvar(:,4),ud.Xminvar(:,2));xlabel('min');ylabel('max');
axis tight;axis equal; ax=axis;grid on;
axes(ud.h(4))
plot(ud.Xminvar(:,3),ud.Xminvar(:,2));xlabel('interm');ylabel('max');
axis equal; grid on;
l_str=['L1=' num2str(l(1),3) ' L2=' num2str(l(2),3) ' L3=' num2str(l(3),3) '\newline'];
lratio_str=['L1/L2=' num2str(l(1)/l(2),2) ' L2/L3=' num2str(l(2)/l(3),2) '\newline'];
v1_str=['v1=[' num2str(v(1,:),'%6.2f') '] \newline'];
v2_str=['v2=[' num2str(v(2,:),'%6.2f') '] \newline'];
v3_str=['v3=[' num2str(v(3,:),'%6.2f') '] \newline'];
v_str=[v1_str v2_str v3_str];
set(ud.result_text,'string',[l_str lratio_str v_str],'verticalalignment','top');

