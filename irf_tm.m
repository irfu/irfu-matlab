function status = irf_tm(h)
% IRF_TM  change the time intervals for plots with several subplots
%   h is the array with all subplot graphics handles
%

if nargin == 0, h=gca; end

tlim = [];

% use the first subplot to estimate available time interval, also check the
% user_data.t_start_epoch of the figure
  hh=h(1,1);  
  xl=get(hh,'XLim');
  hc=get(hh,'Children');
  xd=get(hc(end),'XData');
  avail=[min([xl xd]) max([xl xd])];
  presel=xl;
  
  % read t_start_epoch from figures userdata field and update time
  % intervals
  user_data=get(gcf,'userdata');
  if isfield(user_data,'t_start_epoch')
      avail=avail+user_data.t_start_epoch;
      presel=presel+user_data.t_start_epoch;
  end

  dt = 0.02*diff(avail);
  xlim = [avail(1)-dt avail(2)+dt];
  ttics = timeaxis(xlim);

% initialize dgud variable that has all the properties of time manager
dgud=get(gcf,'userdata'); % take existing values in user_data and add necessary
dgud.figure=get(h(1),'parent'); % add the number of figure which is controled by irf_tm
dgud.autoY=1; % default is always automatically fix Y component
dgud.tlim = avail;
dgud.from = 1;
dgud.cancel = 0;
dgud.h=h;

screensize = get(0, 'screensize');
xsize = 460;
ysize = 150;
left = 0.5*screensize(3)-xsize/2;
bottom = 0.5*screensize(4)-ysize/2;
dgh = figure('name', 'time limits', ...
  'position', [left bottom xsize ysize], ...
  'pointer', 'left', ...
  'defaultaxestickdir', 'out', ...
  'defaultaxestickdirmode', 'manual', ...
  'handlevisibility', 'on');
storecf = get(0, 'currentfigure');
set(0, 'currentfigure', dgh);
axh = axes('units', 'pixel', ...
    'position', [5 ysize-15 xsize-10 10], ...
    'box', 'on', ...
    'ticklength', [0.005 0.025], ...
    'xticklabelmode', 'manual', ...
    'xtick', ttics{1}, ...
    'xticklabel', ttics{2}, ...
    'xlimmode', 'manual', ...
    'xlim', xlim, ...
    'ytick', [], ...
    'ylim', [-5 5], ...
    'buttondownfcn', 'irf_fromto(''ax'')');
xlabel('UT');
dgud.lnh = line('color', [1 0 0], ...
    'linewidth', 6, ...
    'xdata', presel, ...
    'ydata', [0 0], ...
    'buttondownfcn', 'irf_fromto(''ax'')');

dgud.helph = uicontrol('style', 'text', ...
    'string', 'Click on axis selects ''From'' time');
ext = get(dgud.helph, 'extent');
set(dgud.helph, 'position', [5 ysize-65 ext(3:4)]);

uch = uicontrol('style', 'pushbutton', ...
    'string', 'toggle', ...
    'callback', 'irf_fromto(''toggle'')');
ext = get(uch, 'extent');
set(uch, 'position', [250-ext(3)/2 ysize-78 ext(3:4)]);

uch1 = uicontrol('style', 'text', 'string', 'From:');
ext1 = get(uch1, 'extent');
set(uch1, 'position', [5 ysize-100 ext1(3:4)]);
x0 = 5 + ext1(3);
dgud.fromh = uicontrol('style', 'edit', ...
      'string', epoch2iso(presel(1)), ...
    'callback', 'irf_fromto(''from'')', ...
    'backgroundcolor','white');
ext2 = get(dgud.fromh, 'extent');
set(dgud.fromh, 'position', [x0 ysize-98 ext2(3)+50 ext(4)]);

uch1 = uicontrol('style', 'text', 'string', 'Step:');
extstep=get(uch1,'extent');
dgud.step = uicontrol('style', 'edit', ...
    'string', datestr(datenum(fromepoch(diff(presel))), 13), ...
    'callback', 'irf_fromto(''step'')','backgroundcolor','white');
set(uch1, 'position', [x0+ext2(3)+55 ysize-100 extstep(3) ext(4)]);
set(dgud.step, 'position', [x0+ext2(3)+extstep(3)+55 ysize-100 ext2(3) ext(4)]);

uch1 = uicontrol('style', 'text', 'string', 'To:');
dgud.toh = uicontrol('style', 'edit', ...
    'string', epoch2iso(presel(2)), ...
    'callback', 'irf_fromto(''to'')','backgroundcolor','white');
set(uch1, 'position', [ext1(1) 30 ext1(3:4)]);
set(dgud.toh, 'position', [x0 30 ext2(3)+50 ext(4)]);

ypos=30;xpos=x0+ext2(3)+55;
uch1 = uicontrol('style', 'pushbutton','string', ' < ');ext=get(uch1,'extent');
set(uch1, 'position', [xpos ypos ext(3) ext(4)],'callback','irf_fromto(''prev'')');xpos=xpos+ext(3);
uch1 = uicontrol('style', 'pushbutton','string', ' Update ');ext=get(uch1,'extent');
set(uch1, 'position', [xpos ypos ext(3) ext(4)],'callback','irf_fromto(''update'')');xpos=xpos+ext(3);
uch1 = uicontrol('style', 'pushbutton','string', ' > ');ext=get(uch1,'extent');
set(uch1, 'position', [xpos ypos ext(3) ext(4)],'callback','irf_fromto(''next'')');xpos=xpos+ext(3);
uch1 = uicontrol('style', 'pushbutton','string', ' << >> ');ext=get(uch1,'extent');
set(uch1, 'position', [xpos+10 ypos ext(3) ext(4)],'callback','irf_fromto(''all_interval'')');xpos=xpos+ext(3)+10;
uch1 = uicontrol('style', 'pushbutton','string', ' 1 ');ext=get(uch1,'extent');
set(uch1, 'position', [xpos+10 ypos ext(3) ext(4)],'callback','irf_fromto(''zoom_to_1'')');xpos=xpos+ext(3)+10;

uimenu('label','&Update','accelerator','u','callback','irf_fromto(''update'')');
uimenu('label','Auto &YLim','accelerator','y','callback','irf_fromto(''autoY'')');

ypos=ypos-20;xpos=20;
dgud.ylab = uicontrol(...
    'Style','popupmenu',...
    'BackgroundColor','white',...
    'String',{'enter','Vps [V]','NVps [cc]',...
        'B [nT] GSE;B_X;B_Y;B_Z;|B|','E [mV/m];E_X;E_Y;E_Z','|B| [nT]',...
               'Ex [mV/m] DSI','Ey [mV/m] DSI',...
               });
set(dgud.ylab, 'position', [xpos ypos 100 15]);xpos=xpos+100;

pop{1}='Ylabel panel'; for j=1:max(size(h)); pop{j+1}=num2str(j); end
dgud.ylabpanel = uicontrol(...
    'Style','popupmenu',...
    'BackgroundColor','white',...
    'String',pop,...
    'callback','irf_fromto(''ylabel'')');
set(dgud.ylabpanel, 'position', [xpos ypos 100 15]);xpos=xpos+110;

uch2 = uicontrol('style', 'pushbutton', ...
    'string', 'Cancel', ...
    'callback', 'irf_fromto(''cancel'')');
ext2 = get(uch2, 'extent');
set(uch2, 'position', [290 5 ext2(3:4)]);

set(dgh, 'userdata', dgud);

%close(dgh);
