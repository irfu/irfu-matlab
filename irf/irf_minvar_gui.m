function h=irf_minvar_gui(x,column)
%IRF_MINVAR_GUI interactively do the minimum variance analysis
%
% IRF_MINVAR_GUI(X,COLUMN)
%  X - vector to use, [x(:,column(1)) x(:,column(2)) x(:,column(3))]
%  COLUMN - which columns to use, if not given use 2,3,4
%
% You can access the results through variable 'ud' that is defined as global
% or in the figures data 'userdata'.
% data=get(gcf,'userdata');data.ud
%
% ud.l - eigenvalues  ud.l(1), ud.l(2),ud.l(3)
% ud.v - eigenvectors (ud.v(1,:), ..), also ud.v1, ud.v2. ud.v3
% ud.Xminvar - data in minimum variance coordinates
%
% See also IRF_MINVAR

global ud
persistent message;

if isempty(message) % run only the first time during the session
  message='You can anytime access all the results from the variable "ud" or from get(gcf,''userdata'').';
  disp(message);
end

if      nargin < 1, help irf_minvar_gui;return;
elseif  (nargin==1 && ischar(x)), action=x;%disp(['action=' action]);
elseif  isnumeric(x)
  if size(x,2)<3, disp('Vector has too few components');return;end
  if nargin < 2
    if size(x,2)==3, column=[1 2 3];end
    if size(x,2)>3, column=[2 3 4];end
  end
  action='initialize';
elseif isa(x,'TSeries') && x.tensorOrder == 1 && size(x.data,2)==3
  action='initialize';
  xNew = [x.time.epochUnix double(x.data)];
  column = [2 3 4];
  x = xNew;
end

switch action
  case 'initialize'
    % X is used for minimum variance estimates
    evalin('base','clear ud; global ud;');

    if min(column)==1, time_vector=1:size(x,1);
    elseif min(column)>1, time_vector=x(:,1);
    end

    X=[time_vector x(:,column)];X=irf_abs(X);
    irf_plot(1,'newfigure');
    h(1)=subplot(4,1,1);
    set(h(1),'outerposition',[0 0.75 1 0.25]);
    irf_plot(h(1),X);
    axis(h(1),'tight');
    zoom(h(1),'off');
    ud=get(gcf,'userdata');
    if isfield(ud,'t_start_epoch'), ud.t0=ud.t_start_epoch;else, ud.t0=0; end

    ud.X=X;
    ud.from = 1; % first click with mouse is 'from', second is 'to'
    ud.cancel = 0;
    ud.tlim = [min(X(:,1)) max(X(:,1))];
    ud.tlim_mva=ud.tlim+[-1 1]; % default tlim_mva includes all interval, add 1s to help later in program

    %  irf_pl_info(h(1),['irf\_minvar\_gui() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))]); % add information to the plot
    set(h(1),'layer','top');
    grid(h(1),'on');
    ax=axis(h(1));
    ud.patch_mvar_intervals=patch([ud.tlim(1) ud.tlim(2) ud.tlim(2) ud.tlim(1)]-ud.t0,[ax(3) ax(3) ax(4) ax(4)],[-1 -1 -1 -1],'y','parent',h(1));

    h(2)=subplot(4,1,2);set(h(2),'outerposition',[0 0.5 1 0.25]);
    irf_plot(h(2),X);
    axis(h(2),'tight');
    zoom(h(2),'off');

    h(3)=subplot(4,2,5);

    h(4)=subplot(4,2,6);

    ud.h=h;

    xp=0.2;yp=0.2;
    ud.fromtext=uicontrol('style', 'text', 'string', 'From:','units','normalized', 'position', [xp yp 0.1 0.04],'backgroundcolor','red');
    ud.fromh = uicontrol('style', 'edit', ...
      'string', irf_time(ud.tlim(1),'utc'), ...
      'callback', 'irf_minvar_gui(''from'')', ...
      'backgroundcolor','white','units','normalized','position', [xp+0.11 yp 0.25 0.05]);

    yp=0.15;
    ud.totext=uicontrol('style', 'text', 'string', 'To:','units','normalized', 'position', [xp yp 0.1 0.04],'backgroundcolor','white');
    ud.toh=uicontrol('style', 'edit', ...
      'string', irf_time(ud.tlim(2),'utc'), ...
      'callback', 'irf_minvar_gui(''from'')','backgroundcolor','white','units','normalized', 'position', [xp+0.11 yp 0.25 0.05]);


    xp=0.1;yp=0.1;
    uicontrol('style', 'text', 'string', 'Low pass filter f/Fs = ','units','normalized','position', [xp yp 0.2 0.04],'backgroundcolor','white');
    ud.filter = uicontrol('style', 'edit', ...
      'string', '1', ...
      'backgroundcolor','white','units','normalized','position', [xp+0.21 yp 0.1 0.05]);

    xp=0.1;yp=0.05;
    uicontrol('style', 'text', 'string', 'MVAR method','units','normalized','position', [xp yp 0.2 0.04],'backgroundcolor','white');
    ud.mvar_method_handle=uicontrol('Style', 'popup',...
      'String', 'Unconstrained|Constrained min(<Bn^2>)|Constrained <Bn>=0',...
      'backgroundcolor','white','units','normalized',...
      'Position', [xp+0.21 yp 0.25 0.04],...
      'Callback', @setmethod);
    ud.mvar_method='mvar'; % default method

    uimenu('label','&Recalculate','accelerator','r','callback','irf_minvar_gui(''mva'')');

    h(5)=subplot(4,2,8);
    axis(h(5),'off');
    irf_legend(0,['irf\_minvar\_gui() ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss"))],[0.02 0.02],'fontsize',8); % add information to the plot
    ud.result_text=text(0,0.8,'result','parent',h(5));

    set(gcf,'userdata',ud);

    irf_minvar_gui('from');
    fix_legends;
    fix_hittest(ud.h(1:2));

  case 'ax'
    ud=get(gcf,'userdata');
    tlim = get(ud.patch_mvar_intervals, 'xdata'); tlim=tlim(:)';tlim(3:4)=[];
    tlim=tlim+ud.t0;
    p = get(gca, 'currentpoint')+ud.t0;
    tlim_interval=get(gca,'xlim')+ud.t0;
    if ud.from
      tlim(1) = max(tlim_interval(1), p(1));
      tlim(2) = max(p(1),tlim(2));
      set(ud.fromtext,'backgroundcolor','w');
      set(ud.totext,'backgroundcolor','r');
      ud.from = 0;
    else
      tlim(2) = min(tlim_interval(2), p(1));
      tlim(1) = min(tlim(1), p(1));
      set(ud.totext,'backgroundcolor','w');
      set(ud.fromtext,'backgroundcolor','r');
      ud.from = 1;
    end
    set(ud.fromh, 'string', epoch2iso(tlim(1),1));
    set(ud.toh, 'string', epoch2iso(tlim(2),1));
    set(ud.patch_mvar_intervals,'xdata',[tlim(1) tlim(2) tlim(2) tlim(1)]-ud.t0);
    ud.tlim=tlim;
    set(gcf,'userdata',ud);
    irf_minvar_gui('update_mva_axis');
  case 'from'
    ud=get(gcf,'userdata');
    tlim(1) = irf_time(get(ud.fromh,'string'),'utc>epoch');
    tlim(2) = irf_time(get(ud.toh,'string'),'utc>epoch');
    set(ud.patch_mvar_intervals,'xdata',[tlim(1) tlim(2) tlim(2) tlim(1)]-ud.t0);
    ud.tlim=tlim;
    set(gcf,'userdata',ud);
    irf_minvar_gui('update_mva_axis');
  case 'update_mva_axis'
    ud=get(gcf,'userdata');
    if ud.tlim==ud.tlim_mva % plot first time after 'mva'
      irf_plot(ud.h(2),ud.Xminvar);
      axis(ud.h(2),'fill');
      axis(ud.h(2),'tight');
      irf_timeaxis(ud.h(2),'date');
      plot(ud.h(3),ud.Xminvar(:,4),ud.Xminvar(:,2));
      hold(ud.h(3),'on');plot(ud.h(3),ud.Xminvar(1,4),ud.Xminvar(1,2),'k+');hold(ud.h(3),'off');
      xlabel(ud.h(3),'min');ylabel(ud.h(3),'max');
      axis(ud.h(3),'equal');
      grid(ud.h(3),'on');
      plot(ud.h(4),ud.Xminvar(:,3),ud.Xminvar(:,2),'handlevisibility','off');
      hold(ud.h(4),'on');plot(ud.h(4),ud.Xminvar(1,3),ud.Xminvar(1,2),'k+');hold(ud.h(4),'off');
      xlabel(ud.h(4),'interm');
      ylabel(ud.h(4),'max');
      legend(ud.h(4),'start','location','best');
      axis(ud.h(4),'equal');
      grid(ud.h(4),'on');
    elseif (ud.tlim(1)>=ud.tlim_mva(1) && ud.tlim(2)<=ud.tlim_mva(2)) % zoom to something within tlim_mva
      irf_zoom(ud.h(2),'x',ud.tlim);
    else                   % zoom to interval outside mva
      X=irf_tlim(ud.X,ud.tlim);
      clear ud.Xminvar;
      ud.Xminvar=irf_newxyz(X,ud.v1,ud.v2,ud.v3);
      irf_plot(ud.h(2),ud.Xminvar);
      axis(ud.h(2),'tight');
      irf_timeaxis(ud.h(2),'date');
    end
    if (ud.tlim(1)<ud.tlim_mva(1) || ud.tlim(2)>ud.tlim_mva(2)) % if zooming outside tlim_mva mark mva interval
      set(ud.h(2),'layer','top');
      ax=axis(ud.h(2));
      grid(ud.h(2),'on');
      ud.mvar_interval_2nd=patch([ud.tlim_mva(1) ud.tlim_mva(2) ud.tlim_mva(2) ud.tlim_mva(1)],[ax(3) ax(3) ax(4) ax(4)],[-1 -1 -1 -1],'y','buttondownfcn', {@click_ax},'parent',ud.h(2));
    end
    set(gcf,'userdata',ud);
    fix_legends;
    fix_hittest(ud.h(1:2));

  case 'mva'
    ud=get(gcf,'userdata');
    ud.tlim_mva=ud.tlim;
    X = ud.X;
    if eval(get(ud.filter,'string'))<1
      Fs = 1/(X(2,1)-X(1,1));
      flim = Fs*eval(get(ud.filter,'string'));
      X = irf_tlim(X, ud.tlim + [-20/Fs 20/Fs]);
      X = irf_filt(X,0,flim,Fs,5);
    else
      if eval(get(ud.filter,'string'))>1, disp('f/Fs must be <1!!!'), end
      set(ud.filter,'string','1')
    end
    X = irf_tlim(X,ud.tlim);
    clear ud.Xminvar;
    [ud.Xminvar, l, v]=irf_minvar(X,ud.mvar_method);
    ud.l=l;ud.v=v;ud.v1=v(1,:);ud.v2=v(2,:);ud.v3=v(3,:);
    l_str=['L1=' num2str(l(1),3) ' L2=' num2str(l(2),3) ' L3=' num2str(l(3),3) '\newline'];
    lratio_str=['L1/L2=' num2str(l(1)/l(2),2) ' L2/L3=' num2str(l(2)/l(3),2) '\newline'];
    v1_str=['v1=[' num2str(v(1,:),'%6.2f') '] \newline'];
    v2_str=['v2=[' num2str(v(2,:),'%6.2f') '] \newline'];
    v3_str=['v3=[' num2str(v(3,:),'%6.2f') '] \newline'];
    v_str=[v1_str v2_str v3_str];
    set(ud.result_text,'string',[l_str lratio_str v_str],'verticalalignment','top');
    % disp(l_str);disp(lratio_str);disp(v1_str);disp(v2_str);disp(v3_str);
    set(gcf,'userdata',ud);
    irf_minvar_gui('update_mva_axis');
end
ud=get(gcf,'userdata'); % assign ud that can be accessed because it is global
if nargout == 0, h=[]; end
end


function fix_legends
ud=get(gcf,'userdata');

switch size(ud.X,2)-1 % how many components
  case 3 %
    legend(ud.h(1),'x','y','z','Location','EastOutside');
    legend(ud.h(2),'max','interm','min','Location','EastOutside');
  case 4 %
    legend(ud.h(1),'x','y','z','abs','Location','EastOutside');
    legend(ud.h(2),'max','interm','min','abs','Location','EastOutside');
end
end

function fix_hittest(h)
% function call when marking with mouse
set(h,    'buttondownfcn', {@click_ax});
% fixes that buttondownfcn of axes is called instead of children
hChildren = get(h,'children');
if ~iscell(hChildren)
  hChildren = {hChildren};
end
for iH = 1:numel(hChildren)
  set(hChildren{iH},'hittest','off');
end
  function click_ax(varargin)
    irf_minvar_gui('ax');
  end
end

function setmethod(hObj,event) %#ok<INUSD>
% Called when user activates popup menu of minvar method
val = get(hObj,'Value');
data=get(gcf,'userdata');
if val ==1
  data.mvar_method='mvar';
elseif val == 2
  data.mvar_method='td';
elseif val == 3
  data.mvar_method='<Bn>=0';
end
set(gcf,'userdata',data);
end
