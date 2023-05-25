function irf_zoom(varargin)
%IRF_ZOOM   Zoom in to x and y axes
%  Zooms to specified interval, avoids labels at ends for y zooming
%
%   IRF_ZOOM(AX,...) zooms in specified axes
%
%   IRF_ZOOM('x',xlim) zooms X axis
%       X axis are usually time. in this case xlim can be in different form
%       xlim=[tlim1 tlim2] - time interval specified in EPOCH
%       xlim={[yyyy mm dd hh mm ss] [yyyy mm dd hh mm ss]}
%           left side of vector can be skipped, then uses it from axes
%
%   IRF_ZOOM('x',xlim,'tref',tref) tref is EPOCH of time=0 point
%
%   IRF_ZOOM('y') zooms Y axis (avoiding labels at top and bottom)
%       useful when having many subpanels
%
%   IRF_ZOOM('y',ylim) zooms Y axis (avoiding labels at top and bottom)
%

% Old syntax
% irf_zoom(interval,c,axis_handles,t_ref)
%  irf_zoom([min max],c,axis_handles,t_ref)
%  irf_zoom([min max],c,axis_handles)
%  irf_zoom({[yyyy mm dd hh mm ss] [yyyy mm dd hh mm ss]},c,axis_handles)
%  irf_zoom({[yyyy mm dd hh mm ss] [yyyy mm dd hh mm ss]},c,axis_handles,t_ref)
%  left side of the date vectors can be skipped, then one uses the
%  values from axis
%  c='x' for x-axis, 'y' for y-axis
%  t_ref is isdat_epoch of time=0 point
%

flag_use_t_start_epoch=0; % if 1 use userdata.t_start_epoch as tref
if nargin==0, help irf_zoom, return; end
t_ref=0;flag_tref=0; % default value
flag_old_syntax=0;

%% check axes
[ax,args,nargs] = irf.axescheck(varargin{:});
if isempty(ax)
  if any(ishandle(args{1})) % first argument is axis handles
    ax=args{1};
    args=args(2:end);
    nargs=nargs-1;
  else
    if nargs >= 3 % check the OLD syntax
      if any(ishandle(args{3})) % OLD syntax
        disp('WARNING!!!!!!!!!!!!!!!!!!!!!!!!')
        disp('you use old syntax of IRF_ZOOM!')
        disp('will be disabled soon! see help')
        flag_old_syntax=1;
        ax=args{3};
        c=args{2};
        interval=args{1};
        if nargs==5 % tref in OLD syntax
          t_ref=args{5};
          flag_tref=1;
        end
      else
        ax=gca;
      end
    else
      ax=gca;
    end
  end
end

%% NEW syntax case
if ~flag_old_syntax
  c=args{1};
  if nargs == 1 && c=='y' % auto y zooming
    interval=[]; % empty interval is auto y-zooming
  else
    interval=args{2};
  end
  if nargs==4 % check if tref
    if strcmpi(args{3},'tref')
      t_ref=args{4};
    end
  end
end

%% Set tref
if ~flag_tref && nargs <4 % no tref specified
  % Try to read the reference time from figures user_data variable
  user_data=get(gcf,'userdata');
  if isfield(user_data,'t_start_epoch')
    t_ref=user_data.t_start_epoch;
    flag_use_t_start_epoch=1;
  else
    t_ref=0;
  end
end

axis_handles = reshape(ax,1,numel(ax));

if strcmpi(c,'x')
  if iscell(interval)  % Simplified time zooming
    ax=get(axis_handles(1),'xlim');
    if ( ax(1)+t_ref>1e8 && ax(1)+t_ref<1e10 )
      int_min=fromepoch(ax(1)+t_ref);
      int_max=fromepoch(ax(2)+t_ref);
      int_min(7-size(interval{1},2):6)=interval{1};
      int_max(7-size(interval{2},2):6)=interval{2};
      clear interval;
      interval=[toepoch(int_min) toepoch(int_max)];
    end
  elseif isnumeric(interval) && size(interval,2) == 2 % Interval must be vector with two values
  elseif ischar(interval) % assume interval is specified in ISO format
    interval = irf_time(interval,'utc>tint');
  elseif isa(interval,'GenericTimeArray')
    interval = [interval.start.epochUnix interval.stop.epochUnix];
  else
    errStr = 'zooming interval in wrong format';
    irf.log('critical',errStr);
    error('irf_zoom:time_zoom:wrong_format',errStr);
  end
  if flag_use_t_start_epoch % Account for reference time from userdata.t_start_epoch
    interval=interval-t_ref;
  end
end

% Make interval finite if it has only one point
if isnumeric(interval)
  if diff(interval)==0, interval(2)=interval(1)+1e-10; end
end

% Remove XTickLabel and XLabel from all panels but the last one
if strcmpi(c,'x') && numel(axis_handles)>1
  if isgraphics(axis_handles( 1 ),'line')
    parent_handles = cell2mat(get(axis_handles,'Parent'));
    p = cell2mat(get(parent_handles,'Position'));
  else %axis
    p = cell2mat(get(axis_handles,'Position'));
  end
  pymin = min(p(:,2));
end

for hii=axis_handles
  if isgraphics(hii,'line')
    h = get(hii,'Parent');
  else
    h = hii;
  end
  switch lower(c)
    case 'x'
      set(h,'XLim',interval);
      if t_ref>1e8 && t_ref<1e10
        if flag_use_t_start_epoch % Read t_ref from userdata.t_start_epoch
          p = get(h,'position');
          if numel(axis_handles)>1 % in case of multiple handles only last handle gets date label
            if p(2)==pymin
              irf_timeaxis(h);
            else
              irf_timeaxis(h,'nolabels');
            end
          else
            irf_timeaxis(h);
          end
        else
          irf_timeaxis(h,t_ref);
          if length(axis_handles)>1
            p = get(h,'position');
            if p(2)>pymin, xlabel(h,''), set(h,'XTickLabel',''), end
          end
        end
      end
      ud=get(h,'userdata');
      ud.zoom_x=interval+t_ref;
      set(h,'userdata',ud);
    case 'y'
      ud=get(h,'userdata');
      if isempty(interval) % auto y zooming
        if isfield(ud,'zoom_y')
          interval_to_use=ud.zoom_y;
        else
          zoom_y_auto(h)
          interval_to_use=get(h,'ylim');
        end
      else
        ud.zoom_y=interval;
        interval_to_use=interval;
        set(h,'userdata',ud);
      end
      if interval_to_use(1)>0
        interval_to_use(1)=interval_to_use(1)*(1+1e-9);
      elseif interval_to_use(1)==0
        interval_to_use(1)=interval_to_use(1)+1e-9*diff(interval_to_use(1:2));
      else
        interval_to_use(1)=interval_to_use(1)*(1-1e-9);
      end
      if interval_to_use(2)>0
        interval_to_use(2)=interval_to_use(2)*(1-1e-9);
      elseif interval_to_use(2)==0
        interval_to_use(2)=interval_to_use(2)-1e-9*diff(interval_to_use(1:2));
      else
        interval_to_use(2)=interval_to_use(2)*(1+1e-9);
      end
      if interval_to_use(1) > interval_to_use(2) && interval_to_use(1) < 0
        interval_to_use(2)=interval_to_use(2)+.001;
      end
      set(h,'Ylim',interval_to_use);
  end
end

function zoom_y_auto(h)
% make more space related auto zoom than Matlab
hlines=findall(h,'Type','line');
hlines = findobj(hlines,'-not','Tag','irf_pl_mark');

ud=get(h,'userdata');
uf=get(get(h,'parent'),'userdata');
xzero=0; % reference point
if isfield(uf,'t_start_epoch'), xzero=uf.t_start_epoch;end

ylims=[];
for ih=1:numel(hlines)
  hh=hlines(ih);
  xd=get(hh,'XData')+xzero;
  yd=get(hh,'YData');
  if isfield(ud,'zoom_x') % use zoom x values
    xlims=ud.zoom_x;
    ydlim=yd(xd>xlims(1) & xd<xlims(2));
  else
    ydlim=yd;
  end
  ydlim=ydlim(isfinite(ydlim)); % remove NaN and Inf points
  if numel(ydlim)<2
    ylimd=ylims; % don't change if zooming to 1 or less points
  else
    ylimd=[min(ydlim) max(ydlim)];
  end
  if isempty(ylims)
    ylims=ylimd;
  else
    if ylimd(1)<ylims(1)
      ylims(1)=ylimd(1);
    end
    if ylimd(2)>ylims(2)
      ylims(2)=ylimd(2);
    end
  end
end
if isempty(ylims) % has been too few data points to estimate limits
  ylims=get(h,'ylim');
end
if isa(ylims, 'integer'), ylims = double(ylims); end

yscale=get(h,'yscale');

switch lower(yscale)
  case 'linear'
    diffy=diff(ylims);
    if diffy==0
      if ylims(1)==0
        ymin=-1;ymax=1;
      else
        ymin=ylims(1)-abs(ylims(1))/10;
        ymax=ylims(1)+abs(ylims(1))/10;
      end
    else
      dy=double(diffy)/4; % 1st approx
      dy10power=10^(floor(log10(dy)));
      dy1stcipher=round(dy/dy10power);
      if dy1stcipher>5
        dy = 5*dy10power;
      else
        dy=dy1stcipher*dy10power;
      end
      ymin=dy*floor(ylims(1)/dy);
      ymax=dy*ceil(ylims(2)/dy);
    end
    set(h,'ylim',[ymin ymax]);
  case 'log'
    set(h,'YLimMode','auto');
end
