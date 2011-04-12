function irf_zoom(varargin)
%IRF_ZOOM   Zoom in to x and y axes
%  Zooms to specified interval, avoids labels at ends for y zooming
%
%   IRF_ZOOM(AX,...) zooms in specified axes
%
%   IRF_ZOOM('x',xlim) zooms X axis
%       X axis are usualy time. in this case xlim can be in different form
%       xlim=[tlim1 tlim2] - time interval specified in EPOCH
%       xlim={[yyyy mm dd hh mm ss] [yyyy mm dd hh mm ss]}
%           left side of vector can be skipped, then uses it from axes
%
%   IRF_ZOOM('x',xlim,'tref',tref) tref is EPOCH of time=0 point
%
%   IRF_ZOOM('y',ylim) zooms Y axis (avoiding labels at top and bottom)
%       useful when having many subpanels
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

% $Id$

flag_use_t_start_epoch=0; % if 1 use userdata.t_start_epoch as tref
if nargin==0, help irf_zoom, return; end
t_ref=0;flag_tref=0; % default value
flag_old_syntax=0;

%% check axes
[ax,args,nargs] = axescheck(varargin{:});
if isempty(ax),
    if any(ishandle(args{1})), % first argument is axis handles
        ax=args{1};
        args=args(2:end);
        nargs=nargs-1;
    else
        if nargs >= 3, % check the OLD syntax
            if any(ishandle(args{3})), % OLD syntax
                disp('WARNING!!!!!!!!!!!!!!!!!!!!!!!!')
                disp('you use old syntax of IRF_ZOOM!')
                disp('will be disabled soon! see help')
                flag_old_syntax=1;
                ax=args{3};
                c=args{2};
                interval=args{1};
                if nargs==5, % tref in OLD syntax
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
if ~flag_old_syntax,
    c=args{1};
    interval=args{2};
    if nargs==4, % check if tref
        if strcmpi(args{3},'tref'),
            t_ref=args{4};
        end
    end
end

%% Set tref
if ~flag_tref && nargs <4, % no tref specified
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

if size(interval,2) ~= 2, % check if interval input is ok
    disp('zooming interval in wrong format');
    return;
end

if strcmpi(c,'x'),
    if iscell(interval),  % Simplified time zooming
        ax=get(axis_handle(1),'xlim');
        if ( ax(1)+t_ref>1e8 && ax(1)+t_ref<1e10 ),
            int_min=fromepoch(ax(1)+t_ref);
            int_max=fromepoch(ax(2)+t_ref);
            int_min(7-size(interval{1},2):6)=interval{1};
            int_max(7-size(interval{2},2):6)=interval{2};
            clear interval;
            interval=[toepoch(int_min) toepoch(int_max)]-t_ref;
        end
    else % Interval must be vector with two values
        if flag_use_t_start_epoch, % Account for reference time from userdata.t_start_epoch
            interval=interval-t_ref;
        end
    end
end

% Make interval finite if it has only one point
if isnumeric(interval),
    if diff(interval)==0, interval(2)=interval(1)+1e-10; end
end

% Remove XTickLabel and XLabel from all panels but the last one
if strcmpi(c,'x') && numel(axis_handles)>1
    p = cell2mat(get(axis_handles,'Position'));
    pymin = min(p(:,2));
end

for h=axis_handles
    switch lower(c)
        case 'x'
            set(h,'XLim',interval);
            if ax(1)+t_ref>1e8 && ax(1)+t_ref<1e10
                if flag_use_t_start_epoch % Read t_ref from userdata.t_start_epoch
                    p = get(h,'position');
                    if numel(axis_handles)>1, % in case of multiple handles only last handle gets date label
                        if p(2)==pymin,
                            add_timeaxis(h);
                        else
                            add_timeaxis(h,'nolabels');
                        end
                    else
                        add_timeaxis(h);
                    end
                else
                    add_timeaxis(h,t_ref);
                    if length(axis_handles)>1
                        p = get(h,'position');
                        if p(2)>pymin, xlabel(h,''), set(h,'XTickLabel',''), end
                    end
                end
            end
        case 'y'
            set(h,'Ylim',interval);
            yzoom=mean(get(h,'YLim'))+diff(get(h,'YLim'))*[-.499999 .499999];
            set(h,'Ylim',yzoom);
    end
end
