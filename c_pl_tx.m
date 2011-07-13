function out=c_pl_tx(varargin)
%C_PL_TX   Plot data from all four Cluster spacecraft in the same plot
%
% c=c_pl_tx('x?',column,[linestyle],[dt1 dt2 dt3 dt4])
% c=c_pl_tx(x1,x2,x3,x4,[column],[linestyle],[dt1 dt2 dt3 dt4])
%   plot all 4 cluster values, time is in the first column
% c=c_pl_tx(AX,...) plot in the specified axis
%
%   column - gives which column to plot. All columns will be plotted
%            if set to empty string or ommited.
%   linestyle - string or cell (size 4) in format accepted by plot.
%            Usefull to set line style and marker (but not color).
%   dt1 dt2 dt3 dt4 - timeshifts array
%
%   Example:
%      c_pl_tx('irf_abs(B?)')
%      % plot 3 components + magnitude of B1:4.
%      c_pl_tx('B?',3:4)
%      % plot 3th and 4th components of B1:4.
%      c_pl_tx('B?',3:4,[0 2 3 .5])
%      % plot 3th and 4th components of B1:4 with timeshifts
%      c_pl_tx('B?','',[0 2 3 .5])
%      % plot all components of B1:4 with timeshifts
%      c_pl_tx('B?','.-')
%      % plot all components of B1:4 using '.-' (line with dot markers)
%      c_pl_tx('B?','',[0 2 3 .5],{'.-','*','+','-'})
%      % plot all components of B1:4 using timeshifts and individual
%      % linestyles for each sc
%
% See also IRF_PLOT, PLOT
%
% $Id$

[ax,args,nargs] = axescheck(varargin{:});
if isempty(ax), % if empty axes
    ax=gca;
end
hcf=get(ax,'parent'); % get figure handle

if nargs == 0, % show help if no input parameters
    help c_pl_tx;
    return
end

sc_list=1:4; % default plot all s/c data
error(nargchk(1,8,nargs))

if ischar(args{1})
    % We have variables defines in style B?
    flag_get_variables_from_caller=1;
    varable_name_to_get_from_caller=args{1};
    %     for cl_id=1:4
    % 		ttt = evalin('caller',irf_ssub(args{1},cl_id),'[]');
    % 		eval(irf_ssub('x? =ttt;',cl_id)); clear ttt
    % 	end
    if length(args) > 1, args = args(2:end);
    else args = ''; end
else
    % We have four variables as input
    if length(args)<4, error('use c_pl_tx(x1,x2,x3,x4) or c_pl_tx(''x?'')'), end
    % We have x1,x2..x4
    c_eval('x? = args{?};');
    if length(args) > 4, args = args(5:end);
    else args = ''; end
    flag_get_variables_from_caller=0;
end

% Check for deprecated 1
if length(args)>1
    if length(args{2})==1 && args{2}==1
        irf_log('fcal','this usage of c_pl_tx is deprecated. see help c_pl_tx')
        args(2) = [];
    end
end

column = [];
if ~isempty(args)
    if isnumeric(args{1})
        column = args{1};
        args = args(2:end);
    elseif ischar(args{1})
        % empty string means default matrix size
        if isempty(args{1}), args = args(2:end); end
    end
end

delta_t = [];
l_style = {};

if ~isempty(args), have_args = 1;
else have_args = 0;
end

while have_args
    if ischar(args{1})
        if strcmp(args{1},'sc_list')
            args = args(2:end);
            sc_list=args{1};
        else
            % assume that argument defines Linestyle
            if isempty(l_style), c_eval('l_style(?)={args{1}};')
            else irf_log('fcal','L_STYLE is already set')
            end
        end
        args = args(2:end);
    elseif iscell(args{1}) && length(args{1})==4
        % Individual linestyles for each sc
        if isempty(l_style), l_style = args{1};
        else irf_log('fcal','L_STYLE is already set')
        end
        args = args(2:end);
    elseif isnumeric(args{1}) && length(args{1})==4
        % Dt
        if isempty(delta_t), delta_t = args{1};
        else irf_log('fcal','DELTA_T is already set')
        end
        args = args(2:end);
    else
        irf_log('fcal','ignoring input argument')
        args = args(2:end);
    end
    
    if ~isempty(args), have_args = 1;
    else have_args = 0;
    end
end
if flag_get_variables_from_caller==1,
    for cl_id=sc_list,
        ttt = evalin('caller',irf_ssub(varable_name_to_get_from_caller,cl_id),'[]');
        c_eval('x? =ttt;',cl_id); clear ttt
    end
end


if isempty(column) && ~isempty(x1)
    % try to guess the size of the matrix
    column = size(x1,2);
    if column > 2, column = 2:column; end
elseif isempty(column) && ~isempty(x2)
    column = size(x2,2);
    if column > 2, column = 2:column; end
elseif isempty(column) && ~isempty(x3)
    column = size(x3,2);
    if column > 2, column = 2:column; end
elseif isempty(column) && ~isempty(x4)
    column = size(x4,2);
    if column > 2, column = 2:column; end
elseif isempty(column)
    irf_log('fcal','all inputs are empty')
    return
end


if isempty(delta_t), delta_t = [0 0 0 0]; end

% define Cluster colors
cluster_colors={'[0 0 0]';'[1 0 0]';'[0 0.5 0]';'[0 0 1]'};
if isempty(l_style),
    l_style=cell(1,4);
    for ic=1:4, l_style(ic)={['''color'','  cluster_colors{ic}]};end
else
    for ic=1:4, l_style(ic)={['''' l_style{ic} ''',''color'','  cluster_colors{ic}]};end
end

% if isempty(l_style), l_style= {'k','r','g','b'};
% else
%     cls = 'krgb';
%     c_eval('l_style(?)={[cls(?) l_style{?}]};')
%     clear cls
% end

% t_start_epoch is saved in figures user_data variable
% check first if it exist otherwise assume zero
ud=get(hcf,'userdata');
if isfield(ud,'t_start_epoch'),
    t_start_epoch = double(ud.t_start_epoch);
elseif (~isempty(x1) && x1(1,1)>1e8) || (~isempty(x1) && x2(1,1)>1e8) || ...
        (~isempty(x3) && x3(1,1)>1e8) || (~isempty(x4) && x4(1,1)>1e8)
    % Set start_epoch if time is in isdat epoch,
    % warn about changing t_start_epoch
    tt = [];
    c_eval('if ~isempty(x?), tt=[tt; x?(1,1)]; end')
    t_start_epoch = double(min(tt)); clear tt
    ud.t_start_epoch = t_start_epoch; set(hcf,'userdata',ud);
    irf_log('proc',['user_data.t_start_epoch is set to ' ...
        epoch2iso(t_start_epoch)]);
else
    t_start_epoch = double(0);
end

c_eval('ts?=t_start_epoch+delta_t(?);')
% check which data are available
sc_list=[];for ic=1:4,c_eval('if ~isempty(x?), sc_list=[sc_list ?];end',ic);end


if length(column) == 1
    for ic=1:numel(sc_list)
        c_eval(['h=irf_plot(ax, [x?(:,1)-delta_t(?) x?(:,column)],' l_style{sc_list(ic)} ');'],sc_list(ic));
        hold(ax,'on');
    end
    hold(ax,'off');
    irf_zoom(ax,'y'); % optimize Y zoom to skip labels at top and bottom
    grid(ax,'on');
    irf_timeaxis(ax);
    irf_figmenu;
else
    clf
    c = irf_plot(length(column));
    for j=1:length(column)
        %		c(j)=irf_subplot(length(column),1,-j);
        pl = '';
        for jj=sc_list
            if eval(irf_ssub('~isempty(x?)',jj))
                c_eval(['s_s=''(x?(:,1)-ts?),x?(:,column(j)),''''' l_style{jj} ''''''';'],jj);
                if isempty(pl), pl = s_s; else pl = [pl ',' s_s]; end
                clear s_s
            end
        end
        eval(['plot(c(j),' pl ')'])
        irf_zoom(ax,'y'); % optimize Y zoom to skip labels at top and bottom
        grid(c(j),'on');
        
        ud = get(hcf,'userdata'); ud.subplot_handles = c; set(hcf,'userdata',ud);
    end
    irf_timeaxis(c);
    irf_figmenu;
end

if nargout > 0, out = c; end
