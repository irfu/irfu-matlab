function c=irf_plot(varargin)
%IRF_PLOT   Flexible plotting routine for time series
%
% c=irf_plot([H],X,[arguments...]);
%   H axes handle
%   X is one of:
%      - matrix in AV Cluster format
%      - cell array data, each of cells containing a matrix in AV Cluster
%      format
%      - string defining variable (can be CAA variable)
%      - number (initialize figure with so many subplots)
%
%   arguments can be:
%     'subplot' - plot all x values in separate subplots
%     'comp'    - plot vector component in separate subplots
%     ['dt', [dt1, dt2, dt3, dt4]] - specify time shifts, new time = old time - dt
%     ['yy',factor_to_multiply] - add second axis on right, miltiply by factor_to_multiply
%     ['linestyle',LineStyle] - define line style. Simple LineStyle can be be
%     given as last argument, 'linestyle' keyword is not necessary in this
%     case. LineStyle can be given as cell array to specify style for different variables/subplots.
%
% irf_plot, to improve zooming, will sometimes set t_start_epoch within the
% 'userdata' field of the figure and internally use it as origo but in most cases you should not care about this.
%
% Examples:
%   irf_plot(B1) - plot variable B1 (all components), assuming that the
%                    first column is time
%   irf_plot('B1') - plot variable B1, if it does not exist try to load it
%                    with c_load('B1') and try to put ylabel from c_desc('B1')
%   irf_plot('B?') - Cluser oriented, plot B1.. B4 in separate subplots
%   irf_plot({B1,B2}) - plot B1 and B2 in separate subplots
%   irf_plot('B1 B2') - -"- but if B1,B2 do not exist try to load them and
%                    put labels according to c_desc
%   irf_plot({B1,B2},'comp') - plot in 1. subplot B1_X and B2_X, in second
%                    subplot B1_Y and B2_Y etc.
%   irf_plot({B1,B2},'dt',[dt1 dt2]) - separate subplots with B1 and B2,
%                    but in addition B1 and B2 time axis are shifted by dt1
%                    and dt2 correspondingly
%
% See also C_PL_TX, C_DESC

%
% $Id$

% flag_subplot 0 - one plot
%              1 - separate subplots for every component
%              2 - separate subplots for all variables in the cell array
%              3 - components of vectors in separate panels

[ax,args,nargs] = axescheck(varargin{:});
x=args{1};
if isnumeric(x), % check if single number argument, to initialize only subplots
    if numel(x)==1, % only one number
        if x>1 && x<20,
            c=initialize_figure(x);
        else
            disp('Only 1-20 number of subplots supported.;)');
        end
        return
    end
end
if isempty(ax),
    ax=gca;
end
args=args(2:end);
original_args=args;

var_desc{1} = '';
flag_subplot = 0;
have_options = 0;
caa_dataobject={[]}; % by default assume we are not working with CAA variables

if nargs > 1, have_options = 1; end

% Default values that can be override by options
dt = 0;
flag_yy = 0;
scaleyy = 1;
plot_type = '';
marker = '-';

while have_options
    l = 1;
    switch(lower(args{1}))
        case 'subplot'
            plot_type = 'subplot';
        case 'comp'
            plot_type = 'comp';
        case 'dt'
            if nargs>1
                if isnumeric(args{2})
                    dt = args{2};
                    l = 2;
                else irf_log('fcal,','wrongArgType : dt must be numeric')
                end
            else irf_log('fcal,','wrongArgType : dt value is missing')
            end
        case 'yy'
            if nargs>1
                if isnumeric(args{2})
                    flag_yy = 1;
                    scaleyy = args{2};
                    l = 2;
                else irf_log('fcal,','wrongArgType : yy must be numeric')
                end
            else irf_log('fcal,','wrongArgType : yy value is missing')
            end
        case 'linestyle'
            marker = args{2};
            l = 2;
        otherwise
            %irf_log('fcal',['Assuming ''' args{1} ''' is a LineStyle'])
            marker = args{1};
            args = args(2:end);
            break
    end
    args = args(l+1:end);
    if isempty(args), break, end
end

% Plot separate subplots for all x components
if strcmp(plot_type,'subplot') && isnumeric(x), flag_subplot = 1; end
if ischar(x), % Try to get variable labels etc.
    var_nam = tokenize(x); % White space separates variables
    jj = 1;
    for ii=1:length(var_nam), % construct varibale names var_names
        if regexp(var_nam{ii},'?'),
            c_eval(['var_names{jj}=''' var_nam{ii} ''';jj=jj+1;']);
        else
            var_names{jj} = var_nam{ii}; jj=jj+1;
        end
    end
    x = {}; ix = 1;
    for ii=1:length(var_names) % get variables
        try % Try to get variable from calling workspace
            x{ix} = evalin('caller',var_names{ii});
        catch
            try % If there is none try to load variable
                if strfind(var_names{ii},'__') % CAA variable
                    caa_varname{ix}=var_names{ii};
                    [~,caa_dataobject{ix},x{ix}]=c_caa_var_get(var_names{ii});
                else
                    c_load(var_names{ii});eval(['x{ix}=' var_names{ii} ';']);
                end
            catch % If nothing works give up
                irf_log('load',...
                    ['skipping, do not know where to get variable >'...
                    var_names{ii}]);
            end
        end
        if length(x)==ix,
            try
                var_desc{ix} = c_desc(var_names{ii});
            catch
                var_desc{ix} = {};
            end
            ix = ix +1;
        end
    end
end
if iscell(x), % Plot several variables
    
    % No ylabels are given
    % But no way to now the name of variables
    if size(var_desc,2)<size(x,2), var_desc = cell(1,length(x)); end
    
    if dt==0, dt(1:size(x,2)) = double(0); end
    
    switch plot_type
        case ''
            flag_subplot = 2;
            if length(x)==1, x = x{1}; flag_subplot = 0; end
        case 'comp'
            flag_subplot = 3;
        case 'subplot'
            flag_subplot = 2;
    end
else
    try % try to obtain variable description
        var_desc{1} = c_desc(inputname(1));
    catch %#ok<CTCH>
        var_desc{1} = {};
    end
end

% For zooming to work even in cases of wide band it is important that time
% axis is not big number. Isdat epoch is too big. Therefore if time is
% isdat epoch we choose reference time the first point of first variable
% (in practices it does not matter).

if ~isempty(caa_dataobject{1}) % plot CAA variable
    plot(ax,caa_dataobject{1},caa_varname{1},original_args{:});
    if isstruct(x), tt=x.t(1);
    elseif iscell(x), tt=x{1}(1,1);
    else tt=x(1,1);
    end
    c=ax; % axis to which apply add_timeaxis
    flag_subplot=-1; % dont make more plots
end

if flag_subplot==0,  % One subplot
    if isstruct(x)
        % Plot a spectrogram
        irf_spectrogram(ax,x);
        hcbar = colorbar('peer',ax);
        if ~isempty(var_desc{1})
            lab = cell(1,length(var_desc{1}.size));
            for v = 1:length(var_desc{1}.size)
                lab{v} = [var_desc{1}.labels{v} '[' var_desc{1}.units{v} ...
                    '] sc' var_desc{1}.cl_id];
            end
            ylabel(hcbar, lab);
        end
        
        tt = x.t(~isnan(x.t),1);
        tt = tt(1);
    else % x is matrix
        ts = t_start_epoch(x(:,1)); % t_start_epoch is saved in figures user_data variable
        ii = 2:length(x(1,:));
        if flag_yy == 0,
            h = plot(ax,(x(:,1)-ts-dt),x(:,ii),marker,args{:});
        else
            h = plotyy(ax,(x(:,1)-ts),x(:,ii),(x(:,1)-ts),x(:,ii).*scaleyy);
        end
        grid(ax,'on');
        
        % Put YLimits so that no labels are at the end (disturbing in
        % multipanel plots)
        yl = get(ax,'YLim');
        if ~(any(any(x(:,2:end) == yl(1))) || any(any(x(:,2:end) == yl(2))))
            set(ax,'YLim', mean(yl) + diff(yl)*[-.499999 .499999])
        end
        
        if ~isempty(var_desc{1}) && isfield(var_desc{1},'size')
            lab = cell(1,length(var_desc{1}.size));
            for v = 1:length(var_desc{1}.size)
                lab{v} = [var_desc{1}.labels{v} '[' var_desc{1}.units{v} ...
                    '] sc' var_desc{1}.cl_id];
            end
            ylabel(ax,lab);
        end
        
        c = get(h(1),'Parent');
        
        tt = x(~isnan(x(:,1)),1);
        tt = tt(1);
    end
    
elseif flag_subplot==1, % Separate subplot for each component
    if isstruct(x), error('cannot plot spectra in COMP mode'), end
    
    % t_start_epoch is saved in figures user_data variable
    ts = t_start_epoch(x(:,1));
    
    npl = size(x,2) -1;
    c = zeros(1,npl);
    for ipl=1:npl
        c(ipl) = subplot(npl,1,ipl);
        
        if iscell(marker)
            if length(marker)==npl, marker_cur = marker{ipl};
            else marker_cur = marker{1};
            end
        else marker_cur = marker;
        end
        
        plot((x(:,1)-ts-dt),x(:,ipl+1),marker_cur,args{:}); grid on;
        
        % Put YLimits so that no labels are at the end (disturbing in
        % multipanel plots)
        set(gca,'YLim', ...
            mean(get(gca,'YLim'))+diff(get(gca,'YLim'))*[-.499999 .499999])
        
        if ~isempty(var_desc) && ~isempty(var_desc{1})
            scu = cumsum(var_desc{1}.size);
            isz = find( scu == min(scu(ipl<=scu)) );
            sz = var_desc{1}.size(isz); % Size of a data vector
            if sz == 1 % Scalar data
                lab = [var_desc{1}.labels{isz} ' ['...
                    var_desc{1}.units{isz} '] sc' var_desc{1}.cl_id];
            else % Vector data
                % Vector component
                if isz==1, comp = ipl;
                else comp = ipl -scu(isz-1);
                end
                lab = [var_desc{1}.labels{isz} ...
                    '_{' var_desc{1}.col_labels{isz}{comp} '} ['...
                    var_desc{1}.units{isz} '] sc' var_desc{1}.cl_id ];
            end
            ylabel(lab);
        end
    end
    
    tt = x(~isnan(x(:,1)),1);
    tt = tt(1);
    
elseif flag_subplot==2, % Separate subplot for each variable
    if isempty(x), return, end
    
    %   t_start_epoch is saved in figures user_data variable
    if isstruct(x{1}), ts = t_start_epoch(x{1}.t);
    else ts = t_start_epoch(x{1}(:,1));
    end
    
    t_st = []; t_end = [];
    xlen = [];
    
    npl = size(x,2);
    c=initialize_figure(npl);
    for ipl=1:npl
        y = x{ipl};
        if isstruct(y), t_tmp = double(y.t);
        else t_tmp = double(y(:,1));
        end
        t_tmp = t_tmp -double(ts) -double(dt(ipl));
        tt = t_tmp(~isnan(t_tmp));
        if isempty(t_st), t_st = tt(1);
        else if tt(1)<t_st, t_st = tt(1); end
        end
        if isempty(t_end), t_end = tt(end);
        else if tt(end)>t_end, t_end = tt(end); end
        end
        clear tt
        
        if isstruct(y)
            irf_spectrogram(c(ipl),y.t-dt(ipl), y.p, y.f);
            hcbar = colorbar;
            if ~isempty(var_desc{ipl})
                lab = cell(1,length(var_desc{ipl}.size));
                for v = 1:length(var_desc{ipl}.size)
                    lab{v} = [var_desc{ipl}.labels{v} '[' var_desc{ipl}.units{v} ...
                        '] sc' var_desc{ipl}.cl_id];
                end
                ylabel(hcbar, lab);
                disp(lab)
            end
            tt = y.t(~isnan(y.t),1);
            % Save panel width to resize the rest of the panels accordingly
            if isempty(xlen)
                xlen = get(c(ipl),'Position');
                xlen = xlen(3);
            end
        else
            if iscell(marker)
                if length(marker)==npl, marker_cur = marker{ipl};
                else marker_cur = marker{1};
                end
            else marker_cur = marker;
            end
            plot(c(ipl),t_tmp,y(:,2:end),marker_cur); 
            grid(c(ipl),'on');
            
            % Put YLimits so that no labels are at the end (disturbing in
            % multipanel plots)
            set(c(ipl),'YLim',...
                mean(get(c(ipl),'YLim'))+diff(get(c(ipl),'YLim'))*[-.499999 .499999])
            
            if ~isempty(var_desc) && ~isempty(var_desc{ipl})
                for v = 1:length(var_desc{ipl}.size)
                    lab{v} = [var_desc{ipl}.labels{v} '[' ...
                        var_desc{ipl}.units{v} '] sc' var_desc{ipl}.cl_id];
                end
                ylabel(c(ipl),lab); clear lab
            end
            tt = y(~isnan(y(:,1)),1);
        end
    end
    % Set common XLim
    for ipl=1:npl
        set(c(ipl),'XLim',[t_st t_end])
        if ~isempty(xlen)
            p = get(c(ipl),'Position');
            set(c(ipl),'Position',[p(1) p(2) xlen p(4)])
        end
    end
    clear t_st t_end
    
    tt = tt(1);
    
elseif flag_subplot==3,  % components of vectors in separate panels
    if isstruct(x), error('cannot plot spectra in COMP mode'), end
    % t_start_epoch is saved in figures user_data variable
    ts = t_start_epoch(x{1}(:,1));
    
    npl = size(x{1},2) -1;
    if npl==1,     % We make new figure with subplots only if more than 1 component to plot
        c = gca;
    else
        c=initialize_figure(npl);
    end
    for ipl=1:npl
        
        line_colors=get(c(ipl),'ColorOrder');
        for jj=1:size(x,2)
            use_color = 1;
            if iscell(marker)
                if length(marker)==size(x,2), marker_cur = marker{jj};  use_color = 0;
                else marker_cur = marker{1};
                end
            else marker_cur = marker;
            end
            
            if size(x{jj},2)>=ipl+1
                y = x{jj};
                if use_color
                    plot(c(ipl),(y(:,1)-ts-dt(jj)), y(:,ipl+1),...
                        'Color', line_colors(jj,:), 'LineStyle',marker_cur)
                else
                    plot(c(ipl),(y(:,1)-ts-dt(jj)), y(:,ipl+1),marker_cur)
                end
                hold(c(ipl),'on');
            end
        end
        grid(c(ipl),'on');
        
        % Put YLimits so that no labels are at the end (disturbing in
        % multipanel plots)
        set(c(ipl),'YLim',...
            mean(get(c(ipl),'YLim'))+diff(get(c(ipl),'YLim'))*[-.499999 .499999])
        
    end
    tt = y(~isnan(y(:,1)),1);
    tt = tt(1);
end

irf_figmenu;

% Add information about subplot handles to userdata of figure
user_data = get(gcf,'userdata');
if flag_subplot>0, user_data.subplot_handles = c; end
set(gcf,'userdata',user_data);


% In case time is in isdat_epoch add time_axis
if ((tt > 1e8) && (tt < 1e10))
    if flag_subplot == 0, add_timeaxis(ax);
    else add_timeaxis(c);
    end
end

% Do not give axis handle as answer if not asked for
if nargout==0, clear c; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_st_e = t_start_epoch(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gives back the value of t_start_epoch of the figure
% if not  set, sets t_start_epoch of the figure
ud = get(gcf,'userdata');
ii = find(~isnan(t));
if ~isempty(ii), valid_time_stamp = t(ii(1)); else valid_time_stamp = []; end

if isfield(ud,'t_start_epoch')
    t_st_e = double(ud.t_start_epoch);
elseif ~isempty(valid_time_stamp)
    if valid_time_stamp > 1e8
        % Set start_epoch if time is in isdat epoch
        % Warn about changing t_start_epoch
        t_st_e = double(valid_time_stamp);
        ud.t_start_epoch = t_st_e;
        set(gcf,'userdata',ud);
        irf_log('proc',['user_data.t_start_epoch is set to ' ...
            epoch2iso(t_st_e,1)]);
    else
        t_st_e = double(0);
    end
else
    t_st_e = double(0);
end

end
function c=initialize_figure(number_of_subplots,flag)
% flag = "newfigure" % [optional] if to open a new figure
if nargin==1, flag='';end
if isempty(get(0,'CurrentFigure')) % no current figures opened
    flag='newfigure';
elseif isempty(get(gcf,'children')) % current figure is empty
    flag='newfigure';    
end
if number_of_subplots>1 && number_of_subplots<20,
    number_of_subplots=floor(number_of_subplots);
    c=zeros(1,number_of_subplots);
    if strcmpi(flag,'newfigure'), % if to open new figure
        set(gcf,'color','white'); % white background for figures (default is grey)
        set(gcf,'renderer','zbuffer'); % opengl has problems on Mac (no log scale in spectrograms)
        set(gcf,'PaperUnits','centimeters');
        set(gcf,'defaultlinelinewidth',1.0);
        xSize = 10;
        ySize = 5+5*sqrt(number_of_subplots);
        xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
        sz=get(0,'screensize'); 
        xx=min(min(600,sz(3))/xSize,min(900,sz(4))/ySize); % figure at least 600 wide or 900 height but not outside screen
        set(gcf,'Position',[10 10 xSize*xx ySize*xx])
        clear xSize sLeft ySize yTop
    end
    for j=1:number_of_subplots,
        c(j)=irf_subplot(number_of_subplots,1,-j);
        cla(c(j));
    end
    figure(gcf); % bring figure to front
end
end