function c=irf_plot(x,option,varargin);
%IRF_PLOT   Flexible plotting routine
%
% c=irf_plot(x,option,varargin);
% c=irf_plot(x,'subplot') to plot in separate subplots all x values
% c=irf_plot(x,'','yy',factor_to_multiply) add second axis on right
% c=irf_plot({x, y, z},'','dt',[dt1 dt2 ..]) 1st subplot with x, 
%   2nd with y etc. Subtract time shifts if given
% c=irf_plot({x, y, z},'comp','dt',[dt1 dt2 ...]) 1.subplot is x(:,2),y(:,2)...
%   2nd subplot is x(:,3) y(:,3 .. etc.
% c=irf_plot({p1, p2, p3, p4 ...},[dt1 dt2 dt3 dt4 ...]) to plot subplots with 
%   x y z in them with given time shifts
% c=irf_plot(x,'',varargin) to pass different options to plot routines 
%   within irf_plot
%
% varargin is pairs 
%   option, value
% where option can be any option for plot routine or any of the below
%  'dt' - specify time shifts value=[dt1, dt2, dt3, dt4, ...]
%  'yy' - add second axis on right, value=factor_to_multiply
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
%                    subplot B1_Y and B2_Y etc . 
%   irf_plot({B1,B2},'','dt',[dt1 dt2]) - separate subplots with B1 and B2, 
%                    but in addition B1 and B2 time axis are shifted by dt1 
%                    and dt2 correspondingly
%
% $Id$


% flag_subplot 0 - one plot
%              1 - separate subplots for every component
%              2 - separate subplots for all variables in the cell array
%              3 - components of vectors in separate panels

ylabels{1}='';
flag_subplot=0;
if nargin==1, option='';end


% default values that can be override by options
dt=0;
flag_yy=0;scaleyy=1;

ind_varargin_clean=[]; % which indexes to remove after checking varargin
for j=1:2:length(varargin),
    switch varargin{j},
        case 'dt'
            dt=varargin{j+1};
        case 'yy'
            flag_yy=1;
            scaleyy=varargin{j+2};
        otherwise
    end
    ind_varargin_clean=[ind_varargin_clean j j+1];
end
varargin(ind_varargin_clean)=[];

switch option
    case ''
    case 'subplot'
        if isnumeric(x),flag_subplot=1;end    % plot separate subplots for all x components
    case 'comp'
    otherwise
        error('option not allowed');
end

if ischar(x), % try to get variable labels etc.
    var_nam=tokenize(x); % white space separates variables
    jj=1;
    for ii=1:length(var_nam),
        if regexp(var_nam{ii},'?'),
            c_eval(['var_names{jj}=''' var_nam{ii} ''';jj=jj+1;']);
        else
            var_names{jj}=var_nam{ii};jj=jj+1;
        end
    end
    x={};
    for ii=1:length(var_names)
        try % try to get variable from calling workspace
            x{ii}=evalin('caller',var_names{ii});
        catch
            try % if there is none try to load variable
                c_load(var_names{ii});eval(['x{ii}=' var_names{ii} ';']);
            catch % if nothing works give up
                warning(['do not know where to get variable >' var_names{ii}]);
                return;
            end
        end
        try
            var_desc=c_desc(var_names{ii});
            ylabels{ii}=[var_desc.labels{1} '[' var_desc.units{1} '] sc' var_desc.cl_id];
        catch
            ylabels{ii}='';
        end
    end
end


if iscell(x), % plot several variables
    if size(ylabels,2)<size(x,2), % no ylabels are given
        for ii=1:length(x);ylabels{ii}='';end % no way to now the name of variables
    end
    switch option
        case ''
            option='subplot';
            flag_subplot=2;
        case 'comp'
            flag_subplot=3;
        case 'subplot'
            flag_subplot=2;
    end
    if dt==0, dt(1:size(x,2))=0; end
else
    try
        var_desc=c_desc(inputname(1));
        ylabels{1}=[var_desc.labels{1} '[' var_desc.units{1} '] sc' var_desc.cl_id];
    catch
        ylabels{1}='';
    end
end

% for zooming to work even in cases of wide band it is important that time
% axis is not big number. isdat epoch is too big. therefore if time is
% isdat epoch we choose reference time the first point of first variable
% (in practices it does not matter).
%  


if flag_subplot==0,  % one subplot
    %   t_start_epoch is saved in figures user_data variable
    % check first if it exist otherwise assume zero
    ud=get(gcf,'userdata');
    if isfield(ud,'t_start_epoch'), 
        t_start_epoch=ud.t_start_epoch;
    else
        t_start_epoch=0;
    end
    ts=t_start_epoch;
    
    i=2:length(x(1,:));
    if flag_yy == 0, h=plot((x(:,1)-ts),x(:,i),varargin{:});grid on;
    else, h=plotyy((x(:,1)-ts),x(:,i),(x(:,1)-ts),x(:,i).*scaleyy);grid on;
    end
    ylabel(ylabels{1});
    c=get(h(1),'Parent');
    tt=x(1,1);
    
    
elseif flag_subplot==1, % separate subplot for each component 
    %   t_start_epoch is saved in figures user_data variable
    if x(1,1)> 1e8, ts=x(1,1);t_start_epoch=ts;else t_start_epoch=0;end

    npl=size(x,2)-1;
    for ipl=1:npl,
        c(ipl)=subplot(npl,1,ipl);
        i=ipl+1;
        plot((x(:,1)-ts),x(:,i),varargin{:});grid on;
    end
    tt=x(1,1);
    
elseif flag_subplot==2, % separate subplot for each variable
    %   t_start_epoch is saved in figures user_data variable
    qq=x{1};ts=qq(1,1);clear qq; if ts > 1e8, t_start_epoch=ts;else ts=0;t_start_epoch=0;end

    npl=size(x,2);
    for ipl=1:npl,
        c(ipl)=irf_subplot(npl,1,-ipl);
        y=x{ipl};
        i=2:length(y(1,:));
        plot((y(:,1)-ts-dt(ipl)),y(:,i),varargin{:});grid on;
        ylabel(ylabels{ipl});
    end
    tt=y(1,1);
    
elseif flag_subplot==3,  % components of vectors in separate panels
    %   t_start_epoch is saved in figures user_data variable
    qq=x{1};ts=qq(1,1);clear qq; if ts > 1e8, t_start_epoch=ts;else ts=0;t_start_epoch=0;end

    npl=size(x{1},2)-1;
    for ipl=1:npl,
        c(ipl)=irf_subplot(npl,1,-ipl);
%        line_colors={'b','g','r','c','m','y','k'};
         line_colors=get(gca,'ColorOrder');
        for j=1:size(x,2),
            y=x{j};
            plot((y(:,1)-ts-dt(j)),y(:,ipl+1),'color',line_colors(j,:),varargin{:});grid on;hold on;
        end
    end
    tt=y(1,1);
end

% execute irf_figmenu if there is no such menu
user_data=get(gcf,'userdata');
if isstruct(user_data)
    if isfield(user_data,'irf_figmenu')
    else
        irf_figmenu;user_data.irf_figmenu=1;
    end
else
    irf_figmenu;user_data.irf_figmenu=1;
end
set(gcf,'userdata',user_data);


% add t_start_epoch, used by add_timeaxis and subplot handles
user_data=get(gcf,'userdata');
if flag_subplot>0, 
    user_data.subplot_handles=c;
    user_data.t_start_epoch=t_start_epoch;
end % add information about subplot handles to userdata of figure
set(gcf,'userdata',user_data);


% in case time is in isdat_epoch add time_axis 
if ((tt > 1e8) & (tt < 1e10))
    if flag_subplot == 0, add_timeaxis(gca);
    else, add_timeaxis(c);
    end
end

