function c=av_tplot(x,t_unit_in_original_units,t_origo_in_original_units,varargin);
% function c=av_tplot(x,t_unit_in_original_units,t_origo_in_original_units,varargin);
% function c=av_tplot(x,'subplot') to plot in separate subplots all x values
% function c=av_tplot(x,'yy',factor_to_multiply) add second axis on right
% function c=av_tplot({x, y, z},[dt1 dt2 ..]) 1st subplot with x, 2nd with y etc. Subtract time shifts if given
% function c=av_tplot({x, y, z},[dt1 dt2 ...],'comp') 1.subplot is x(:,2),y(:,2)... 2nd subplot is x(:,3) y(:,3 .. etc.
% function c=av_tplot({p1, p2, p3, p4 ...},[dt1 dt2 dt3 dt4 ...]) to plot subplots with x y z in them with given time shifts
% function c=av_tplot(x,1,0,varargin) to pass different options to plot routines within av_tplot

% flag_subplot 0 - one plot
%              1 - separate subplots for every component
%              2 - separate subplots for all variables in the cell array

if now<datenum(2004,12,07)
    disp('********************************************************************');
    disp('new features, e.g. try >av_tplot(''B?'') or av_tplot(''B1 R2'');');
    disp('everything you find in c_desc works in av_tplot, even if not loaded!');
    disp('***************** this message only until 2004-12-07 **********************');
end

flag_subplot=0;
flag_yy=0;
if nargin >=2, inp2=t_unit_in_original_units;end
if nargin >=3, inp3=t_origo_in_original_units;end

if ((nargin >= 2) & isstr(t_unit_in_original_units)),
 q=t_unit_in_original_units;
 t_unit_in_original_units=1;
 if strcmp(q,'subplot'),
  if isnumeric(x),flag_subplot=1;end    % plot separate subplots for all x components
 elseif strcmp(q,'yy'),
  if isnumeric(x),flag_yy=1;end    % add second yy axis
  if nargin > 2, scaleyy=t_origo_in_original_units; else scaleyy=1;end
 end
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
                c_load(var_names{ii});eval(['x{ii}=' var_names{ii}]);
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
 if nargin == 1,
    flag='subplot';
    dt(1:size(x,2))=0;
    t_unit_in_original_units=1;t_origo_in_original_units=0;
 elseif nargin == 2,
   if ischar(inp2),
    flag=inp2;
    dt(1:size(x,2))=0;
    t_unit_in_original_units=1;t_origo_in_original_units=0;
   else
    dt=t_unit_in_original_units;
    t_unit_in_original_units=1;t_origo_in_original_units=0;
    flag='subplot';
   end
 elseif nargin>=3,
   dt=t_unit_in_original_units;
   t_unit_in_original_units=1;
   flag=t_origo_in_original_units;
   t_origo_in_original_units=0;
 end
 switch flag
   case 'comp'
     flag_subplot=3;
   case 'subplot'
     flag_subplot=2;
   otherwise
     error('input not allowed');
 end
end

if ((nargin <2))
 t_unit_in_original_units=1;
end
if (nargin <3)
 t_origo_in_original_units=0;
end

tu=t_unit_in_original_units;
ts=t_origo_in_original_units;

if flag_subplot==0,
  i=2:length(x(1,:));
  if flag_yy == 0, h=plot((x(:,1)-ts)/tu,x(:,i),varargin{:});grid on;
  else, h=plotyy((x(:,1)-ts)/tu,x(:,i),(x(:,1)-ts)/tu,x(:,i).*scaleyy);grid on;
  end
  ylabel(ylabels{1});
  c=get(h(1),'Parent');
  tt=x(1,1);
elseif flag_subplot==1,
  npl=size(x,2)-1;
  for ipl=1:npl,
    c(ipl)=subplot(npl,1,ipl);
    i=ipl+1;
    plot((x(:,1)-ts)/tu,x(:,i));grid on;
  end
  tt=x(1,1);
elseif flag_subplot==2,
  npl=size(x,2);
  for ipl=1:npl,
    c(ipl)=av_subplot(npl,1,-ipl);
    y=x{ipl};
    i=2:length(y(1,:));
    plot((y(:,1)-ts-dt(ipl))/tu,y(:,i),varargin{:});grid on;
    ylabel(ylabels{ipl});
  end
  tt=y(1,1);
elseif flag_subplot==3,  % components of vectors in separate panels
  npl=size(x{1},2)-1;
  for ipl=1:npl,
    c(ipl)=av_subplot(npl,1,-ipl);
    line_colors={'b','g','r','c','m','y','k'}; 
    for j=1:size(x,2),
      y=x{j};
      plot((y(:,1)-ts-dt(j))/tu,y(:,ipl+1),varargin{:},line_colors{j});grid on;hold on;
    end
  end
  tt=y(1,1);
end

  % in case time is in isdat_epoch add time_axis
  if ((tt > 1e8) & (tt < 1e10))
    if flag_subplot == 0, add_timeaxis(gca);
    else, add_timeaxis(c);
    end
  end

% execute av_figmenu if there is no such menu
user_data=get(gcf,'userdata');
if isstruct(user_data)
  if isfield(user_data,'av_figmenu')
  else 
    av_figmenu;user_data.av_figmenu=1;
  end
else
  av_figmenu;user_data.av_figmenu=1;
end
set(gcf,'userdata',user_data);
