function c=av_tplot(x,t_unit_in_original_units,t_origo_in_original_units,varargin);
% function c=av_tplot(x,t_unit_in_original_units,t_origo_in_original_units,varargin);
% function c=av_tplot(x,'subplot') to plot in separate subplots all x values
% function c=av_tplot(x,'yy',factor_to_multiply) to plot in separate subplots all x values
% function c=av_tplot({x y z}) to plot subplots with x y z in them
% function c=av_tplot({p1 p2 p3 p4 ...},[dt1 dt2 dt3 dt4 ...]) to plot subplots with x y z in them with given time shifts
% function c=av_tplot(x,1,0,varargin) to pass different options to plot routines within av_tplot

flag_subplot=0;flag_yy=0;
if ((nargin >= 2) & isstr(t_unit_in_original_units)),
 q=t_unit_in_original_units;t_unit_in_original_units=1;
 if strcmp(q,'subplot'),
  if isnumeric(x),flag_subplot=1;end    % plot separate subplots for all x components
 elseif strcmp(q,'yy'),
  if isnumeric(x),flag_yy=1;end    % add second yy axis
  if nargin > 2, scaleyy=t_origo_in_original_units; else scaleyy=1;end
 else,
  return;
 end
end

if iscell(x), % plot separate subplots for all x cell arrays
 flag_subplot=2;
 if nargin == 2,
  dt=t_unit_in_original_units;
  t_unit_in_original_units=1;t_origo_in_original_units=0;
 else,
  dt(1:size(x,2))=0;
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
  end
  tt=y(1,1);
end

  % in case time is in isdat_epoch add time_axis
  if ((tt > 1e8) & (tt < 1e10))
    if flag_subplot == 0, add_timeaxis(gca);
    else, add_timeaxis(c);
    end
  end

