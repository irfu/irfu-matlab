function [z]=irf_cross(x,y,flag)
%IRF_CROSS   Cross product
%
% [z]=irf_cross(x,y)
% [z]=irf_cross(x,y,flag)
%
% returns cross product of x and y
% x,y, are vectors - [time],x,y,z
% if necessary linearly interpolates (and/or extrapolates) y to x
% if x == [x1 x2 x3] (x is single vector with 3 components) then use y time axis
%
% if given flag=1 then return only values without time column
%

time_axis=[];
if nargin<3, flag=0;end

if size(x,2) > 3 % if more than 3 columns assume that the first one is time
  xx=x(:,[2  3 4]); % columns 2 3 4 are x, y, z
  time_axis=x(:,1);
elseif size(x,2) == 3 % check if first is not time otherwise assume there is no time
  if x(1,1) > 1e8 % first column most probably time (epoch > 1e8)
    if min(diff(x(:,1)))>0 % all values are monotonically increasing, definitely time
      xx=x(:,2:3); % use last two columns
      xx(:,3)=xx(:,1)*0; % put third column to zeros
    else % assume there is no time as first column is not monotonically increasing
      xx=x;
    end
  else
    xx=x;
  end
else
  error('not enough components for x vector');
end

if size(y,2) > 3 % if more than 3 columns assume that the first one is time
  yy=y(:,[2  3 4]);
  if isempty(time_axis) % if x does not have time axis use the time axis of y (if exist)
    time_axis=y(:,1);
  end
elseif size(y,2) == 3
  yy=y;
else
  error('not enough components for y vector');
end

if size(xx,1) == 1 && size(x,2)==3 % if x is single vector WITHOUT TIME use it for all y vectors
  xx=ones(size(yy,1),1)*xx;
end
if size(yy,1) == 1 % if y is single vector CAN BE WITH TIME use it for all x vectors
  yy=ones(size(xx,1),1)*yy;
end

if size(xx,1) ~= size(yy,1)  % input vectors not of the same length (and none of them is single vector)
  if size(x,2)>3 && size(y,2)>3 % then both vectors should have time axis
    %    if debug ==1, disp('interpolating y to x in irf_cross(x,y)');end
    yy=interp1(y(:,1),yy,x(:,1),'linear','extrap');
  else
    error('do not know how to interpret input');
  end
end


zout=[	xx(:,2).*yy(:,3)-xx(:,3).*yy(:,2), - (xx(:,1).*yy(:,3)-xx(:,3).*yy(:,1)), xx(:,1).*yy(:,2)-xx(:,2).*yy(:,1)];

if flag==1
  z=zout;
else
  z=[time_axis zout];
end
