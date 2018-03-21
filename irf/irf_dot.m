function [z]=irf_dot(x,y,flag_output)
%IRF_DOT   Calculate dot product between vectors in 3D space
%
% [z]=IRF_DOT(x,y,flag) calculate dot product of  vectors x and y 
%    having 3 components
%
% x,y,z  are column vectors  [time],x,y,z,[r]
%
% If x or y has only one row then this vector value is used
% in dot product for all y(or x) data points.
% If time axis of x and y are different then y is interpolated to x time axis
%
% if flag==1 then y returns only values of dot product without time axis

z=x;

if size(x,2) < 3
  errS = 'not enough components for y vector.';
  irf.log('critical',errS), error(errS)
elseif size(x,2) == 3, xx=x;
else, xx=x(:,2:4);
end

if size(y,2) > 3, yy=y(:,2:4);
elseif size(y,2) == 3,  yy=y;
else
  errS = 'not enough components for y vector.';
  irf.log('critical',errS), error(errS)
end

if size(x,1) ~= size(y,1)
  if (size(x,1) == 1)
    qq=yy;qq(:,1)=xx(1,1);qq(:,2)=xx(1,2);qq(:,3)=xx(1,3);xx=qq;
  elseif size(y,1) == 1
    qq=xx;qq(:,1)=yy(1,1);qq(:,2)=yy(1,2);qq(:,3)=yy(1,3);yy=qq;
  else
    irf.log('warning','interpolating y to x, assuming that first column is time');
    qq=irf_resamp(y,x);yy=qq(:,2:end);
  end
end

zout=xx(:,1).*yy(:,1)+xx(:,2).*yy(:,2)+xx(:,3).*yy(:,3);

z(:,end)=[];z(:,end)=[];z(:,end)=zout;
if size(z,2)>2  % if input is vector [t x y z r ...], result should be [t dotproduct]
  z(:,2:end-1)=[];
end

% if flag=1 only abs(y) should be returned
if exist('flag_output','var')
  if (flag_output == 1) && size(z,2)>1
    z=z(:,end);  
  end
end

