function [y]=irf_abs(x,flag)
%IRF_ABS   Absolute value
%
% [y]=irf_abs(x,flag)
% x is given as a vector [[t,] x1, x2, [x3]]
%  if 3 components then assumed [x1,x2,x3]
% y returns x adding column with abs(x)
%
% if flag==1 then y returns only abs(x)
%
% If input is of type TSeries and flag==1 then the returned value y is a
% scalar time series corresponding to abs(x).

if isempty(x), y=[];return;end % empty output for empty input

if isa(x,'TSeries') % Time series
  Ts = x.abs();
  if nargin == 2 && flag == 1, y = Ts.data;
  else, y = Ts;
  end
else
  lx = size(x,2); % the number of vector components

  y=[x x(:,1)*0];
  if lx == 2
    y(:,lx+1)=sqrt(x(:,1).^2+x(:,2).^2);
  elseif lx == 3
    y(:,lx+1)=sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
  elseif lx > 3
    y(:,lx+1)=sqrt(x(:,2).^2+x(:,3).^2+x(:,4).^2);
  else
    disp('Not enough vector components in irf_abs()')
  end

  % if flag=1 only abs(y) should be returned
  if nargin == 2
    if flag == 1
      yy=y(:,lx+1);clear y;y=yy;
    end
  end
end

end