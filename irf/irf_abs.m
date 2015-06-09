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
%
% $Id$

if isempty(x), y=[];return;end % empty output for empty input

if isa(x,'TSeries')
  % Time series
  if nargin == 2 && x.tensorOrder==1 && flag == 1
    if(regexp(x.tensorBasis,'xy')) % cartesian xy (2d) or xyz (3d)
      data = sqrt( sum(abs(x.data).^2, 2) );
      y = irf.ts_scalar(x.time, data);
    elseif(regexp(x.tensorBasis,'r[tl]p')) % spherical, colatitude/latitude
      data = abs(x.r.data);
      y = irf.ts_scalar(x.time, data);
    elseif(regexp(x.tensorBasis,'rpz')) % cylindrical
      data = sqrt(abs(x.r.data).^2 + abs(x.z.data).^2);
      y = irf.ts_scalar(x.time, data);
    else
      error('Not yet implemented!');
    end
  else
    error('Not yet implemented!');
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
  if nargin == 2,
    if flag == 1,
      yy=y(:,lx+1);clear y;y=yy;
    end
  end
end

end