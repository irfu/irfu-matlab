function [y]=av_norm(x)
%function [y]=av_norm(x)
%x is given as a vector [[t,] x1, x2, [x3]]
%  if 3 components then assumed [x1,x2,x3]
%y returns normalized x 

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_norm')

lx= size(x,2); % the number of vector components

y=x;
if lx == 2
  norm=sqrt(x(:,1).^2+x(:,2).^2);
  y=x./[norm norm];
elseif lx == 3
  norm=sqrt(x(:,1).^2+x(:,2).^2+x(:,3).^2);
  y=x./[norm norm norm];
elseif lx > 3
  norm=sqrt(x(:,2).^2+x(:,3).^2+x(:,4).^2);
  n=zeros(size(x)-[0 1]);for ii=1:lx-1,n(:,ii)=norm;end
  y(:,2:lx)=x(:,2:lx)./n;
else
  disp('Not enough vector components in av_abs()')
end
