function [xm,xw]=irf_mf(x,pol_order)
%function [xm,xw]=irf_mf(x,pol_order)
% estimate mean field xm and wave field xw using polynomial fit of order pol_order
% for the number of columns larger than 3 assume that first column is time
%
xm=x; xw=x;
nc=size(x,2);
nr=size(x,1);
if nc>=4
  iref=2;
  t=(x(:,1)-x(1,1))/(x(nr,1)-x(1,1));
else
  iref=1;
  t=1:size(x,1)';
end
for ii=iref:nc
 [p,~] = polyfit(t,x(:,ii),pol_order);
 xm(:,ii)=polyval(p,t);
 xw(:,ii)=x(:,ii)-xm(:,ii);
end

end
