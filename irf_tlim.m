function y=irf_tlim(x,t1,t2)
%IRF_TLIM   iFind a subinterval defined by time limits
%
% y=irf_tlim(x,t1,t2)
% y=irf_tlim(x,[t1 t2])
% y is part of the x that is within interval t1 < x(:,1) < t2
%
% $Id$

if nargin == 2
 t2=t1(2)+1e-7;t1=t1(1)-1e-7;
end

in=find((x(:,1) > t1) & (x(:,1) < t2));
y=x(in,:);
