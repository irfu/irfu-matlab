function [y,in]=irf_tlim(x,t1,t2,mode)
%IRF_TLIM   Find a subinterval defined by time limits
%
% [Y, IN] = IRF_TLIM(X,T1,T2,[MODE])
% [Y, IN] = IRF_TLIM(X,[T1 T2],[MODE])
% Y is part of the X that is within interval T1 <= X(:,1) < T2
% if MODE is nonzero, part of X outside the interval :
% X(:,1) < T1 & X(:,1) > T2
% 
% [y,ind]=irf_tlim(x,t1,t2)
% ind - returns inndexes of rows that were within interval
%
% $Id$

error(nargchk(2,4,nargin))

if nargin==2, t2 = t1(2) + 1e-7; t1 = t1(1) - 1e-7; end
if nargin<4, mode = 0; end

[t1,dt] = irf_stdt(t1,t2);
t2 = t1 + dt;

if mode==0, in = find((x(:,1) >= t1) & (x(:,1) < t2));
else in = find((x(:,1) < t1) | (x(:,1) > t2));
end

y = x(in,:);
