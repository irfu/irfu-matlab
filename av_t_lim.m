function y=av_t_lim(x,t1,t2)
%function y=av_t_lim(x,t1,t2)
%function y=av_t_lim(x,[t1 t2])
% y is part of the x that is within interval t1 < x(:,1) < t2

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_tlim')

if nargin == 2
 t2=t1(2)+1e-7;t1=t1(1)-1e-7;
end

in=find((x(:,1) > t1) & (x(:,1) < t2));
y=x(in,:);
