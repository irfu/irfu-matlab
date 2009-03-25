function [x,in]=irf_tlim(x,t1,t2,mode)
%IRF_TLIM   Find a subinterval defined by time limits
%
% [Y, IN] = IRF_TLIM(X,T1,T2,[MODE])
% [Y, IN] = IRF_TLIM(X,[T1 T2])
% Y is part of the X that is within interval T1 <= X(:,1) < T2
% if MODE is nonzero, part of X outside the interval :
% X(:,1) < T1 & X(:,1) > T2
% 
% [y,ind]=irf_tlim(x,t1,t2)
% ind - returns inndexes of rows that were within interval
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

error(nargchk(2,4,nargin))

if isempty(x), in = []; return, end
if nargin==2, t2 = t1(2) + 1e-7; t1 = t1(1) - 1e-7; end
if nargin<4, mode = 0; end

[t1,dt] = irf_stdt(t1,t2);
t2 = t1 + dt;

if isstruct(x)
	if ~isfield(x,'t'), error('struct X must have field ''t'''), end
	t = x.t;
else
	t = x(:,1);
end

if mode==0, in = find((t >= t1) & (t < t2));
else in = find((t < t1) | (t > t2));
end

if isstruct(x)
	fn = fieldnames(x);
	ndata = length(t);
	for fi=1:length(fn)
		if iscell(x.(fn{fi}))
			for i = 1:length(x.(fn{fi}))
				if size(x.(fn{fi}){i},1) == ndata
					x.(fn{fi}){i} = x.(fn{fi}){i}(in,:);
				end
			end
		else
			if size(x.(fn{fi}),1) == ndata
				x.(fn{fi}) = x.(fn{fi})(in,:);
			end
		end
	end
else
	x = x(in,:);
end
