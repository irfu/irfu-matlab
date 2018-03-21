function [x,in]=irf_tlim(x,tStart,tEnd,mode)
%IRF_TLIM   Returns data within specified time interval
%
% Y = IRF_TLIM(X,tStart,tEnd,[MODE])
% Y = IRF_TLIM(X,[tStart tEnd],[MODE])
%
% Y is part of the X that is within interval tStart <= X(:,1) < tEnd
% if MODE is nonzero, Y is part of X outside the interval :
% X(:,1) < tStart & X(:,1) > tEnd
% 
% [Y,IND]=IRF_TLIM(...)
%	IND contains indexes of rows that were returned 

if nargin == 0 
	help irf_tlim;
	return;
end

if isempty(x), in = []; return, end

if ischar(tStart) % time interval specified as string
	tStart=irf_time(tStart,'utc>tint');
end

if nargin==2
	mode = 0;
	tEnd = tStart(2) + 1e-7; tStart = tStart(1) - 1e-7;
elseif nargin==3
	if length(tStart)==2 && length(tEnd)==1
		mode = tEnd;
		tEnd = tStart(2) + 1e-7; tStart = tStart(1) - 1e-7;
	else, mode = 0;
	end
elseif nargin == 1 && nargin > 4
	irf.log('critical','incorrect number of input arguments');
	error('incorrect number of input arguments');
end

[tStart,dt] = irf_stdt(tStart,tEnd);
tEnd = tStart + dt;

if isstruct(x)
	if ~isfield(x,'t'), error('struct X must have field ''t'''), end
	t = x.t;
else
	t = x(:,1);
end

if mode==0, in = find((t >= tStart) & (t < tEnd));
else, in = find((t < tStart) | (t > tEnd));
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
