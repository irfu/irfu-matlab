function out=makeFName(st,fmt)
%makeFName construct a filename string from time
% fname = makeFName(st, [FORMAT])
% Input:
% st - ISDAT epoch
%
% Output:
% fname - string formatted accordint to FORMAT:
%	0: YYYYMMDD_hhmm (default)
%	1: YYMMDDhhmmss
%	2: YYYYMMDD_hhmmss_hhmmss (CAA)
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

if nargin < 2, fmt = 0; end
if fmt==2 & length(st)<=1, error('ST must have two elements for FORMAT=2'), end

d = fromepoch(st(1));
s{1} = num2str(d(1));

for k=2:6
	s{k} = num2str(round(d(k)));
	if d(k)<10, s{k} = ['0' s{k}]; end
end

if fmt==0
	out = [s{1} s{2} s{3} '_' s{4} s{5}];
elseif fmt==1
	out = [s{1}(3:4) s{2} s{3} s{4} s{5} s{6}(1:2)];
elseif fmt==2
	d = fromepoch(st(end));
	for k=4:6
		se{k-3} = num2str(round(d(k)));
		if d(k)<10, se{k-3} = ['0' se{k-3}]; end
	end
	out = [s{1} s{2} s{3} '_' s{4} s{5} s{6} '_' se{1} se{2} se{3}];
else
	error('unknown format')
end
