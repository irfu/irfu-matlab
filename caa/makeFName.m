function out=makeFName(st,fmt)
%makeFName construc a filename string from time
% fname = makeFName(st, [FORMAT])
% Input:
% st - ISDAT epoch
%
% Output:
% fname - string formatted accordint to FORMAT:
%	0: YYYYMMDD_hhmm (default)
%	1: YYMMDDhhmmss
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

if nargin < 2, fmt = 0; end

d = fromepoch(st);
s{1} = num2str(d(1));

for k=2:6
	s{k} = num2str(d(k));
	if d(k)<10, s{k} = ['0' s{k}]; end
end

if fmt==0
	out = [s{1} s{2} s{3} '_' s{4} s{5}];
elseif fmt==1
	out = [s{1}(3:4) s{2} s{3} s{4} s{5} s{6}(1:2)];
else
	error('unknown format')
end
