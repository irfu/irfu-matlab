function out=makeFName(st)
%makeFName construc a filename string from time
% fname = makeFName(st)
% Input:
% st - ISDAT epoch
%
% Output:
% fname - string YYYYMMDD_HHMM
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

d = fromepoch(st);

s{1} = num2str(d(1));

for k=2:5
	s{k} = num2str(d(k));
	if length(s{k})<2, s{k} = ['0' s{k}]; end
end

out = [s{1} s{2} s{3} '_' s{4} s{5}];
