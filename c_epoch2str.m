function res=c_epoch2str(t)
%c_epoch2str convert ISDAT epoch to yyyy-mm-ddTHH:MM:SS.SSSZ format
% res=c_epoch2str(t)
%
% $Id$
%
% see also FROMEPOCH

% Copyright 2004 Yuri Khotyaintsev

error(nargchk(1,1,nargin))

s = fromepoch(t);
for j=1:5
	ss{j} = num2str(s(j+1));
	if s(j+1)<10,	ss{j} = ['0' ss{j}]; end
end
res = sprintf('%4d-%s-%sT%s:%s:%sZ',s(1),ss{1},ss{2},ss{3},ss{4},ss{5});
