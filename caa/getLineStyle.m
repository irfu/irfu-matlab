function res = getLineStyle(n)
%res = getLineStyle(n)
% n - integer starts from 0
warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',
...
mfilename,'irf_lstyle')

colors = {'k' 'r' 'g' 'b' 'm' 'c' 'y'};
markers = {'-' '--' '.-' '-.-'};

mr = floor(n/length(colors));
cl = n - mr*length(colors);

res = sprintf('%s%s',colors{cl+1},markers{mr+1});
