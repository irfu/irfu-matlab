function res = irf_lstyle(n)
%IRF_LSTYLE    Generate line style
%
% res = irf_lstyle(n)
% n - integer starts from 0
%
% $Id$

colors = {'k' 'r' 'g' 'b' 'm' 'c' 'y'};
markers = {'-' '--' '.-' '-.-'};

mr = floor(n/length(colors));
cl = n - mr*length(colors);

res = sprintf('%s%s',colors{cl+1},markers{mr+1});
