function res = getLineStyle(n)
%res = getLineStyle(n)
% n - integer starts from 0

colors = {'k' 'r' 'g' 'b' 'm' 'c' 'y'};
markers = {'-' '--' '.-' '-.-'};

mr = floor(n/length(colors));
cl = n - mr*length(colors);

res = sprintf('%s%s',colors{cl+1},markers{mr+1});
