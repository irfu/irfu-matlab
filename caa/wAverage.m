function out = wAverage(data,fsample)
%wAverage compute weinghted average
% exportAscii(var,fsample)
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev
warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_waverage')

if nargin<2
	fsample = 1/(data(2,1)-data(1,1));
	c_log('proc',['Sampling frequency ' num2str(fsample)])
end

ndata = round((data(end,1)-data(1,1))*fsample);
fsample = ndata/(data(end,1)-data(1,1));
dt = 1/fsample;
c_log('proc',['Sampling period ' num2str(dt)])

out = zeros(ndata+1,2);
out(:,1) = linspace(data(1,1),data(end,1),ndata+1);
ind = round((data(:,1)-data(1,1))/dt + 1);
out(ind,2) = data(:,2);
dtmp =  [0 0 0 out(:,2)' 0 0 0];
for j=1:ndata+1
	out(j,2) = w_ave(dtmp(j:j+6));
end

function av = w_ave(x)
%m = [.1 .25 .3 .25 .1];
m = [.07 .15 .18 .2 .18 .15 .07];
cor = sum(m(find(x==0))); % find missing points==0
if cor==1, av = 0;
else, av = sum(x.*m)/(1-cor);
end
