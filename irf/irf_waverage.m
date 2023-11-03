function out = irf_waverage(data,fsample,nAveragePoints)
%IRF_WAVERAGE  Compute weighted average
%
% out=irf_waverage(var,fsample,[nAveragePoints])
%
% nAveragePoints = 7 (default) or 5 - number of point to average over

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin <3, nAveragePoints = 7; end
if ~ismember(nAveragePoints,[5 7])
  error('nAveragePoints must be 5 or 7')
end
if size(data,1)<=1
  irf_log('proc',['Not enough points (' num2str(size(data,1)) ') to average'])
  out = data;
  return
end
if nargin<2
  fsample = 1/(data(2,1)-data(1,1));
  irf_log('proc',['Sampling frequency ' num2str(fsample)])
end

ndata = round((data(end,1)-data(1,1))*fsample);
ncol = size(data,2);
fsample = ndata/(data(end,1)-data(1,1));
dt = 1/fsample;
irf_log('proc',['Sampling period ' num2str(dt)])

out = zeros(ndata+1,ncol);
out(:,1) = linspace(data(1,1),data(end,1),ndata+1);
ind = round((data(:,1)-data(1,1))/dt + 1);
out(ind,2:end) = data(:,2:end);
out(isnan(out)) = 0; % set NaNs to zeros
for col=2:ncol
  switch nAveragePoints
    case 5
      dtmp =  [0 0 out(:,col)' 0 0];
    case 7
      dtmp =  [0 0 0 out(:,col)' 0 0 0];
    otherwise
      error('nAveragePoints must be 5 or 7')
  end

  for j=1:ndata+1
    out(j,col) = w_ave(dtmp(j:j+nAveragePoints-1),nAveragePoints);
  end
end
% Make sure we do return matrix of the same size
out = out(ind,:);

function av = w_ave(x,nAveragePoints)
switch nAveragePoints
  case 5
    m = [.1 .25 .3 .25 .1];
  case 7
    m = [.07 .15 .18 .2 .18 .15 .07];
  otherwise
    error('nAveragePoints must be 5 or 7')
end
cor = sum(m(x==0)); % find missing points==0
if cor==1, av = 0;
else, av = sum(x.*m)/(1-cor);
end
