function [bout,ii,ttGap] = th_clean_eb(bs,mode)
%TH_CLEAN_EB  Basic leaning of data (spikes)
%
%  [bCleaned,idxGap,ttGap] = th_clean_eb(bs,mode)
%
%  Input: 
%         bs - time series
%       mode - mode for filling the spikes LINEAR(default), ZERO, NaN
%
%  Output:
%      bCleaned - data with filled spikes
%        idxGap - indeces of the spikes
%         ttGap - time table corresponding to the spikes

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<2, mode = 'linear'; end
if nargout==3, ttGap=irf.TimeTable(); end
if isempty(bs), bout = []; ii = []; return, end

DT2=60;
DEV_MAX = 60;
bMedian = bs;
t = bs(:,1); nData = length(t);
for i=2:nData-1
  tMin = t(i)-DT2; tMax = t(i)+DT2;
  if tMin<t(1), tMin = t(1); tMax = 2*t(i)-t(1);
  elseif tMax>t(end),
    tMax = t(end); tMin=2*t(i)-t(end);
  end
  bMedian(i,2:4) = median(bs(t>=tMin & t<=tMax,2:4));
end

% Find Range change between 200 and 500 nT using polyfit
bTmp = irf_abs(bs);
ii = find(bTmp(:,5)>200 & bTmp(:,5)<500);
t0 = bTmp(ii(1),1); X = bs(ii,1)-t0; idxModeCh = [];
for iComp=2:4
  Y = bs(ii,iComp);
  P = polyfit(X,Y,3);
  F = polyval(P,X); E=abs(Y-F); idxModeCh = [idxModeCh; find(E>5*median(E))];
end
if ~isempty(idxModeCh)
  idxModeCh = sort(idxModeCh);
  % Leave only the jumps ocurring in at least two components
  idxModeCh = unique(idxModeCh(diff(idxModeCh)==0));
  idxJ = find(diff(idxModeCh)>1);
  idxModeCh = ii(idxModeCh);
  % Display info
  idxEnd = unique([idxModeCh(idxJ+1); idxModeCh(end)]);
  idxSt  = unique([idxModeCh(1); idxModeCh(idxJ)]);
  for iChunk=1:length(idxSt)
    irf.log('warning',...
      sprintf('Range change jump at %s (%.1f nT, %d points)',...
      epoch2iso(bs(idxSt(iChunk),1)),bTmp(idxSt(iChunk),5),idxEnd(iChunk)-idxSt(iChunk)))
  end
end

% Clean spikes
db = bMedian(:,2:4)-bs(:,2:4);
normb=db(:,1).^2+db(:,2).^2+db(:,3).^2;
ii = find(normb>DEV_MAX);

if ~isempty(ii)
  ii = sort(unique([ii; ii-1; ii-2; ii+1; ii+2; ii+3; ii+4; ii+5]));
  ii(ii<1) = []; ii(ii>nData) = [];
end
ii = sort(unique([ii; idxModeCh]));
if ~isempty(ii)  
  switch lower(mode)
    case 'linear'
      % Fill gaps
      x = bs; x(ii,:) = [];
      bout = [t interp1(x(:,1),x(:,2:end),t)];
    case 'zero'
      bout = bs; bout(ii,2:end) = 0;
    case 'nan'
      bout = bs; bout(ii,2:end) = NaN;
    otherwise
      error('Unknown methos. must be one of LINEAR, ZERO, NaN')
  end
  
  if nargout<3, return, end
  idxGap=[1; find(diff(ii)>1)+1];
  for i=1:length(idxGap)
    iiStart = ii(idxGap(i));
    if i==length(idxGap), iiEnd = ii(end);
    else iiEnd = ii(idxGap(i+1)-1);
    end
    ttGap = add(ttGap,[t(iiStart) t(iiEnd)],{},...
      sprintf('gap #%d : %d points',i,iiEnd-iiStart+1));
  end
else
  bout = bs;
end
