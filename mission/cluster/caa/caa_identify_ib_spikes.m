function [pdata,idx] = caa_identify_ib_spikes(idata)
% CAA_IDENTIFY_IB_SPIKES  Identify spikes in the internal burst data
%
% [PDATA,IDX] = CAA_IDENTIFY_IB_SPIKES(PDATA)
%
% Input:
% PDATA - input in TM units (from mEFWburstTM1.mat)
%
% Output:
% PDATA - data with spikes removed and interpolated by nearest good values
% IDX - indeces if the points with spikes
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

pdata = idata;
pdata(:,isnan(pdata(1,:))) = [];

ndata = size(pdata,1);

PL = 8192; % page length
if ceil(ndata/PL) < 12, PL = PL/2; end

SDTTHRESH = 30; % Threshold for spikes

pos = ndata; idx = []; prevSpike = [];

%figure(73)

while pos-PL >= 1
  pos = pos - PL;
  di = sum( abs(diff( pdata(pos:pos+PL-1,2:end) ,1,1)) ,2);

  %      clf, h(1) = subplot(2,1,1); h(2) = subplot(2,1,2);
  %      plot(h(1), pdata(:,2:end),'.-'), plot(h(2),di,'.-')
  %      hold(h(1),'on'),hold(h(2),'on'),set(h(1),'XLim',pos + [0 PL-1],'YLimMode','auto')
  %      ylabel(h(1),'pdata'), ylabel(h(2),'diff')

  nSpike = 0;
  while nSpike<=3
    OFF_P=2; OFF_M=2;
    if nSpike == 0 && ~isempty(prevSpike) % Spike position indicated by spike on the prev page
      spikeProxy = prevSpike - PL - pos;
      % Check if spikeProxy is too close to the boundary
      if spikeProxy <= OFF_M, spikeProxy = OFF_M + 1;
      elseif spikeProxy >= PL - OFF_P -1, spikeProxy = PL - OFF_P -2;
      end

      if spikeProxy < PL
        if any(di(spikeProxy + (-OFF_M:1:OFF_P)) > 0.25*SDTTHRESH*std(di))
          ii = find( di(spikeProxy + (-OFF_M:1:OFF_P)) == ...
            max(di(spikeProxy + (-OFF_M:1:OFF_P))) ) +spikeProxy -OFF_M -1;
          spike = ii + pos;
          %fprintf('(1) spike at %d\n',spike)
        else
          prevSpike = [];
          continue
        end
      else
        prevSpike = [];
        continue
      end
    else % Other spikes
      mm = max(di);
      if mm > SDTTHRESH*std(di)
        ii = find(di == mm);
        ii = round(mean(ii));
        spike = ii + pos;
        %fprintf('(2) spike at %d\n',spike)
      else
        break
      end
    end

    if isempty(prevSpike), prevSpike = spike; end
    if size(spike,1)==2
      idx_add = [spike(1) + ((-OFF_M+1):1:0) spike(2) + (0:1:(OFF_P-1))]
    else
      idx_add = spike + ((-OFF_M+1):1:(OFF_P-1));
    end
    if size(idx_add,1)>1
      idx_add = idx_add(:)';
    end
    idx = [idx idx_add]; %#ok<AGROW>
    nSpike = nSpike + 1;

    %plot(h(1),idx_add,pdata(idx_add,2:end),'ko')
    %plot(h(2),ii,di(ii),'ko')

    pdata(idx_add,2:end) = ones(length(idx_add),1)*(...
      pdata(min(idx_add)-1,2:end) + pdata(max(idx_add)+1,2:end) )*0.5;
    di = sum( abs(diff( pdata(pos:pos+PL-1,2:end) ,1,1)) ,2);

    %plot(h(1),idx_add,pdata(idx_add,2:end),'m*'), set(h(1),'XLim',[min(idx_add)-5 max(idx_add)+5])
    %plot(h(2),ii,di(ii),'m*'), set(h(2),'XLim',ii + [-10 10])
    %keyboard
  end
end

% N_CONST sets the minimum number of points of constant TM value
% which we consider bad
N_CONST = 4;
% Check for constant TM, means probe is in a strange state
di = sum( diff( pdata(:,2:end) ,1,1) ,2);
ii = find(di==0);
if length(ii)>1
  % at least three consequetive points are the same
  kk = find(ii(1:end-1)-ii(2:end)==-1);
  if ~isempty(kk)
    jj = 1;
    while jj<=length(kk)
      bad_i = find(ii-ii(kk(jj))-(1:length(ii))'+kk(jj)==0);
      if length(bad_i)>N_CONST
        idx = [idx ((ii(bad_i)-1):(ii(bad_i(end))+1))]; %#ok<AGROW>
      end
      if isempty(bad_i), jj = jj + 1;
      else
        ll = find(kk>bad_i(end));
        if isempty(ll), break
        else, jj = ll(1);
        end
      end
    end
  end
end

% First points are usually bad
idx = unique([1:6 sort(idx)]);
idx(idx<1) = [];
idx(idx>ndata) = [];
pdata = idata;
pdata(1:6,2:end) = ones(6,1)*pdata(7,2:end);

% Fill the gaps in the data
di = diff(idx);
ii = find(di > 1);
if ~isempty(ii)
  for i=1:length(ii)
    first = idx(ii(i)+1);
    if i==length(ii)
      last = idx(end);
    else
      last = idx(ii(i+1));
    end
    idx_add = first:1:last;
    if last+1<ndata
      pdata(idx_add,2:end) = ones(length(idx_add),1)*(...
        pdata(first-1,2:end) + pdata(last+1,2:end) )*0.5;
    else
      pdata(idx_add,2:end) = ones(length(idx_add),1)*pdata(first-1,2:end);
    end
  end
end