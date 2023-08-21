function res = caa_filter_e(data,wind)
%CAA_FILTER_E  simple high pass filter EFW E
%
% res = caa_filter_e(data,[window])
%	simple high pass filter, will removemoving average over WINDOW
%	WINDOW in seconds, delault 0.5 sec (2 Hz)
%
% See also: IRF_FILT
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------


MAXDATA = 100000;

if nargin<2, wind=.5; end

ii = find( ~isnan(data(:,2)) );

if length(data(ii,1))<=1
  irf_log('proc','interval too short')
  res=[];
  return
end

sf = c_efw_fsample(data(ii,1),'hx');
if ~sf, error('no sampling frequency'),end
nw2 = ceil(sf*wind/2);

% Order data
ndata = round((data(end,1)-data(1,1))*sf+1);

nkomp = size(data,2)-1;
value = 1e10;

ind = round((data(:,1)-data(1,1))*sf+1);

if 0
  E = ones(ndata,nkomp)*value; %#ok<UNRCH>
  E(ind,:) = data(:,2:end);
  E( E(:,1)==value ,:) = NaN;

  ttt = ones(nw2,nkomp)*NaN;
  E_tmp = [ttt; E; ttt];
  clear ttt

  n_start = 1;
  if ndata>MAXDATA; n_end = MAXDATA;
  else, n_end = ndata;
  end

  EE = zeros(n_end,2*nw2+1);
  while n_start<=ndata
    for komp = 1:nkomp
      for k=n_start:n_end, EE(k-n_start+1,:) = E_tmp(k:k+nw2*2,komp); end
      m = mean(EE(1:n_end-n_start+1,:),2);
      E(n_start:n_end,komp) = E(n_start:n_end,komp) - m;
    end
    if n_end==ndata, break, end
    n_start = n_end + 1;
    n_end = n_end + MAXDATA;
    if n_end>ndata, n_end = ndata; end
  end
  clear EE
  res = data;
  res(:,2:end) = E(ind,:);
else
  E = ones(ndata,nkomp+1)*value;
  E(:,1) = linspace(data(1,1),data(1,1)+(ndata-1)/sf,ndata)';
  E(ind,2:end) = data(:,2:end);
  E( E(:,2)==value ,2:end) = NaN;

  ttt = E(:,2);
  ttt(~isnan(ttt)) = 1;
  ttt(isnan(ttt)) = 0;
  ii = irf_find_diff(ttt);

  if isempty(ii) && ttt(1)==0
    irf_log('proc','data is NaN')
    res = [];
    return
  end

  if isempty(ii)
    % No gaps
    clear E
    res = irf_filt(data,1/wind,0,sf,3);
    return
  end

  if ttt(ii(1))==0, ii = [1; ii]; end
  if ttt(ii(end))==1, ii = [ii; length(ttt)]; end
  clear ttt

  for in=1:length(ii)-1
    if (ii(in+1)-1 - ii(in))/sf<wind, E(ii(in):ii(in+1)-1,2:end) = NaN;
    else
      E(ii(in):ii(in+1)-1,:) = irf_filt(E(ii(in):ii(in+1)-1,:),1/wind,0,sf,3);
      % Put edges to NaNs
      E(ii(in):ii(in)+nw2,2:end) = NaN;
      E(ii(in+1)-1-nw2:ii(in+1)-1,2:end) = NaN;
    end
  end
  res = E(ind,:);
end

