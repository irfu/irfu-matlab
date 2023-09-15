function [res] = sdp_wash_spin(data,method,fillGap)
%MMS.SDP_WASH_SPIN  Wash spin tones
%
% RES = C_EFW_WASH_SCP(DATA,[METHOD,GAPFILL])
%   returns RES that is washed satellite potential DATA
%
% Washing requires filling the data gaps defined by METHOD:
%           'PREVSPIN' - copy data from previous spin [default]
%           'LINEAR'   - linear interpolation through the gap
%
% RES contains washed data with original gaps preserved unless GAPFILL=1
% specified.
%
% MMS.SDP_WASH_SPIN(DATA,[METHOD,GAPFILL])
%   when nargout==0 plot the data

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if ~isa(data,'TSeries'), error('TSeries input expected'), end
if nargin < 2, method = 'prevspin'; end
if nargin < 3, fillGapFlag = false;
else
  fillGapFlag = any(fillGap);
end

SPIN_PERIOD = 60/3.1;

res = data; epoch = data.time.epochUnix;
fSample = 1/median(diff(epoch));
jj = [];
ii = find( isnan(data.data(:,1)) );

% Fill gaps before spectral cleaning
if ~isempty(ii)
  if length(ii) == length(data), return, end

  % First point is already a gap
  if ii(1) == 1
    jj = find(~isnan(data.data(:,1)));
    res.data(1:jj(1)-1,:) = res.data(jj(1),:);
    ii(1:jj(1)-1) = [];
  end
end

% Data ends by a gap
if ~isempty(ii) && ii(end) == size(data,1)
  if isempty(jj), jj = find(~isnan(data.data(:,1))); end
  res.data(jj(end)+1:end,:) = res.data(jj(end),:);
  ii(jj(end)+1:end) = [];
end

if ~isempty(ii)
  ju = find(diff(ii)-1);
  begs = [ii(1)-1; ii(ju+1)-1];
  ends = [ii(ju)+1; ii(end)+1];
  switch lower(method)
    case 'linear'
      for idx=1:length(begs)
        res.data(begs(idx)+1:ends(idx)-1,:) = interp1(...
          [epoch(begs(idx)); epoch(ends(idx))],...
          [res(begs(idx),2); res(ends(idx),2)],...
          epoch(begs(idx)+1:ends(idx)-1),'linear');
      end
    case 'prevspin'
      for idx=1:length(begs)
        prevs = res.data( epoch >= epoch(begs(idx)+1)-SPIN_PERIOD, :);
        res.data(begs(idx)+1:ends(idx)-1,:) = ...
          prevs(1:ends(idx)-begs(idx)-1,:);
      end
    otherwise
      error('METHOD can be : LINEAR, PREVSPIN')
  end

end

for iComp = 1:3
  res.data(:,iComp) = clean_spec(res.data(:,iComp),4,0.75,nargout);
  res.data(:,iComp) = clean_spec(res.data(:,iComp),4,0.2,nargout);
  res.data(:,iComp) = clean_spec(res.data(:,iComp),4,0.2,nargout);
end

if ~fillGapFlag
  res.data(isnan(data.data(:,1)),:) = NaN;
end

  function res = clean_spec(res,k1,k2,nn)
    % nn if nn=0, do not plot, if nn=1 plot (nn is nargout of calling c_efw_wash_scp)
    m = mean(res);
    if isnan(m), error('data still contains NaNs'), end
    f = fft(res - m);
    nff = length(f);
    if(rem(nff,2)==0), kfft=nff/2 +1;
    else, kfft=(nff+1)/2;
    end
    freq = fSample*(1:(kfft-1))'/nff;
    a = abs(f(1:kfft-1));
    if nn==0
      clf,plot(freq,a,'k.-'), hold on
    end
    cls = 'brgmkbrgm';
    for n = 1:8
      fref = 1/SPIN_PERIOD * n;
      ip = find(freq>fref*4/5 & freq<fref*6/5); lip8 = fix(length(ip)/8);
      iiii = find( a > mean(a(ip))+(k1+(n-1)*k2)*std(a(ip)) );
      iiii = iiii( iiii>ip(1)+lip8 & iiii<ip(end)-lip8);
      if isempty(iiii), continue, end
      if nn==0
        plot(freq(iiii),a(iiii),[cls(n) '*'])
      end
      f(iiii) = mean(f(ip));
    end
    if nn==0
      hold off
    end
    res = ifft(f,'symmetric') + m;
  end

end




