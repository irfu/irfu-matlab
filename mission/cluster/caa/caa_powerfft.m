function [outspecrec,outPxx,outF] = caa_powerfft(data,nfft,sfreq,overlap)
%CAA_POWERFFT  compute power spectrum
%
% [t,power,f] = caa_powerfft(data,nfft,sfreq,[overlap])
% [specrec] = caa_powerfft(data,nfft,sfreq,[overlap])
%	SPECREC is a structure:
%		SPECREC.T - time
%		SPECREC.P - spectrum
%		SPECREC.F - frequency
%
% See also FFT, CAA_SPECTROGRAM
%

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

disp('');
disp('WARNING!!!!!')
disp('');
disp('caa_powerfft is replaced with irf_powerfft');
disp('caa_powerfft will be removed in the near future');
disp('');


narginchk(3,4)
if nargin<4, overlap = 0; end
if overlap<0 || overlap>100, error('OVERLAP must be in a range 0..99'), end

ii = find(~isnan(data(:,1)));
if isempty(ii), error('time is NaN')
else, ts = data(ii(1),1);
end

% Number of intervals must be computed from time
nint = fix(((data(ii(end),1)-ts)*sfreq+1)*(1+overlap*.01)/nfft);
ncomp = size(data,2) - 1;

% Check if there is enough data
if nint<1
  outF = [];
  outPxx = [];
  outspecrec = [];
  return
end

if nfft/2==fix(nfft/2), nf = nfft/2;
else, nf = nfft/2 + 1;
end
specrec.f = sfreq*((1:nf) -1)'/nfft;
for jj=1:ncomp, specrec.p(jj) = {zeros(nint,nf)}; end
specrec.t = zeros(nint,1);

w = hanning(nfft);
wnorm = sum(w.^2)/nfft;	% normalization factor from windowing

tcur = ts;
nnorm = 2.0/nfft/sfreq/wnorm;

for jj=1:nint
  specrec.t(jj) = tcur+(nfft-1)/sfreq*.5;
  X = order_data(irf_tlim(data,tcur, tcur+(nfft-1)/sfreq),nfft,sfreq,tcur);
  if isempty(X), for comp=1:ncomp,specrec.p{comp}(jj,:) = NaN;end
  else
    for comp=2:ncomp+1
      if any(~isnan(X(:,comp)))
        ff = fft(detrend(X(:,comp)) .* w,nfft);
        pf = ff .*conj(ff) *nnorm;
        specrec.p{comp-1}(jj,:) = pf(1:nf);
      else, specrec.p{comp-1}(jj,:) = NaN;
      end
    end
  end
  tcur = tcur + (1-overlap*.01)*(nfft-1)/sfreq;
end

if nargout==1, outspecrec = specrec;
else
  outspecrec = specrec.t;
  outPxx = specrec.p;
  outF = specrec.f;
end

% Help function to clear datagaps
% We throw away intervals with less than 90% of data
function out = order_data(in,ndata,sfreq,ts)
if isempty(in), out = []; return, end
ncomp = size(in,2);
out = ones(ndata,ncomp)*NaN;
out(:,1) = linspace(ts,ts+(ndata-1)/sfreq,ndata)';
ind = round((in(:,1)-ts)*sfreq+1);
out(ind,2:ncomp) = in(:,2:ncomp);
for comp=2:size(in,2)
  ii = find(isnan(out(:,comp)));
  if length(ii)>ndata*.1
    out(:,comp) = NaN;
  else
    m = mean(out(~isnan(out(:,comp)),comp));
    out(ii,comp) = m;
  end
end

