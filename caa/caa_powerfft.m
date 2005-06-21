function [Pxx,F] = caa_powerfft(data,nfft,sfreq,overlap)
%CAA_POWERFFT  compute fft
%
% [Pxx,F] = caa_powerfft(data,nfft,sfreq,[overlap])
%
% See also FFT, CAA_SPECTROGRAM
%
% $Id$

% Copyright 2005 Yuri Khotyaintsev

error(nargchk(3,4,nargin))
if nargin<4, overlap = 0; end
if overlap<0 | overlap>100, error('OVERLAP must be in a range 0..99'), end

nint = fix(size(data,1)*(1+overlap*.01)/nfft);
ncomp = size(data,2) - 1;

% Check if there is enough data
if nint<1,
	F = [];
	Pxx = [];
	return
end

if nfft/2==fix(nfft/2), nf = nfft/2;
else, nf = nfft/2 + 1;
end
F = sfreq*((1:nf) -1)'/nfft;
Pxx = zeros(nint,ncomp+1,nf);

w = hanning(nfft);
wnorm = sum(w.^2)/nfft;	% normalization factor from windowing

ii = find(~isnan(data(:,1)));
if isempty(ii), error('time is NaN')
else, ts = data(ii(1),1);
end
tcur = ts;
nnorm = 2.0/nfft/sfreq/wnorm;

for jj=1:nint
	Pxx(jj,1,:) = tcur+(nfft-1)/sfreq*.5;
	X = order_data(irf_tlim(data,tcur, tcur+(nfft-1)/sfreq),nfft,sfreq,tcur);
	if isempty(X), Pxx(jj,comp,:) = NaN;
	else
		for comp=2:ncomp+1
			if ~isempty(find(~isnan(X(:,comp))))
				ff = fft(detrend(X(:,comp)) .* w,nfft);
				pf = ff .*conj(ff) *nnorm;
				Pxx(jj,comp,:) = pf(1:nf);
			else, Pxx(jj,comp,:) = NaN;
			end
		end
	end
	tcur = tcur + (1-overlap*.01)*(nfft-1)/sfreq;
end

% Help function to clear datagaps
% We throw away intervals with less than 20% of data
function out = order_data(in,ndata,sfreq,ts)
	if isempty(in), out = []; return, end
	ncomp = size(in,2);
	out = ones(ndata,ncomp)*NaN;
	out(:,1) = linspace(ts,ts+(ndata-1)/sfreq,ndata)';
	ind = round((in(:,1)-ts)*sfreq+1);
	out(ind,2:ncomp) = in(:,2:ncomp);
	for comp=2:size(in,2)
		ii = find(isnan(out(:,comp)));
		if length(ii)>ndata*.2
			out(:,comp) = NaN;
		else
			m = mean(out(find(~isnan(out(:,comp))),comp));
			out(ii,comp) = m;
		end
	end
	
