function [outSpecrecOrT,outPxx,outF] = irf_powerfft(data,nFft,samplFreqHz,overlap,smoothWidth)
%IRF_POWERFFT  compute power spectrum
%
% [Specrec]          = irf_powerfft(data,nFft,samplFreqHz,[overlap,smoothWidth])
% [t,power,f]        = irf_powerfft(data,nFft,samplFreqHz,[overlap,smoothWidth])
%       Return the same results as either one struct, or as separate values.
%
% (no return value)    irf_powerfft(data,nFft,samplFreqHz,[overlap])
%       Uses results to call irf_spectrogram.
%   
%   data        - Either (1) a vector (first column is time; other columns are
%                 data), or (2) a TSeries object.
%   nFft        - Length of time covered by every spectrum, when the samples
%                 are sampled at frequency samplFreqHz. This is also the number
%                 of samples (real or reconstructed) used for every spectrum.
%   samplFreqHz - The sampling frequency [samples/s].
%   overlap     - Overlap between spectras in percent (0..99).
%   smoothWidth - 
%	Specrec     - Structure with below fields:
%		.t - Column array. Time of each individual spectrum [seconds].
%            If using TSeries, then epoch unix.
%            If not using TSeries, then same as data time.
%		.p - Spectrum values. NaN if too little underlying data for
%            the particular spectrum.
%		.f - Frequencies [Hz] for spectrum. One array for all spectras.
%       The structure can be passed to IRF_SPECTROGRAM with/without
%       modifications.
%
% Note: Input data does not have to be sampled at samplFreqHz, and sampling
% frequency may vary. Missing samples, needed for FFT, are reconstructed by
% using the mean value is used, assuming not too many samples are missing. If
% the actual sampling rate is higher than samplFreqHz, then samples may be
% dropped.
%
% See also FFT, IRF_SPECTROGRAM

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

MIN_FRACTION_NAN = 0.9;



narginchk(3,5);
if nargin<4
  overlap     = 0;
  smoothWidth = 0;
elseif nargin<5
  smoothWidth = 0; 
end

if overlap<0 || overlap>=100
    error('Illegal OVERLAP. Must be in the range 0..100-.')
end



usingTSeries = isa(data,'TSeries');

if(usingTSeries)
  % CASE: "data" is TSeries.
  nSpectras      = fix(((data.time.stop-data.time.start)*samplFreqHz+1) * (1+overlap*.01) / nFft);
  nComp          = size(data.data, 2);    % Number of components of data.
  tIntervalStart = data.time.start;
else
  % CASE: "data" is NOT TSeries.
  ii = find(~isnan(data(:,1)));
  if isempty(ii)
    error('All timestamps are NaN.');
  end
  ts = data(ii(1),1);
  
  % Number of intervals must be computed from time
  nSpectras      = fix(((data(ii(end),1)-ts)*samplFreqHz+1) * (1+overlap*.01)/nFft);
  nComp          = size(data,2) - 1;    % Number of components of data. First column is time.
  tIntervalStart = ts;
end

% Check if there is enough data
if( nSpectras<1 )
  outF          = [];
  outPxx        = [];
  outSpecrecOrT = [];
  return
end

intervalLengthSec = (nFft-1)/samplFreqHz;    % Length of time interval used for spectrum.

% if nFft/2==fix(nFft/2), nf = nFft/2;
% else, nf = (nFft+1)/2;
% end
nf = ceil(nFft/2);

Specrec.f = samplFreqHz*((1:nf) -1)'/nFft;
for iComp=1:nComp
  Specrec.p(iComp) = {zeros(nSpectras,nf)};
end
Specrec.t = zeros(nSpectras,1);

w     = hanning(nFft);
wnorm = sum(w.^2)/nFft;	      % Normalization factor from windowing
nnorm = 2.0/nFft/samplFreqHz/wnorm;



%==================================
% Iterate over individual spectras
%==================================
for jj = 1:nSpectras
  
  if(usingTSeries)
    % Possibly FIXME (Specrec.t to GenericTime once irf_spectrogram can handle TSeries)
    Specrec.t(jj) = tIntervalStart.epochUnix + intervalLengthSec*0.5;    % Time interval midpoint. Very slow.
    
    TintTT = irf.tint(tIntervalStart, intervalLengthSec);    % Spectrum time interval (start+stop). Very slow?
    dataIntervalRaw     = tlim(data, TintTT);      % Select data points within spectrum time interval.   
    dataIntervalPreproc = preprocess_data(...      % "Preprocessed" data, suitable for doing FFT on.
      dataIntervalRaw, ...
      nFft, samplFreqHz, tIntervalStart);
  else
    Specrec.t(jj) = tIntervalStart + intervalLengthSec*0.5;
    dataIntervalPreproc = preprocess_data(...
      irf_tlim(data, tIntervalStart, tIntervalStart+intervalLengthSec), ...
      nFft, samplFreqHz, tIntervalStart);
  end

  for iComp = 1:nComp
    if( ~isempty(dataIntervalPreproc) && all(~isnan(dataIntervalPreproc(:,iComp))) )
        
      %==============
      % FFT + window
      %==============
      ff = fft(detrend(dataIntervalPreproc(:,iComp)) .* w, nFft);
      pf = ff .*conj(ff) * nnorm;
      
      Specrec.p{iComp}(jj,:) = pf(1:nf);
    else
      Specrec.p{iComp}(jj,:) = NaN;
    end
  end

  % Derive next start time (GenericTimeArray adds "double" as seconds)
  tIntervalStart = tIntervalStart + (1-overlap*.01)*intervalLengthSec;
end



if smoothWidth
  if (smoothWidth > 2/samplFreqHz), smoothSpectrum(); 
  else, irf.log('warn','smoothing not done - smoothWidth too small')
  end
end

if nargout==1, outSpecrecOrT = Specrec;
elseif nargout==3
  outSpecrecOrT = Specrec.t;
  outPxx        = Specrec.p;
  outF          = Specrec.f;
elseif nargout==0
  irf_spectrogram(Specrec)
else
  error('irf_powerfft: unknown number of output parameters');
end



  % Help function to handle datagaps. Given a vector/TSeries of samples to base
  % a spectrum on (possibly lacking samples, jumps in timestamps), construct a
  % "complete" 1D vector that can be used for FFT. Throws away time intervals
  % with less than some percentage (90%) of data.
  % 
  % inData   : Data (samples) to construct final vector from. Must only contain
  %            samples for the final spectrum. May contain jumps in timestamps
  %            (if TSeries).
  % nOutData : Size of final vector.
  % ts       : Start time of final vector.
  % outData  : Vector of size nOutData x N. Each dimension 1 index corresponds
  %            to constant time increments.
  %
  function outData = preprocess_data(inData, nOutData, samplFreqHz, tIntervalStart2)
    if(usingTSeries)
      if isempty(inData.data)
        outData = [];
        return
      end
      nComp2 = size(inData.data, 2);

      % Construct "outData" vector/matrix with NaN, except where there is actual data.
      % ind = Indices into final vector that can be filled with actual samples.
      ind = round((inData.time.epochUnix - tIntervalStart2.epochUnix) * samplFreqHz + 1);   
      outData         = NaN(nOutData, nComp2);
      outData(ind, :) = inData.data;
    else
      if isempty(inData)
          outData = [];
          return
      end

      nComp2          = size(inData,2) - 1;       % Exclude time column
      ind             = round((inData(:,1)-tIntervalStart2)*samplFreqHz + 1);
      outData         = NaN(nOutData, nComp2);
      outData(ind, :) = inData(:, 2:end);         % Exclude time column
    end

    % If data has less than minimum allowed fraction of NaN, then set all
    % components to NaN. Otherwise replace NaN with mean (mean calculated by
    % excl. NaN).
    indNaN = isnan(outData);
    if( any(indNaN) )
      m = irf.nanmean(outData, 1, MIN_FRACTION_NAN);
      for col = 1:nComp2
        if(any(indNaN(:,col)))
          outData(indNaN(:,col),col) = m(col);
        end
      end
    end
    
  end

  
  
  function smoothSpectrum
    % Basic smoothing procedure, borrowed from mms_fft().
    % Smoothes the spectras, not the input data.
    
    nc = floor(smoothWidth*samplFreqHz);
    % Make sure nc is even.
    if (mod(nc,2) == 1), nc = nc-1; end
    
    idx = nc/2:nc:nf-nc/2; nFreq = length(idx);
    freqs = zeros(nFreq,1);
    for ij = 1:nFreq
      freqs(ij) = mean(Specrec.f(idx(ij)-nc/2+1:idx(ij)+nc/2-1));
    end
    Specrec.f = freqs;
    
    powers = zeros(nSpectras,nFreq);
    for iComp2=1:nComp
      for ij = 1:nFreq
        powers(:,ij) = mean(Specrec.p{iComp2}(:,idx(ij)-nc/2+1:idx(ij)+nc/2-1),2);
      end
      Specrec.p{iComp2} = powers;
    end 
  end



end    % irf_powerfft
