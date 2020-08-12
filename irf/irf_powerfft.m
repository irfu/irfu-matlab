function [outSpecrecOrT,outPxx,outF] = irf_powerfft(data,nFft,samplFreqHz,overlapPercent,smoothWidth)
%IRF_POWERFFT  compute power spectrum
%
% [Specrec]        = irf_powerfft(data,nFft,samplFreqHz,[overlapPercent,smoothWidth])
% [t,power,f]      = irf_powerfft(data,nFft,samplFreqHz,[overlapPercent,smoothWidth])
%       Return the same results as either one struct, or as separate return
%       values.
%
% (no return value)  irf_powerfft(data,nFft,samplFreqHz,[overlapPercent, smoothWidth])
%       Uses results to call irf_spectrogram.
%
%   data           - Either
%                     (1) a vector (first column is time in seconds; remaining
%                         columns are data), or
%                     (2) a TSeries object.
%   nFft           - Length of time covered by every spectrum, when the samples
%                    are sampled at frequency samplFreqHz. This is also the
%                    number of samples (real or reconstructed) used for every
%                    spectrum.
%   samplFreqHz    - The sampling frequency [samples/s].
%   overlapPercent - Overlap between spectras in percent (0..100-).
%   smoothWidth    -
%
%	Specrec        - Structure with below fields:
%		.t(iTime)
%          - Numeric column array. Time of each individual spectrum [seconds].
%            If using TSeries, then epoch unix.
%            If not using TSeries, then same time format as data time.
%		.p{iComponent}(iTime, iFreq)
%          - Cell array of numeric arrays of spectrum values.
%            NaN if too little underlying data for the particular spectrum.
%            iComponent = Which type of samples in data (of multiple samples
%            per timestamp).
%		.f(iFreq)
%          - Numeric column array. Frequencies [Hz] for spectrum. One array
%            for all spectras.
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



narginchk(3,5);
if nargin<4
  overlapPercent = 0;
  smoothWidth    = 0;
elseif nargin<5
  smoothWidth    = 0;
end

if overlapPercent<0 || overlapPercent>=100
  error('Illegal OVERLAP. Must be in the range 0..100-.')
end
overlap = overlapPercent * 0.01;



usingTSeries = isa(data,'TSeries');



if(usingTSeries)
  % CASE: "data" is TSeries.
  nSpectras      = fix(((data.time.stop-data.time.start)*samplFreqHz+1) * (1+overlap) / nFft);
  nComp          = size(data.data, 2);    % Number of components of data.
  tIntervalStart = EpochTT(data.time.start - EpochTT(0.5/samplFreqHz));
  % NOTE: Have the time interval used for a given spectrum begin "0.5 samples"
  % before the first sample, and end "0.5 samples" after the last sample so that
  % the spectrum time interval is proportional to the number of samples
  % (assuming the nominal sampling frequency).
else
  % CASE: "data" is NOT TSeries.
  iFinite = find(~isnan(data(:,1)));
  if isempty(iFinite)
    error('All timestamps are NaN.');
  end
  dataTime1 = data(iFinite(1),   1);
  dataTime2 = data(iFinite(end), 1);
  
  % Number of intervals must be computed from time
  nSpectras      = fix(((dataTime2-dataTime1)*samplFreqHz+1) * (1+overlap)/nFft);
  nComp          = size(data,2) - 1;    % Number of components of data. First column is assumed to be time.
  tIntervalStart = dataTime1 - 0.5/samplFreqHz;
end

% Check if there is enough data. If not, then EXIT early.
if( nSpectras<1 )
  outF          = [];
  outPxx        = [];
  outSpecrecOrT = [];
  return
end

intervalLengthSec = nFft/samplFreqHz;    % Length of time interval used for ONE spectrum.

nFreqs = ceil(nFft/2);   % Number of frequencies in spectrum. ~nFft/2 due to only considering real-valued samples.

Specrec.f = samplFreqHz*((1:nFreqs) -1)'/nFft;
for iComp = 1:nComp
  Specrec.p(iComp) = {zeros(nSpectras, nFreqs)};   % Pre-allocate
end
Specrec.t = zeros(nSpectras,1);    % Pre-allocate



% NOTE: Using "hanning" function and not "hann" or "hamming".
w     = hanning(nFft);
wnorm = sum(w.^2)/nFft;	      % Normalization factor from windowing
nnorm = 2.0/nFft/samplFreqHz/wnorm;



%==================================
% Iterate over individual spectras
%==================================
for iSpectrum = 1:nSpectras
  
  if(usingTSeries)
    % NOTE: Assigning Specrec.t using .epochUnix.
    % Possibly FIXME (Specrec.t to GenericTime once irf_spectrogram can handle TSeries)
    Specrec.t(iSpectrum) = tIntervalStart.epochUnix + intervalLengthSec*0.5;    % Time interval midpoint. Very slow.
    
    TintTT = irf.tint(tIntervalStart, intervalLengthSec);    % Spectrum time interval (start+stop). Very slow?
    dataIntervalRaw     = tlim(data, TintTT);      % Select data points within spectrum time interval.
    dataIntervalPreproc = preprocess_data(...      % "Preprocessed" data, suitable for doing FFT on.
      dataIntervalRaw, ...
      nFft, samplFreqHz, tIntervalStart);
  else
    Specrec.t(iSpectrum) = tIntervalStart + intervalLengthSec*0.5;   % Center of time interval
    
    % irf_tlim seems to misinterpret tEnd and therefore select the wrong time interval (at least for test cases).
%     dataIntervalPreproc = preprocess_data(...
%       irf_tlim(data, tIntervalStart, tIntervalStart+intervalLengthSec), ...
%       nFft, samplFreqHz, tIntervalStart);
    dataIntervalPreproc = preprocess_data(...
      select_data(data, tIntervalStart, tIntervalStart+intervalLengthSec), ...
      nFft, samplFreqHz, tIntervalStart);
  end

  for iComp = 1:nComp
    if( ~isempty(dataIntervalPreproc) && all(~isnan(dataIntervalPreproc(:,iComp))) )
      
      %==============
      % FFT + window
      %==============
      ff = fft(detrend(dataIntervalPreproc(:,iComp)) .* w, nFft);
      pf = ff .*conj(ff) * nnorm;
      
      Specrec.p{iComp}(iSpectrum,:) = pf(1:nFreqs);
    else
      % There is at least one NaN in data underlying spectrum. ==> Use NaN for entire spectrum.
      Specrec.p{iComp}(iSpectrum,:) = NaN;
    end
  end
  
  % Derive next start time (GenericTimeArray adds "double" as seconds)
  tIntervalStart = tIntervalStart + (1-overlap)*intervalLengthSec;
end



if smoothWidth
  if (smoothWidth > 2/samplFreqHz)
      Specrec = smoothSpectrum(Specrec, smoothWidth, samplFreqHz, nFreqs, nSpectras, nComp);
  else
      irf.log('warn','smoothing not done - smoothWidth is too small')
  end
end



%=================================================================================
% Return result, depending on the number of return values specified by the caller
%=================================================================================
if nargout==0
  irf_spectrogram(Specrec)
elseif nargout==1
  outSpecrecOrT = Specrec;
elseif nargout==3
  outSpecrecOrT = Specrec.t;
  outPxx        = Specrec.p;
  outF          = Specrec.f;
else
  error('irf_powerfft: illegal number of output parameters');
end



end    % irf_powerfft



% Help function to handle datagaps. Given a vector/TSeries of samples to base
% a spectrum on (possibly lacking samples, jumps in timestamps), construct a
% "complete" 1D vector that can be used for FFT. Throws away time intervals
% with less than some percentage (90%) of data.
%
% inData   : Data (samples) to construct final vector from. Must only contain
%            samples for the final spectrum. May contain jumps in timestamps
%            (if TSeries).
% nOutData : Size of final vector.
% outData  : Vector of size nOutData x N. Each dimension 1 index corresponds
%            to constant time increments.
%
function outData = preprocess_data(inData, nSamplesOut, samplFreqHz, tIntervalStart)
    
  MAX_FRACTION_NAN = 0.1;

  if isa(inData,'TSeries')
      
    if isempty(inData.data)
      outData = [];
      return
    end
    nComp = size(inData.data, 2);

    % Construct "outData" vector/matrix with NaN, except where there is actual data.
    % iOut = Indices into final vector that can be filled with actual samples.
    % NOTE: Using .tts instead of .epochUnix, since it is considerably faster.
    
    iIn  = [1:length(inData.time)]';
    iOut = round((inData.time.tts - tIntervalStart.tts) * samplFreqHz + 0.5);
    
    % SPECULATIVE BUGFIX: Remove indices representing timestamps outside designated time interval.
    %b = (1<=iOut) & (iOut <= nSamplesOut);
    %iIn( ~b) = [];
    %iOut(~b) = [];
    
    outData          = NaN(nSamplesOut, nComp);
    outData(iOut, :) = inData.data(iIn, :);
    
  else
      
    if isempty(inData)
      outData = [];
      return
    end
      
    nComp         = size(inData,2) - 1;       % Exclude time column
    iIn           = 1:size(inData, 1);
    iOut          = round((inData(:,1)-tIntervalStart)*samplFreqHz + 0.5);
    
    % SPECULATIVE BUGFIX: Remove indices representing timestamps outside designated time interval.
    % Should not be necessary depending on boundary handling in select_data().
%     b = (1<=iOut) & (iOut <= nSamplesOut);
%     iIn( ~b) = [];
%     iOut(~b) = [];
    
    outData          = NaN(nSamplesOut, nComp);
    outData(iOut, :) = inData(iIn, 2:end);         % Exclude time column
    
  end
  assert(size(outData, 1) == nSamplesOut)
    
  
  
  % Replace NaN with the mean of non-NaN values, unless there are too many NaN
  % In which case NaN will be used anyway.
  bNan = isnan(outData);
  if any(bNan)
    m = irf.nanmean(outData, 1, 1-MAX_FRACTION_NAN);
    for iCol = 1:nComp
      if(any(bNan(:,iCol)))
        outData(bNan(:,iCol),iCol) = m(iCol);
      end
    end
  end
  
end    % preprocess_data



% Basic smoothing procedure, borrowed from mms_fft().
% Smoothes the spectras, not the input data.
function Specrec = smoothSpectrum(Specrec, smoothWidth, samplFreqHz, nf, nSpectras, nComp)
    
  nc = floor(smoothWidth*samplFreqHz);
  % Make sure nc is even.
  if (mod(nc,2) == 1)
    nc = nc-1;
  end
    
  idx   = nc/2 : nc : nf-nc/2;
  nFreq = length(idx);
  freqArray = zeros(nFreq,1);
  for ij = 1:nFreq
    freqArray(ij) = mean(Specrec.f(idx(ij)-nc/2+1 : idx(ij)+nc/2-1));
  end
  Specrec.f = freqArray;
    
  powers = zeros(nSpectras,nFreq);
  for iComp2=1:nComp
    for ij = 1:nFreq
      powers(:,ij) = mean(Specrec.p{iComp2}(:, idx(ij)-nc/2+1 : idx(ij)+nc/2-1), 2);
    end
    Specrec.p{iComp2} = powers;
  end
end



% Select subset of data which fits into time window.
%
% NOTE: Older implementation used irf_tlim which is overkill and calls irf_stdt(?) which may parse the specified time
% interval wrong.
%
% data : [timeColumn, ...]
%
function data = select_data(data, t1, t2)
  t = data(:, 1);
  %b = (t1 <= t) & (t <= t2);
  b = (t1 <= t) & (t < t2);
  data = data(b, :);
end