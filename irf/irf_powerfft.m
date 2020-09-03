function [outSpecrecOrT,outPxx,outF] = irf_powerfft(data,nFft,samplFreqHz,overlapPercent,smoothWidth)
%IRF_POWERFFT  compute power spectrum
%
% [Specrec]        = irf_powerfft(data,nFft,samplFreqHz,[overlapPercent,smoothWidth])
% [t,p,f]          = irf_powerfft(data,nFft,samplFreqHz,[overlapPercent,smoothWidth])
%       Return the same results as either one struct, or as separate return
%       values. See "Return value".
%
% (no return value)  irf_powerfft(data,nFft,samplFreqHz,[overlapPercent, smoothWidth])
%       Uses results to call irf_spectrogram.
%
%
% ARGUMENTS
% =========
%   data           - Either
%                     (1) a vector (first column is time in seconds; remaining
%                         columns are data), or
%                     (2) a TSeries object. May contain one sample per
%                         timestamp, or a 1D vector of samples per timestamp.
%   nFft           - Length of time covered by every spectrum, as measured by
%                    the number of samples when sampled at exactly sampling
%                    frequency samplFreqHz. This is also the number of samples
%                    (real or reconstructed) used for every spectrum.
%   samplFreqHz    - The sampling frequency [samples/s].
%   overlapPercent - Overlap between spectras in percent (0..100-).
%   smoothWidth    -
%
%
% RETURN VALUE(S)
% ===============
%	Specrec - Structure with below fields:
%		.t(iTime)
%          - Numeric column array. Time of each individual spectrum [seconds].
%            If "data" is a numeric matrix, then same time format as data time.
%            If "data" is a TSeries, then time format epoch unix.
%		.p{iDataComp}(iTime, iFreq)
%          - Cell array of same-sized numeric 2D arrays of spectrum values.
%            NaN if too little underlying data for the particular spectrum.
%            iDataComp = Which data component out of multiple samples per
%            timestamp.
%		.f(iFreq)
%          - Numeric column array. Frequencies [Hz] for spectrum. One array
%            for all spectras.
%
%       Note: The structure can be passed to IRF_SPECTROGRAM with/without
%       modifications.
%       Note: IRF_SPECTROGRAM interprets .t as epoch unix, implying that if
%       "data" is a numeric array, then its time format must also be epoch unix.
%
%
% Note: If a sample=NaN is used for deriving a specific spectrum, then the
% entire spectrum is set to NaN.
%
% Note: Input data does not have to be sampled at samplFreqHz, and sampling
% frequency may vary over time. To derive spectra, shorter sequences of samples
% are sent to FFT. These sequences of samples are at exactly sampling frequency
% samplFreqHz, and are constructed from the submitted samples, by e.g. adjusting
% the timestamps.
% * When the actual sampling rate is HIGHER than samplFreqHz, then some samples
%   are dropped from (not used for) calculating the given spectrum.
% * When the actual sampling rate is LOWER than samplFreqHz, then missing
%   values are reconstructed by using the mean value. If the percentage of
%   missing samples for a given spectrum is over a set (hardcoded)
%   treeshold, then only NaN is used. ==> That entire spectrum is set to NaN.
% 
%   WARNING: When the actual sampling rate is sufficiently much lower than
%   samplFreqHz, then there will be so many data gaps that all NaN is used for
%   all samples, rendering the affected spectras NaN, despite the presence of
%   valid data. To resolve this problem, the caller can
%   (1) split up the samples into sequences of approximately constant sampling
%       frequency and call IRF_POWERFFT (and optionally IRF_SPECTROGRAM) once
%       for each such sequence, or
%   (2) use the lowest sampling frequency in the data, sacrificing the
%       higher-frequency parts of the spectras.
%
% See also FFT, IRF_SPECTROGRAM

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% ~BUG: When multiple samples map to the same timestamp in sample sequences used
% for FFT, only one (the last one) will be used. Could instead be e.g. average,
% or nearest sample in time.
% BUG: Special case (early exit) for nSpectras<1 does not consider the actual
% number of expected return values.



narginchk(3,5);
if nargin<4
  overlapPercent = 0;
  smoothWidth    = 0;
elseif nargin<5
  smoothWidth    = 0;
end

assert(isscalar(overlapPercent))
assert( 0<=overlapPercent || overlapPercent<100, ...
    'Argument overlapPercent is not in the range 0..100-.')
overlap = overlapPercent * 0.01;
assert(isscalar(samplFreqHz), 'Argument samplFreqHz is not scalar.')

% IMPLEMENTATION NOTE: If samplFreqHz is single-precision (or integer,
% probably), rather than double-precision, then this can lead to numerical
% problems. Historically, this had lead to spectrumTimeSec1 not incrementing.
samplFreqHz = double(samplFreqHz); 



%=======================================================================
% Normalize the input data to ONE data format used throughout the code
% --------------------------------------------------------------------
% timeSecArray : Nx1. Time stamps in seconds, using "arbitrary" epoch.
% samples      : NxM.
%=======================================================================
if isa(data,'TSeries')
  timeSecArray = data.time.tts;   % TT2000, but in seconds.
  samples      = data.data;
  usingTSeries = true;
elseif ismatrix(data) && isnumeric(data) && (size(data, 2)>=1)
  timeSecArray = data(:, 1);
  samples      = data(:, 2:end);
  usingTSeries = false;
else
  error('Argument "data" is neither (1) TSeries, nor (2) numeric 2D array with at least one column.')
end
% NOTE: samples(iTime, i)
assert(ndims(samples) <= 2, ...
    ['Argument "data" must be at most 1D per timestamp but is not.', ...
    ' This function is not yet designed to handle this case (yet).'])
clear data



% Remove samples with timestamp=NaN.
iKeep = find(~isnan(timeSecArray));
if isempty(iKeep)
  error('All timestamps are NaN.');
end
timeSecArray = timeSecArray(iKeep);
% NOTE: Limits the dimensionality of data samples. Therefore uses some extra
% ":". However, other aspects of this code also limits the dimensionality even
% more.
samples      = samples(iKeep, :, :, :, :, :);



timeSec1  = timeSecArray(1,   1);
timeSec2  = timeSecArray(end, 1);
nComp     = size(samples, 2);
% fix() : Round toward zero.
% nSpectras = Number of time intervals for which spectras should be made.
% STI = Spectrum Time Interval
nSti         = fix(((timeSec2-timeSec1)*samplFreqHz+1) * (1+overlap)/nFft);
stiFirstSec1 = timeSec1 - 0.5/samplFreqHz;   % Time of beginning of first STI.
stiLengthSec = nFft/samplFreqHz;             % Length of one STI.

% Check if there is enough data. If not, then EXIT early.
if( nSti<1 )
  outF          = [];
  outPxx        = [];
  outSpecrecOrT = [];
  return
end

% Number of frequencies in spectrum. ~nFft/2 due to only considering real-valued
% samples.
nFreqs = ceil(nFft/2);



%==========================================================
% Initialize Specrec: Set frequencies + pre-allocate other
%==========================================================
Specrec.f = samplFreqHz*((1:nFreqs) -1)'/nFft;
for iComp = 1:nComp
  Specrec.p(iComp) = {zeros(nSti, nFreqs)};   % Pre-allocate
end
Specrec.t = zeros(nSti,1);    % Pre-allocate



% NOTE: Using "hanning" function and not "hann" or "hamming".
w     = hanning(nFft);
wnorm = sum(w.^2)/nFft;	      % Normalization factor from windowing
nnorm = 2.0/nFft/samplFreqHz/wnorm;



%==============================================
% Iterate over STIs (spectrum time intervals)
%==============================================
% IMPLEMENTATION NOTE: Empirically, calling select_preprocess_data in a separate
% loop in advance speeds up the code (cuts execution time by ~11% in SolO CWF
% test, ~43 s --> ~37 s; ~12% in SolO SWF test that calls for every snapshot).
% Unknown why. RAM caching?
% tTt = tic;
samplesStiCa = cell(nSti, 1);
for iSti = 1:nSti
  stiSec1 = stiFirstSec1 + (1-overlap)*stiLengthSec*(iSti-1);
  samplesStiCa{iSti} = select_preprocess_data(timeSecArray, samples, stiSec1, samplFreqHz, nFft);
end
% toc(tTt)
%==============================================
% Iterate over STIs (spectrum time intervals)
%==============================================
% tTt = tic;
for iSti = 1:nSti
  stiSec1 = stiFirstSec1 + (1-overlap)*stiLengthSec*(iSti-1);
  
  Specrec.t(iSti) = stiSec1 + stiLengthSec*0.5;   % Center of time interval
  if usingTSeries
      Specrec.t(iSti) = EpochUnix.from_ttns(int64(Specrec.t(iSti) * 1e9));
  end
  
  %samplesSti = select_preprocess_data(timeSecArray, samples, stiSec1, samplFreqHz, nFft);
  samplesSti = samplesStiCa{iSti};
  
  for iComp = 1:nComp
    if( ~isempty(samplesSti) && all(~isnan(samplesSti(:,iComp))) )
      
      %==============
      % FFT + window
      %==============
      ff = fft(detrend(samplesSti(:,iComp)) .* w, nFft);
      pf = ff .*conj(ff) * nnorm;
      
      Specrec.p{iComp}(iSti,:) = pf(1:nFreqs);
    else
      % There is at least one NaN in the data that underlies the spectrum.
      % ==> Use NaN for entire spectrum.
      Specrec.p{iComp}(iSti,:) = NaN;
    end
  end
  
end    % for iSti
% toc(tTt)



if smoothWidth
  if (smoothWidth > 2/samplFreqHz)
    Specrec = smoothSpectrum(Specrec, smoothWidth, samplFreqHz, nFreqs, nSti, nComp);
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



% Extract a sequence of samples, representing samples at equidistant points in
% time (sampling frequency) that FFT can be applied to. Handle data gaps and
% varying sampling rate (according to time stamps).
%
% NOTE: Code assumes that "samples" is at most 2D.
%
% ARGUMENTS
% =========
% timeSecArray
% samples
% intervalSec1 : Timestamp of beginning of time interval for which to extract a
%                sequence of samples.
%                NOTE: End of time interval is derived from this argument
%                combined with nSamplesOut and samplFreqHz.
% samplFreqHz
% nSamplesOut  : Number of samples (per component) in output.
%
% RETURN VALUE
% ============
% samplesOut  : Vector of size N x nSamplesOut. Each dimension 1 index corresponds
%               to constant time increments, according to samplFreqHz.
%
function samplesOut = select_preprocess_data(timeSecArray, samples, intervalSec1, samplFreqHz, nSamplesOut)
  % Only keep samples based on timestamps
  % -------------------------------------
  % IMPLEMENTATION NOTE: Only keeping samples based on timestamps first speeds
  % up the code substantially since rounding actually seems to take a
  % considerable amount of time (it would otherwise be done repeatedly for the
  % same samples, when extracting data for different spectrums).
  % IMPLEMENTATION NOTE: Indexing using iKeep = find(bKeep) (non-logical
  % indexing), instead of bKeep (logical indexing), also seems to speed up code
  % somewhat, maybe.
  intervalSec2 = intervalSec1 + nSamplesOut*samplFreqHz;
  bKeep        = (intervalSec1 <= timeSecArray) & (timeSecArray <= intervalSec2);
  %iKeep        = find(bKeep);
  %iKeep        = find(bKeep, 1, 'first') : find(bKeep, 1, 'last');
  timeSecArray = timeSecArray(bKeep);
  samples      = samples(     bKeep, :);
  
  nTimestamps = size(samples, 1);
  nComp       = size(samples, 2);
  
  % Find mapping between (original) sample indices and out sample indices.
  % NOTE: This is NOT a 1-to-1 mapping. Zero, one, or multiple input indices
  %       may be mapped to the same output index. Mapping multiple samples to
  %       one output sample and then arbitrarily picking only one of those input
  %       samples is in principle suboptimal (averaging is better?) but is
  %       probably OK for most applications. Improve?!
  iIn  = [1:nTimestamps]';
  iOut = round((timeSecArray-intervalSec1)*samplFreqHz + 0.5);
    
  % Ensure that code only extracts the desired data, and only assigns
  % the desired (and legal) indices.
  % NOTE: Historically, many bugs have been associated with getting this wrong.
  b = (1 <= iOut) & (iOut <= nSamplesOut);
  iIn( ~b) = [];
  iOut(~b) = [];
  
  samplesOut          = NaN(nSamplesOut, nComp);   % Pre-allocate as NaN, to be on the safe side.
  samplesOut(iOut, :) = samples(iIn, :);
    
  samplesOut          = replace_NaN(samplesOut);
end



% Replace NaN with the mean of non-NaN values, unless there are too many NaN in
% which case NaN will be used anyway for all values (all timestamps for the same
% component).
function samples = replace_NaN(samples)
  MAX_FRACTION_NAN = 0.1;
    
  bNan = isnan(samples);
  if any(bNan)
    % NOTE: One mean per component/channel.
    meanArray = irf.nanmean(samples, 1, 1-MAX_FRACTION_NAN);
    
    for iComp = 1:size(samples, 2)
      if(any(bNan(:,iComp)))
        samples(bNan(:,iComp), iComp) = meanArray(iComp);
      end
    end
  end
end



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



% Help function to handle data gaps. Given a matrix/TSeries of samples to base
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
% function outData = preprocess_data(inData, nSamplesOut, samplFreqHz, tIntervalStart)
%     
%   MAX_FRACTION_NAN = 0.1;
% 
%   if isa(inData,'TSeries')
%       
%     if isempty(inData.data)
%       outData = [];
%       return
%     end
%     nComp = size(inData.data, 2);
% 
%     % Construct "outData" vector/matrix with NaN, except where there is actual data.
%     % iOut = Indices into final vector that can be filled with actual samples.
%     % NOTE: Using .tts instead of .epochUnix, since it is considerably faster.
%     
%     iIn  = [1:length(inData.time)]';
%     iOut = round((inData.time.tts - tIntervalStart.tts) * samplFreqHz + 0.5);
%     
%     % SPECULATIVE BUGFIX: Remove indices representing timestamps outside designated time interval.
%     b = (1<=iOut) & (iOut <= nSamplesOut);
%     iIn( ~b) = [];
%     iOut(~b) = [];
%     
%     outData          = NaN(nSamplesOut, nComp);
%     outData(iOut, :) = inData.data(iIn, :);
%     
%   else
%       
%     if isempty(inData)
%       outData = [];
%       return
%     end
%       
%     nComp         = size(inData,2) - 1;       % Exclude time column
%     iIn           = 1:size(inData, 1);
%     iOut          = round((inData(:,1)-tIntervalStart)*samplFreqHz + 0.5);
%     
%     % SPECULATIVE BUGFIX: Remove indices representing timestamps outside designated time interval.
%     % Should not be necessary depending on boundary handling in select_data().
% %     b = (1<=iOut) & (iOut <= nSamplesOut);
% %     iIn( ~b) = [];
% %     iOut(~b) = [];
%     
%     outData          = NaN(nSamplesOut, nComp);
%     outData(iOut, :) = inData(iIn, 2:end);         % Exclude time column
%     
%   end
%   assert(size(outData, 1) == nSamplesOut)
%     
%   
%   
%   % Replace NaN with the mean of non-NaN values, unless there are too many NaN
%   % In which case NaN will be used anyway.
%   bNan = isnan(outData);
%   if any(bNan)
%     m = irf.nanmean(outData, 1, 1-MAX_FRACTION_NAN);
%     for iCol = 1:nComp
%       if(any(bNan(:,iCol)))
%         outData(bNan(:,iCol),iCol) = m(iCol);
%       end
%     end
%   end
%   
% end    % preprocess_data



% Select subset of data which fits into time window.
%
% NOTE: Older implementation used irf_tlim which is overkill and calls irf_stdt(?) which may parse the specified time
% interval wrong.
%
% data : [timeColumn, ...]
%
% function data = select_data(data, t1, t2)
%   t = data(:, 1);
%   %b = (t1 <= t) & (t <= t2);
%   b = (t1 <= t) & (t < t2);
%   data = data(b, :);
% end