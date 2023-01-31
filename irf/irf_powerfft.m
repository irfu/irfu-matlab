function [outSpecrecOrT,outPxx,outF] = irf_powerfft(...
  data, nFft, samplFreqHz, overlapPercent, smoothWidth)
%IRF_POWERFFT  compute power spectrum
%
% [Specrec]       = irf_powerfft(data, nFft, samplFreqHz, [overlapPercent, smoothWidth])
% [t,p,f]         = irf_powerfft(data, nFft, samplFreqHz, [overlapPercent, smoothWidth])
%       Return the same results as either one struct, or as separate return
%       values. See "Return value".
%
% (no return value) irf_powerfft(data, nFft, samplFreqHz, [overlapPercent, smoothWidth])
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
%                         Note: (2) is more likely to correctly handle leap
%                         seconds.
%   nFft           - Length of time covered by every spectrum, as measured by
%                    the number of samples when sampled at exactly sampling
%                    frequency samplFreqHz. This is also the number of samples
%                    (real or reconstructed) used for every spectrum.
%   samplFreqHz    - The sampling frequency [samples/s].
%   overlapPercent - Overlap between spectras in percent (0..100-).
%   smoothWidth    - (Unit: seconds?)
%
%
% RETURN VALUE(S)
% ===============
%	Specrec - Structure with below fields:
%		.t(iTime)
%          - Numeric column array. Time of each individual spectrum [seconds].
%            If "data" is a numeric matrix, then same time format as data time.
%            If "data" is a TSeries, then time format is epoch unix.
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
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

% ~BUG: When multiple samples map to the same timestamp in sample sequences used
% for FFT, only one (the last one) will be used. Could instead be e.g. average,
% or nearest sample in time.
% BUG?: Might not handle smoothWidth > 0 for zero spectras.
%
% INTERNAL NAMING CONVENTIONS
% ===========================
% Spec    = A single spectrum (one FFT). May also refer to the time interval
%           from whih it is derived.
% AllData = ~All data (timestamps+samples) underlying all spectras.



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

% IMPLEMENTATION NOTE: If nFft is single-precision (or integer, probably),
% rather than double-precision, then this can lead to numerical problems.
% Historically, this had lead to i1Spec being larger than tAllDataSec indices.
nFft = double(nFft);



%======================================================================
% Normalize the input data (one of two formats) to ONE data format.
% --------------------------------------------------------------------
% tAllDataSec : Nx1. Time stamps in seconds, using "arbitrary" epoch.
% samples      : NxM.
%======================================================================
if isa(data,'TSeries')
  tAllDataSec    = data.time.tts;   % TT2000, but in seconds.
  samplesAllData = data.data;
  usingTSeries   = true;
elseif ismatrix(data) && isnumeric(data) && (size(data, 2)>=1)
  tAllDataSec    = data(:, 1);
  samplesAllData = data(:, 2:end);
  usingTSeries   = false;
else
  error(['Argument "data" is neither (1) TSeries, nor', ...
    ' (2) numeric 2D array with at least one column.'])
end
% NOTE: samplesAllData(iTime, i)
% NOTE: Assertion only relevant when argument "data" is TSeries.
assert(ndims(samplesAllData) <= 2, ...
  ['Argument "data" must be at most 1D per timestamp but is not.', ...
  ' This function is not yet designed to handle this case (yet).'])

% Just to convince the reader of this source code that "data" is not used
% herafter.
clear data



%======================================================
% Remove samples with TIMESTAMP==NaN (not sample==NaN)
%======================================================
iKeep = find(~isnan(tAllDataSec));
if isempty(iKeep)
  error('All timestamps are NaN.');
end
tAllDataSec = tAllDataSec(iKeep);
% NOTE: This limits the dimensionality of data samples. Therefore uses some
% extra ":". However, other aspects of this code also limits the dimensionality
% even more.
samplesAllData = samplesAllData(iKeep, :, :, :, :, :);



assert(issorted(tAllDataSec, 'strictascend'), ...
  ['Timestamps are not monotonically increasing. This may e.g. be due to', ...
  ' leap seconds if the caller is not using TSeries for argument "data".'])



t1AllDataSec = tAllDataSec(1,   1);
t2AllDataSec = tAllDataSec(end, 1);
nTimeAllData = numel(tAllDataSec);
nComp        = size(samplesAllData, 2);
t1FirstSpecSec = t1AllDataSec - 0.5/samplFreqHz;   % Beginning of the FIRST spectrum.
t2LastSpecSec  = t2AllDataSec + 0.5/samplFreqHz;   % End       of the LAST  spectrum.
lenSpecSec     = nFft/samplFreqHz;                 % Length (in time) of one spectrum.
% Approximate (ideal) time between beginnings of successive spectra (help
% variable).
specDistSec = (1-overlap)*lenSpecSec;

% IMPLEMENTATION NOTE: Timestamps might not follow the stated sampling rate
% (the argument) exactly, but only approximately. This is important for
% snapshots where one might want to have exactly one spectrum per snapshot.
% Not rounding may lead to bad choices of number of spectra, e.g. two
% spectras based on exactly the same samples (empirical), or no spectrum
% (potentially) when the snapshot is just slightly too short (timestamps too
% close together; correct number of samples). Therefore, do not use e.g. "ceil".
% Ex: solo_L2_rpw-lfr-surv-swf-e-cdag_20210102_V04.cdf
nSpec       = round((t2LastSpecSec - t1FirstSpecSec - lenSpecSec) / specDistSec + 1);



% Check if there is enough data for at least one spectrum. If not, then EXIT
% early.
% if nSpec<1
%   outF          = [];
%   outPxx        = [];
%   outSpecrecOrT = [];
%   return
% end

% Number of frequencies in a spectrum. ~nFft/2 due to only considering
% real-valued samples.
nFreqs = ceil(nFft/2);



% NOTE: Using "hanning" function (not "hann" or "hamming").
w     = hanning(nFft);
wnorm = sum(w.^2)/nFft;	      % Normalization factor from windowing
nnorm = 2.0/nFft/samplFreqHz/wnorm;



%==========================================================
% Initialize Specrec: Set frequencies + pre-allocate other
%==========================================================
Specrec.f = samplFreqHz*((1:nFreqs) -1)'/nFft;
for iComp = 1:nComp
  Specrec.p(iComp) = {zeros(nSpec, nFreqs)};   % Pre-allocate
end
Specrec.t = zeros(nSpec,1);    % Pre-allocate



%=========================================================
% Iterate over spectra: Extract samples for each spectrum
%=========================================================
% NOTE: Calling select_preprocess_data() is by far the slowest part of this
% function. /Erik P G Johansson 2020-09-04
%
% IMPLEMENTATION NOTE: Empirically, calling select_preprocess_data() in a
% separate loop in advance speeds up the code (cuts total execution time by ~11%
% in SolO CWF test, ~43 s --> ~37 s; ~12% in SolO SWF test that calls for every
% snapshot). Unknown why. RAM caching? /Erik P G Johansson 2020-09-03/04
%
% IMPLEMENTATION NOTE: Changing "for"-->"parfor" speeds up the LOOP:
% 32.596355 s --> 13.654587 s; irony, SolO CWF test. (Numbers assumes that the
% "parallel pool" has already started. /Erik P G Johansson 2020-09-04
%
% IMPLEMENTATION NOTE: Changing "for"-->"parfor" makes code SLOWER for SolO SWF
% quicklook plotting with one IRF_powerfft call for every snapshot.
% 30.592659 s --> 89.429904 s (solo.ql.plot_SWF; not just irf_powerfft).
% /Erik P G Johansson 2020-09-04
% ==> NOT USING "parfor".
%
samplesSpecCa = cell(nSpec, 1);
i1Spec        = int32(1);   % Index to first sample in spectrum iSpec.
%-----------------------------------------------------------------------------
% Pre-calculate beginning and end (in time) of all spectra.
% IMPLEMENTATION NOTE: Storing pre-calculated values, to ensure that both loops
% over spectra use the same values.
%-----------------------------------------------------------------------------
iSpecArray = 1:nSpec;
t1SpecSec = linspace(t1FirstSpecSec, t2LastSpecSec-lenSpecSec, nSpec);
t2SpecSec = t1SpecSec + lenSpecSec;
for iSpec = iSpecArray     % NOT USING "parfor"
  
  %=============================================================================
  % Find index interval i1:i2 for the samples of the current spectrum
  % -----------------------------------------------------------------
  % IMPLEMENTATION NOTE: Doing this outside of select_preprocess_data() while
  % making use of tAllDataSec being sorted (increasing) speeds up this
  % function. If not, then one has to use "find" inside select_preprocess_data()
  % which is slower.
  % /Erik P G Johansson 2020-09-14.
  %=============================================================================
  while tAllDataSec(i1Spec) < t1SpecSec(iSpec)
    i1Spec = i1Spec + 1;
  end
  % NOTE: May have that tAllDataSec(i1Spec) > t2SpecSec.
  % ==> Must begin with t2SpecSec = i1Spec-1.
  i2Spec = i1Spec-1;
  while (i2Spec+1 <= nTimeAllData) && (tAllDataSec(i2Spec+1) <= t2SpecSec(iSpec))
    i2Spec = i2Spec + 1;
  end
  
  % Only use below assertions when debugging/testing (code is slow).
  % NOTE: Even if disabled, the assertions still illustrate conditions that must
  % be satisfied.
  %assert(  find(tAllDataSec >= t1SpecSec, 1, 'first') == i1 )
  %assert( (find(tAllDataSec <= t2SpecSec, 1, 'last')  == i2) || (i2 == nTimeAllData))
  
  %==============================================
  % Extract samples to use for deriving spectrum
  %==============================================
  %samplesSpecCa{iSpec} = select_preprocess_data(tAllDataSec, samplesAllData, t1SpecSec, samplFreqHz, nFft);
  samplesSpecCa{iSpec} = select_preprocess_data3(tAllDataSec, samplesAllData, t1SpecSec(iSpec), samplFreqHz, nFft, i1Spec, i2Spec);
  %samplesSpecCa{iSpec} = select_preprocess_data2(tAllDataSec, samplesAllData, t1SpecSec, samplFreqHz, nFft);
end



%==========================================
% Iterate over spectra: Set Specrec .p, .t
%==========================================
for iSpec = 1:nSpec
  
  % Center of time interval
  Specrec.t(iSpec) = t1SpecSec(iSpec) + lenSpecSec*0.5;
  if usingTSeries
    % IMPLEMENTATION NOTE: Converts time
    % TTNS (TT2000 in seconds; with leap seconds)
    % --> Epoch UNIX (seconds without(?) leap seconds).
    % Could do this time conversion in the normalization stage at beginning of
    % function. Doing it here has historical reasons, plus the slight
    % advantage of correctly(?) handling leap seconds, if the caller supplies
    % a TSeries object. Converting earlier should lead to timestamps during
    % positive leap seconds being represented as NOT incrementing
    % monotonically, which could screw up this function's the algorithms. This
    % might be in vain though since the final time format still does not
    % handle leap seconds correctly. Should lead to graphic error in
    % irf_spectrogram().
    Specrec.t(iSpec) = EpochUnix.from_ttns(int64(Specrec.t(iSpec) * 1e9));
  end
  
  %samplesSpec = select_preprocess_data(tAllDataSec, samplesAllData, t1SpecSec, samplFreqHz, nFft);
  samplesSpec = samplesSpecCa{iSpec};
  
  for iComp = 1:nComp
    if( ~isempty(samplesSpec) && all(~isnan(samplesSpec(:,iComp))) )
      
      %==============
      % FFT + window
      %==============
      ff = fft(detrend(samplesSpec(:,iComp)) .* w, nFft);
      pf = ff .*conj(ff) * nnorm;
      
      Specrec.p{iComp}(iSpec,:) = pf(1:nFreqs);
    else
      % There is at least one NaN in the data that underlies the spectrum.
      % ==> Use NaN for entire spectrum.
      Specrec.p{iComp}(iSpec,:) = NaN;
    end
  end
  
end    % for iSpec



if smoothWidth
  if (smoothWidth > 2/samplFreqHz)
    Specrec = smoothSpectrum(Specrec, smoothWidth, samplFreqHz, nFreqs, nSpec, nComp);
  else
    irf.log('warn','smoothing not done - smoothWidth is too small')
  end
end



%==========================================================================
% Return result, depending on the number of return values specified by the
% caller
%==========================================================================
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
% tAllDataSec  : 1D array. Length nTime.
% samples      : 2D array. Size nTime x nComp.
% intervalSec1 : Timestamp of beginning of time interval for which to extract a
%                sequence of samples.
%                NOTE: End of time interval is derived from this argument
%                combined with nSamplesOut and samplFreqHz.
% samplFreqHz
% nSamplesOut  : Number of samples (per component) in output.
%
% RETURN VALUE
% ============
% samplesOut  : Numeric array.
%       Size nTimeAllData x nSamplesOut. Each dimension 1 index corresponds to
%       constant time increments, according to samplFreqHz.
%
function samplesOut = select_preprocess_data(tSec, samples, intervalSec1, samplFreqHz, nSamplesOut)
% Only keep samples within specified time interval
% ------------------------------------------------
% IMPLEMENTATION NOTE: Only keeping samples based on timestamps before doing
% anything else speeds up the code substantially since rounding actually seems
% to take a considerable amount of time (it would otherwise be done repeatedly
% for the same samples, when extracting data for different spectrums).
%
% IMPLEMENTATION NOTE: Indexing using iKeep = find(bKeep) ("non-logical"
% indexing), instead of bKeep (logical indexing), might speed up code
% somewhat, maybe.

% Using "find" to find indices (instead of bKeep).
% 40-46 s --> 35-42 s (irony; SolO RPW CWF test file)

intervalSec2 = intervalSec1 + nSamplesOut/samplFreqHz;

%=============================================================================
% Implementation 1
%     bKeep   = (intervalSec1 <= tSec) & (tSec <= intervalSec2);
%     tSec    = tSec(bKeep);
%     samples = samples(     bKeep, :);
%=============================================================================
bKeep        = (intervalSec1 <= tSec) & (tSec <= intervalSec2);
tSec = tSec(bKeep);
%samples      = samples(     bKeep, :);

%assert(numel(samples) == numel(tSec))
nTimeSel = size(tSec, 1);
nComp    = size(samples,      2);

%=============================================================================
% Find mapping between (~original) sample indices and out sample indices.
% NOTE: This is NOT a 1-to-1 mapping. Zero, one, or multiple input indices
%       may be mapped to the same output index. Mapping multiple samples to
%       one output sample and then arbitrarily picking only one of those input
%       samples is in principle suboptimal (averaging is better?) but is
%       probably OK for most applications. Improve?!
%=============================================================================

%iIn  = [1:nTimestamps]';
% EXPERIMENTAL. ==> Makes samples=samples(bKeep) unnecessary.
% ASSUMPTION: Argument tSec is sorted (incrementing).
i1   = find(bKeep, 1, 'first'); % Separate for profiling purposes.
iIn  = i1 + [0:nTimeSel-1];

iOut = round((tSec-intervalSec1)*samplFreqHz + 0.5);

% Ensure that code only extracts the desired data, and only assigns
% the desired (and legal) indices.
% NOTE: Historically, many bugs have been associated with getting this wrong.
b = (1 <= iOut) & (iOut <= nSamplesOut);
iIn( ~b) = [];
iOut(~b) = [];

% Pre-allocate as NaN, to be on the safe side.
samplesOut          = NaN(nSamplesOut, nComp);
samplesOut(iOut, :) = samples(iIn, :);

samplesOut          = replace_NaN(samplesOut);
end



% EXPERIMENTAL: For speeding up.
%
% ARGUMENTS
% =========
% i1, i2: Indices into "tSec, samples" representing first and last
%         sample in a time interval.
%
function samplesOut = select_preprocess_data3(...
  tSec, samples, intervalSec1, samplFreqHz, nSamplesOut, i1, i2)
nComp = size(samples, 2);
iIn   = i1:i2;
iOut  = round((tSec(i1:i2)-intervalSec1)*samplFreqHz + 0.5);

b = (1 <= iOut) & (iOut <= nSamplesOut);
iIn( ~b) = [];
iOut(~b) = [];

% Pre-allocate as NaN, to be on the safe side.
samplesOut          = NaN(nSamplesOut, nComp);
samplesOut(iOut, :) = samples(iIn, :);

samplesOut          = replace_NaN(samplesOut);
end



% EXPERIMENTAL: Use interp1 because maybe it is faster (maybe implemented using
% e.g. C).
% PROBLEM: interp1 does NOT return NaN when nearest value is farther than half
% an integration width.
% Slower using interp1 for every spectrum separately ==> One may guess probably
% also slower if used globally). /Erik P G Johansson 2020-09-04
% function specSamples = select_preprocess_data2(tSec, samples, intervalSec1, samplFreqHz, nSamplesOut)
%   intervalSec2 = intervalSec1 + nSamplesOut/samplFreqHz;
%
%   tAllDataSec = linspace(...
%       intervalSec1 + 0.5/samplFreqHz, ...
%       intervalSec2 - 0.5/samplFreqHz, ...
%       nSamplesOut);
%
%   % NOTE: interp1(X,V,Xq, ...) works for length(X) == size(V,1)
%   [specSamples] = interp1(tSec, samples, tAllDataSec, 'nearest');
%
%   %specSamples  = replace_NaN(specSamples);
% end



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
function Specrec = smoothSpectrum(...
  Specrec, smoothWidth, samplFreqHz, nf, nSpectras, nComp)

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