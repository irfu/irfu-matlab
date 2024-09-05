%
% Apply transfer function to signal according to chosen algorithm. Potentially
% modify TF and data before and after application of TF.
%
%
% TERMINOLOGY
% ===========
% De-trending : REMOVING fit on data BEFORE applying the TF. It does NOT
%               automatically imply RE-trending.
% Re-trending : ADDING BACK a scaled version of the previously removed fit when
%               de-trending, AFTER applying the TF.
% SNF         : Split by Non-Finite values. Splits time series into separate
%               smaller time series before applying de-trending and (modified)
%               TF.
%
%
% De-trending / Re-trending
% =========================
% NOTE: Detrending makes it impossible to modify the amplitude & phase for the
% frequency components in the trend itself (the fit), e.g. to delay the signal
% in the trend. If the input signal is interpreted as N-periodic (e.g. when
% applying TF using FFT), then de-trending affects the jump between the
% beginning and end of the signal (reduces it in the case of linear or higher
% order de-trending), which should reduce erroneous high-frequency content
% (counterexample: sine wave, one cycle). The implementation scales the "trend"
% (polynomial fit) by tfZ(omega==0).
% --
% NOTE: Retrending is bad for non-lowpass filters (TFs) since the retrending
% requires scaling the fit by tfZ(omega=0) which is only meaningful for lowpass
% filters.
% --
% ** Code has the ability to enable/disable de-trending:
%       -- To handle both DC and AC signals.
%       -- Make testing easier.
% ** Code has the ability to make TF zero above frequency cutoff. This cut-off
%    is naturally sampling frequency-dependent and is therefore NOT a natural
%    part of the TF itself.
%
%
% ARGUMENTS
% =========
% dt
%       Scalar. Time difference between samples.
% y1
%       Column vector of samples. May be modified by this function before
%       actually applying the TF.
% tf
%       Function handle. Transfer function. Same as for bicas.tf.apply_TF_freq().
% varargin
%       Optional settings arguments as interpreted by
%       irf.utils.interpret_settings_args().
%       Available settings:
%         * 'detrendingDegreeOf'
%               Integer.
%               >=0 : Degree of the polynomical fit used for de-trending.
%                <0 : No de-trending.
%               Default = -1.
%         * 'retrendingEnabled'
%         * 'tfHighFreqLimitFraction'
%               Fraction of Nyquist frequency (1/dt). TF is regarded as zero
%               above this frequency. Can be Inf.
%         +more. See implementation.
%
%
% RETURN VALUES
% =============
% y2
%       y1 after applying modifications and TF. Always column vector.
% Debug
%       Struct. For debugging and automatic tests.
%       .y1ModifCa
%           Potentially modified y1 on which the TF is applied.
%       .y2ModifCa
%           Data after applying TF and before it is potentially modified
%           (again).
%       .tfModif
%           The actual (potentially modified) TF used.
%       .i1Array, .i2Array
%           1D column arrays with indices into y1 for the first and last index
%           for the respective time intervals separated by non-finite y1 values.
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-11-04.
%
function [y2, Debug] = apply_TF(dt, y1, tf, varargin)
% PROPOSAL: Move bicas.tf to bicas.proc.L1L2.tf.
%   PRO: Code only used for L1/L1R-->L2.
%   CON: Implies that code is less generic.
%
% PROPOSAL: Check that data is finite. Only call bicas.tf.apply_TF_freq
%           if all data is non-finite.
%   PRO: bicas.tf.apply_TF_freq() can assume (needs to be updated) that
%        always Z<>NaN and thereby detect if TF can not be evaluated via NaN.
%       PRO: Can construct TFs in steps/parts where each part does not have
%            to be evaluated for all omega (return NaN if can not be
%            evaluated).
%           CON: Not necessarily best solution. TFs could give error when
%                not being able to return value.
%
% PROPOSAL: Somehow make function/code reusable for case of applying TF
%           using arbitrary method and algorithm modifications.
%   Ex: Freq-domain application of TF (FFT+multiplication)
%   Ex: Time-domain application of TF (convolution+kernel).
%   Ex: Time & freq.:  De- & re-trending
%   Ex: Freq. only(?): Cutting freq. TF at Nyquist frequency
%   Ex: Time only:     Hann windows on kernel, time-domain edge handling.
%   --
%   PROPOSAL: Integrate all different applications of TF into same function
%             together with all algorithm modifications.
%   PROPOSAL: Refactor functionality into separate reusable code.
%   PROPOSAL: Accept function handle to equivalent of
%             bicas.tf.apply_TF_freq().
%       PRO: Forces application of both de- & re-trending.
%       CON: Less clear code.
%           PRO: Kludgy if having to do logic between retrending & retrending.
%               Ex: Select how to apply TF (frequency domain, time domain).
%               CON: Can still put that logic into one function that is always being called.
%           PRO: Kludgy if using same principle (function accepting function
%                handle) for multiple functionalities consisting of before &
%                after operations.
%               Ex: ?!!
%   PROPOSAL: Replace with pair of functions: before and after call to bicas.tf.apply_TF_freq().
%   PROPOSAL: Replace with class with two non-constructor methods
%             y=c.method(y) before & after.
%
% TODO-DEC: Format for RVs of modified signals split up?
%   PROPOSAL: Cell arrays signals
%       PRO: Minimizes amount of data.
%       PRO: Explicitly enumerates/isolates time intervals.
%       PROPOSAL: Cell array for first, last indices into initial array.
%
%   PROPOSAL: Analogous with actual returned signal: One long array with
%             non-finite values between intervals of data.
%       PRO: Good for plotting.
%
% PROPOSAL: Better abbbreviation for SNF.
%   ~split, values, samples, time series
%   FV
%     Reflects the sample values in the actual datasets, i.e. what the user
%     sees. Relevant if using term for settings.
%     FVs do not represent +-inf.
%   Non-finite
%     Reflects the sample values in the input to this function, what the
%     function actually sees.
%   --
%   SNF  = Split by Non-Finite -- IMPLEMENTED
%   SNFS = Split by Non-Finite Samples
%   SFV = Split by FV
%   NFS = Non-Finite Splitting


DEFAULT_SETTINGS.detrendingDegreeOf         = -1;
DEFAULT_SETTINGS.retrendingEnabled          = false;
DEFAULT_SETTINGS.tfHighFreqLimitFraction    = Inf;
DEFAULT_SETTINGS.method                     = 'FFT';
DEFAULT_SETTINGS.kernelEdgePolicy           = 'mirror';
DEFAULT_SETTINGS.kernelHannWindow           = false;
DEFAULT_SETTINGS.snfEnabled                 = false;
DEFAULT_SETTINGS.snfSubseqMinSamples        = 1;

S = irf.utils.interpret_settings_args(...
  DEFAULT_SETTINGS, varargin);
irf.assert.struct(S, fieldnames(DEFAULT_SETTINGS), {})
clear DEFAULT_SETTINGS

assert(...
  isscalar( S.snfEnabled) && ...
  islogical(S.snfEnabled))
assert(...
  isscalar(S.snfSubseqMinSamples) && ...
  S.snfSubseqMinSamples >= 1)

% ASSERTION: Arguments
assert(iscolumn(y1), 'Argument y1 is not a column vector.')



%=========================================================================
% Create modified version of TF which is set to zero for high frequencies
%=========================================================================
% NOTE: Permit Settings.tfHighFreqLimitFraction to be +Inf.
assert(...
  isnumeric(  S.tfHighFreqLimitFraction) ...
  && isscalar(S.tfHighFreqLimitFraction) ...
  && ~isnan(  S.tfHighFreqLimitFraction) ...
  && (        S.tfHighFreqLimitFraction >= 0))
% Nyquist frequency [rad/s] =
% = 2*pi [rad/sample] * (1/2 * 1/dt [samples/s])
% = pi/dt
nyquistFreqRps     = pi/dt;
tfHighFreqLimitRps = S.tfHighFreqLimitFraction * nyquistFreqRps;
tfModif = @(omegaRps) (tf(omegaRps) .* (omegaRps < tfHighFreqLimitRps));



if S.snfEnabled
  %===================================================================
  % Split up time interval into sub-intervals separated by non-finite
  % samples (fill values)
  %===================================================================
  % SS = SubSequence
  [i1Array, i2Array] = irf.utils.split_by_false(isfinite(y1));
  nSs = numel(i1Array);
else
  i1Array = 1;
  i2Array = numel(y1);
  nSs     = 1;
end

% Pre-allocate and initialize values that will not later be overwritten.
y2 = NaN(size(y1));

Debug = struct();
Debug.y1ModifCa = cell(nSs, 1);   % Pre-allocate
Debug.y2ModifCa = cell(nSs, 1);   % Pre-allocate
Debug.i1Array   = i1Array;
Debug.i2Array   = i2Array;
Debug.tfModif   = tfModif;

for iSs = 1:nSs
  i1 = i1Array(iSs);
  i2 = i2Array(iSs);

  y1ss = y1(i1:i2);

  if numel(y1ss) >= S.snfSubseqMinSamples
    [y2ss, D] = apply_TF_with_DRT(dt, y1ss, tfModif, S);
    Debug.y1ModifCa{iSs} = D.y1Modif;
    Debug.y2ModifCa{iSs} = D.y2Modif;
  else
    y2ss = NaN(size(y1ss));
    Debug.y1ModifCa{iSs} = [];
    Debug.y2ModifCa{iSs} = [];
  end

  y2(i1:i2) = y2ss;
end
end



function [y2, Debug] = apply_TF_with_DRT(dt, y1, tf, Settings)


%#####################
% Optionally DE-trend
%#####################
Drt = bicas.tf.drt(...
  Settings.detrendingDegreeOf, ...
  Settings.retrendingEnabled);
y1Modif = Drt.detrend(y1);



%#########################
% APPLY TRANSFER FUNCTION
%#########################
switch(Settings.method)

  case 'FFT'
    y2Modif = bicas.tf.apply_TF_freq(dt, y1Modif, tf);
    %[y2B, tfOmegaLookups, tfZLookups] = bicas.tf.apply_TF_freq(dt, y1B, tfB);

  case 'kernel'
    % TODO-NI: Kernel length == Signal length
    %          ==> Bad for very long time series? E.g. CWF?
    % NOTE: Length also affects amount of allocated memory (kernel,
    % padding).
    lenKernel = length(y1);
    %lenKernelMax = ceil(10 / dt);
    %lenKernel = min(lenKernel, lenKernelMax);

    % NOTE: The called function applies the Hann window instead of
    % current function since it only applies to kernel method (as
    % opposed to de-trending & re-trending).
    y2Modif = bicas.tf.apply_TF_time(...
      dt, y1Modif, tf, lenKernel, Settings.kernelEdgePolicy, ...
      'hannWindow', Settings.kernelHannWindow);

  otherwise
    error('BICAS:Assertion:IllegalArgument', ...
      'Illegal setting "method" value.')
end



%#####################
% Optionally RE-trend
%#####################
y2 = Drt.retrend(y2Modif, tf(0));



Debug = struct();
Debug.y1Modif = y1Modif;
Debug.y2Modif = y2Modif;
end
