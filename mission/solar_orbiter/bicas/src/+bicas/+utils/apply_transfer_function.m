function y2 = apply_transfer_function(dt, y1, tfOmega, tfZ, varargin)
% y2 = apply_transfer_function_in_freq(dt, y1, tfOmega, tfZ, varargin)
% Generic general-purpose function for applying a spectrum TF to a sequence of samples
% (real-valued, time domain).
%
% ALGORITHM
% =========
% (1) de-trend (if enabled)
% (2) DFT
% (3) Interpret DFT component frequencies as pairs of positive and negative frequencies (lower and higher half
%     of DFT components. (Interpret TF as symmetric function covering positive & negative frequencies. Zero at
%     zero frequency.)
% (4) Multiply DFT coefficients with complex TF values.
% (5) Inverse DFT
% (6) Re-trend (if de-trending enabled)
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-13
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% NOTE: All arguments/return value vectors are column vectors. TF = transfer function.
% dt       : Time between each sample. Unit: seconds
% y1       : Samples. Must be real-valued.
% tfOmega  : TF frequencies: Unit: radians/s (RPS = Radians per second)
% tfZ      : Complex TF values (multiplication factors) at frequencies tfOmega. No unit.
% y2       : y1 after the application of the TF.
%            If y1 contains at least one NaN, then all components in y2 will be NaN. No error will be thrown.
% varargin : Sequence of optional options as parameter/value pairs.
%   'EnableDetrending', enableDetrending : Override the default on whether de-trending is used.
%
%
% NOTES
% =====
% NOTE: This function effectively implements an approximate convolution. For an inverse application of a TF
% (de-convolution), the caller has to invert the TF first.
% NOTE: irfu-matlab contains at least two other functions for applying transfer functions to data but which are
% not general-purpose:
% 1) c_efw_invert_tf.m      (extensive; in both time domain and frequency domain; multiple ways of handling edges)
% 2) c_efw_burst_bsc_tf.m   (short & simple)
%
%
% IMPLEMENTATION NOTES
% ====================
% -- The only reason for that this function is public is to make it possible for external test code to access it.
% -- Added ability to enable/disable de-trending to make testing easier.
% -- Conversion of transfer functions to fit the input format should be done by wrapper functions and NOT by
%    this function.
% -- This function is only represents the pure mathematical algorithm and therefore only
% works with "mathematically pure" units: radians, amplitudes (no dB!). This is useful since it
% (1) separates (a) the core processing code from (b) related but simple processing of data (changing units,
% different ways of representing transfer functions, checking for constant sampling rate)
% (2) makes the potentially tricky TF-code easier to understand and check (due to (1)),
% (3) makes a better unit for testing,
% (4) makes it easier to simultaneously support different forms of input data (in wrapper functions),
% (5) it is easy to combine multiple TF:s on the TF format that this function accepts,
% (6) easier to use it for mathematically calculated transfer functions, e.g. due to RPW's parasitic capacitance
% (although that should not be done in isolation, but rather by combining it with other TF:s.

% ------------------------------------------------------------------------------
% PROPOSAL: Process multiple sequences functions in one call.
%   TODO-DECISION: Handle sequences of different length?
%       Ex: Different-length snapshots (LFR)
%       Ex: Longer CWF (SWF?!) sequences divided by missing data into different-length sequences.
% PROPOSAL: Option for using inverse TF? Can easily be implemented in the actual call to the function though
%           (dangerous?).
% PROPOSAL: Always use inverse TF since that is what is going to be used anyway?!
% PROPOSAL: Option for error on NaN/Inf.
% PROPOSAL: Assert tfOmega sorted?
% PROPOSAL: Read TF from a function (pointer) instead.
%   PRO: Caller can control the inter-/extrapolation of table (linear, spline etc).
%   PRO: Useful when combining table TFs with other table TFs (other set of frequencies), or with analytical TFs.
%   CON: Slower?
% ------------------------------------------------------------------------------



% Set the type of polynomial that should be used for detrending.
N_POLYNOMIAL_COEFFS_TREND_FIT = 1;    % 1 = Linear function.

% ASSERTIONS
if ~iscolumn(y1) || ~iscolumn(tfOmega) || ~iscolumn(tfZ)
    error('BICAS:calibration:Assertion', 'Argument y1, tfOmega, or tfZ is not a column vector.')
elseif ~isscalar(dt)
    error('BICAS:calibration:Assertion', 'dt is not scalar.')
elseif (numel(tfOmega) ~= numel(tfZ))
    error('BICAS:calibration:Assertion', 'tfOmega, tfZ have different sizes.')
elseif ~isreal(y1)
    error('BICAS:calibration:Assertion', 'y1 is not real.')
    % NOTE: The algorithm itself does not make sense for non-real functions.
elseif any(tfOmega<=0)
    error('BICAS:calibration:Assertion', 'tfOmega contains non-positive frequencies.')
    % NOTE: Requires TF to be undefined for 0 Hz since such a value should not be used since it represents
    % adding a constant offset to the in signal (sort of).
end



% Extract values of options.
% IMPLEMENTATION NOTE: Written to make it possible to have many options. I just forgot which other ones I
% wanted to add...
enableDetrending = 1;   % Default value
while length(varargin) >= 1
    if length(varargin) >= 2
        if strcmp(varargin{1}, 'EnableDetrending')
            enableDetrending = varargin{2};
        else
            error('BICAS:calibration:Assertion', 'Can not interpret extra option.')
        end
        varargin(1:2) = [];
    else
        error('BICAS:calibration:Assertion', 'Can not interpret options.')
    end
end



N = length(y1);

if enableDetrending
    % De-trend
    trendFitsCoeffs = polyfit((1:N)', y1, N_POLYNOMIAL_COEFFS_TREND_FIT);
    yTrendFit       = polyval(trendFitsCoeffs, (1:N)');
    y1              = y1 - yTrendFit;
end

% Derive DFT
yDft1 = fft(y1);

%===========================================================================================================
% Define the frequencies we use to interpret the DFT
% --------------------------------------------------
% IMPLEMENTATION NOTE:
% We work only with REAL-valued signals. Therefore,
% (1) We want to interpret the signal as consisting of pairs of positive and negative frequencies (pairs of
% complex bases).
% (2) We want to interpret the TF as being a symmetric function, defined for both positive and negative
% frequencies, Z(-omega)=Z(omega).
%
% The DFT components k=1..N can be thought of as representing different frequencies
%    omega_k = 2*pi*(k-1) / (N*dt) .
% Since
%    exp(i*2*pi*omega_k*t_n) = exp(i*2*pi*omega_(k+N)*t_n),
%    where t_n = (n-1)*dt ,
% the exact frequencies are however subject to a choice/interpretation, where
%    omega_k <--> omega_(k+N) .
% Since we only work with real-valued signals, we want to interpret the DFT components as having frequencies
%    omega_1, ..., omega_ceil(N/2), omega_ceil(-N/2+1), ..., omega_0
% but to look up values in the TF, we have to use the absolute values of the above frequencies.
%
% NOTE: The above must work for both even & odd N. For even N, the DFT component k=N/2 should be zero due to
% being a real-valued signal. Whether one interprets k=N/2 as a positive or negative frequency (using
% floor/ceil) should therefore be unimportant.
%===========================================================================================================
k1 = ( 1            : ceil(N/2) )';
k2 = ( ceil(-N/2+1) : 0         )';
omegaDft1 = 2*pi * (k1 - 1) / (N * dt);          % Non-negative frequencies.
omegaDft2 = 2*pi * (k2 - 1) / (N * dt);          % Negative     frequencies.
omegaDftLookupHighest = omegaDft1(end);

% Frequencies that should be USED for finding TF values to use for the respective DFT components.
% Assumes that the TF can be interpreted as a symmetric functions.
tfOmegaLookups = [omegaDft1; - omegaDft2];

% ASSERTION
if omegaDftLookupHighest > max(tfOmega)
    error('BICAS:calibration:Assertion', 'Transfer function does not cover the highest frequency %d [rad/s] in the data samples.', omegaDftLookupHighest)
end



% Find complex TF values, i.e. complex factors to be multiplied with every DFT component
% --------------------------------------------------------------------------------------
% NOTE: interp1 does NOT seem to require that submitted table of (x,y) values is sorted in x.
% NOTE: interp1 requires that submitted table (x,y) has unique x values.
% NOTE: "extrap" ==> Extrapolate TF to frequencies below/above the lowest/highest stated frequencies. ==>
% Non-zero TF value for 0 Hz. ==> Applying TF for 0 Hz-component means adding a constant offset (sort of)
% which one does not want. ==> Force the TF 0 Hz values to be one.
% Note that de-trending (if enabled) should already have removed the zero-frequency component from the in
% signal.
tfZLookups = interp1(tfOmega, tfZ, tfOmegaLookups, 'linear', 'extrap' );
tfZLookups(tfOmegaLookups == 0) = 1;    % Should only be needed when de-trending is disabled.

% Apply TF to data.
% NOTE: For real input signal and even N, this should produce complex yDft2(N/2+1) values.
yDft2 = yDft1 .* tfZLookups;

% IDFT
% ----
% Forces yDft1 to be (interpreted as) conjugate symmetric due to (1) possible rounding errors, and (2) a
% complex yDft2(N/2+1) for even N.
%
% ifft options:
%     "ifft(..., 'symmetric') causes ifft to treat X as conjugate symmetric
%     along the active dimension.  This option is useful when X is not exactly
%     conjugate symmetric merely because of round-off error.  See the
%     reference page for the specific mathematical definition of this
%     symmetry."
y2 = ifft(yDft2, 'symmetric');

% Re-trend
if enableDetrending
    y2 = y2 + yTrendFit;
end

end
