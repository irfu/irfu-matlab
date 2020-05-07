function [y2] = apply_transfer_function(dt, y1, tf, varargin)
%
% Generic general-purpose function for applying a spectrum TF(transfer function) to a sequence of (real-valued, time
% domain) samples.
%
%
% ALGORITHM
% =========
% (1) De-trend (if enabled)
% (2) Compute DFT using MATLAB's "fft" function.
% (3) Interpret DFT component frequencies as pairs of positive and negative frequencies (lower and higher half
%     of DFT components. (Interpret TF as symmetric function, Z(omega) = Z*(-omega), *=conjugate, covering positive &
%     negative frequencies.)
% (4) Multiply DFT coefficients with complex TF values.
% (5) Compute inverse DFT using MATLAB's "ifft(... , 'symmetric')" function.
% (6) Re-trend (if de-trending enabled)
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% NOTE: All arguments/return value vectors are column vectors. TF = transfer function.
% dt       : Time between each sample. Unit: seconds
% y1       : Samples. Must be real-valued (assertion). May contain NaN.
% tf       : Function handle to function z=tf(omega). z is a complex value (amplitude+phase) and has not unit.
%            omega unit: rad/s.
%            Will only be called for omega>=0. tf(0) must be real.
%            NOTE: If the caller wants to use a tabulated TF, then s/he should construct an anonymous function using
%            "interp1" and submit it as argument.
% y2       : y1 after the application of the TF.
%            If y1 contains at least one NaN, then all components in y2 will be NaN. No error will be thrown.
% varargin : Optional settings arguments as interpreted by EJ_library.utils.interpret_settings_args.
%   Possible settings:
%       enableDetrending : Override the default on whether de-trending is used. Default=0.
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
% NOTE/BUG: Detrending makes it impossible to modify the amplitude & phase for the frequency components in the trend,
% e.g. to delay the signal. If the input signal is interpreted as N-periodic, then de-trending affects the jump between
% the beginning and end of the signal (reduces it in the case of linear de-trending), which affects the high-frequency
% content(?) but probably in a good way. The implementation scales the "trend" (polynomial fit) by tfZ(omega==0).
% --
% NOTE: Presently not sure if MATLAB has standard functions for applying a transfer function in the
% frequency domain and that is tabulated or function handle.   /Erik P G Johansson 2019-09-11
%
%
% IMPLEMENTATION NOTES
% ====================
% -- Has the ability to enable/disable de-trending to make testing easier.
% -- Conversion of transfer functions to fit the input format should be done by wrapper functions and NOT by this
%    function.
% -- This function only represents the pure mathematical algorithm and therefore only works with "mathematically pure"
%    variables and units: radians, complex amplitudes (no dB, no volt^2, no amplitude+phase). This is useful since it
% (1) separates (a) the core processing code from (b) related but simple processing of data (changing units,
% different ways of representing transfer functions, checking for constant sampling rate)
% (2) makes the potentially tricky TF-code easier to understand and check (due to (1)),
% (3) makes a better code unit for code testing,
% (4) makes it easier to simultaneously support different forms of input data (in wrapper functions),
% (5) it is easy to combine multiple TF:s on the TF format that this function accepts,
% (6) easier to use it for mathematically calculated transfer functions, e.g. due to RPW's parasitic capacitance
% (although that should not be done in isolation, but rather by combining it with other TF:s.
%
%
% TERMINOLOGY
% ===========
% DFT = Discrete Fourier Transform
% TF  = Transfer function, ("spectrum") transfer functions, i.e. transfer functions which modify the spectrum content of
% a signal, represented in the pure mathematical form as Z=Z(omega), where Z is a complex number (practically,
% multiply frequency component of the signal in volt; not volt^2) and omega is a frequency (radians/s).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-13



% ------------------------------------------------------------------------------
% TODO-NEED-INFO: WHY DOES THIS FUNCTION EXIST? DOES NOT MATLAB HAVE THIS FUNCTIONALITY?
% PROPOSAL: Option for using inverse TF? Can easily be implemented in the actual call to the function though
%           (dangerous?).
% PROPOSAL: Option for error on NaN/Inf.
% PROPOSAL: Eliminate dt from function. Only needed for interpreting tfOmega. Add in wrapper.
% PROPOSAL: Eliminate de-trending. Add in wrapper.
%   CON/NOTE: Might not be compatible with future functionality (Hann Windows etc).
%
% PROPOSAL: If slow to call function handle for transfer function tf, permit caller to submit table with implicit frequencies.
%   PROPOSAL: Return the Z values actually used, so that caller can call back using them.
%   PROPOSAL: Separate function for generating such vector.
%
% TODO-NEED-INFO: How does algorithm handle X_(N/2+1) (which has no frequency twin)? Seems like implemention should
%   multiply it by a complex Z (generic situation) ==> Complex y2. Still, no such example has been found yet.
%   Should be multiplied by abs(Z)?! Z-imag(z)?! Keep as is?!
%
% PROPOSAL: Not require column vectors. Only require vectors.
% ------------------------------------------------------------------------------



% Set the order of the polynomial that should be used for detrending.
N_POLYNOMIAL_COEFFS_TREND_FIT = 1;    % 1 = Linear function.

%============
% ASSERTIONS
%============
if ~iscolumn(y1)
    error('BICAS:apply_transfer_function:Assertion:IllegalArgument', 'Argument y1 is not a column vector.')
elseif ~isnumeric(y1)
    error('BICAS:apply_transfer_function:Assertion:IllegalArgument', 'Argument y1 is not numeric.')
elseif ~isreal(y1)
    error('BICAS:apply_transfer_function:Assertion:IllegalArgument', 'y1 is not real.')
    % NOTE: The algorithm itself does not make sense for non-real functions.
elseif ~isscalar(dt)
    error('BICAS:apply_transfer_function:Assertion:IllegalArgument', 'dt is not scalar.')
elseif ~(dt>0)
    error('BICAS:apply_transfer_function:Assertion:IllegalArgument', 'dt is not positive.')
elseif ~isa(tf, 'function_handle')
    % EJ_library.assert.func does not seem to handle return values correctly.
    error('BICAS:apply_transfer_function:Assertion:IllegalArgument', 'tf is not a function.')
elseif ~isreal(tf(0))
    error('BICAS:apply_transfer_function:Assertion:IllegalArgument', 'tf(0) is not real.')    
end



DEFAULT_SETTINGS.enableDetrending = 0;
Settings = EJ_library.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
EJ_library.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})



N = length(y1);

if Settings.enableDetrending
    %##########
    % De-trend
    %##########
    trendFitsCoeffs1 = polyfit((1:N)', y1, N_POLYNOMIAL_COEFFS_TREND_FIT);
    yTrend1          = polyval(trendFitsCoeffs1, (1:N)');
    y1               = y1 - yTrend1;
end



%#############
% Compute DFT
%#############
yDft1 = fft(y1);



%================================================================================================================
% Define the frequencies used to interpret the DFT components X_k (yDft1)
% -----------------------------------------------------------------------
% IMPLEMENTATION NOTE:
% The code only works with REAL-valued time-domain signals. Therefore,
% (1) We want to interpret the signal as consisting of pairs of positive and negative frequencies (pairs of
% complex bases).
% (2) We want to interpret the TF as being a symmetric function, defined for both positive and negative
% frequencies, Z(omega)=Z*(-omega), *=conjugate.
%
% The DFT components X_k, k=1..N can be thought of as representing different frequencies
%    omega_k = 2*pi*(k-1) / (N*dt) .
% Since
%    exp(i*2*pi*omega_k*t_n) = exp(i*2*pi*omega_(k+N)*t_n),
%    where t_n = (n-1)*dt ,
% the exact frequencies associated with DFT components X_k are however subject to a choice/interpretation, where
%    omega_k <--> omega_(k+N) .
% Since we only work with real-valued signals, we want to interpret the DFT components as having frequencies
%    omega_1, ..., omega_ceil(N/2), omega_[ceil(N/2)+1-N], ..., omega_0
% but to look up values in the TF, we have to use the absolute values of the above frequencies and conjugate Z when
% necessary.
%
% NOTE: omega_0 = 0.
% NOTE: The above must work for both even & odd N. For even N, the DFT component X_N/2+1 (which does not have a
% frequency twin) should be real for real signals.
%================================================================================================================
%tfOmegaLookups     = 2*pi * ((1:N) - 1) / (N*dt);
%i = (tfOmegaLookups >= pi/dt);    % Indicies for which omega_k should be replaced by omega_(k-N).
%tfOmegaLookups(i) = abs(tfOmegaLookups(i)  - 2*pi/dt);
kOmegaLookup   = [1:ceil(N/2), (ceil(N/2)+1-N):0 ];   % Modified k values used to calculate omega_k for every X_k.
tfOmegaLookups = 2*pi * (kOmegaLookup - 1) / double(N*dt);



%=============================================================================================================
% Find complex TF values, i.e. complex factors to multiply every DFT component with
% ---------------------------------------------------------------------------------
% NOTE: De-trending (if enabled) should already have removed the zero-frequency component from the in signal.
%=============================================================================================================
tfZLookups    = tf(abs(tfOmegaLookups));
i             = tfOmegaLookups < 0;
tfZLookups(i) = conj(tfZLookups(i));
% % ASSERTION:
if ~all(isfinite(tfZLookups) | isnan(tfZLookups))
    error('BICAS:apply_transfer_function:Assertion', 'Transfer function tf returned non-finite value for at least one frequency.')
end



%##################
% Apply TF to data
%##################
% NOTE: For real input signal and even N, this should produce complex yDft2(N/2+1) values.
% IMPORTANT NOTE: Must transpose complex vector in a way that does not negate the imaginary part.
%                 Transposing with ' (apostrophe) negates the imaginary part.
yDft2 = yDft1 .* transpose(tfZLookups);



%##############
% Compute IDFT
%##############
% IMPLEMENTATION NOTE: Uses ifft options to force yDft2 to be (interpreted as) conjugate symmetric due to possible
% rounding errors.
%
% ifft options:
%     "ifft(..., 'symmetric') causes ifft to treat X as conjugate symmetric
%     along the active dimension.  This option is useful when X is not exactly
%     conjugate symmetric merely because of round-off error.  See the
%     reference page for the specific mathematical definition of this
%     symmetry."
y2 = ifft(yDft2, 'symmetric');
y2p = ifft(yDft2);    % TEST



%##########
% Re-trend
%##########
if Settings.enableDetrending
    % Use Z(omega=0) to scale trend, including higher order polynomial components.
    trendFitsCoeffs2 = trendFitsCoeffs1 * tfZLookups(1);
    
    yTrend2 = polyval(trendFitsCoeffs2, (1:N)');
    y2      = y2 + yTrend2;
end



% ASSERTION: Real output.
% IMPLEMENTATION NOTE: Will react sometimes if "ifft" with 'symmetric' is not used.
if ~isreal(y2)
    maxAbsImag = max(abs(imag(y2)))
    error('BICAS:apply_transfer_function:Assertion', 'y2 is not real. Bug.')
end


end
