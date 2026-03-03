%
% Code for applying a TF in the frequency domain. In practice, one main
% function.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef freq



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Generic general-purpose function for applying a TF (linear
    % frequency-dependent transfer function) to a sequence of (real-valued,
    % time domain) samples. The operation takes place in the frequency domain.
    %
    %
    % ALGORITHM
    % =========
    % (1)  Compute DFT of samples using MATLAB's "fft" function.
    % (2a) Interpret the sample DFT component frequencies as pairs of positive
    %      and negative frequencies (lower and higher half of DFT components).
    % (2b) Interpret the specified TF as a ~symmetric function,
    %      Z(omega) = Z*(-omega),
    %      which covers both positive & negative frequencies, and where
    %      *=complex conjugate.
    % (3)  Multiply the sample DFT coefficients with the corresponding complex
    %      TF values.
    % (4)  Compute inverse DFT using MATLAB's "ifft(... , 'symmetric')"
    %      function.
    %
    %
    % NOTES
    % =====
    % NOTE: This function effectively implements an approximate convolution.
    % For an inverse application of a TF (de-convolution), the caller has to
    % invert the TF first.
    % --
    % NOTE: irfu-matlab contains at least two other functions for applying
    % transfer functions to data. These two other functions are however not
    % general-purpose (can not easily be reused):
    % 1) c_efw_invert_tf.m      (extensive; in both time domain and frequency
    %                            domain; multiple ways of handling edges)
    % 2) c_efw_burst_bsc_tf.m   (short & simple)
    % Above functions do not *APPEAR TO* make use of any particular standard
    % MATLAB functions for applying transfer functions to data more than
    % general functions like fft(), ifft(), conv(), cconv(). (One needs to
    % check the code more to be sure though.)
    % --
    % NOTE: I am presently not sure if MATLAB has standard functions for
    % applying a transfer function in the frequency domain and that is
    % tabulated or function handle.   /Erik P G Johansson 2019-09-11
    % --
    % NOTE: This function does not support any kind of special edge handling.
    %
    %
    % EXPONENT SIGN CONVENTION IN TRANSFER FUNCTIONS
    % ==============================================
    % This function/algorithm uses the following:
    %   y1(t)     ~ e^(i*omega*t)                # Exponent sign convention
    %                                            # used by MATLAB's fft()
    %                                            # & ifft().
    %   tf(omega) ~ e^(i*omega*(-tau(omega)))    # Transfer function supplied
    %                                            # to this function.
    %   y2(t)     ~ e^(i*omega*t) * tf(omega)
    %             = e^(i*omega*(t-tau(omega)))
    % (Weighted summing/integration over exponentials is implicit.) Therefore,
    % a TF component with a positive tau represents a positive phase delay of
    % tau for that frequency, i.e.
    %   y2(t) == y1(t-tau)
    % if e.g. y1(t) only has one frequency component.
    % NOTE: This should be the same convention as used by the Laplace
    % transform.
    %
    %
    % IMPLEMENTATION NOTES, DESIGN INTENT
    % ===================================
    % -- Modifications of the TRANSFER FUNCTION should be made by wrapper
    %    functions and NOT by this function.
    %       Ex: Modifications to fit the input format.
    %           Ex: Turn a given tabulated TF into an actual MATLAB function
    %               (handle).
    %       Ex: Remove high frequency components for inverted lowpass filter.
    %       Ex: Remove/dampen low frequencies for inverted highpass filter.
    % -- Modification of input/output SAMPLES should be done by wrapper
    %    functions and NOT in this function.
    %       Ex: De-trending, re-trending
    %       Ex: Splitting by non-finite values (e.g. data gaps represented by
    %           NaN).
    % -- This function only represents the pure mathematical algorithm and
    %    therefore only works with "mathematically pure" variables and units:
    %    radians, complex amplitudes (no dB, no volt^2, no amplitude+phase).
    %    This is useful since it
    % (1) separates
    %       (a) the core processing code from
    %       (b) related but simple processing of data (changing units,
    %           different ways of representing transfer functions, checking for
    %           constant sampling rate),
    % (2) makes the potentially tricky TF-code easier to understand and check
    %     (due to (1)),
    % (3) makes a better code unit for code testing,
    % (4) makes it easier to simultaneously support different forms of input
    %     data (in wrapper functions),
    % (5) it is easy to combine multiple TFs on the TF format that this
    %     function accepts,
    % (6) easier to use it for mathematically calculated transfer functions,
    %     e.g. due to RPW's parasitic capacitance (although that should not be
    %     done in isolation, but rather by combining it with other TFs.
    %
    %
    % ARGUMENTS
    % =========
    % dt
    %       Time between each sample. Unit: seconds
    % y1
    %       Nx1. Samples. Must be real-valued (assertion).
    %       NOTE: May contain non-finite values.
    % tf
    %       Function handle to function z=tf(omega). z is a complex value
    %       (amplitude+phase) and has no intrinsic unit. omega unit: rad/s.
    %       Will only be called for omega>=0. tf(0) must be real.
    %       NOTE: Permitted to return NaN, but not infinity.
    %       NOTE: If the caller wants to use a tabulated TF, then s/he should
    %       construct an anonymous function that interpolates the tabulated TF
    %       (e.g. using "interp1") and submit it as argument.
    %       NOTE: The transfer function is permitted to return NaN. This will
    %             set y2 to NaN.
    %
    %
    % RETURN VALUES
    % =============
    % y2
    %       y1 after the application of the TF.
    %       If y1 contains at least one NaN, then all components in y2 will be
    %       NaN. No error will be thrown if that is the case.
    %
    %
    % Author: Erik P G Johansson, IRF, Uppsala, Sweden
    % First created 2017-02-13
    %
    function [y2] = apply_TF(dt, y1, tf)
      % TODO-NI: WHY DOES THIS FUNCTION NEED TO EXIST? DOES NOT MATLAB HAVE
      %          THIS FUNCTIONALITY?
      %
      % PROPOSAL: Assert not inf in signal.
      %   NOTE: NaN is deliberately permitted in signal.
      %   NOTE: NaN but not infinity is deliberately permitted in TF.
      %
      % TODO-NI: How does algorithm handle X_(N/2+1) (which has no frequency
      %          twin)? Seems like the implemention should multiply it by a
      %          complex Z (generic situation) ==> Complex y2. Still, no such
      %          example has been found yet.
      %   Should be multiplied by abs(Z)?! Z-imag(z)?! Keep as is?!
      %
      % PROPOSAL: Permit submitting multiple y1 at the same time (same length,
      %           same dt, same tf).
      %   PRO: Faster?
      %       PRO: Can call tf once.
      %   CON: Not useful for CWF data which tends to have unique lengths.
      %   CON: Snapshots tend to have the same length but rotates between
      %        sampling rates, i.e. dt values.



      % EMID = Error Message ID
      EMID_ARG = 'BICAS:Assertion:IllegalArgument';

      %============
      % ASSERTIONS
      %============
      % NOTE: The number of arguments has changed historically.
      assert(nargin == 3)
      if ~iscolumn(y1)
        error(EMID_ARG, 'Argument y1 is not a column vector.')
      elseif ~isnumeric(y1)
        error(EMID_ARG, 'Argument y1 is not numeric.')
      elseif ~isreal(y1)
        error(EMID_ARG, 'y1 is not real.')
        % NOTE: The algorithm itself does not make sense for non-real
        %       functions.
      elseif ~isnumeric(dt)
        error(EMID_ARG, 'dt is not numeric..')
      elseif ~isscalar(dt)
        error(EMID_ARG, 'dt is not scalar.')
      elseif ~(dt > 0)
        error(EMID_ARG, 'dt is not positive.')
      elseif ~isa(tf, 'function_handle')
        % irf.assert.func does not seem to handle return values
        % correctly.
        error(EMID_ARG, 'tf is not a function.')
      elseif ~isreal(tf(0))
        error(EMID_ARG, 'tf(0) is not real.')
      end



      nSamples = length(y1);



      %##################################################
      % Compute TF Z(omega) values to later apply to DFT
      %##################################################
      %========================================================================
      % Define the frequencies used to interpret the DFT components X_k
      % (yDft1)
      % ---------------------------------------------------------------
      % IMPLEMENTATION NOTE:
      % The code only works with REAL-valued time-domain signals. Therefore,
      % (1) We want to interpret the signal as consisting of pairs of positive
      %     and negative frequencies (pairs of complex bases), and therefore
      %     with complex-conjugated weights.
      % (2) We want to interpret the TF as being a symmetric function, defined
      %     for both positive and negative frequencies,
      %       Z(omega) = Z*(-omega),
      %     where *=conjugate.
      %
      % Excerpt from MATLAB's help text for fft():
      % """"
      %   For length N input vector x, the DFT is a length N vector X,
      %   with elements
      %                    N
      %      X(k) =       sum  x(n)*exp(-j*2*pi*(k-1)*(n-1)/N), 1 <= k <= N.
      %                   n=1
      %   The inverse DFT (computed by IFFT) is given by
      %                    N
      %      x(n) = (1/N) sum  X(k)*exp( j*2*pi*(k-1)*(n-1)/N), 1 <= n <= N.
      %                   k=1
      % """"
      %
      % The DFT components X_k, k=1..N can be thought of as representing
      % different frequencies. The exponents above can be interpreted as
      %   j*2*pi*(k-1)*(n-1)/N = j*omega_k*t
      % where
      %   omega_k = 2*pi*(k-1) / (N*dt)
      %   t_n     = (n-1)*dt
      % .
      % Since
      %    e^(i * omega_k * t_n) = e^(i * omega_(k+m*N) * t_n),
      % where
      %    m = any integer,
      % the exact frequencies associated with DFT components X_k are however
      % subject to a choice/interpretation, where
      %    omega_k <--> omega_(k+m*N) .
      % Since we only work with real-valued signals, we want to interpret the
      % DFT components as having frequencies
      %    omega_1, ..., omega_ceil(N/2), omega_[ceil(N/2)-N+1], ..., omega_0
      % but to look up values in the TF, we have to use the absolute values of
      % the above frequencies and conjugate Z when necessary.
      %
      % NOTE: omega_0 = 0.
      % NOTE: The above must work for both even & odd N. For even N, the DFT
      %       component X_N/2 (which does not have a frequency twin) should be
      %       real for real signals.
      %========================================================================
      % Modified k values (~indices) used to calculate omega_k for every X_k.
      kOmegaSamplesAr = [1:ceil(nSamples/2), (ceil(nSamples/2)-nSamples+1):0];
      % At least roughly equivalent. Not tested.
      % kOmegaSamplesAr2 = [1:nSamples] - [(1:nSamples)>(nSamples/2+1)] * nSamples;

      omegaSamplesAr  = 2*pi * (kOmegaSamplesAr - 1) / double(nSamples*dt);



      %======================================================================
      % Find complex TF values, i.e. complex factors to multiply every DFT
      % component with
      % ------------------------------------------------------------------
      % NOTE: De-trending (outside function) should already have removed the
      %       zero-frequency component from the in signal.
      %======================================================================
      tfZAr                  = tf(abs(omegaSamplesAr));
      bNegativeFreqAr        = omegaSamplesAr < 0;
      tfZAr(bNegativeFreqAr) = conj(tfZAr(bNegativeFreqAr));   % Modify some indices.
      % ASSERTION
      %if ~all(isfinite(tfZLookups) | isnan(tfZLookups))
      if ~all(~isinf(tfZAr))
        % NOTE: Deliberately permits Z=NaN (but not infinity) since
        % bicas.proc.L1L2.cal.VoltageCalibration is designed to create TFs that
        % return Z=NaN for impossible combinations where it does not matter
        % anyway.
        % /EJ 2020-11-05
        error(...
          'BICAS:Assertion', ...
          ['Transfer function "tf" returned non-finite value (not NaN)', ...
          ' for at least one frequency.'])
      end



      %#############
      % Compute DFT
      %#############
      yDft1 = fft(y1);



      %##################
      % Apply TF to data
      %##################
      % NOTE: For real input signal and even N, this should produce complex
      % yDft2(N/2+1) values.
      % IMPORTANT NOTE: Must transpose complex vector in a way that does not
      % negate the imaginary part. Transposing with ' (apostrophe) negates the
      % imaginary part. transpose() does not negate the imaginary part.
      yDft2 = yDft1 .* transpose(tfZAr);



      %#####################
      % Compute inverse DFT
      %#####################
      % IMPLEMENTATION NOTE: Uses ifft() options to force yDft2 to be
      % (interpreted as) conjugate symmetric due to possible rounding errors.
      %
      % ifft options:
      %     "ifft(..., 'symmetric') causes ifft to treat X as conjugate
      %     symmetric along the active dimension.  This option is useful when X
      %     is not exactly conjugate symmetric merely because of round-off
      %     error.  See the reference page for the specific mathematical
      %     definition of this symmetry."
      y2 = ifft(yDft2, 'symmetric');



      % ASSERTION: Real (numbers) output.
      % IMPLEMENTATION NOTE: Will react sometimes if "ifft" with 'symmetric' is
      % not used.
      if ~isreal(y2)
        maxAbsImag = max(abs(imag(y2)));
        error('BICAS:Assertion', ...
          'y2 is not real (non-complex). Bug. maxAbsImag=%g.', maxAbsImag)
      end

    end



  end    % methods(Static)



end
