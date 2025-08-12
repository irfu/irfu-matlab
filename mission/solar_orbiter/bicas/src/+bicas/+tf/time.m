%
% Code for applying TF in the time domain. In practice, one main function plus
% helper functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef time



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Generic general-purpose function for applying a TF (linear
    % frequency-dependent transfer function) to a sequence of real-valued (time
    % domain) samples.
    %
    % Wrapper around bicas.tf.kernel.apply_kernel().
    %
    %
    % ARGUMENTS
    % =========
    % dt
    %       Scalar, numeric. Seconds. Time between samples.
    % y1
    %       Column vector. Signal.
    % tf
    %       Function handle. Z = tf(omegaRps). Transfer function.
    % lenKernel
    %       Length of kernel that shall be generated from tf.
    % edgePolicy
    % hannWindowEnabled
    %
    %
    % RETURN VALUES
    % =============
    % y2
    %       Column vector. Transformed signal.
    % yKernelB
    %       Column vector.
    %       Kernel used internally, after the Hann window was applied to it (if
    %       enabled).
    %       NOTE: Only returned for testing purposes.
    %
    %
    % Author: Erik P G Johansson, Uppsala, Sweden
    % First created 2020-08-21.
    %
    function [y2, yKernel] = apply_TF(...
        dt, y1, tf, lenKernel, edgePolicy, hannWindowEnabled)
      %
      % PROPOSAL: MTEST code.
      %   PROPOSAL: Use data from L1R datasets. Compare freq. and time domain
      %             applications.
      %     PROPOSAL: Load directly from CDF.
      %       PROBLEM: Needs HK:DLR, mux mode, ACHG/ACLG
      %   PROPOSAL: Good enough plots for meetings.

      %=============
      % ~ASSERTIONS
      %=============
      % TODO-DEC: Which argument assertions should one bother to have?
      % bicas.tf.freq.apply_TF() and bicas.tf.kernel.apply_kernel() check most
      % arguments.
      assert(isnumeric(lenKernel), 'lenKernel is not numeric.')
      assert(isscalar(lenKernel),  'lenKernel is not scalar.')
      assert(lenKernel>0,          'lenKernel is not positive.')
      assert(islogical(hannWindowEnabled) & isscalar(hannWindowEnabled))



      %===========================
      % Obtain time domain kernel
      %===========================
      [yKernel, iKernelOrigin] = bicas.tf.time.get_impulse_kernel(...
        lenKernel, dt, tf);



      if hannWindowEnabled
        %=============================
        % Apply Hann window to kernel
        %=============================
        shiftedHannWin = bicas.tf.time.get_shifted_Hann_window(...
          iKernelOrigin, lenKernel);
        yKernel        = yKernel .* shiftedHannWin;
      end



      %================
      % Process signal
      %================
      y2 = bicas.tf.kernel.apply_kernel(y1, yKernel, iKernelOrigin, edgePolicy);
    end



    % NOTE: Set "kernel origin" (see bicas.tf.kernel.apply_kernel()) to middle,
    % and rounded down for even-length kernels. The algorithm is designed so
    % that this can be set quite arbitrarily but in reality one probably wants
    % to set it around the middle index.
    %
    % NOTE: Uses bicas.tf.freq.apply_TF(), BICAS' other main function for
    % applying transfer functions to signals, using FFT. Here it is only used
    % for obtaining an impulse response in the time domain, i.e. kernel.
    function [yKernel, iKernelOrigin] = get_impulse_kernel(lenKernel, dt, tf)
      iKernelOrigin           = floor(1 + (lenKernel-1)/2);
      yImpulse                = zeros(lenKernel, 1);
      yImpulse(iKernelOrigin) = 1;
      yKernel                 = bicas.tf.freq.apply_TF(dt, yImpulse, tf);
    end



    % Obtain a Hann window, with the maximum at iMax.
    function shiftedHannWin = get_shifted_Hann_window(iMax, n)
      unshiftedHannWin = hann(n, 'periodic');

      % (Periodic) Hann window center index, i.e. where the Hann window=max=1,
      % for the initial Hann window produced by hann().
      % --
      % NOTE: Index is a half-integer for ODD-numbered-length PERIODIC Hann
      % windows (i.e. ODD-numbered-length kernels). Rounding is therefore
      % ~arbitrary.
      iInitialHannWinCenter = 1 + ceil((n-1)/2);   % Round up.

      % Disable?
      assert(...
        (mod(n, 2) == 1) || ...
        unshiftedHannWin(iInitialHannWinCenter) == 1)

      % Circularly shift Hann window so that the Hann window max is at the
      % kernel origin.
      shiftedHannWin = circshift(...
        unshiftedHannWin, ...
        iMax - iInitialHannWinCenter);
    end



  end    % methods(Static)



end
