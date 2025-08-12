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
    % varargin
    %       Optional settings arguments as interpreted by
    %       irf.utils.interpret_settings_args().
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
    function [y2, yKernelB] = apply_TF(dt, y1, tf, lenKernel, edgePolicy, varargin)
      %
      % PROPOSAL: Function for generating kernels.
      % PROPOSAL: Separate function for Hann Window modification of kernel.
      %
      % PROPOSAL: edgePolicy as setting, not named argument.
      %
      % PROPOSAL: MTEST code.
      %   PROPOSAL: Use data from L1R datasets. Compare freq. and time domain
      %             applications.
      %     PROPOSAL: Load directly from CDF.
      %       PROBLEM: Needs HK:DLR, mux mode, ACHG/ACLG
      %   PROPOSAL: Good enough plots for meetings.
      %
      % BUG?! ~Applies detrending+retrending, but so does calling function
      %       bicas.tf.apply_TF().

      EMID = 'BICAS:Assertion:IllegalArgument';

      DEFAULT_SETTINGS.detrendingDegreeOf = -1;
      DEFAULT_SETTINGS.retrendingEnabled  = false;
      DEFAULT_SETTINGS.hannWindow         = false;

      Settings = irf.utils.interpret_settings_args(...
        DEFAULT_SETTINGS, varargin);
      irf.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})
      clear DEFAULT_SETTINGS



      %=============
      % ~ASSERTIONS
      %=============
      % TODO-DEC: Which argument assertions should one bother to have?
      % bicas.tf.freq.apply_TF() and bicas.tf.kernel.apply_kernel() check most
      % arguments.
      if ~isnumeric(lenKernel)
        error(EMID, 'lenKernel is not numeric.')
      elseif ~isscalar(lenKernel)
        error(EMID, 'lenKernel is not scalar.')
      elseif ~(lenKernel>0)
        error(EMID, 'lenKernel is not positive.')
      end
      assert(islogical(Settings.hannWindow))



      %===========================
      % Obtain time domain kernel
      %===========================
      % NOTE: Set "kernel origin" (see bicas.tf.kernel.apply_kernel()) to middle, and
      % rounded down for even-length kernels. The algorithm is designed so that this
      % can be set quite arbitrarily but in reality one probably wants to set it
      % around the middle index.
      iKernelOrigin           = floor(1 + (lenKernel-1)/2);
      yImpulse                = zeros(lenKernel, 1);
      yImpulse(iKernelOrigin) = 1;
      yKernel                 = bicas.tf.freq.apply_TF(dt, yImpulse, tf);
      % NOTE: Uses bicas.tf.freq.apply_TF(), BICAS' other main function for applying
      % transfer functions to signals, using FFT. Here it is only used for obtaining
      % an impulse response in the time domain, i.e. kernel.



      if Settings.hannWindow
        %=============================
        % Apply Hann window to kernel
        %=============================
        unshiftedHannWin = hann(lenKernel, 'periodic');

        % (Periodic) Hann window center index, i.e. where the Hann window=max=1,
        % for the initial Hann window produced by hann().
        % --
        % NOTE: Index is a half-integer for ODD-numbered-length PERIODIC Hann
        % windows (i.e. ODD-numbered-length kernels). Rounding is therefore
        % ~arbitrary.
        iInitialHannWinCenter = 1 + ceil((lenKernel-1)/2);   % Round up.

        % Disable?
        assert(...
          (mod(lenKernel, 2) == 1) || ...
          unshiftedHannWin(iInitialHannWinCenter) == 1)

        % Circularly shift Hann Window so that the Hann window max is at the kernel
        % origin.
        shiftedHannWin = circshift(...
          unshiftedHannWin, ...
          iKernelOrigin - iInitialHannWinCenter);

        yKernelB = yKernel .* shiftedHannWin;
      else
        % Do nothing.
        yKernelB = yKernel;
      end



      %================
      % Process signal
      %================
      Drt = bicas.tf.Deretrending(Settings.detrendingDegreeOf, Settings.retrendingEnabled);
      y1b = Drt.detrend(y1);

      y2b = bicas.tf.kernel.apply_kernel(y1b, yKernelB, iKernelOrigin, edgePolicy);

      % NOTE: Using frequency domain-TF for scaling.
      y2 = Drt.retrend(y2b, tf(0));
    end



  end    % methods(Static)



end
