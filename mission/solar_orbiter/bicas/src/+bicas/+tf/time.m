%
% Code for applying a TF in the time domain. In practice, one main function
% plus helper functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef time
  % PROPOSAL: MTEST code.
  %   PROPOSAL: Use data from L1R datasets. Compare freq. and time domain
  %             applications.
  %     PROPOSAL: Load directly from CDF.
  %       PROBLEM: Needs HK:DLR, mux mode, ACHG/ACLG
  %   PROPOSAL: Good enough plots for meetings.
  %
  % PROPOSAL: Separate argument for Hann window size ("radius"/"diameter")
  %           independently of kernel size.
  %   CON: Function has argument lenKernel. Hann window sets values outside
  %        window to zero.
  %        ==> Always implicit that Hann window should have the same "radius"
  %            as the kernel. Otherwise the caller would have shortened the
  %            kernel.
  % PROPOSAL: Separate arguments for kernel length before and after origin.
  %
  % PROPOSAL: Separate public function for obtaining the exact kernel used,
  %           including Hann window and any other possible modification.
  %   PRO: Can be used for manually inspecting the kernel.
  %
  % PROBLEM(?): Kernels contain non-zero samples on the wrong side of the
  %             kernel origin.
  %   THEORY: Should be because
  %           applying bicas.tf.freq.apply_TF(dtSec, yImpulse, tf) to the
  %           impulse, the kernel is treated as if it was periodic, meaning
  %           that the impulse response "tail" wraps around to the other side.
  %     THEORY: The "wrong side" of the kernel origin should change more than
  %             the other side when the kernel length is changed. In the limit
  %             of long kernel lengths, the right side should not change at
  %             all.
  %     TODO: Investigate with bicas.tools.tfkernel_playground.
  %     TODO: Compare impulse response with different kernel lengths.
  %     PROPOSAL: Force wrong side of kernel to have zero-valued samples.
  %       PROPOSAL: Set kernel origin at the edge of the kernel samples array.
  %         PROBLEM: Currently has no arguments for configuring this.
  %           PROPOSAL: Argument for whether kernel origin should be at the
  %                     beginning or end of kernel samples.
  %           PROPOSAL: Two kernel length arguments: Kernel length before &
  %                     after kernel origin.
  %             PROPOSAL: Exactly one of the two kernel arguments has to be
  %                       non-zero.
  %
  % TODO: Research whether MATLAB contains library functions for this
  %       functionality.
  %
  % PROPOSAL: Permit submitting multiple y1 at the same time (same dt, same
  %           tf).
  %   PRO: Faster?
  %     PRO: Only needs to obtain kernel (bicas.tf.freq.apply_TF() and Hann
  %          window) once.



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
    % dtSec
    %       Scalar, numeric. Seconds. Time between samples.
    % y1
    %       Column vector. Signal.
    % tf
    %       Function handle. Z = tf(omegaRps). Transfer function.
    % lenKernel
    %       Length of kernel that shall be generated from tf.
    %       NOTE: The kernel covers both delayed and advanced (future) samples,
    %       i.e. the kernel covers time both before and after an impulse.
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
        dtSec, y1, tf, lenKernel, edgePolicy, hannWindowEnabled)
      %
      %=============
      % ~ASSERTIONS
      %=============
      % TODO-DEC: Which argument assertions should one bother to have?
      % bicas.tf.freq.apply_TF() and bicas.tf.kernel.apply_kernel() check most
      % arguments.
      assert(isnumeric(lenKernel),   'lenKernel is not numeric.')
      assert(isscalar(lenKernel),    'lenKernel is not scalar.')
      assert(mod(lenKernel, 1) == 0, 'lenKernel is not an integer.')
      assert(lenKernel>0,            'lenKernel is not positive.')
      assert(islogical(hannWindowEnabled) & isscalar(hannWindowEnabled))



      %===========================
      % Obtain time domain kernel
      %===========================
      [yKernel, iKernelOrigin] = bicas.tf.time.get_raw_kernel(...
        lenKernel, dtSec, tf);



      if hannWindowEnabled
        %=============================
        % Apply Hann window to kernel
        %=============================
        truncatedHannWin = bicas.tf.time.get_truncated_Hann_window(...
          iKernelOrigin, lenKernel);
        yKernel          = yKernel .* truncatedHannWin;
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
    %
    function [yKernel, iKernelOrigin] = get_raw_kernel(lenKernel, dtSec, tf)
      % PROPOSAL: Better name
      %   kernel
      %   impulse, impulse response
      %   TF
      %   NOTE: Without modifications, e.g. Hann window.
      %   raw
      %   --
      %   get_TF_kernel()
      %   get_impulse_kernel()
      %   get_impulse_response_kernel()

      iKernelOrigin           = floor(1 + (lenKernel-1)/2);
      yImpulse                = zeros(lenKernel, 1);
      yImpulse(iKernelOrigin) = 1;
      yKernel                 = bicas.tf.freq.apply_TF(dtSec, yImpulse, tf);

      if 0
        %=============
        % Plot kernel
        %=============
        % NOTE: See bicas.tools.tfkernel_playground.

        % t = 0 <=> iKernelOrigin
        t = [(-iKernelOrigin+1):(lenKernel-iKernelOrigin)] * dtSec;
        titleStr = sprintf(...
          "lenKernel=%d; iKernelOrigin=%d; dtSec=%d; lenKernel*dtSec=%g", ...
          lenKernel, iKernelOrigin, dtSec, lenKernel*dtSec);

        figure
        plot(t, yKernel, '.-')
        title(titleStr)
        xlabel('t [s]')
      end
    end



    % Obtain a Hann window that can be applied to the kernel, with the maximum
    % value=1 at iMax. Only one element has the maximum value.
    %
    function truncatedHannWin = get_truncated_Hann_window(iMax, n)
      % IMPLEMENTATION NOTE: Does not use MATLAB's hann() function.
      % PRO: hann() max value varies depending on parameters.
      %   (aperiodic), n=even : max<1; has two max elements
      %   (aperiodic), n=odd  : max=1; has one max element
      %   'periodic',  n=even : max=1; has one max element
      %   'periodic',  n=odd  : max<1; has two max elements
      %   It is however useful to (1) always have max value=1 and (2) better
      %   (probably?; or at least more elegant/symmetric) to have the max value
      %   in exactly one element.
      % PRO: One can not set the center and "radius" of the Hann window using
      %      hann().
      % It is therefore easier to just implement the Hann window function
      % oneself.

      if n == 1
        truncatedHannWin = 1;
      else
        % "Radius" of the hann window. Sets the boundary on the first element
        % *OUTSIDE* the kernel.
        R = (max(iMax-1, n-iMax) + 1);

        % Normalized locations within the Hann window:
        %   center:     x=0
        %   boundaries: |x|=1
        %     Where Hann weight=0
        xAr = ([1:n]' - iMax) / R;

        truncatedHannWin = cos(pi/2*xAr).^2;
      end
    end



  end    % methods(Static)



end
