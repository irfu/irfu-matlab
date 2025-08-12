%
% matlab.unittest automatic test code for bicas.tf.time.
%
%
% LOCAL NAMING CONVENTIONS
% ========================
% HW : Hann Window
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-11
%
classdef time___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)
    % TODO-DEC: Non-trivial tests.
    %   PROPOSAL: Hann window with unity TF.
    %   PROPOSAL: Hann window with max delay TF. ==> Kernel=0.
    %   PROPOSAL: Hann window with delay=max_delay/2
    %             ==> Kernel=half amplitude.
    %   PROPOSAL: Delay TF with delay longer than kernel?
    %
    % PROPOSAL: Odd-length kernel+Hann window
    %
    % PROPOSAL: General delay TF, both odd- and even-length kernel/HW.
    %           Iterate over kernel lengths and delays (integer samples).
    %
    % PROPOSAL: Detrending/retrending.



    % No HW, no detrending.
    %
    function test_apply_TF___0(T)
      % PROPOSITION: Almost unnecessary since bicas.tf.kernel.apply_kernel() tests for the same.
      %   CON: Tests conversion TF-->kernel.

      dt          = 0.1;
      N           = 100;
      nDelaySmpls = 10;
      tf = bicas.tf.utest_utils.get_TF_delay(nDelaySmpls*dt);

      t  = [0:N-1]' * dt;
      % NOTE: Arbitrary NON-TRIVIAL input signal with magnitude of order ~1.
      y1 = 3 + 2*t - exp(-t);

      %=====================
      % edgePolicy == ZEROS
      %=====================
      % Manually shift/delay signal.
      expY2 = [zeros(nDelaySmpls, 1); y1(1:end-nDelaySmpls)];

      actY2 = bicas.tf.time.apply_TF(dt, y1, tf, N, 'ZEROS', false);

      T.verifyEqual(actY2, expY2, 'AbsTol', 1e-13)

      %======================
      % edgePolicy == MIRROR
      %======================
      % Manually pad with mirrored samples and shift.
      expY2 = [y1(nDelaySmpls:-1:1); y1(1:end-nDelaySmpls)];

      actY2 = bicas.tf.time.apply_TF(dt, y1, tf, N, 'MIRROR', false);

      T.verifyEqual(actY2, expY2, 'AbsTol', 1e-13)
    end



    % Constant delay+HW. ==> Constant scaling due to HW.
    %
    function test_apply_TF___HW_delay_0(T)
      dt     = 0.1;
      nSmpls = 100;

      t      = [0:nSmpls-1]' * dt;
      % NOTE: Arbitrary NON-TRIVIAL input signal with magnitude of order ~1.
      y1     = 2 + 3*t - exp(-t);

      % NOTE: Use both even and odd numbers.
      for lenKernel = [7, 17, 32]

        nds0 = 3;
        nds1 = floor(lenKernel / 8);
        % NOTE: Yields hwFactor=0 for even lenKernel!
        nds2 = floor(lenKernel / 2);

        for nDelaySmpls = [nds0, nds1, nds2]
          T.test_apply_TF___HW_delay( ...
            y1=y1, dt=dt, ...
            nDelaySmpls=nDelaySmpls, ...
            lenKernel  =lenKernel)
        end
      end

    end    % function



    % NOTE: Expected values have not been verified, but only seem plausible.
    %
    function test_get_shifted_Hann_window(T)
      function test(iMax, n, expShiftedHw)
        actShiftedHw = bicas.tf.time.get_shifted_Hann_window(iMax, n);
        T.assertEqual(actShiftedHw, expShiftedHw, "AbsTol", 1e-15)
      end

      % iMax=1
      test(1, 1, [1])
      test(1, 2, [1    0]')
      test(1, 3, [0.75 0.75 0.00]')
      test(1, 4, [1.00 0.50 0.00 0.50]')

      % iMax=end
      test(1, 1, [1])
      test(2, 2, [0    1]')
      test(3, 3, [0.75 0.00 0.75]')
      test(4, 4, [0.50 0.00 0.50 1.00]')

      % iMax = middle
      test(2, 4, [0.50 1.00 0.50 0.00]')
      test(3, 5, [...
        0; ...
        0.345491502812526; ...
        0.904508497187474; ...
        0.904508497187474; ...
        0.345491502812526])
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Test delaying arbitrary signal by an arbitrary amount. Assume slight
    % numerical errors.
    %
    function test_apply_TF___HW_delay(T, A)
      arguments
        T
        A.y1
        A.lenKernel
        A.nDelaySmpls
        A.dt
      end

      %===========================================================
      % Derive factor by which the signal will be (very slightly)
      % weakened/multiplied due to the Hann window.
      %===========================================================
      % IMPLEMENTATION NOTE: The effect of the Hann window can be boiled down
      % to multiplication with a scalar number since the kernel (but not the
      % frequency-domain TF) has only one non-zero value.
      % NOTE: This requires knowledge of the implementation.
      if mod(A.lenKernel, 2) == 0
        % CASE: EVEN-numbered length kernel
        % Normalized delay: -1 <= x <= 1
        xHw = A.nDelaySmpls/A.lenKernel;
      else
        % CASE: ODD-numbered length kernel
        % ==> Has to compensate for rounding to samples.
        xHw = (A.nDelaySmpls-0.5)/A.lenKernel;
      end
      hwFactor = 0.5 * (1 + cos(xHw*2*pi));

      tf = bicas.tf.utest_utils.get_TF_delay(A.nDelaySmpls*A.dt);



      %=====================
      % edgePolicy == ZEROS
      %=====================
      % Manually delay, and multiply with Hann window factor.
      expY2 = [zeros(A.nDelaySmpls, 1); A.y1(1:end-A.nDelaySmpls)] * hwFactor;
      actY2 = bicas.tf.time.apply_TF(...
        A.dt, A.y1, tf, A.lenKernel, 'ZEROS', true);

      T.assertEqual(actY2, expY2, 'AbsTol', 1e-13)

      %======================
      % edgePolicy == MIRROR
      %======================
      % Manually pad with mirrored samples, delay, and multiply by Hann window
      % factor.
      expY2 = [A.y1(A.nDelaySmpls:-1:1); A.y1(1:end-A.nDelaySmpls)] * hwFactor;
      actY2 = bicas.tf.time.apply_TF(...
        A.dt, A.y1, tf, A.lenKernel, 'MIRROR', true);

      T.assertEqual(actY2, expY2, 'AbsTol', 1e-13)
    end



  end    % methods(Access=private)



end
