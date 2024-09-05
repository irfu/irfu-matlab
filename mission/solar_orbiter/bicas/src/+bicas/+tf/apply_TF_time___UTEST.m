%
% matlab.unittest automatic test code for bicas.tf.apply_TF_time().
%
% 2021-08-12: Does not test detrending/retrending. Low priority since tested
% elsewhere.
%
%
% NAMING CONVENTION
% =================
% HW : Hann Window
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-11
%
classdef apply_TF_time___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Odd-length kernels/HWs.



  properties(TestParameter)
  end



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
    function test0(testCase)
      % PROPOSITION: Almost unnecessary since bicas.tf.apply_TF_kernel() tests for same.
      %   CON: Tests conversion TF-->kernel.

      dt          = 0.1;
      N           = 100;
      nDelaySmpls = 10;
      tf = bicas.tf.utest_utils.get_TF_delay(nDelaySmpls*dt);

      t  = [0:N-1]' * dt;
      y1 = 3 + 2*t - exp(-t);    % Arbitrary input signal.

      %=====================
      % edgePolicy == ZEROS
      %=====================
      % Manually shift/delay signal.
      y2_exp = [zeros(nDelaySmpls, 1); y1(1:end-nDelaySmpls)];

      y2 = bicas.tf.apply_TF_time(dt, y1, tf, N, 'ZEROS');

      testCase.verifyEqual(y2, y2_exp, 'AbsTol', 1e-13)

      %======================
      % edgePolicy == MIRROR
      %======================
      % Manually pad with mirrored samples and shift.
      y2_exp = [y1(nDelaySmpls:-1:1); y1(1:end-nDelaySmpls)];

      y2 = bicas.tf.apply_TF_time(dt, y1, tf, N, 'MIRROR');

      testCase.verifyEqual(y2, y2_exp, 'AbsTol', 1e-13)
    end



    % Constant delay+HW. ==> Constant scaling due to HW.
    % Even-length kernel/HW.
    %
    function test_HW0(testCase)
      dt     = 0.1;
      nSmpls = 100;

      t      = [0:nSmpls-1]' * dt;
      y1     = 2 + 3*t - exp(-t);    % Arbitrary input signal.

      % NOTE: Use both even and odd numbers.
      for lenKernel = [7, 17, 32]

        nds0 = 3;
        nds1 = floor(lenKernel / 8);
        % NOTE: Yields hwFactor=0 for even lenKernel!
        nds2 = floor(lenKernel / 2);

        for nDelaySmpls = [nds0, nds1, nds2]

          %============================================================
          % Factor by which the signal will be weakened/multiplied due
          % to the Hann window.
          %============================================================
          % IMPLEMENTATION NOTE: The effect of the Hann window can be
          % boiled down to multiplication with a scalar number since
          % the kernel only has one non-zero value.
          % --
          if mod(lenKernel, 2) == 0
            % CASE: EVEN-numbered length kernel
            % Normalized delay: -1 <= x <= 1
            xHw      = nDelaySmpls/lenKernel;
            hwFactor = 0.5 * (1 + cos(xHw*2*pi));
          else
            % CASE: ODD-numbered length kernel
            % ==> Has to compensate for rounding to samples.
            xHw      = (nDelaySmpls-0.5)/lenKernel;
            hwFactor = 0.5 * (1 + cos(xHw*2*pi));
          end

          tf = bicas.tf.utest_utils.get_TF_delay(nDelaySmpls*dt);

          %=====================
          % edgePolicy == ZEROS
          %=====================
          % Manually delay, and multiply Hann window factor.
          y2exp = [zeros(nDelaySmpls, 1); y1(1:end-nDelaySmpls)] * hwFactor;
          y2act = bicas.tf.apply_TF_time(...
            dt, y1, tf, lenKernel, ...
            'ZEROS', 'hannWindow', true);

          testCase.verifyEqual(y2act, y2exp, 'AbsTol', 1e-13)

          %======================
          % edgePolicy == MIRROR
          %======================
          % Manually pad with mirrored samples, delay, and multiply
          % Hann window factor.
          y2exp = [y1(nDelaySmpls:-1:1); y1(1:end-nDelaySmpls)] * hwFactor;
          y2act = bicas.tf.apply_TF_time(...
            dt, y1, tf, lenKernel, ...
            'MIRROR', 'hannWindow', true);

          testCase.verifyEqual(y2act, y2exp, 'AbsTol', 1e-13)

        end    % for
      end    % for

    end    % function



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)
  end    % methods(Static, Access=private)



end
