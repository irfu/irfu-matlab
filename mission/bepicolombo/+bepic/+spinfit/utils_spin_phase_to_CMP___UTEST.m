%
% matlab.unittest automatic test code for
% bepic.spinfit.utils.spin_phase_to_cumulative_spin_phase().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils_spin_phase_to_CMP___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty(T)
      T.test(...
        double.empty(0, 1), ...
        double.empty(0, 1))
    end



    function test_one_value(T)
      T.test(...
        [3], ...
        [3])
    end



    function test_0(T)
      T.test(...
        [0.2; 0.3] * 2*pi, ...
        [0.2; 0.3] * 2*pi)
    end



    function test_1(T)
      % Long jump (<~2*pi) but not wrapping.
      T.test(...
        [0.1; 0.9] * 2*pi, ...
        [0.1; 0.9] * 2*pi)
    end



    function test_2(T)
      % Long jump (<~2*pi) but wrapping.
      T.test(...
        [0.5; 0.4] * 2*pi, ...
        [0.5; 1.4] * 2*pi)
    end



    function test_3(T)
      actCumulSpinPhaseRadAr = bepic.spinfit.utils.spin_phase_to_cumulative_spin_phase(...
        [0.2, 0.3, 0.7, 0.1, 0.5, 0.4]' * 2*pi);

      T.assertEqual(...
        actCumulSpinPhaseRadAr, ...
        [0.2, 0.3, 0.7, 1.1, 1.5, 2.4]' * 2*pi, ...
        RelTol=1e-14)
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Internal helper function to shorten some tests.
    %
    % NOTE: Can not set T.assertEqual() arguments AbsTol or RelTol.
    %
    function test(T, spinPhaseRadAr, expCumulSpinPhaseRadAr)

      actCumulSpinPhaseRadAr = bepic.spinfit.utils.spin_phase_to_cumulative_spin_phase(...
        spinPhaseRadAr);

      T.assertEqual(actCumulSpinPhaseRadAr, expCumulSpinPhaseRadAr)
    end



  end    % methods(Access=private)



end
