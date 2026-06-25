%
% matlab.unittest automatic test code for
% bepic.spinfit.spin_phase_to_cumulative_spin_phase().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_spin_phase_to_cumulative_spin_phase___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_spin_phase_to_cumulative_spin_phase_empty(T)
      T.test_spin_phase_to_cumulative_spin_phase(...
        double.empty(0, 1), ...
        double.empty(0, 1))
    end



    function test_spin_phase_to_cumulative_spin_phase_one_value(T)
      T.test_spin_phase_to_cumulative_spin_phase(...
        [3], ...
        [3])
    end



    function test_spin_phase_to_cumulative_spin_phase_0(T)
      T.test_spin_phase_to_cumulative_spin_phase(...
        [0.2; 0.3] * 2*pi, ...
        [0.2; 0.3] * 2*pi)
    end



    function test_spin_phase_to_cumulative_spin_phase_1(T)
      % Long jump (<~2*pi) but not wrapping.
      T.test_spin_phase_to_cumulative_spin_phase(...
        [0.1; 0.9] * 2*pi, ...
        [0.1; 0.9] * 2*pi)
    end



    function test_spin_phase_to_cumulative_spin_phase_2(T)
      % Long jump (<~2*pi) but wrapping.
      T.test_spin_phase_to_cumulative_spin_phase(...
        [0.5; 0.4] * 2*pi, ...
        [0.5; 1.4] * 2*pi)
    end



    function test_spin_phase_to_cumulative_spin_phase_3(T)
      actCumulSpinPhaseRadAr = bepic.spinfit.spin_phase_to_cumulative_spin_phase(...
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
    function test_spin_phase_to_cumulative_spin_phase( ...
        T, spinPhaseRadAr, expCumulSpinPhaseRadAr)

      actCumulSpinPhaseRadAr = bepic.spinfit.spin_phase_to_cumulative_spin_phase(...
        spinPhaseRadAr);

      T.assertEqual(actCumulSpinPhaseRadAr, expCumulSpinPhaseRadAr)
    end



  end    % methods(Access=private)



end
