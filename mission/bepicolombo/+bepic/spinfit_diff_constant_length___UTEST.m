%
% matlab.unittest automatic test code for bepic.spinfit.diff_constant_length().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_diff_constant_length___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_arrays(T)
      actR = bepic.spinfit.diff_constant_length( ...
        tt2000Ar           = int64.empty( 0, 1), ...
        spinPhaseRadiansAr = double.empty(0, 1), ...
        samplesAr          = double.empty(0, 1) ...
      );

      T.assertSize(actR.offsetAr,          [0, 1])
      T.assertSize(actR.coefficientCos1Ar, [0, 1])
      T.assertSize(actR.coefficientSin1Ar, [0, 1])
      T.assertSize(actR.coefficientCos2Ar, [0, 1])
      T.assertSize(actR.coefficientSin2Ar, [0, 1])
    end



    function test_NaN_samples(T)
      % PROPOSAL: Vary N_IN.
      N_IN = 30;

      actR = bepic.spinfit.diff_constant_length( ...
        tt2000Ar           = int64(linspace(1, 7e9, N_IN))', ...
        spinPhaseRadiansAr = linspace(pi, 3*pi, N_IN)', ...
        samplesAr          = NaN(N_IN, 1) ...
      );

      T.assertEqual(actR.offsetAr,          [nan; nan])
      T.assertEqual(actR.coefficientCos1Ar, [nan; nan])
      T.assertEqual(actR.coefficientSin1Ar, [nan; nan])
      T.assertEqual(actR.coefficientCos2Ar, [nan; nan])
      T.assertEqual(actR.coefficientSin2Ar, [nan; nan])
    end



    function test_samples_constant(T)
      N_IN = 30;

      % Data in two time windows.
      actR = bepic.spinfit.diff_constant_length( ...
        tt2000Ar           = int64(linspace(1e9, 7e9, N_IN))', ...
        spinPhaseRadiansAr = linspace(pi, 3*pi, N_IN)', ...
        samplesAr          = 3*ones(N_IN, 1) ...
      );

      % NOTE: Testing the offset return value.
      T.verifyEqual(actR.offsetAr,          [3; 3], RelTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [0; 0], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [0; 0], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [0; 0], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [0; 0], AbsTol=1e-13)
    end



    function test_spin_samples_perfect_fit(T)
      N_IN = 3*16;

      % Data in two time windows.
      tt2000Ar       = int64(linspace(1e9, 7e9, N_IN))';
      spinPhaseRadAr = linspace( pi, 5*pi, N_IN)';
      samplesAr      = 2 + 3*cos(spinPhaseRadAr) + 4*sin(spinPhaseRadAr) + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

      actR = bepic.spinfit.diff_constant_length( ...
        tt2000Ar           = tt2000Ar, ...
        spinPhaseRadiansAr = spinPhaseRadAr, ...
        samplesAr          = samplesAr ...
      );

      T.verifyEqual(actR.offsetAr,          [2; 2], AbsTol=1e-15)
      T.verifyEqual(actR.coefficientCos1Ar, [3; 3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4; 4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5; 5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6; 6], AbsTol=1e-13)
    end



  end    % methods(Test)



end
