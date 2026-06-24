%
% UNFINISHED
%
% matlab.unittest automatic test code for bepic.spinfit.spin_aligned().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_spin_aligned___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty(T)
      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = int64.empty( 0, 1), ...
        spinPhaseRadAr      = double.empty(0, 1), ...
        samplesAr           = double.empty(0, 1), ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = 0, ...
        dataGapMinNs        = int64(1));

      T.assertEqual(actR.tt2000Ar, int64.empty(0, 1))
    end



    function test_0(T)
      N_IN = 32;

      tt2000Ar       = int64(100e9 + linspace(1.5e9, 3.5e9, N_IN))';
      spinPhaseRadAr = wrapTo2Pi(    linspace(  -pi,    pi, N_IN))';
      samplesAr      = 2 ...
        + 3*cos(  spinPhaseRadAr) + 4*sin(  spinPhaseRadAr) ...
        + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = tt2000Ar, ...
        spinPhaseRadAr      = spinPhaseRadAr, ...
        samplesAr           = samplesAr, ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = 0, ...
        dataGapMinNs        = int64(10e9));

      T.assertEqual(actR.tt2000Ar,          int64([102.5e9]))
      T.verifyEqual(actR.offsetAr,          [2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6], AbsTol=1e-13)
    end



  end    % methods(Test)



end
