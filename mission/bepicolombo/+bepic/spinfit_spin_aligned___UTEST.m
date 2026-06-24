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



    function test_data_gap(T)
      N_IN = 2*32;

      tt2000Ar       = int64(100e9 + [...
        linspace( 1.5e9,  3.5e9, N_IN/2), ...
        linspace(11.5e9, 13.5e9, N_IN/2)...
        ])';
      % IMPLEMENTATION NOTE: Must not have any identical spin phase values,
      % which will then be translated to identical fake TT2000 values. ==> Crash
      % This indicates a bug.
      spinPhaseRadAr = wrapTo2Pi([...
        linspace(  -pi+0.01, pi-0.01, N_IN / 2) ...
        linspace(  -pi,      pi,      N_IN / 2) ...
        ])';
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
        dataGapMinNs        = int64(1e9));

      T.assertEqual(actR.tt2000Ar,          int64([102.5e9; 112.5e9]))
      T.verifyEqual(actR.offsetAr,          [2; 2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3; 3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4; 4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5; 5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6; 6], AbsTol=1e-13)
    end



  end    % methods(Test)



end
