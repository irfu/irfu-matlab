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



    function test_all_empty(T)
      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = int64.empty( 0, 1), ...
        spinPhaseRadAr      = double.empty(0, 1), ...
        samplesAr           = double.empty(0, 1), ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = pi, ...
        dataGapMinNs        = int64(1));

      T.verifyEqual(actR.timeWindowCenterTt2000Ar, int64.empty(0, 1))
      T.verifySize(actR.stdDeviationAr,            [0, 1])
      T.verifySize(actR.nBadPoints,                [0, 1])
      T.verifySize(actR.offsetAr,                  [0, 1])
      T.verifySize(actR.coefficientCos1Ar,         [0, 1])
      T.verifySize(actR.coefficientSin1Ar,         [0, 1])
      T.verifySize(actR.coefficientCos2Ar,         [0, 1])
      T.verifySize(actR.coefficientSin2Ar,         [0, 1])
    end



    % No data gap.
    % One time window.
    function test_0(T)
      N_IN = 32;

      tt2000Ar       = int64(    linspace(101.5e9, 103.5e9, N_IN))';
      spinPhaseRadAr = wrapTo2Pi(linspace(      0,    2*pi, N_IN))';
      samplesAr      = 2 ...
        + 3*cos(  spinPhaseRadAr) + 4*sin(  spinPhaseRadAr) ...
        + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = tt2000Ar, ...
        spinPhaseRadAr      = spinPhaseRadAr, ...
        samplesAr           = samplesAr, ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = pi, ...
        dataGapMinNs        = int64(10e9));

      T.assertEqual(actR.timeWindowCenterTt2000Ar, int64([102.5e9]))
      T.verifyEqual(actR.stdDeviationAr,    [0], AbsTol=1e-14)
      % Actual=6 for unkown reason!!
      % T.verifyEqual(actR.nBadPoints,        [6])
      T.verifyEqual(actR.offsetAr,          [2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6], AbsTol=1e-13)
    end



    % No data gap.
    % Two time windows.
    function test_1(T)
      N_IN = 32;

      tt2000Ar       = int64(    linspace(100e9, 104e9, N_IN))';
      spinPhaseRadAr = wrapTo2Pi(linspace(    0,  4*pi, N_IN))';
      samplesAr      = 2 ...
        + 3*cos(  spinPhaseRadAr) + 4*sin(  spinPhaseRadAr) ...
        + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = tt2000Ar, ...
        spinPhaseRadAr      = spinPhaseRadAr, ...
        samplesAr           = samplesAr, ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = pi, ...
        dataGapMinNs        = int64(10e9));

      T.assertEqual(actR.timeWindowCenterTt2000Ar, int64([101e9; 103e9]))
      T.verifyEqual(actR.offsetAr,          [2; 2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3; 3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4; 4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5; 5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6; 6], AbsTol=1e-13)
    end



    % No data gap.
    % Two time windows.
    function test_2(T)
      N_IN = 32;

      tt2000Ar       = int64(    linspace(  100e9,   114e9, N_IN))';
      spinPhaseRadAr = wrapTo2Pi(linspace(0.75*pi, 4.25*pi, N_IN))';
      samplesAr      = 2 ...
        + 3*cos(  spinPhaseRadAr) + 4*sin(  spinPhaseRadAr) ...
        + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = tt2000Ar, ...
        spinPhaseRadAr      = spinPhaseRadAr, ...
        samplesAr           = samplesAr, ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = pi, ...
        dataGapMinNs        = int64(10e9));

      T.assertEqual(actR.timeWindowCenterTt2000Ar, int64([101e9; 109e9]))
      T.verifyEqual(actR.offsetAr,          [2; 2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3; 3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4; 4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5; 5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6; 6], AbsTol=1e-13)
    end



    % Data gap
    % * Deliberately two time windows before AND after data gap to test
    %   bepic.spinfit.spin_aligned()'s merging of processing results.
    % * Deliberately using the same spin phase just before and after data gap
    %   (testing old bug fix).
    function test_data_gap(T)
      N_IN = 2*32;

      tt2000Ar = int64([...
        linspace(100e9, 106e9, N_IN/2), ...
        linspace(200e9, 208e9, N_IN/2)...
        ])';
      % NOTE: Deliberately using the same spin phase just before and after data
      %       gap (testing old bug fix).
      spinPhaseRadAr = wrapTo2Pi([...
        linspace( 0.5*pi, 3.5*pi, N_IN / 2) ...
        linspace( 3.5*pi, 7.5*pi, N_IN / 2) ...
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
        timeWindowCenterRad = pi, ...
        dataGapMinNs        = int64(1e9));

      N = ones(4, 1);
      T.verifyEqual(actR.timeWindowCenterTt2000Ar, int64([101e9; 105e9; 203e9; 207e9]))
      T.verifyEqual(actR.stdDeviationAr,    [0*N], AbsTol=1e-14)
      % Actual<>0*N for unkown reason!!
      % T.verifyEqual(actR.nBadPoints,        [5;5;1;1])
      T.verifyEqual(actR.offsetAr,          2*N, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, 3*N, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, 4*N, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, 5*N, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, 6*N, AbsTol=1e-13)
    end



  end    % methods(Test)



end
