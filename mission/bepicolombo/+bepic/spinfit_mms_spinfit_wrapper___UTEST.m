%
% matlab.unittest automatic test code for bepic.spinfit.mms_spinfit_wrapper().
%
% NOTE: Does not currently test varying parameters which are expected to be
% (very) constant.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_mms_spinfit_wrapper___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Test for one sample, two samples, extrapolation of center
  % timestamps.



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    % Constants for when using constant-length time windows. These values may or
    % may not be the ones used outside the test code.
    TIME_WINDOW_PERIOD_NS     = int64(4e9);
    TIME_WINDOW_LENGTH_NS     = int64(4e9);
    TIME_WINDOW_CENTER_TT2000 = int64(2e9);
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_arrays(T)
      actR = bepic.spinfit.mms_spinfit_wrapper( ...
        tt2000Ar       = int64.empty( 0, 1), ...
        spinPhaseRadAr = double.empty(0, 1), ...
        samplesAr      = double.empty(0, 1), ...
        timeWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        timeWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000);

      T.verifyEqual(actR.tt2000Ar,         int64.empty(0, 1))
      T.verifySize(actR.offsetAr,          [0, 1])
      T.verifySize(actR.coefficientCos1Ar, [0, 1])
      T.verifySize(actR.coefficientSin1Ar, [0, 1])
      T.verifySize(actR.coefficientCos2Ar, [0, 1])
      T.verifySize(actR.coefficientSin2Ar, [0, 1])
    end



    % All samples = NaN
    % Samples in two time windows.
    function test_NaN_samples(T)
      % PROPOSAL: Vary N_IN.
      N_IN = 30;

      actR = bepic.spinfit.mms_spinfit_wrapper( ...
        tt2000Ar       = int64(linspace(1, 7e9, N_IN))', ...
        spinPhaseRadAr = wrapTo2Pi(linspace(pi, 3*pi, N_IN))', ...
        samplesAr      = NaN(N_IN, 1), ...
        timeWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        timeWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000);

      T.verifyEqual(actR.tt2000Ar,          int64([2e9; 6e9]))
      T.verifyEqual(actR.offsetAr,          [nan; nan])
      T.verifyEqual(actR.coefficientCos1Ar, [nan; nan])
      T.verifyEqual(actR.coefficientSin1Ar, [nan; nan])
      T.verifyEqual(actR.coefficientCos2Ar, [nan; nan])
      T.verifyEqual(actR.coefficientSin2Ar, [nan; nan])
    end



    % Samples = constant
    % Samples in two time windows.
    function test_samples_constant(T)
      N_IN = 30;

      actR = bepic.spinfit.mms_spinfit_wrapper( ...
        tt2000Ar       = int64(linspace(1e9, 7e9, N_IN))', ...
        spinPhaseRadAr = wrapTo2Pi(linspace(pi, 3*pi, N_IN))', ...
        samplesAr      = 3*ones(N_IN, 1), ...
        timeWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        timeWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000);

      % NOTE: Testing the offset return value.
      T.verifyEqual(actR.tt2000Ar,          int64([2e9; 6e9]))
      T.verifyEqual(actR.offsetAr,          [3; 3], RelTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [0; 0], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [0; 0], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [0; 0], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [0; 0], AbsTol=1e-13)
    end



    % Samples can be fit perfectly with nonzero values for all coefficients.
    % Samples in two time windows.
    function test_spin_samples_perfect_fit(T)
      N_IN = 3*16;

      tt2000Ar       = int64(linspace(2e9, 6e9, N_IN))';
      spinPhaseRadAr = wrapTo2Pi(linspace( pi, 3*pi, N_IN))';
      samplesAr      = 2 ...
        + 3*cos(  spinPhaseRadAr) + 4*sin(  spinPhaseRadAr) ...
        + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

      actR = bepic.spinfit.mms_spinfit_wrapper( ...
        tt2000Ar       = tt2000Ar, ...
        spinPhaseRadAr = spinPhaseRadAr, ...
        samplesAr      = samplesAr, ...
        timeWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        timeWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000);

      T.verifyEqual(actR.tt2000Ar,          int64([2e9; 6e9]))
      T.verifyEqual(actR.offsetAr,          [2; 2], AbsTol=1e-15)
      T.verifyEqual(actR.coefficientCos1Ar, [3; 3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4; 4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5; 5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6; 6], AbsTol=1e-13)
    end



  end    % methods(Test)



end
