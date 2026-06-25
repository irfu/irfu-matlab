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
    % These values may or may not be the ones used outside the test code.
    TIME_WINDOW_PERIOD_NS     = int64(4e9);
    TIME_WINDOW_LENGTH_NS     = int64(4e9);
    TIME_WINDOW_CENTER_TT2000 = int64(2e9);
    N_MIN_FIT_SAMPLES         = 5+3;
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
        timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
        nMinFitSamples         = T.N_MIN_FIT_SAMPLES);

      T.verifySize(actR.timeWindowCenterTt2000Ar, [0, 1])
      T.verifySize(actR.stdDeviationAr,           [0, 1])
      T.verifySize(actR.nBadPoints,               [0, 1])
      T.verifySize(actR.offsetAr,                 [0, 1])
      T.verifySize(actR.coefficientCos1Ar,        [0, 1])
      T.verifySize(actR.coefficientSin1Ar,        [0, 1])
      T.verifySize(actR.coefficientCos2Ar,        [0, 1])
      T.verifySize(actR.coefficientSin2Ar,        [0, 1])
    end



    % All samples = NaN
    % Samples in two time windows.
    function test_NaN_samples(T)
      % PROPOSAL: Vary N_IN.
      N_IN = 30;

      actR = bepic.spinfit.mms_spinfit_wrapper( ...
        tt2000Ar       = int64(linspace(101e9, 107e9, N_IN))', ...
        spinPhaseRadAr = wrapTo2Pi(linspace(pi, 3*pi, N_IN))', ...
        samplesAr      = NaN(N_IN, 1), ...
        timeWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        timeWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
        nMinFitSamples         = T.N_MIN_FIT_SAMPLES);

      NAN_AR = NaN(2, 1);
      T.verifyEqual(actR.timeWindowCenterTt2000Ar, int64([102e9; 106e9]))
      T.verifyEqual(actR.stdDeviationAr,           NAN_AR)
      T.verifyEqual(actR.nBadPoints,               [0; 0])
      T.verifyEqual(actR.offsetAr,                 NAN_AR)
      T.verifyEqual(actR.coefficientCos1Ar,        NAN_AR)
      T.verifyEqual(actR.coefficientSin1Ar,        NAN_AR)
      T.verifyEqual(actR.coefficientCos2Ar,        NAN_AR)
      T.verifyEqual(actR.coefficientSin2Ar,        NAN_AR)
    end



    % Samples can be fit perfectly with nonzero values for all coefficients.
    % Samples in two time windows.
    function test_spin_samples_perfect_fit(T)
      % PROPOSAL: Split in three separate test functions.
      %   PROPOSAL: Shared function for calling
      %     bepic.spinfit.mms_spinfit_wrapper() with standard values.
      %   PROPOSAL: Shared function for generating samples.

      N_IN = 64;

      function actR = test(beginTt2000, endTt2000)
        tt2000Ar       = int64(linspace(beginTt2000, endTt2000, N_IN))';
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
          timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
          nMinFitSamples         = T.N_MIN_FIT_SAMPLES);
      end

      % Timestamps to return two time windows.
      actR = test(101e9, 107e9);

      T.verifyEqual(actR.timeWindowCenterTt2000Ar, int64([102e9; 106e9]))
      T.verifyEqual(actR.offsetAr,          [2; 2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3; 3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4; 4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5; 5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6; 6], AbsTol=1e-13)

      % Timestamps to return only the second time window.
      actR = test(103e9, 107e9);

      T.verifyEqual(actR.timeWindowCenterTt2000Ar, int64([106e9]))
      T.verifyEqual(actR.offsetAr,          [2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6], AbsTol=1e-13)

      % Timestamps to return only the first time window.
      actR = test(101e9, 105e9);

      T.verifyEqual(actR.timeWindowCenterTt2000Ar, int64([102e9]))
      T.verifyEqual(actR.offsetAr,          [2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1Ar, [3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, [4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, [5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, [6], AbsTol=1e-13)
    end



    % Samples can be fit perfectly with nonzero values for all coefficients.
    % Samples in two time windows.
    function test_spin_samples_too_few(T)
      N_IN = 6+6;

      tt2000Ar       = int64(linspace(1e9, 7e9, N_IN))';
      spinPhaseRadAr = wrapTo2Pi(linspace( pi, 5*pi, N_IN))';
      samplesAr      = 2 + 3*cos(  spinPhaseRadAr);

      actR = bepic.spinfit.mms_spinfit_wrapper( ...
        tt2000Ar       = tt2000Ar, ...
        spinPhaseRadAr = spinPhaseRadAr, ...
        samplesAr      = samplesAr, ...
        timeWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        timeWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        timeWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
        nMinFitSamples         = T.N_MIN_FIT_SAMPLES);

      NAN_AR = NaN(2, 1);
      T.verifyEqual(actR.offsetAr,          NAN_AR, AbsTol=1e-15)
      T.verifyEqual(actR.coefficientCos1Ar, NAN_AR, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1Ar, NAN_AR, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2Ar, NAN_AR, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2Ar, NAN_AR, AbsTol=1e-13)
    end



  end    % methods(Test)



end
