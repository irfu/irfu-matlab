%
% matlab.unittest automatic test code for bepic.spinfit.fit_TAFW().
%
% NOTE: Does not currently test varying parameters which are expected to be
% (very) constant.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_fit_TAFW___UTEST < matlab.unittest.TestCase
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
    N_MIN_FIT_SAMPLES         = 5+1;
  end



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  % Technically, additional properties of testCase objects with cell array
  % default values. Test methods with arguments with the same name will be
  % called once for every element in the cell arrays.
  properties(TestParameter)
    N_FIT_COEFF = {3, 5}
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_arrays(T, N_FIT_COEFF)

      actR = bepic.spinfit.fit_TAFW( ...
        tt2000Ar       = int64.empty( 0, 1), ...
        spinPhaseRadAr = double.empty(0, 1), ...
        samplesAr      = double.empty(0, 1), ...
        fitWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        fitWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        fitWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
        nMinFitSamples         = T.N_MIN_FIT_SAMPLES, ...
        nFitCoefficients       = N_FIT_COEFF);

      T.verifySize(actR.fitWindowCenterTt2000, [0, 1])
      T.verifySize(actR.stdDeviation,          [0, 1])
      T.verifySize(actR.nBadPoints,            [0, 1])
      T.verifySize(actR.offset,                [0, 1])
      T.verifySize(actR.coefficientCos1,       [0, 1])
      T.verifySize(actR.coefficientSin1,       [0, 1])

      switch(N_FIT_COEFF)
        case 3
          T.verifyTrue(~ismember("coefficientCos2Ar", actR.Properties.VariableNames))
          T.verifyTrue(~ismember("coefficientSin2Ar", actR.Properties.VariableNames))
        case 5
          T.verifySize(actR.coefficientCos2,   [0, 1])
          T.verifySize(actR.coefficientSin2,   [0, 1])
      end
    end



    % All samples = NaN
    % Samples in two fit windows.
    function test_NaN_samples(T, N_FIT_COEFF)
      % PROPOSAL: Vary N_IN.
      N_IN = 30;

      actR = bepic.spinfit.fit_TAFW( ...
        tt2000Ar       = int64(linspace(101e9, 107e9, N_IN))', ...
        spinPhaseRadAr = wrapTo2Pi(linspace(pi, 3*pi, N_IN))', ...
        samplesAr      = NaN(N_IN, 1), ...
        fitWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        fitWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        fitWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
        nMinFitSamples         = T.N_MIN_FIT_SAMPLES, ...
        nFitCoefficients       = N_FIT_COEFF);

      NAN_AR = NaN(2, 1);
      T.verifyEqual(actR.fitWindowCenterTt2000, int64([102e9; 106e9]))
      T.verifyEqual(actR.stdDeviation,           NAN_AR)
      T.verifyEqual(actR.nBadPoints,             [0; 0])
      T.verifyEqual(actR.offset,                 NAN_AR)
      T.verifyEqual(actR.coefficientCos1,        NAN_AR)
      T.verifyEqual(actR.coefficientSin1,        NAN_AR)

      switch(N_FIT_COEFF)
        case 3
          T.verifyTrue(~ismember("coefficientCos2Ar", actR.Properties.VariableNames))
          T.verifyTrue(~ismember("coefficientSin2Ar", actR.Properties.VariableNames))
        case 5
          T.verifyEqual(actR.coefficientCos2,    NAN_AR)
          T.verifyEqual(actR.coefficientSin2,    NAN_AR)
      end
    end



    % Samples can be fit perfectly with nonzero values for all coefficients.
    % Samples in two fit windows.
    function test_spin_samples_perfect_fit(T)
      % PROPOSAL: Split in three separate test functions.
      %   PROPOSAL: Shared function for calling
      %     bepic.spinfit.fit_TAFW() with standard values.
      %   PROPOSAL: Shared function for generating samples.

      N_IN = 64;

      function actR = test(beginTt2000, endTt2000)
        tt2000Ar       = int64(linspace(beginTt2000, endTt2000, N_IN))';
        spinPhaseRadAr = wrapTo2Pi(linspace( pi, 3*pi, N_IN))';
        samplesAr      = 2 ...
          + 3*cos(  spinPhaseRadAr) + 4*sin(  spinPhaseRadAr) ...
          + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

        actR = bepic.spinfit.fit_TAFW( ...
          tt2000Ar       = tt2000Ar, ...
          spinPhaseRadAr = spinPhaseRadAr, ...
          samplesAr      = samplesAr, ...
          fitWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
          fitWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
          fitWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
          nMinFitSamples         = T.N_MIN_FIT_SAMPLES, ...
          nFitCoefficients       = 5);
      end

      % Timestamps to return two fit windows.
      actR = test(101e9, 107e9);

      T.verifyEqual(actR.fitWindowCenterTt2000, int64([102e9; 106e9]))
      T.verifyEqual(actR.offset,          [2; 2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1, [3; 3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1, [4; 4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2, [5; 5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2, [6; 6], AbsTol=1e-13)

      % Timestamps to return only the second fit window.
      actR = test(103e9, 107e9);

      T.verifyEqual(actR.fitWindowCenterTt2000, int64([106e9]))
      T.verifyEqual(actR.offset,          [2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1, [3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1, [4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2, [5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2, [6], AbsTol=1e-13)

      % Timestamps to return only the first fit window.
      actR = test(101e9, 105e9);

      T.verifyEqual(actR.fitWindowCenterTt2000, int64([102e9]))
      T.verifyEqual(actR.offset,          [2], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos1, [3], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1, [4], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientCos2, [5], AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin2, [6], AbsTol=1e-13)
    end



    % Samples can be fit perfectly with nonzero values for all coefficients.
    % Samples in two fit windows.
    function test_spin_samples_too_few(T, N_FIT_COEFF)
      N_IN = 6+6;

      tt2000Ar       = int64(linspace(1e9, 7e9, N_IN))';
      spinPhaseRadAr = wrapTo2Pi(linspace( pi, 5*pi, N_IN))';
      samplesAr      = 2 + 3*cos(  spinPhaseRadAr);

      actR = bepic.spinfit.fit_TAFW( ...
        tt2000Ar       = tt2000Ar, ...
        spinPhaseRadAr = spinPhaseRadAr, ...
        samplesAr      = samplesAr, ...
        fitWindowPeriodNs     = T.TIME_WINDOW_PERIOD_NS, ...
        fitWindowLengthNs     = T.TIME_WINDOW_LENGTH_NS, ...
        fitWindowCenterTt2000 = T.TIME_WINDOW_CENTER_TT2000, ...
        nMinFitSamples         = 5+3, ...    % Required for this test.
        nFitCoefficients       = N_FIT_COEFF);

      NAN_AR = NaN(2, 1);
      T.verifyEqual(actR.offset,          NAN_AR, AbsTol=1e-15)
      T.verifyEqual(actR.coefficientCos1, NAN_AR, AbsTol=1e-13)
      T.verifyEqual(actR.coefficientSin1, NAN_AR, AbsTol=1e-13)
      switch(N_FIT_COEFF)
        case 3
          T.verifyTrue(~ismember("coefficientCos2Ar", actR.Properties.VariableNames))
          T.verifyTrue(~ismember("coefficientSin2Ar", actR.Properties.VariableNames))
        case 5
          T.verifyEqual(actR.coefficientCos2, NAN_AR, AbsTol=1e-13)
          T.verifyEqual(actR.coefficientSin2, NAN_AR, AbsTol=1e-13)
      end

    end



  end    % methods(Test)



end
