%
% matlab.unittest automatic test code for
% bepic.spinfit.fit_pot_SAFW().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_fit_pot_SAFW___UTEST < matlab.unittest.TestCase



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    N_COEFF = 2
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_timestamps(T)
      actT = bepic.spinfit.fit_pot_SAFW( ...
        tt2000Ar                = int64.empty(0, 1), ...
        spinPhaseRadAr          = double.empty(0, 1), ...
        samplesAr               = double.empty(0, 2), ...
        fitWindowPeriodRad      = 2*pi, ...
        fitWindowLengthRad      = 2*pi, ...
        fitWindowCenterRefRad   = pi, ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF);

      T.verifyEqual(actT.fitWindowCenterTt2000, int64.empty( 0, 1))
      T.verifyEqual(actT.A,                     double.empty(0, 3, 2))
    end



    function test_one_timestamp(T)
      actT = bepic.spinfit.fit_pot_SAFW( ...
        tt2000Ar                = int64([10]'), ...
        spinPhaseRadAr          = [0.5]' * 2*pi, ...
        samplesAr               = [3, 4], ...
        fitWindowPeriodRad      = 2*pi, ...
        fitWindowLengthRad      = 2*pi, ...
        fitWindowCenterRefRad   = pi, ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF);

      T.verifyEqual(actT.fitWindowCenterTt2000, int64.empty( 0, 1))
      T.verifyEqual(actT.A,                     double.empty(0, 3, 2))
    end



    % Exists for historical reasons. Is maybe overkill now.
    function test_2_fit_windows(T)
      samples1Ar = [1 2 1 2 3 4 3 4]';
      samples2Ar = [6 7 6 7 8 9 8 9]';
      actT = bepic.spinfit.fit_pot_SAFW( ...
        tt2000Ar                = int64(    [ 01  02  06  07    11  12  16  17]'), ...
        spinPhaseRadAr          = wrapTo2Pi([0.1 0.2 0.6 0.7   1.1 1.2 1.6 1.7]' * 2*pi), ...
        samplesAr               =           [samples1Ar, samples2Ar], ...
        fitWindowPeriodRad      = 2*pi*1, ...
        fitWindowLengthRad      = 2*pi*1, ...
        fitWindowCenterRefRad   = 2*pi*0.5, ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF);

      T.verifyEqual(actT.fitWindowCenterTt2000, int64([  5;  15]))
      T.verifySize( actT.A,                     [2, 3, 2])
      T.verifyEqual(actT.A(:, 3, 1),            [1.5; 3.5], "RelTol", 1e-15)
      T.verifyEqual(actT.A(:, 3, 2),            [6.5; 8.5])
    end



  end    % methods(Test)



end
