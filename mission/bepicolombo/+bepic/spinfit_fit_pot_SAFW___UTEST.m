%
% matlab.unittest automatic test code for
% bepic.spinfit.fit_pot_SAFW().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_fit_pot_SAFW___UTEST < matlab.unittest.TestCase
  % PROPOSAL: One test for testing the fit. Keep it simple: One FTWN etc.



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    N_COEFF_OFFICIAL = bepic.spinfit.N_COEFF_OFFICIAL_POT_FIT
    N_COEFF_IRFU     = bepic.spinfit.N_COEFF_IRFU_POT_FIT
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
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF_OFFICIAL);

      T.verifyEqual(actT.fitWindowCenterTt2000, int64.empty( 0, 1))
      T.verifyEqual(actT.c,                     double.empty(0, T.N_COEFF_OFFICIAL, 2))
      T.verifyEqual(actT.cIrfu,                 double.empty(0, T.N_COEFF_IRFU,     2))
    end



    function test_one_timestamp(T)
      actT = bepic.spinfit.fit_pot_SAFW( ...
        tt2000Ar                = int64([10]'), ...
        spinPhaseRadAr          = [0.5]' * 2*pi, ...
        samplesAr               = [3, 4], ...
        dataGapMinNs            = int64(1000), ...
        nMinSamplesPerFitWindow = T.N_COEFF_OFFICIAL);

      T.verifyEqual(actT.fitWindowCenterTt2000, int64.empty( 0, 1))
      T.verifyEqual(actT.c,                     double.empty(0, T.N_COEFF_OFFICIAL, 2))
      T.verifyEqual(actT.cIrfu,                 double.empty(0, T.N_COEFF_IRFU,     2))
    end



    % Test the content of the fits.
    %
    % NOTE: Exists for historical reasons. Is maybe overkill now (hard to
    % maintain).
    % function test_2xFTWN(T)
    %   samples1Ar = [1 2 1 2 1 2   3 4 3 4 3 4]';
    %   samples2Ar = [6 7 6 7 6 7   8 9 8 9 8 9]';
    %   actT = bepic.spinfit.fit_pot_SAFW( ...
    %     tt2000Ar                = int64(    [ 01  02  03  04  06  07    11  12  13  14  16  17]'), ...
    %     spinPhaseRadAr          = wrapTo2Pi([0.1 0.2 0.3 0.4 0.6 0.7   1.1 1.2 1.3 1.4 1.6 1.7]' * 2*pi), ...
    %     samplesAr               =           [samples1Ar, samples2Ar], ...
    %     dataGapMinNs            = int64(1000), ...
    %     nMinSamplesPerFitWindow = T.N_COEFF_OFFICIAL);
    %
    %   T.verifyEqual(actT.fitWindowCenterTt2000, int64([  5;  15]))
    %   T.verifySize( actT.A,                     [2, 3, 2])
    %   T.verifyEqual(actT.A(:, 3, 1),            [1.5; 3.5], "RelTol", 1e-15)
    %   T.verifyEqual(actT.A(:, 3, 2),            [6.5; 8.5])
    % end



  end    % methods(Test)



end
