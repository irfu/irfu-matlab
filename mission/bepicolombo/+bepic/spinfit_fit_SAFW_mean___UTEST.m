%
% matlab.unittest automatic test code for
% bepic.spinfit.fit_SAFW_mean().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_fit_SAFW_mean___UTEST < matlab.unittest.TestCase
  % PROPOSAL: More tests.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_samples(T)
      actT = bepic.spinfit.fit_SAFW_mean( ...
        tt2000Ar              = int64.empty(0, 1), ...
        spinPhaseRadAr        = double.empty(0, 1), ...
        samplesAr             = double.empty(0, 1), ...
        fitWindowPeriodRad    = 2*pi, ...
        fitWindowLengthRad    = 2*pi, ...
        fitWindowCenterRefRad = pi, ...
        dataGapMinNs          = int64(1000));

      T.verifyEqual(actT.fitWindowCenterTt2000, int64.empty( 0, 1))
      T.verifyEqual(actT.mean,                  double.empty(0, 1))
    end



    function test_one_sample(T)
      actT = bepic.spinfit.fit_SAFW_mean( ...
        tt2000Ar              = int64([10]'), ...
        spinPhaseRadAr        = [0.5]' * 2*pi, ...
        samplesAr             = [  3]', ...
        fitWindowPeriodRad    = 2*pi, ...
        fitWindowLengthRad    = 2*pi, ...
        fitWindowCenterRefRad = pi, ...
        dataGapMinNs          = int64(1000));

      T.verifyEqual(actT.fitWindowCenterTt2000, int64.empty( 0, 1))
      T.verifyEqual(actT.mean,                  double.empty(0, 1))
    end




    % No data gap.
    function test_0(T)
      actT = bepic.spinfit.fit_SAFW_mean( ...
        tt2000Ar              = int64(    [ 01   06   11   16]'), ...
        spinPhaseRadAr        = wrapTo2Pi([0.1, 0.6, 0.1, 0.6]' * 2*pi), ...
        samplesAr             =           [  1,   2,   3,   4]', ...
        fitWindowPeriodRad    = 2*pi, ...
        fitWindowLengthRad    = 2*pi, ...
        fitWindowCenterRefRad = 2*pi*0.5, ...
        dataGapMinNs          = int64(1000));

      T.verifyEqual(actT.fitWindowCenterTt2000, int64([  5,  15]'))
      T.verifyEqual(actT.mean,                        [1.5, 3.5]')
    end



  end    % methods(Test)



end
