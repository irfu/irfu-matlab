%
% matlab.unittest automatic test code for
% bepic.spinfit.utils.get_segment_SAFWs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils_get_segment_SAFWs___UTEST < matlab.unittest.TestCase
  % PROBLEM: Difficult to manually derive the return value timestamps.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_zero_samples(T)
      ActFitWindowTable = bepic.spinfit.utils.get_segment_SAFWs( ...
        tt2000Ar           = int64.empty( 0, 1), ...
        spinPhaseRadAr     = double.empty(0, 1), ...
        fitWindowPeriodRad = 2*pi, ...
        fitWindowLengthRad = 2*pi, ...
        fitWindowBeginRad  = 0*pi);

      ExpFitWindowTable = table();
      ExpFitWindowTable.beginTt2000 = int64.empty( 0, 1);
      ExpFitWindowTable.endTt2000   = int64.empty( 0, 1);

      T.verifyEqual(ActFitWindowTable, ExpFitWindowTable)
    end



    function test_empty_one_sample(T)
      ActFitWindowTable = bepic.spinfit.utils.get_segment_SAFWs( ...
        tt2000Ar           = int64([30]), ...
        spinPhaseRadAr     = [3], ...
        fitWindowPeriodRad = 2*pi, ...
        fitWindowLengthRad = 2*pi, ...
        fitWindowBeginRad  = 0*pi);

      ExpFitWindowTable = table();
      ExpFitWindowTable.beginTt2000 = int64.empty( 0, 1);
      ExpFitWindowTable.endTt2000   = int64.empty( 0, 1);

      T.verifyEqual(ActFitWindowTable, ExpFitWindowTable)
    end



    function test_complex(T)
      % Non-linear time-to-spin phase function.

      ActFitWindowTable = bepic.spinfit.utils.get_segment_SAFWs( ...
        tt2000Ar           = int64(    [  20   30   35   45   50]'), ...
        spinPhaseRadAr     = wrapTo2Pi([ 0.5, 1.0, 1.5, 2.0, 2.5]' * 2*pi), ...
        fitWindowPeriodRad = 1  *2*pi, ...
        fitWindowLengthRad = 2  *2*pi, ...
        fitWindowBeginRad  = 0.5*2*pi);

      %                                                    [0.5, 1.5, 2.5]' * 2*pi);
      %                                                    [2.5, 3.5, 4.5]' * 2*pi);
      T.verifyEqual(ActFitWindowTable.beginTt2000, int64([20,  35,  50]'));
      T.verifyEqual(ActFitWindowTable.endTt2000,   int64([50,  60,  70]'));
    end



  end    % methods(Test)



end
