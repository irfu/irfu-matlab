%
% matlab.unittest automatic test code for
% bepic.spinfit.ftwn.get_SAFWs_when_no_data_gaps().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ftwn_get_SAFWs_when_no_data_gaps___UTEST < matlab.unittest.TestCase
  % PROBLEM: Difficult to manually derive the return value timestamps.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_zero_samples(T)
      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_when_no_data_gaps( ...
        tt2000Ar       = int64.empty( 0, 1), ...
        spinPhaseRadAr = double.empty(0, 1), ...
        ftwnPeriodRad  = 2*pi, ...
        ftwnLengthRad  = 2*pi, ...
        ftwnBeginRad   = 0*pi);

      ExpFtwnTable = table();
      ExpFtwnTable.beginTt2000  = int64.empty( 0, 1);
      ExpFtwnTable.endTt2000    = int64.empty( 0, 1);
      ExpFtwnTable.centerTt2000 = int64.empty( 0, 1);

      T.verifyEqual(ActFtwnTable, ExpFtwnTable)
    end



    function test_empty_one_sample(T)
      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_when_no_data_gaps( ...
        tt2000Ar       = int64([30]), ...
        spinPhaseRadAr = [3], ...
        ftwnPeriodRad  = 2*pi, ...
        ftwnLengthRad  = 2*pi, ...
        ftwnBeginRad   = 0*pi);

      ExpFtwnTable = table();
      ExpFtwnTable.beginTt2000  = int64.empty( 0, 1);
      ExpFtwnTable.endTt2000    = int64.empty( 0, 1);
      ExpFtwnTable.centerTt2000 = int64.empty( 0, 1);

      T.verifyEqual(ActFtwnTable, ExpFtwnTable)
    end



    function test_complex(T)
      % Non-linear time-to-spin phase function.

      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_when_no_data_gaps( ...
        tt2000Ar       = int64(    [  20   30   35   45   50]'), ...
        spinPhaseRadAr = wrapTo2Pi([ 0.5, 1.0, 1.5, 2.0, 2.5]' * 2*pi), ...
        ftwnPeriodRad  = 1  *2*pi, ...
        ftwnLengthRad  = 2  *2*pi, ...
        ftwnBeginRad   = 0.5*2*pi);

      ExpFtwnTable = table();
      ExpFtwnTable.beginTt2000  = int64([20,  35,  50]');
      ExpFtwnTable.endTt2000    = int64([50,  60,  70]');
      ExpFtwnTable.centerTt2000 = int64([35,  50,  60]');   % Empirical.

      T.verifyEqual(ActFtwnTable, ExpFtwnTable)
    end



  end    % methods(Test)



end
