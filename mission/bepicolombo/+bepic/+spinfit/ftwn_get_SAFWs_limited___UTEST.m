%
% matlab.unittest automatic test code for
% bepic.spinfit.ftwn.get_SAFW_limited().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ftwn_get_SAFWs_limited___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_0xSamples(T)
      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_limited( ...
        tt2000Ar           = int64.empty(0, 1), ...
        spinPhaseRadAr     = double.empty(0, 1), ...
        dataGapMinNs       = int64(2e9));

      ExpFtwnTable = table();
      ExpFtwnTable.iBegin       = double.empty(0, 1);
      ExpFtwnTable.iEnd         = double.empty(0, 1);

      T.assertEqual(ActFtwnTable, ExpFtwnTable)
    end



    function test_1xFTWN_nominal(T)
      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_limited( ...
        tt2000Ar           = int64([3.00; 3.25; 3.50; 3.75] * 1e9), ...
        spinPhaseRadAr     =       [0.00; 0.25; 0.50; 0.75] * 2*pi, ...
        dataGapMinNs       = int64(2e9));

      ExpFtwnTable = table();
      ExpFtwnTable.iBegin       = 1;
      ExpFtwnTable.iEnd         = 4;

      T.assertEqual(ActFtwnTable, ExpFtwnTable)
    end



    function test_2xFTWN_nominal(T)
      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_limited( ...
        tt2000Ar           = int64([3.25; 3.50; 3.75;   4.00; 4.60] * 1e9), ...
        spinPhaseRadAr     =       [0.25; 0.50; 0.75;   0.00; 0.60] * 2*pi, ...
        dataGapMinNs       = int64(2e9));

      ExpFtwnTable = table();
      ExpFtwnTable.iBegin       = [1; 4];
      ExpFtwnTable.iEnd         = [3; 5];

      T.assertEqual(ActFtwnTable, ExpFtwnTable)
    end



    function test_1xSample__minSamplesPerFtwn(T)
      % NOTE: Having spins/(rejected) FTWNs with one and two samples
      % respectively in the test data. Must still have spin wraps (decreasing
      % spin) so that the algorithm could even in principle (by mistake)
      % detect one sample as an FTWN, which is somewhat artificial but needed
      % for the test to be meaningful.

      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_limited( ...
        tt2000Ar           = int64([4.20] * 1e9), ...
        spinPhaseRadAr     =       [0.20] * 2*pi, ...
        dataGapMinNs       = int64(2e9));

      ExpFtwnTable = table();
      ExpFtwnTable.iBegin       = double.empty(0, 1);
      ExpFtwnTable.iEnd         = double.empty(0, 1);

      T.assertEqual(ActFtwnTable, ExpFtwnTable)
    end



    function test_3xFTWN__minSamplesPerFtwn(T)
      % NOTE: Having spins/(rejected) FTWNs with one and two samples
      % respectively in the test data. Must still have spin wraps (decreasing
      % spin) so that the algorithm could even in principle (by mistake)
      % detect one sample as an FTWN, which is somewhat artificial but needed
      % for the test to be meaningful.

      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_limited( ...
        tt2000Ar           = int64([3.25; 3.50; 3.75;   4.20;   5.00; 5.30;   6.00; 6.20; 6.40; 6.60; ] * 1e9), ...
        spinPhaseRadAr     =       [0.25; 0.50; 0.75;   0.20;   0.00; 0.30;   0.00; 0.20; 0.40; 0.60; ] * 2*pi, ...
        dataGapMinNs       = int64(2e9));

      ExpFtwnTable = table();
      ExpFtwnTable.iBegin       = [1;  5;  7];
      ExpFtwnTable.iEnd         = [3;  6; 10];

      T.assertEqual(ActFtwnTable, ExpFtwnTable)
    end



  end
  properties(Constant)
    % The spin phase increments but the data gap should (not) split the
    % samples into separate FTWNs anyway.
    TT2000_AR         = int64([3.10; 3.20;   5.80; 5.90;   6.00; 6.20; 6.40; ] * 1e9), ...
    SPIN_PHASE_RAD_AR =       [0.10; 0.20;   0.80; 0.90;   0.00; 0.20; 0.40; ] * 2*pi, ...
  end
  methods(Test)



    % Lower dataGapMinNs. ==> One identified data gap
    function test_3xFTWN__dataGap_0(T)
      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_limited( ...
        tt2000Ar           = T.TT2000_AR, ...
        spinPhaseRadAr     = T.SPIN_PHASE_RAD_AR, ...
        dataGapMinNs       = int64(2e9));    % Lower data gap threshold.

      ExpFtwnTable = table();
      ExpFtwnTable.iBegin       = [1;  3;  5];
      ExpFtwnTable.iEnd         = [2;  4;  7];
      ExpFtwnTable = ExpFtwnTable(3, :);    % Remove old test data.

      T.assertEqual(ActFtwnTable, ExpFtwnTable)
    end



    % Higher dataGapMinNs. ==> No identified data gap.
    function test_3xFTWN__dataGap_1(T)
      ActFtwnTable = bepic.spinfit.ftwn.get_SAFWs_limited( ...
        tt2000Ar           = T.TT2000_AR, ...
        spinPhaseRadAr     = T.SPIN_PHASE_RAD_AR, ...
        dataGapMinNs       = int64(3e9));

      ExpFtwnTable = table();
      ExpFtwnTable.iBegin       = [1;  5];
      ExpFtwnTable.iEnd         = [4;  7];

      T.assertEqual(ActFtwnTable, ExpFtwnTable)
    end



  end    % methods(Test)



end
