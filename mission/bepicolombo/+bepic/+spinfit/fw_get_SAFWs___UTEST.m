%
% matlab.unittest automatic test code for bepic.spinfit.fw.get_SAFWs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef fw_get_SAFWs___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_zero_samples(T)
      ActFitWindowTable = bepic.spinfit.fw.get_SAFWs( ...
        tt2000Ar             = int64.empty( 0, 1), ...
        spinPhaseRadAr       = double.empty(0, 1), ...
        fitWindowPeriodRad   = 2*pi, ...
        fitWindowLengthRad   = 2*pi, ...
        fitWindowBeginRefRad = 0*pi, ...
        dataGapMinNs         = int64(100));

      ExpFitWindowTable = table();
      ExpFitWindowTable.beginTt2000  = int64.empty(0, 1);
      ExpFitWindowTable.endTt2000    = int64.empty(0, 1);
      ExpFitWindowTable.centerTt2000 = int64.empty(0, 1);

      T.verifyEqual(ActFitWindowTable, ExpFitWindowTable)
    end



    function test_empty_one_sample(T)
      ActFitWindowTable = bepic.spinfit.fw.get_SAFWs( ...
        tt2000Ar             = int64([30]), ...
        spinPhaseRadAr       = [3], ...
        fitWindowPeriodRad   = 2*pi, ...
        fitWindowLengthRad   = 2*pi, ...
        fitWindowBeginRefRad = 0*pi, ...
        dataGapMinNs         = int64(100));

      ExpFitWindowTable = table();
      ExpFitWindowTable.beginTt2000  = int64.empty(0, 1);
      ExpFitWindowTable.endTt2000    = int64.empty(0, 1);
      ExpFitWindowTable.centerTt2000 = int64.empty(0, 1);

      T.verifyEqual(ActFitWindowTable, ExpFitWindowTable)
    end



    % Check against producing duplicated fit windows.
    function test_no_duplicates(T)
      ActFitWindowTable = bepic.spinfit.fw.get_SAFWs( ...
        tt2000Ar             = int64(    [  121, 126, 131, 139, 144]'), ...
        spinPhaseRadAr       = wrapTo2Pi([  2.1, 2.6, 3.1, 3.9, 4.4]' * 2*pi), ...
        fitWindowPeriodRad   = 1  *2*pi, ...
        fitWindowLengthRad   = 1  *2*pi, ...
        fitWindowBeginRefRad = 0.0*2*pi, ...
        dataGapMinNs         = int64(7));

      ExpFitWindowTable = table();
      ExpFitWindowTable.beginTt2000  = int64([120, 130, 140]');
      ExpFitWindowTable.endTt2000    = int64([130, 140, 150]');
      ExpFitWindowTable.centerTt2000 = int64([125, 135, 145]');

      T.verifyEqual(ActFitWindowTable, ExpFitWindowTable)
    end



    % Non-linear time-to-spin phase function.
    % Data gap
    function test_complex(T)
      ActFitWindowTable = bepic.spinfit.fw.get_SAFWs( ...
        tt2000Ar             = int64(    [   20   30  200  210]'), ...
        spinPhaseRadAr       = wrapTo2Pi([ 0.75, 1.0, 2.0, 2.5]' * 2*pi), ...
        fitWindowPeriodRad   = 1  *2*pi, ...
        fitWindowLengthRad   = 2  *2*pi, ...
        fitWindowBeginRefRad = 0.5*2*pi, ...
        dataGapMinNs         = int64(100));

      ExpFitWindowTable = table();
      ExpFitWindowTable.beginTt2000  = int64([ 10, 190, 210]');
      ExpFitWindowTable.endTt2000    = int64([ 90, 230, 250]');
      ExpFitWindowTable.centerTt2000 = int64([ 50, 210, 230]');

      T.verifyEqual(ActFitWindowTable, ExpFitWindowTable)
    end



  end    % methods(Test)



end
