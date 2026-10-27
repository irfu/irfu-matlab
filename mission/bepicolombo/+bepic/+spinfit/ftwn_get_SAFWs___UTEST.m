%
% matlab.unittest automatic test code for bepic.spinfit.ftwn.get_SAFWs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ftwn_get_SAFWs___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty_zero_samples(T)
      ActFtwn = bepic.spinfit.ftwn.get_SAFWs( ...
        tt2000Ar        = int64.empty( 0, 1), ...
        spinPhaseRadAr  = double.empty(0, 1), ...
        ftwnPeriodRad   = 2*pi, ...
        ftwnLengthRad   = 2*pi, ...
        ftwnBeginRefRad = 0*pi, ...
        dataGapMinNs    = int64(100));

      ExpFtwn = table();
      ExpFtwn.beginTt2000  = int64.empty(0, 1);
      ExpFtwn.endTt2000    = int64.empty(0, 1);
      ExpFtwn.centerTt2000 = int64.empty(0, 1);

      T.verifyEqual(ActFtwn, ExpFtwn)
    end



    function test_empty_one_sample(T)
      ActFtwn = bepic.spinfit.ftwn.get_SAFWs( ...
        tt2000Ar        = int64([30]), ...
        spinPhaseRadAr  = [3], ...
        ftwnPeriodRad   = 2*pi, ...
        ftwnLengthRad   = 2*pi, ...
        ftwnBeginRefRad = 0*pi, ...
        dataGapMinNs    = int64(100));

      ExpFtwn = table();
      ExpFtwn.beginTt2000  = int64.empty(0, 1);
      ExpFtwn.endTt2000    = int64.empty(0, 1);
      ExpFtwn.centerTt2000 = int64.empty(0, 1);

      T.verifyEqual(ActFtwn, ExpFtwn)
    end



    % Check against producing duplicated fit windows.
    function test_no_duplicates(T)
      ActFtwn = bepic.spinfit.ftwn.get_SAFWs( ...
        tt2000Ar        = int64(    [  121, 126, 131, 139, 144]'), ...
        spinPhaseRadAr  = wrapTo2Pi([  2.1, 2.6, 3.1, 3.9, 4.4]' * 2*pi), ...
        ftwnPeriodRad   = 1  *2*pi, ...
        ftwnLengthRad   = 1  *2*pi, ...
        ftwnBeginRefRad = 0.0*2*pi, ...
        dataGapMinNs    = int64(7));

      ExpFtwn = table();
      ExpFtwn.beginTt2000  = int64([120, 130, 140]');
      ExpFtwn.endTt2000    = int64([130, 140, 150]');
      ExpFtwn.centerTt2000 = int64([125, 135, 145]');

      T.verifyEqual(ActFtwn, ExpFtwn)
    end



    % Non-linear time-to-spin phase function.
    % Data gap
    function test_complex(T)
      ActFtwn = bepic.spinfit.ftwn.get_SAFWs( ...
        tt2000Ar        = int64(    [   20   30  200  210]'), ...
        spinPhaseRadAr  = wrapTo2Pi([ 0.75, 1.0, 2.0, 2.5]' * 2*pi), ...
        ftwnPeriodRad   = 1  *2*pi, ...
        ftwnLengthRad   = 2  *2*pi, ...
        ftwnBeginRefRad = 0.5*2*pi, ...
        dataGapMinNs    = int64(100));

      ExpFtwn = table();
      ExpFtwn.beginTt2000  = int64([ 10, 190, 210]');
      ExpFtwn.endTt2000    = int64([ 90, 230, 250]');
      ExpFtwn.centerTt2000 = int64([ 50, 210, 230]');

      T.verifyEqual(ActFtwn, ExpFtwn)
    end



  end    % methods(Test)



end
