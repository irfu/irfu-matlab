%
% UNFINISHED
%
% matlab.unittest automatic test code for bepic.spinfit.spin_aligned().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_spin_aligned___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty(T)
      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = int64.empty( 0, 1), ...
        spinPhaseRadAr      = double.empty(0, 1), ...
        samplesAr           = double.empty(0, 1), ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = 0, ...
        dataGapMinNs        = int64(1));

      T.assertEqual(actR.tt2000Ar, int64.empty(0, 1))
    end



    function test_0(T)
      N_IN = 32;

      % Data in two time windows.
      tt2000Ar       = int64(logspace(9, 10, N_IN))';     % Not linearly incrementing!!!
      spinPhaseRadAr = wrapTo2Pi(linspace( pi, 5*pi, N_IN))';
      samplesAr      = 2 + 3*cos(spinPhaseRadAr) + 4*sin(spinPhaseRadAr) + 5*cos(2*spinPhaseRadAr) + 6*sin(2*spinPhaseRadAr);

      actR = bepic.spinfit.spin_aligned( ...
        tt2000Ar            = int64.empty( 0, 1), ...
        spinPhaseRadAr      = double.empty(0, 1), ...
        samplesAr           = double.empty(0, 1), ...
        timeWindowPeriodRad = 2*pi, ...
        timeWindowLengthRad = 2*pi, ...
        timeWindowCenterRad = 0, ...
        dataGapMinNs        = int64(1e10));

      T.assertEqual(actR.tt2000Ar, int64.empty(0, 1))
    end



  end    % methods(Test)



end
