%
% matlab.unittest automatic test code for
% bepic.spinfit.cumulative_spin_phase_to_TT2000().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef spinfit_cumulative_spin_phase_to_TT2000___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_cumulative_spin_phase_to_TT2000_all_empty(T)
      T.test_cumulative_spin_phase_to_TT2000(...
        int64.empty( 0, 1), ...
        double.empty(0, 1), ...
        double.empty(0, 1), ...
        int64.empty(0, 1))
    end



    function test_cumulative_spin_phase_to_TT2000_empty_new_points(T)
      T.test_cumulative_spin_phase_to_TT2000(...
        int64( [10 20 30]'), ...
        double([ 1  2  3]'), ...
        double.empty(0, 1), ...
        int64.empty( 0, 1))
    end



    function test_cumulative_spin_phase_to_TT2000_0(T)
      T.test_cumulative_spin_phase_to_TT2000(...
        int64( [10 20 30]'), ...
        double([ 1  2  3]'), ...
        double([ 1.5 2.5]'), ...
        int64( [  15  25]'))
    end



    function test_cumulative_spin_phase_to_TT2000_1(T)
      % Big jump in TT2000 and cumulative spin phase.
      T.test_cumulative_spin_phase_to_TT2000(...
        int64( [10 20 300]'), ...
        double([ 1  2  30]'), ...
        double([ 1.5  16]'), ...
        int64( [  15 160]'))
    end



    function test_cumulative_spin_phase_to_TT2000_extrapolate(T)
      % Big jump in TT2000 and cumulative spin phase.
      T.test_cumulative_spin_phase_to_TT2000(...
        int64( [10 20 40]'), ...
        double([ 1  2  3]'), ...
        double([ 0  4]'), ...
        int64( [ 0 60]'))
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Internal helper function to shorten some tests.
    %
    function test_cumulative_spin_phase_to_TT2000(T, ...
        dataTt2000Ar, dataCumulSpinPhaseRadAr, inCumulSpinPhaseRadAr, ...
        expTt2000Ar)

      actTt2000Ar = bepic.spinfit.cumulative_spin_phase_to_TT2000( ...
        dataTt2000Ar, dataCumulSpinPhaseRadAr, inCumulSpinPhaseRadAr);

      T.assertEqual(actTt2000Ar, expTt2000Ar)
    end



  end    % methods(Access=private)



end
