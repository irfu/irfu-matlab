%
% matlab.unittest automatic test code for
% bepic.spinfit.fw_fit_window_time_to_indices().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef fw_fit_window_time_to_indices___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_0_samples_0_fit_windows(T)
      actICa = bepic.spinfit.fw.fit_window_time_to_indices(...
        int64.empty(0, 1), ...
        int64.empty(0, 1), ...
        int64.empty(0, 1));

      T.assertEqual(actICa, cell.empty(0, 1))
    end



    function test_0_samples_2_fit_windows(T)
      actICa = bepic.spinfit.fw.fit_window_time_to_indices(...
        int64.empty(0, 1), ...
        int64([10; 20]), ...
        int64([20; 30]));

      T.assertEqual(actICa, {double.empty(0, 1); double.empty(0, 1)})
    end



    function test_1_sample_2_fit_windows(T)
      actICa = bepic.spinfit.fw.fit_window_time_to_indices(...
        int64([13]), ...
        int64([10; 20]), ...
        int64([20; 30]));

      T.assertEqual(actICa, {[1:1]'; double.empty(0, 1)})
    end



    function test_2_samples_2_fit_windows(T)
      actICa = bepic.spinfit.fw.fit_window_time_to_indices(...
        int64([13; 23]), ...
        int64([10; 20]), ...
        int64([20; 30]));

      T.assertEqual(actICa, {[1:1]'; [2:2']})
    end



    function test_7_samples_0_fit_windows(T)
      actICa = bepic.spinfit.fw.fit_window_time_to_indices(...
        int64([3; 8; 13; 18; 23; 28; 33]), ...
        int64.empty(0, 1), ...
        int64.empty(0, 1));

      T.assertEqual(actICa, cell.empty(0, 1))
    end



    function test_samples_not_on_boundaries(T)
      actICa = bepic.spinfit.fw.fit_window_time_to_indices(...
        int64([3; 8; 13; 18; 23; 28; 33]), ...
        int64([10; 20]), ...
        int64([20; 30]));

      T.assertEqual(actICa, {[3:4]'; [5:6]'})
    end



    function test_samples_on_boundaries(T)
      actICa = bepic.spinfit.fw.fit_window_time_to_indices(...
        int64([3; 10; 13; 20; 23; 30; 33]), ...
        int64([10; 20]), ...
        int64([20; 30]));

      T.assertEqual(actICa, {[2:3]'; [4:5]'})
    end



  end    % methods(Test)



end
