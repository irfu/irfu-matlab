%
% matlab.unittest automatic test code for bepic.spinfit.utils.split_by_data_gap().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils_split_by_data_gap___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty(T)
      T.test( ...
        int64.empty( 0, 1), ...
        3, ...
        double.empty(0, 1), ...
        double.empty(0, 1))
    end



    function test_scalar(T)
      T.test( ...
        int64([3]), ...
        4, ...
        [1], ...
        [1])
    end



    % Incrementing jumps in data.
    function test_0(T)
      % NOTE: Indirectly tests for gaps one shorter than specified min data gap.
      T.test( ...
        int64([10 11 13 16 20 25 31]'), ...
        4, ...
        [1 5 6 7]', ...
        [4 5 6 7]')
    end



    % Data gaps at almost ends.
    function test_1(T)
      T.test( ...
        int64([10 16 18 20 22 24 30]'), ...
        6, ...
        [1 2 7]', ...
        [1 6 7]')
    end



    function test_complex(T)
      T.test( ...
        int64([10 11 13 14 20 22 23 30 31]'), ...
        5, ...
        [1 5 8]', ...
        [4 7 9]')
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



      function test(...
          T, tt2000Ar, dataGapNs, expIBeginAr, expIEndAr)
        [actIBeginAr, actIEndAr, actNSegments] = ...
          bepic.spinfit.utils.split_by_data_gap(tt2000Ar, dataGapNs);

        T.assertEqual(actIBeginAr,  expIBeginAr)
        T.assertEqual(actIEndAr,    expIEndAr)
        T.assertEqual(actNSegments, numel(actIBeginAr))
      end



  end    % methods(Access=private)



end
