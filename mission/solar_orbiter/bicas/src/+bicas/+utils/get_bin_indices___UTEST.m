%
% matlab.unittest automatic test code for bicas.utils.get_bin_indices().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-10 from older test code.
%
classdef get_bin_indices___UTEST < matlab.unittest.TestCase



  properties(TestParameter)
    % NOTE: Thresholds chosen to overlap with lengths of
    N_BB_THRESHOLD = {3,4,5,6,10};
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase, N_BB_THRESHOLD)
      % TODO: Data in zero-length bins.
      % TODO: Vary N_BB_MAX



      % NOTE: Implementation might use internal hardcoded constant.
      % Automatic tests must include large enough input data to trigger
      % alternate implementation.

      % NOTE: Changes rows to columns to make typing arguments easier.
      function test(xRow, xBoundariesRow, iInBinCaRowRow)

        % Transpose components in iInBinCaRowRow.
        iInBinCa = cellfun(@(c) (c'), iInBinCaRowRow, 'UniformOutput', false)';

        actOutput = bicas.utils.get_bin_indices(...
          xRow', xBoundariesRow', N_BB_THRESHOLD);
        testCase.verifyEqual(actOutput, iInBinCa)
      end

      function test_exc(varargin)
        testCase.verifyError(...
          @() bicas.utils.get_bin_indices(varargin{:}), ...
          ?MException)
      end

      %===================================================================

      A1x0 = zeros(1,0);
      C1x0 = cell( 1,0);

      % Edge cases: No data, or no boundaries
      test(A1x0,    A1x0,   C1x0)
      test([3],     A1x0,   C1x0)
      test(A1x0,    [5],    C1x0)
      test(A1x0,    [5,10], {A1x0})
      test([3,13],  [5,10], {A1x0})

      test([6,7,8], [5,10], {[1,2,3]})

      test([3,6,7,8,11],  [5,10],      {[2,3,4]})
      test([3,6,7,8,11],  [0,5,10,15], {[1], [2,3,4], [5]})

      % Data on bin boundaries.
      test([5,10], [0,5,10,15], {A1x0, [1], [2]})

      % Data in zero-length bins (on boundaries).
      test([5,10], [0, 5,5, 10,10, 15], {A1x0, A1x0, [1], A1x0, [2]})

      % Zero-length bin boundaries.
      test([3,6,7,8,11],  [0,5,5,10,10,15], {[1], A1x0, [2,3,4], A1x0, [5]})

      % Data not sorted.
      test_exc([11,6,3,8,7],  [0,5,10,15])
    end



  end    % methods(Test)



end
