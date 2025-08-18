%
% matlab.unittest automatic test code for irf.utils.split_by_false().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2022-08-11
%
classdef split_by_false___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_0(testCase)

      % Arbitrary number output variables.
      function test(bArray, expI1Array, expI2Array)
        bArray     = logical(bArray);

        [actI1Array, actI2Array] = irf.utils.split_by_false(bArray);
        testCase.verifyTrue(iscolumn(actI1Array))
        testCase.verifyTrue(iscolumn(actI2Array))
        testCase.verifyEqual(actI1Array, expI1Array)
        testCase.verifyEqual(actI2Array, expI2Array)
      end

      ECA = ones(0, 1);   % Empty Column Array

      test(ECA, ECA, ECA)

      test(0, ECA, ECA)
      test(1, 1,   1  )

      test([0,0,0]', ECA, ECA)
      test([1,1,1]', 1,   3)

      test([1,0,0]', 1, 1);
      test([0,1,0]', 2, 2);
      test([0,0,1]', 3, 3);

      test([0,1,1]', 2, 3);
      test([1,0,1]', [1 3]', [1 3]');
      test([1,1,0]', 1, 2);

      test([1,1, 0,0,0, 1,1,1,1]', [1,6]', [2,9]');
      test([0,0, 1,1,1, 0,0,0,0]', 3,      5     );
      test([1,1,0,0]', 1, 2);
      test([0,0,1,1]', 3, 4);

      % Complex test
      test([0,0,1,1,1,0,0,1]', [3,8]',    [5,8]');

    end



    function test_exc(testCase)
      function test_exc(varargin)
        testCase.verifyError(...
          @() irf.utils.split_by_false(varargin{:}), ...
          ?MException)
      end

      test_exc(true(0,0))
      test_exc(true(2,2))
      test_exc(true(0,2))   % Too many columns (2).
      test_exc(true(2,0))
    end



  end    % methods(Test)



end
