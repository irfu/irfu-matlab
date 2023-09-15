%
% matlab.unittest automatic test code for
% irf.utils.split_by_false().
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



    function test0(testCase)

      % Arbitrary number output variables.
      function test(boolArray, expI1Array, expI2Array)
        % Pre-allocate correct size for later assignment via function
        actOutputs = cell(1, 2);

        [actI1Array, actI2Array] = irf.utils.split_by_false(boolArray);
        testCase.verifyEqual(actI1Array, expI1Array)
        testCase.verifyEqual(actI2Array, expI2Array)

        % Use transposed argument.
        [actOutputs{:}] = irf.utils.split_by_false(boolArray');
        testCase.verifyEqual(actOutputs, {expI1Array, expI2Array})
      end

      function test_exc(varargin)
        testCase.verifyError(...
          @() irf.utils.split_by_false(varargin{:}), ...
          ?MException)
      end

      ECA = ones(0,1);   % Empty Column Array

      test(ECA, ECA, ECA)

      test([0], ECA, ECA)
      test([1], [1], [1])
      test([0], ECA, ECA)

      test([0,0,0], ECA, ECA)
      test([1,1,1], [1], [3])

      test([1,0,0], [1], [1]);
      test([0,1,0], [2], [2]);
      test([0,0,1], [3], [3]);

      test([1,1,0,0,1,1], [1,5]', [2,6]');
      test([0,0,1,1,0,0], [3],    [4]);
      test([1,1,0,0],     [1],    [2]);
      test([0,0,1,1],     [3],    [4]);

      test([0,0,1,1,1,0,0,1], [3,8]',    [5,8]');

      test_exc(ones(2,2))
      test_exc(ones(0,2))
      test_exc(ones(2,0))
    end



  end    % methods(Test)



end
