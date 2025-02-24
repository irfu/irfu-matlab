%
% matlab.unittest automatic test code for solo.hwzv.get_LRX().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_LRX___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_exception(testCase)
      function test_exc(zvR0, zvR1, zvR2, iLsf)
          testCase.assertError(...
              @() solo.hwzv.get_LRX(zvR0, zvR1, zvR2, iLsf), ...
              ?MException)
      end

      % ==================
      % Illegal LSF values
      % ==================
      test_exc([0], [1], [0], [0])
      test_exc([0], [1], [0], [5])

      test_exc([0],   [1],   [0],   [NaN])
      test_exc([0;1], [1;0], [0;1], [1;0])

      % ===================
      % Illegal array sizes
      % ===================
      test_exc([],    [],    [],    []   )
      test_exc([0,1], [1,0], [0,1], [1,2])

      % ==========
      % Output NaN
      % ==========
      test_exc(...
        [nan]', ...
        [0  ]', ...
        [0  ]', ...
        [1  ]')
      test_exc(...
        [0  ]', ...
        [nan]', ...
        [0  ]', ...
        [2  ]')
      test_exc(...
        [0  ]', ...
        [0  ]', ...
        [nan]', ...
        [3  ]')
    end



    function test0(testCase)

      function test(zvR0, zvR1, zvR2, iLsf, expZvLrx)
          actZvLrx = solo.hwzv.get_LRX(zvR0, zvR1, zvR2, iLsf);

          % NOTE: The function return value class depends on the input
          %       arguments. It is NOT guaranteed to be logical.
          % testCase.assertTrue(islogical(actZvLrx))

          testCase.assertEqual(actZvLrx, expZvLrx)
      end

      test(zeros(0, 1), zeros(0, 1), zeros(0,1), zeros(0,1), zeros(0, 1))

      test(...
        [1 0 0 0   0 1 1 1]', ...
        [0 1 0 0   1 0 1 1]', ...
        [0 0 1 0   1 1 0 1]', ...
        [1 2 3 4   1 2 3 4]', ...
        [1 1 1 1   0 0 0 1]')

      % ========
      % Test NaN
      % ========
      % Tests de facto behaviour. Function behaviour should possibly be changed
      % to forbid NaN, and then these tests need to be changed.
      test(...
        [nan 0   0   0]', ...
        [0   nan 0   0]', ...
        [0   0   nan 0]', ...
        [2   3   4   1]', ...
        [0   0   1   0]')

      test(...
        [0   nan nan nan]', ...
        [nan 0   nan nan]', ...
        [nan nan 0   nan]', ...
        [1   2   3   4  ]', ...
        [0   0   0   1  ]')
    end



  end    % methods(Test)



end
