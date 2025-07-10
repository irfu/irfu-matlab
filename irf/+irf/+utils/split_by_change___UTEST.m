%
% matlab.unittest automatic test code for irf.utils.split_by_change().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef split_by_change___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_1(testCase)

      testCase.test_exc({})

      % Unequal number of rows.
      testCase.test_exc({zeros(0, 1), [1]})
      testCase.test_exc({[2;3], [1]})
      testCase.test_exc({[2,3], [1;2]})

      % Empty arrays
      testCase.test({zeros(0, 0)}, zeros(0, 1), zeros(0, 1))
      testCase.test({zeros(1, 0)}, [1], [1])
      testCase.test({zeros(3, 0)}, [1], [3])

      % Non-empty numeric arrays
      testCase.test({[4]},           [1],       [1]);
      testCase.test({[3;3]},         [1],       [2]);
      testCase.test({[7;7;1;1]},     [1;3],     [2;4]);
      testCase.test({[11;12;13;14]}, [1;2;3;4], [1;2;3;4]);
      testCase.test({[NaN;NaN;3]},   [1;3],     [2;3]);
      testCase.test({[ ...
        1 1; ...
        1 1; ...
        1 3; ...
        3 1;...
        ]}, ...
        [1;3;4], ...
        [2;3;4]);

      % 3D arrays
      A(:,:,1) = [1 2; 1 2; 1 3];
      A(:,:,2) = [3 4; 3 4; 3 4];
      testCase.test({A}, [1;3], [2;3])

      % Non-empty char string cell arrays.
      testCase.test({{'abc'; 'DEFG'; 'HIJKL'}}, [1;2;3], [1;2;3]);
      testCase.test({{'abc'; 'abc';  'abc'  }}, [1],     [3]);
      testCase.test({
        {'abc'; 'abc'; 'abc'; 'DEFG'; 'DEFG'; 'DEFG'}, ...
        logical([0;0;1;1;0;0])}, ...
        [1;3;4;5], [2;3;4;6]);

      % Multiple arrays, mixed types.
      DATA_AR = {...
        1, 'ab'; ...
        1, 'ab'; ...
        2, 'ab'; ...
        1, 'ab'; ...
        1, 'cde'; ...
        1, 'cde'};
      % Numeric + char arrays
      testCase.test( ...
        { cell2mat(DATA_AR(:, 1)), DATA_AR(:, 2) }, ...
        [1; 3; 4; 5], ...
        [2; 3; 4; 6])
      % Numeric + string
      testCase.test( ...
        { cell2mat(DATA_AR(:, 1)), string(DATA_AR(:, 2)) }, ...
        [1; 3; 4; 5], ...
        [2; 3; 4; 6])
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function test(testCase, arraysCa, expI1Array, expI2Array)
      expN = numel(expI1Array);

      [actI1Array, actI2Array, actN] = ...
        irf.utils.split_by_change(arraysCa{:});

      testCase.assertEqual(actI1Array, expI1Array)
      testCase.assertEqual(actI2Array, expI2Array)
      testCase.assertEqual(actN,       expN)
    end

    function test_exc(testCase, arraysCa)
      testCase.assertError(...
        @() irf.utils.split_by_change(arraysCa{:}), ...
        ?MException)
    end




  end    % methods(Access=private)



end
