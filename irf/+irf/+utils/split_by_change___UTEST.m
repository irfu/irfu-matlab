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

      function test(vararginCa, expI1Array, expI2Array)
        expN = numel(expI1Array);

        [actI1Array, actI2Array, actN] = ...
          irf.utils.split_by_change(vararginCa{:});

        testCase.assertEqual(actI1Array, expI1Array)
        testCase.assertEqual(actI2Array, expI2Array)
        testCase.assertEqual(actN,       expN)
      end

      function test_exc(vararginCa)
        testCase.assertError(...
          @() irf.utils.split_by_change(vararginCa{:}), ...
          ?MException)
      end

      test_exc({})
      % Unequal number of rows.
      test_exc({zeros(0, 1), [1]})
      test_exc({[2;3], [1]})
      test_exc({[2,3], [1;2]})

      % Empty arrays
      test({zeros(0, 0)}, zeros(0, 1), zeros(0, 1))
      test({zeros(1, 0)}, [1], [1])
      test({zeros(3, 0)}, [1], [3])

      % Non-empty numeric arrays
      test({[4]},           [1],       [1]);
      test({[3;3]},         [1],       [2]);
      test({[7;7;1;1]},     [1;3],     [2;4]);
      test({[11;12;13;14]}, [1;2;3;4], [1;2;3;4]);
      test({[NaN;NaN;3]},   [1;3],     [2;3]);
      test({[ ...
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
      test({A}, [1;3], [2;3])

      % Non-empty char string cell arrays.
      test({{'abc'; 'DEFG'; 'HIJKL'}}, [1;2;3], [1;2;3]);
      test({{'abc'; 'abc';  'abc'  }}, [1],     [3]);
      test({
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
      test( ...
        { cell2mat(DATA_AR(:, 1)), DATA_AR(:, 2) }, ...
        [1; 3; 4; 5], ...
        [2; 3; 4; 6])
      % Numeric + string
      test( ...
        { cell2mat(DATA_AR(:, 1)), string(DATA_AR(:, 2)) }, ...
        [1; 3; 4; 5], ...
        [2; 3; 4; 6])
    end



  end    % methods(Test)



end
