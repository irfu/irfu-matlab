%
% matlab.unittest automatic test code for bicas.utils.group_unique_rows().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef group_unique_rows___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty(testCase)
      % Zero arrays
      testCase.test({}, cell(0, 1))

      % Zero rows
      testCase.test({zeros(0, 1)},                 cell(0, 1))
      testCase.test({zeros(0, 1), zeros(0, 2, 3)}, cell(0, 1))
    end



    function test_nonempty(testCase)
      % One row.
      testCase.test({zeros(1, 0)},              {[1]})
      testCase.test({[3]},                      {[1]})
      testCase.test({[3], [4,5], 'abc', "ABC"}, {[1]})
    end



    function test_1(testCase)
      % Multiple rows
      testCase.test({[3;3;3]}, {[1;2;3]})
      testCase.test({[3;3;3]}, {[1;2;3]})
      testCase.test({[5; 6]}, {[1], [2]})
      testCase.test({[5; 6; 6; 5; 5]}, {[1; 4; 5], [2; 3]})

      % Multiple arrays, MATLAB classes
      testCase.test(...
        {...
        [1; 1; 1; 1; 1; 1; 2; 1], ...
        [3; 3; 9; 3; 3; 3; 3; 3]
        }, ...
        {[1;2;4;5;6;8], [3], [7]})


      testCase.test({["abc"; "cd"; "abc"]}, {[1;3]; [2]})
    end



    function test_high_dim_arrays(testCase)
      M(:, :, 1) = [1, 1, 1; 1, 1, 1; 2, 2, 2];
      M(:, :, 2) = [1, 2, 3; 1, 2, 3; 1, 2, 3];

      testCase.test({M}, {[1:2]', [3]'})
    end



    function test_bug(testCase)
      % Test against old bug: Comparison error due to not extracting the right
      % elements.
      testCase.test({[8, 9; 8, 9]}, {[1:2]'})
    end



    function test_error(testCase)
      testCase.test(   {zeros(1, 1), zeros(1, 1)}, {[1]})
      testCase.test_exc(zeros(0, 1), zeros(1, 1))

      testCase.test(   {zeros(3, 1), zeros(3, 1)}, {[1;2;3]})
      testCase.test_exc(zeros(3, 1), zeros(4, 1))
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function test(testCase, arraysCa, expIGroupCa)
      assert(iscell(arraysCa))
      expIGroupCa = expIGroupCa(:);

      % CALL TESTED FUNCTION
      actIGroupCa = bicas.utils.group_unique_rows(arraysCa{:});

      testCase.assertTrue(iscolumn(actIGroupCa))
      testCase.assertEqual(actIGroupCa, expIGroupCa)

      actAllIndices = sort(cat(1, actIGroupCa{:}));
      if isempty(arraysCa)
        testCase.assertTrue(isempty(actAllIndices))
      else
        nRows = size(arraysCa{1}, 1);
        % NOTE: nRows==0 is possible.
        testCase.assertEqual(actAllIndices(:), [1:nRows]')
      end

      % Check the internal order of elements
      % ------------------------------------
      % IMPLEMENTATION NOTE: The orders checked should in principle not be
      % important for nominal use of the function but is deliberately made
      % well-defined only to make tests rigorous.
      actIFirstInGroup = cellfun(@(x) (x(1)), actIGroupCa);
      testCase.assertTrue(issorted(actIFirstInGroup))

      for i = 1:numel(actIGroupCa)
        groupAr = actIGroupCa{i};
        testCase.assertTrue(issorted(groupAr));
      end
    end



    function test_exc(testCase, varargin)
      testCase.assertError(...
        @() bicas.utils.group_unique_rows(varargin{:}), ...
        ?MException)
    end



  end    % methods(Access=private)



end
