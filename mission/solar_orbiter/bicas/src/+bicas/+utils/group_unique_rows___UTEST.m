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



    function test_empty(T)
      % Zero arrays
      T.test({}, cell(0, 1))

      % Zero rows
      T.test({zeros(0, 1)},                 cell(0, 1))
      T.test({zeros(0, 1), zeros(0, 2, 3)}, cell(0, 1))
    end



    function test_nonempty(T)
      % One row.
      T.test({zeros(1, 0)},              {[1]})
      T.test({[3]},                      {[1]})
      T.test({[3], [4,5], 'abc', "ABC"}, {[1]})
    end



    function test_1(T)
      % Multiple rows
      T.test({[3;3;3]}, {[1;2;3]})
      T.test({[3;3;3]}, {[1;2;3]})
      T.test({[5; 6]}, {[1], [2]})
      T.test({[5; 6; 6; 5; 5]}, {[1; 4; 5], [2; 3]})

      % Multiple arrays, MATLAB classes
      T.test(...
        {...
        [1; 1; 1; 1; 1; 1; 2; 1], ...
        [3; 3; 9; 3; 3; 3; 3; 3]
        }, ...
        {[1;2;4;5;6;8], [3], [7]})


      T.test({["abc"; "cd"; "abc"]}, {[1;3]; [2]})
    end



    function test_high_dim_arrays(T)
      M(:, :, 1) = [1, 1, 1; 1, 1, 1; 2, 2, 2];
      M(:, :, 2) = [1, 2, 3; 1, 2, 3; 1, 2, 3];

      T.test({M}, {[1:2]', [3]'})
    end



    function test_bug(T)
      % Test against old bug: Comparison error due to not extracting the right
      % elements.
      T.test({[8, 9; 8, 9]}, {[1:2]'})
    end



    function test_error(T)
      T.test(   {zeros(1, 1), zeros(1, 1)}, {[1]})
      T.test_exc(zeros(0, 1), zeros(1, 1))

      T.test(   {zeros(3, 1), zeros(3, 1)}, {[1;2;3]})
      T.test_exc(zeros(3, 1), zeros(4, 1))
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function test(T, arraysCa, expIGroupCa)
      assert(iscell(arraysCa))
      expIGroupCa = expIGroupCa(:);

      % CALL TESTED FUNCTION
      actIGroupCa = bicas.utils.group_unique_rows(arraysCa{:});

      T.assertTrue(iscolumn(actIGroupCa))
      T.assertEqual(actIGroupCa, expIGroupCa)

      actAllIndices = sort(cat(1, actIGroupCa{:}));
      if isempty(arraysCa)
        T.assertTrue(isempty(actAllIndices))
      else
        nRows = size(arraysCa{1}, 1);
        % NOTE: nRows==0 is possible.
        T.assertEqual(actAllIndices(:), [1:nRows]')
      end

      % Check the internal order of elements
      % ------------------------------------
      % IMPLEMENTATION NOTE: The orders checked should in principle not be
      % important for nominal use of the function but is deliberately made
      % well-defined only to make tests rigorous.
      actIFirstInGroup = cellfun(@(x) (x(1)), actIGroupCa);
      T.assertTrue(issorted(actIFirstInGroup))

      for i = 1:numel(actIGroupCa)
        groupAr = actIGroupCa{i};
        T.assertTrue(issorted(groupAr));
      end
    end



    function test_exc(T, varargin)
      T.assertError(...
        @() bicas.utils.group_unique_rows(varargin{:}), ...
        ?MException)
    end



  end    % methods(Access=private)



end
