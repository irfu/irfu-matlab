%
% matlab.unittest automatic test code for bicas.utils.group_by_change().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef group_by_change___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % NOTE: Not real thorough testing. Mostly to see that the correct
    % information passes through the wrapper.
    function test_0(testCase)
      testCase.test_exc({})

      % Unequal number of rows.
      testCase.test_exc({zeros(0, 1), [1]})
      testCase.test_exc({[2;3],       [1]})
      testCase.test_exc({[2,3],       [1;2]})

      testCase.test({zeros(0, 1), zeros(0,2,3)}, cell(0, 1))
      testCase.test({ ...
        [1;1;1; 1; 2;2;2], ...
        [3;3;3; 4; 4;4;4], ...
        }, ...
        {[1;2;3]; [4]; [5;6;7]})
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
      actIGroupCa = bicas.utils.group_by_change(arraysCa{:});

      testCase.assertTrue(iscolumn(actIGroupCa))
      testCase.assertEqual(actIGroupCa, expIGroupCa)

      % Further assertions on the order of cells and array elements
      % -----------------------------------------------------------
      % IMPLEMENTATION NOTE: The orders checked should in principle not be
      % important for nominal use of the function but is deliberately made
      % well-defined only to make tests rigorous.

      % ASSERTION: All indices in all groups combined are specified exactly
      %            once.
      actAllIndices = sort(cat(1, actIGroupCa{:}));
      if isempty(arraysCa)
        testCase.assertTrue(isempty(actAllIndices))
      else
        nRows = size(arraysCa{1}, 1);
        % NOTE: nRows==0 is possible.
        testCase.assertEqual(actAllIndices(:), [1:nRows]')
      end

      % ASSERTION: All groups are sorted by the first index within each group.
      actIFirstInGroup = cellfun(@(x) (x(1)), actIGroupCa);
      testCase.assertTrue(issorted(actIFirstInGroup))

      % ASSERTION: Indices are sorted within each group.
      for i = 1:numel(actIGroupCa)
        groupAr = actIGroupCa{i};
        testCase.assertTrue(issorted(groupAr));
      end
    end



    function test_exc(testCase, arraysCa)
      testCase.assertError(...
        @() bicas.utils.group_by_change(arraysCa{:}), ...
        ?MException)
    end



  end    % methods(Access=private)



end
