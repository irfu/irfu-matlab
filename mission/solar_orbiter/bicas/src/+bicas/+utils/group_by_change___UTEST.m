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
    function test_0(T)
      T.test_exc({})

      % Unequal number of rows.
      T.test_exc({zeros(0, 1), [1]})
      T.test_exc({[2;3],       [1]})
      T.test_exc({[2,3],       [1;2]})

      T.test({zeros(0, 1), zeros(0,2,3)}, cell(0, 1))
      T.test({ ...
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



    function test(T, arraysCa, expIGroupCa)
      assert(iscell(arraysCa))
      expIGroupCa = expIGroupCa(:);

      % CALL TESTED FUNCTION
      actIGroupCa = bicas.utils.group_by_change(arraysCa{:});

      T.assertTrue(iscolumn(actIGroupCa))
      T.assertEqual(actIGroupCa, expIGroupCa)

      % Further assertions on the order of cells and array elements
      % -----------------------------------------------------------
      % IMPLEMENTATION NOTE: The orders checked should in principle not be
      % important for nominal use of the function but is deliberately made
      % well-defined only to make tests rigorous.

      % ASSERTION: All indices in all groups combined are specified exactly
      %            once.
      actAllIndices = sort(cat(1, actIGroupCa{:}));
      if isempty(arraysCa)
        T.assertTrue(isempty(actAllIndices))
      else
        nRows = size(arraysCa{1}, 1);
        % NOTE: nRows==0 is possible.
        T.assertEqual(actAllIndices(:), [1:nRows]')
      end

      % ASSERTION: All groups are sorted by the first index within each group.
      actIFirstInGroup = cellfun(@(x) (x(1)), actIGroupCa);
      T.assertTrue(issorted(actIFirstInGroup))

      % ASSERTION: Indices are sorted within each group.
      for i = 1:numel(actIGroupCa)
        groupAr = actIGroupCa{i};
        T.assertTrue(issorted(groupAr));
      end
    end



    function test_exc(T, arraysCa)
      T.assertError(...
        @() bicas.utils.group_by_change(arraysCa{:}), ...
        ?MException)
    end



  end    % methods(Access=private)



end
