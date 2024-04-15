%
% matlab.unittest automatic test code for solo.qli.batch.UmdFmdDictionary.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef UmdFmdDictionary___UTEST < matlab.unittest.TestCase



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    % NOTE: Tests assumes that values increment.
    C = struct(...
      'UMD_1', solo.qli.utils.umdt('2020-01-01'), ...
      'UMD_2', solo.qli.utils.umdt('2020-01-02'), ...
      'FMD_1', datetime('2025-01-01 01:01'), ...
      'FMD_2', datetime('2025-01-01 02:02'), ...
      'FMD_3', datetime('2025-01-01 03:03') ...
      )
  end



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and
    % teardown methods which store/read their own data from the testCase
    % object.
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_constructor_subsasgn_subsref_n(testCase)
      C = testCase.C;

      function assert_DFMDD_1(Dfmdd)
        testCase.assertEqual(Dfmdd(C.UMD_1), C.FMD_1)
        testCase.assertEqual(Dfmdd.n, 1)
        testCase.assertEqual(Dfmdd.DaysDtArray, C.UMD_1)
        testCase.assertEqual(Dfmdd.FmdDtArray,  C.FMD_1)
      end

      function assert_DFMDD_2(Dfmdd)
        testCase.assertEqual(Dfmdd(C.UMD_1), C.FMD_1)
        testCase.assertEqual(Dfmdd(C.UMD_2), C.FMD_2)
        testCase.assertEqual(Dfmdd.n, 2)
        % NOTE: Not assuming order of arrays.
        testCase.assertEqual(sort(Dfmdd.DaysDtArray), sort([C.UMD_1; C.UMD_2 ]))
        testCase.assertEqual(sort(Dfmdd.FmdDtArray),  sort([C.FMD_1; C.FMD_2]))
      end

      Dfmdd = solo.qli.batch.UmdFmdDictionary();
      testCase.assertEqual(Dfmdd.n, 0)
      testCase.assertEqual(Dfmdd.DaysDtArray, solo.qli.const.EMPTY_DT_ARRAY)
      testCase.assertEqual(Dfmdd.FmdDtArray,  datetime.empty(0, 1))

      Dfmdd(C.UMD_1) = C.FMD_1;
      Dfmdd2 = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_1);
      assert_DFMDD_1(Dfmdd)
      assert_DFMDD_1(Dfmdd2)

      Dfmdd(C.UMD_2) = C.FMD_2;
      Dfmdd2 = solo.qli.batch.UmdFmdDictionary([C.UMD_1; C.UMD_2], [C.FMD_1; C.FMD_2]);
      assert_DFMDD_2(Dfmdd)
      assert_DFMDD_2(Dfmdd2)
    end



    function test_isequal(testCase)
      C = testCase.C;

      Dfmdd1 = solo.qli.batch.UmdFmdDictionary();
      Dfmdd2 = solo.qli.batch.UmdFmdDictionary();
      testCase.assertTrue(isequal(Dfmdd1, Dfmdd2))

      Dfmdd1 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_1; C.UMD_2], [C.FMD_1; C.FMD_2]);
      Dfmdd2 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_1; C.UMD_2], [C.FMD_1; C.FMD_2]);
      Dfmdd3 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_2; C.UMD_1], [C.FMD_2; C.FMD_1]);
      Dfmdd4 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_2], [C.FMD_2]);
      testCase.assertTrue(isequal(Dfmdd1, Dfmdd2))
      testCase.assertTrue(isequal(Dfmdd1, Dfmdd3))
      testCase.assertFalse(isequal(Dfmdd1, Dfmdd4))
    end



    function test_set_if_smaller(testCase)
      C = testCase.C;

      Dfmdd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpDfmdd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      Dfmdd = Dfmdd.set_if_smaller(              C.UMD_1, C.FMD_3);
      testCase.assertEqual(Dfmdd, ExpDfmdd)

      Dfmdd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpDfmdd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_1);
      Dfmdd = Dfmdd.set_if_smaller(              C.UMD_1, C.FMD_1);
      testCase.assertEqual(Dfmdd, ExpDfmdd)
    end



    function test_set_if_greater(testCase)
      C = testCase.C;

      Dfmdd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpDfmdd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_3);
      Dfmdd = Dfmdd.set_if_greater(              C.UMD_1, C.FMD_3);
      testCase.assertEqual(Dfmdd, ExpDfmdd)

      Dfmdd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpDfmdd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      Dfmdd = Dfmdd.set_if_greater(              C.UMD_1, C.FMD_1);
      testCase.assertEqual(Dfmdd, ExpDfmdd)
    end



    % NOTE: Mostly tests using double (not datetime) as keys and values.
    function test_merge_max(testCase)

      function test(DfmddCa, ExpDfmdd)
        ActDfmdd = solo.qli.batch.UmdFmdDictionary.merge_max(...
          DfmddCa);

        testCase.assertEqual(ActDfmdd, ExpDfmdd)
      end

      % Call test() for every permutation of input DFMDDs.
      function test_perm(DfmddCa, ExpDfmdd)
        permutations = perms(1:numel(DfmddCa));

        if isempty(permutations)
          permutations = nan(1, 0);    % 1x0 (not 0x1).
        end

        for iPerm = 1:size(permutations, 1)
          DfmddCa2 = DfmddCa(permutations(iPerm, :), 1);

          test(DfmddCa2, ExpDfmdd)
        end
      end

      %=========================================================================

      UMD_1 = solo.qli.utils.umdt('2020-01-01');
      UMD_2 = solo.qli.utils.umdt('2020-01-02');
      FMD_1 = datetime('2030-01-01');
      FMD_2 = datetime('2030-01-02');

      test(...
        cell(0, 1), ...
        solo.qli.batch.UmdFmdDictionary())

      test(...
        {solo.qli.batch.UmdFmdDictionary()}, ...
        solo.qli.batch.UmdFmdDictionary())

      test_perm(...
        { ...
        solo.qli.batch.UmdFmdDictionary(); ...
        solo.qli.batch.UmdFmdDictionary(UMD_1, FMD_1)
        }, ...
        solo.qli.batch.UmdFmdDictionary(UMD_1, FMD_1))

      test_perm(...
        { ...
        solo.qli.batch.UmdFmdDictionary(UMD_1, FMD_1);
        solo.qli.batch.UmdFmdDictionary(UMD_2, FMD_2) ...
        }, ...
        solo.qli.batch.UmdFmdDictionary([UMD_1; UMD_2], [FMD_1; FMD_2]))

      test_perm(...
        { ...
        solo.qli.batch.UmdFmdDictionary(UMD_1, FMD_1);
        solo.qli.batch.UmdFmdDictionary(UMD_1, FMD_2) ...
        }, ...
        solo.qli.batch.UmdFmdDictionary(UMD_1, FMD_2))

      test_perm(...
        { ...
        %solo.qli.batch.UmdFmdDictionary([1, 3, 4], [11, 6, 8]);
        %solo.qli.batch.UmdFmdDictionary([2, 3, 4], [22, 5, 9]) ...
        solo.qli.batch.UmdFmdDictionary(UMD_1+caldays([1; 3; 4]), FMD_1+caldays([11; 6; 8]));
        solo.qli.batch.UmdFmdDictionary(UMD_1+caldays([2; 3; 4]), FMD_1+caldays([22; 5; 9])) ...
        }, ...
        solo.qli.batch.UmdFmdDictionary(UMD_1+caldays([1; 2; 3; 4]), FMD_1+caldays([11; 22; 6; 9])))
        %solo.qli.batch.UmdFmdDictionary([1, 2, 3, 4], [11, 22, 6, 9]))
    end



  end    % methods(Test)



end
