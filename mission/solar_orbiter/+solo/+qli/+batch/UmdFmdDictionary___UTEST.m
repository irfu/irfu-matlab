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
      'UMD_1', irf.dt.um('2020-01-01'), ...
      'UMD_2', irf.dt.um('2020-01-02'), ...
      'FMD_1', datetime('2025-01-01 01:01'), ...
      'FMD_2', datetime('2025-01-01 02:02'), ...
      'FMD_3', datetime('2025-01-01 03:03') ...
      )
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_constructor_subsasgn_subsref_n(testCase)
      C = testCase.C;

      function assert_UFD_1(Ufd)
        testCase.assertEqual(Ufd(C.UMD_1), C.FMD_1)
        testCase.assertEqual(Ufd.n, 1)
        testCase.assertEqual(Ufd.UmdDtArray, C.UMD_1)
        testCase.assertEqual(Ufd.FmdDtArray, C.FMD_1)
      end

      function assert_UFD_2(Ufd)
        testCase.assertEqual(Ufd(C.UMD_1), C.FMD_1)
        testCase.assertEqual(Ufd(C.UMD_2), C.FMD_2)
        testCase.assertEqual(Ufd.n, 2)
        % NOTE: Not assuming order of arrays.
        testCase.assertEqual(sort(Ufd.UmdDtArray), sort([C.UMD_1; C.UMD_2 ]))
        testCase.assertEqual(sort(Ufd.FmdDtArray), sort([C.FMD_1; C.FMD_2]))
      end

      Ufd = solo.qli.batch.UmdFmdDictionary();
      testCase.assertEqual(Ufd.n, 0)
      testCase.assertEqual(Ufd.UmdDtArray, solo.qli.const.EMPTY_DT_ARRAY)
      testCase.assertEqual(Ufd.FmdDtArray, datetime.empty(0, 1))

      Ufd(C.UMD_1) = C.FMD_1;
      Ufd2 = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_1);
      assert_UFD_1(Ufd)
      assert_UFD_1(Ufd2)

      Ufd(C.UMD_2) = C.FMD_2;
      Ufd2 = solo.qli.batch.UmdFmdDictionary([C.UMD_1; C.UMD_2], [C.FMD_1; C.FMD_2]);
      assert_UFD_2(Ufd)
      assert_UFD_2(Ufd2)
    end



    function test_isequal(testCase)
      C = testCase.C;

      Ufd1 = solo.qli.batch.UmdFmdDictionary();
      Ufd2 = solo.qli.batch.UmdFmdDictionary();
      testCase.assertTrue(isequal(Ufd1, Ufd2))

      Ufd1 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_1; C.UMD_2], [C.FMD_1; C.FMD_2]);
      Ufd2 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_1; C.UMD_2], [C.FMD_1; C.FMD_2]);
      Ufd3 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_2; C.UMD_1], [C.FMD_2; C.FMD_1]);
      Ufd4 = solo.qli.batch.UmdFmdDictionary(...
        [C.UMD_2], [C.FMD_2]);
      testCase.assertTrue(isequal(Ufd1, Ufd2))
      testCase.assertTrue(isequal(Ufd1, Ufd3))
      testCase.assertFalse(isequal(Ufd1, Ufd4))
    end



    function test_set_if_smaller(testCase)
      C = testCase.C;

      Ufd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpUfd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      Ufd = Ufd.set_if_smaller(                C.UMD_1, C.FMD_3);
      testCase.assertEqual(Ufd, ExpUfd)

      Ufd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpUfd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_1);
      Ufd = Ufd.set_if_smaller(                C.UMD_1, C.FMD_1);
      testCase.assertEqual(Ufd, ExpUfd)
    end



    function test_set_if_greater(testCase)
      C = testCase.C;

      Ufd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpUfd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_3);
      Ufd = Ufd.set_if_greater(                C.UMD_1, C.FMD_3);
      testCase.assertEqual(Ufd, ExpUfd)

      Ufd    = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      ExpUfd = solo.qli.batch.UmdFmdDictionary(C.UMD_1, C.FMD_2);
      Ufd = Ufd.set_if_greater(                C.UMD_1, C.FMD_1);
      testCase.assertEqual(Ufd, ExpUfd)
    end



    % NOTE: Mostly tests using double (not datetime) as keys and values.
    function test_merge_max(testCase)

      function test(UfdCa, ExpUfd)
        ActUfd = solo.qli.batch.UmdFmdDictionary.merge_max(...
          UfdCa);

        testCase.assertEqual(ActUfd, ExpUfd)
      end

      % Call test() for every permutation of input UFDs.
      function test_perm(UfdCa, ExpUfd)
        permutations = perms(1:numel(UfdCa));

        if isempty(permutations)
          permutations = nan(1, 0);    % 1x0 (not 0x1).
        end

        for iPerm = 1:size(permutations, 1)
          UfdCa2 = UfdCa(permutations(iPerm, :), 1);

          test(UfdCa2, ExpUfd)
        end
      end

      %=========================================================================

      UMD_1 = irf.dt.um('2020-01-01');
      UMD_2 = irf.dt.um('2020-01-02');
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
        solo.qli.batch.UmdFmdDictionary(UMD_1+caldays([1; 3; 4]), FMD_1+caldays([11; 6; 8]));
        solo.qli.batch.UmdFmdDictionary(UMD_1+caldays([2; 3; 4]), FMD_1+caldays([22; 5; 9])) ...
        }, ...
        solo.qli.batch.UmdFmdDictionary(UMD_1+caldays([1; 2; 3; 4]), FMD_1+caldays([11; 22; 6; 9])))
    end



  end    % methods(Test)



end
