%
% matlab.unittest automatic test code for solo.qli.batch.utils.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % NOTE: Mostly tests using double (not datetime) as keys and values.
    function test_merge_dictionaries_max(testCase)

      function test(DictCa, keyType, valueType, ExpDict)
        ActDict = solo.qli.batch.utils.merge_dictionaries_max(...
          DictCa, keyType, valueType);

        testCase.assertEqual(ActDict, ExpDict)
      end

      function test_perm(DictCa, keyType, valueType, ExpDict)
        permutations = perms(1:numel(DictCa));

        if isempty(permutations)
          permutations = nan(1, 0);    % 1x0 (not 0x1).
        end

        for iPerm = 1:size(permutations, 1)
          DictCa2 = DictCa(permutations(iPerm, :), 1);

          test(DictCa2, keyType, valueType, ExpDict)
        end
      end

      %=========================================================================

      test(...
        cell(0, 1), ...
        datetime.empty, datetime.empty, ...
        dictionary(datetime.empty, datetime.empty))

      test(...
        {dictionary(double.empty, double.empty)}, ...
        double.empty, double.empty, ...
        dictionary(double.empty, double.empty))

      test_perm(...
        { ...
        dictionary([], []); ...
        dictionary([3], [33])
        }, ...
        double.empty, double.empty, ...
        dictionary([3], [33]))

      test_perm(...
        { ...
        dictionary([3], [33]);
        dictionary([4], [44]) ...
        }, ...
        double.empty, double.empty, ...
        dictionary([3, 4], [33, 44]))

      test_perm(...
        { ...
        dictionary([3], [8]);
        dictionary([3], [9]) ...
        }, ...
        double.empty, double.empty, ...
        dictionary([3], [9]))

      test_perm(...
        { ...
        dictionary([1, 3, 4], [11, 6, 8]);
        dictionary([2, 3, 4], [22, 5, 9]) ...
        }, ...
        double.empty, double.empty, ...
        dictionary([1, 2, 3, 4], [11, 22, 6, 9]))
    end



  end    % methods(Test)



end
