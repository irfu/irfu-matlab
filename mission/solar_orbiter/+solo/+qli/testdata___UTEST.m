%
% matlab.unittest automatic test code for solo.qli.testdata.
%
% Could add tests for mroe code, but that seems unnecessary. The tested code
% supports automatic and manual tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef testdata___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_2D_array(testCase)

      N_LARGE = 100;

      % NOTE: Does not test the exact return value, only its general
      % properties. The exact return value is not important and can change
      % with future implementations.
      function test(nRows, nCols, aMin, aMax)
        VARIATION_SCALE = 0.8;
        SCALES_CA = {'LIN', 'LOG'};

        for i = 1:numel(SCALES_CA)
          scale = SCALES_CA{i};

          if (aMin < 0) && strcmp(scale, 'LOG')
            continue
          end

          actA = solo.qli.testdata.get_2D_array(nRows, nCols, aMin, aMax, scale);

          actAMin = min(actA, [], 'all');
          actAMax = max(actA, [], 'all');

          testCase.assertEqual(size(actA), [nRows, nCols])
          testCase.assertTrue(actAMin >= aMin)
          testCase.assertTrue(actAMax <= aMax)

          if nRows * nCols >= N_LARGE
            % Somewhat arbitrary condition that is probably true.
            testCase.assertTrue( (actAMax-actAMin) >= (aMax-aMin)*VARIATION_SCALE )
          end
        end
      end

      if 1
        test(1, 1, 10, 20)
        test(3, 1, 10, 20)
        test(1, 4, 10, 20)
        test(3, 4, 10, 20)
      end

      n = ceil(sqrt(N_LARGE));
      test(n, n,       10, 20)
      test(N_LARGE, 1, 10, 20)
      test(1, N_LARGE, 10, 20)
    end



  end    % methods(Test)



end
