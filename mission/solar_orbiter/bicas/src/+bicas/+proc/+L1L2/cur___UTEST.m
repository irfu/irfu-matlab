%
% matlab.unittest automatic test code for bicas.proc.L1L2.cur.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef cur___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_calibrate_bias_currents(testCase)
      % NOTE: Could also check
      %       INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY but does not.

      function [argsCa, expCurrentAampere] = convert_test_arguments(...
          Tt2000CurrentNanoSampereAr1, tt2000CurrentSampereAr2, Ccal, policyId)

        Bso = bicas.create_default_BSO();
        Bso.override_value('PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY', ...
          policyId, 'test')
        Bso.make_read_only();

        L  = bicas.Logger('HUMAN_READABLE', false);

        curTt2000Ar          = int64(Tt2000CurrentNanoSampereAr1(:, 1));
        currentNanoSampereAr =       Tt2000CurrentNanoSampereAr1(:, 2:4);

        sciTt2000Ar       = int64(tt2000CurrentSampereAr2(:, 1));
        expCurrentAampere =       tt2000CurrentSampereAr2(:, 2:4);

        argsCa = {curTt2000Ar, currentNanoSampereAr, sciTt2000Ar, Ccal, Bso, L};
      end

      function test(varargin)
        [argsCa, expCurrentAampere] = convert_test_arguments(varargin{:});

        actCurrentAampere = ...
          bicas.proc.L1L2.cur.calibrate_bias_currents(argsCa{:});

        testCase.assertEqual(actCurrentAampere, expCurrentAampere)
      end


      % Could add tests but it is doubtful if that is needed at this point.
      % Subfunctions are tested more.
      N    = NaN;
      Ccal = bicas.proc.L1L2.cal.CurrentCalibrationTest( ...
        3, ...
        dictionary(...
        int64([20 31 32 33 41 42 43]), ...
        [0 1 1 1  2 2 2]));
      test(...
        [ ...
        31 11e9 N     N; ...
        32 N    12e9  N; ...
        33 N    N    13e9; ...
        41 21e9 N    N; ...
        42 N    22e9 N; ...
        43 N    N    23e9 ...
        ], ...
        [ ...
        20 N         N         N; ...
        31 1+1*3*11  N         N; ...
        32 1+1*3*11  2+1*3*12  N; ...
        33 1+1*3*11  2+1*3*12  3+1*3*13; ...
        41 1+2*3*21  2+2*3*12  3+2*3*13; ...
        42 1+2*3*21  2+2*3*22  3+2*3*13; ...
        43 1+2*3*21  2+2*3*22  3+2*3*23 ...
        ], Ccal, 'WARNING')
    end



    function test_convert_CUR_to_CUR_on_SCI_TIME(testCase)
      % NOTE: Could also check
      %       INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY but does not.

      function [argsCa, expCurrentSampere] = convert_test_arguments(...
          Tt2000CurrentNanoSampereAr1, tt2000CurrentSampereAr2, policyId)

        Bso = bicas.create_default_BSO();
        Bso.override_value('PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY', ...
          policyId, 'test')
        Bso.make_read_only();

        L  = bicas.Logger('HUMAN_READABLE', false);

        curTt2000Ar      = int64(Tt2000CurrentNanoSampereAr1(:, 1));
        currentSampereAr =       Tt2000CurrentNanoSampereAr1(:, 2:4);

        sciTt2000Ar       = int64(tt2000CurrentSampereAr2(:, 1));
        expCurrentSampere =       tt2000CurrentSampereAr2(:, 2:4);

        argsCa = {curTt2000Ar, currentSampereAr, sciTt2000Ar, Bso, L};
      end

      function test(varargin)
        [argsCa, expCurrentSampere] = convert_test_arguments(varargin{:});

        % CALL TESTED CODE
        actCurrentSampere = bicas.proc.L1L2.cur.convert_CUR_to_CUR_on_SCI_TIME(...
          argsCa{:});

        testCase.assertEqual(actCurrentSampere, expCurrentSampere, "RelTol", 1e-15)
      end

      function test_exc(varargin)
        argsCa = convert_test_arguments(varargin{:});

        % CALL TESTED CODE
        testCase.assertError(...
          @() bicas.proc.L1L2.cur.convert_CUR_to_CUR_on_SCI_TIME(argsCa{:}), ...
          ?MException)
      end

      % IOP = Independent Of Policy
      function test_IOP(Tt2000CurrentNanoSampereAr1, tt2000CurrentSampereAr2)
        for policyId = ["WARNING", "ERROR"]
          test(Tt2000CurrentNanoSampereAr1, tt2000CurrentSampereAr2, char(policyId))
        end
      end

      % DOP = Depends on POLICY
      function test_DOP(Tt2000CurrentNanoSampereAr1, tt2000CurrentSampereAr2)
        test(    Tt2000CurrentNanoSampereAr1, tt2000CurrentSampereAr2, 'WARNING')
        test_exc(Tt2000CurrentNanoSampereAr1, tt2000CurrentSampereAr2, 'ERROR')
      end

      N = NaN;

      test_IOP(zeros(0, 4), zeros(0, 4))
      test_IOP( ...
        [33 1    2    3   ], ...
        [34 1e-9 2e-9 3e-9])

      test_DOP( ...
        [33 1  2  3], ...
        [32 N  N  N])

      % Complex test.
      test_DOP( ...
        [ ...
        31 11 N  N; ...
        32 N  12 N; ...
        33 N  N  13; ...
        41 21 N  N; ...
        42 N  22 N; ...
        43 N  N  23 ...
        ], ...
        [ ...
        20 N      N      N; ...
        31 11e-9  N      N; ...
        32 11e-9  12e-9  N; ...
        33 11e-9  12e-9  13e-9; ...
        41 21e-9  12e-9  13e-9; ...
        42 21e-9  22e-9  13e-9; ...
        43 21e-9  22e-9  23e-9 ...
        ])
    end



  end    % methods(Test)



end
