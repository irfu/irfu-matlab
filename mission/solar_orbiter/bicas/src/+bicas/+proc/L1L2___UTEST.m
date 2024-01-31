%
% matlab.unittest automatic test code for bicas.proc.L1L2.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef L1L2___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_autodetect_sweeps(testCase)
      % PROPOSAL: Split up in multiple test functions.
      %   CON: Needs to call static functions in test class which tests multiple
      %        functions (though only one is tested right now).
      %   CON: Long function names to distinguish test functions from test
      %        functions for other tested functions.
      %   PROPOSAL: Redefine as test class for
      %             bicas.proc.L1L2.autodetect_sweeps() only.
      %
      % PROPOSAL: test() accepts one struct. ==> Named arguments.
      % PROPOSAL: Abbreviation for currentMmDiffMinTm.

      % Generic test. All input and expected output is specified in arguments.
      %
      % ARGUMENTS
      % =========
      % expIsSweeping
      %       logical FPA as double/NaN
      % bdm, hkBiasCurrent
      %       FPAs as float/NaN
      %
      function test(...
          tt2000, bdm, hkBiasCurrent, ...
          sbdaEndTt2000, windowLengthPts, currentMinMaxDiffThresholdTm, windowMarginSec, ...
          expIsSweeping)

        assert(issorted(tt2000), 'ascend')
        assert(isa(expIsSweeping, 'double'))
        bdmFpa           = bicas.utils.FPArray.floatNan2int(bdm, 'uint8');
        hkBiasCurrentFpa = bicas.utils.FPArray.floatNan2int(hkBiasCurrent, 'uint16');

        Bso = bicas.create_default_BSO();
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC',                         spdfbreakdowntt2000(sbdaEndTt2000), 'test');
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS',               windowLengthPts,                    'test');
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MINMAX_DIFF_THRESHOLD_TM', currentMinMaxDiffThresholdTm,       'test');
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MARGIN_SEC',               windowMarginSec,                    'test');
        Bso.make_read_only()

        % CALL TESTED FUNCTION
        actIsSweepingFpa = bicas.proc.L1L2.autodetect_sweeps(tt2000, bdmFpa, hkBiasCurrentFpa, Bso);

        actIsSweeping = actIsSweepingFpa.logical2doubleNan();
        %[actIsSweeping, expIsSweeping]
        testCase.assertEqual(actIsSweeping, expIsSweeping)
      end

      % Test using hard-coded data except for specified arguments.
      %
      % Hardcoded arguments: tt2000, hkBiasCurrent, windowLengthPts
      %
      % ARGUMENTS
      % =========
      % bdm, expIsSweeping
      %       Scalar constant values. FPAs as float/NaN.
      %
      function test2(bdm, sbdaEndTt2000, currentMinMaxDiffThresholdTm, expIsSweeping)
        % NOTE: Only test window lengths equal or shorter than length of
        % data.
        %
        assert(isfloat(bdm          ) && isscalar(bdm          ))
        assert(isfloat(expIsSweeping) && isscalar(expIsSweeping))

        % NOTE: Test window length up until exact number of records.
        for windowLengthPts = [1, 3, 6]
          test(...
            [0:1000:5000]', ...       % tt2000
            [1,1,1,1,1,1]' * bdm, ....
            [                         % hkBiasCurrent
            1, 2, 3; ...
            3, 1, 2; ...
            2, 3, 1; ...
            1, 2, 3; ...
            3, 1, 2; ...
            2, 3, 1; ...
            ], ...
            sbdaEndTt2000, ...
            windowLengthPts, ...      % Iterated over
            currentMinMaxDiffThresholdTm, ...
            0, ...                    % windowMarginSec
            [1 1 1 1 1 1]' * expIsSweeping)
        end
      end

      %===================================================================

      ALL_ENABLED = true;

      if ALL_ENABLED
        % SBDA
        test2(4,  10000, 3, 1)
        test2(4,  10000, 1, 1)

        test2(0,  10000, 3, 0)
        test2(0,  10000, 1, 0)

        % SCDA
        test2(4, -10000, 3, 0)
        test2(4, -10000, 1, 1)

        test2(0, -10000, 3, 0)
        test2(0, -10000, 1, 0)
      end

      if ALL_ENABLED
        % Window longer than data. ==> No sweep.

        % function test(...
        %     tt2000, bdm, hkBiasCurrent, ...
        %     sbdaEndTt2000, windowLengthPts, currentMinMaxDiffThresholdTm, windowMarginSec, ...
        %     expIsSweeping)
        test(...
          [0:1000:2000]', ...
          [4 4 4]', ....
          [ ...
          100, 200, 300; ...
          300, 100, 200; ...
          200, 300, 100; ...
          ], ...
          -10000, ...
          4, ...
          100, ...
          0, ...   % windowMarginSec
          [0, 0, 0]')
      end

      if ALL_ENABLED
        DATA = [ ...
          6, 4,   1, 2, 3,   0; ...
          7, 4,   2, 3, 1,   1; ...
          8, 4,   1, 2, 3+2, 1; ...   % BDM=4. Exceed threshold
          9, 4,   3, 1, 2,   1; ...
          10, 4,   2, 3, 1,   0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5);
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, ...
          bdm, ....
          hkBiasCurrent, ...
          6, ...   % Time threshold
          2, ...   % Window length
          3, ...   % Threshold
          0, ...   % windowMarginSec
          expIsSweeping)
      end

      if ALL_ENABLED
        DATA = [ ...
          10, 4,   2, 3, 1,   0; ...
          11, 0,   1, 2, 3+2, 0; ...   % BDM=0. Exceed threshold. ==> No sweep
          12, 4,   3, 1, 2,   0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5);
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, bdm, hkBiasCurrent, ...
          6, ...   % Time threshold
          2, ...   % Window length
          3, ...   % Threshold
          0, ...   % windowMarginSec
          expIsSweeping)
      end

      if ALL_ENABLED
        % Test window margin
        DATA = [ ...
          1, 4,   2, 3, 1,   0; ...
          2, 4,   1, 2, 3,   1; ...   % Set due to window margin
          3, 4,   3, 1, 2,   1; ...   % Set due to window length
          4, 4,   2, 3, 1-2, 1; ...   % Exceed threshold
          5, 4,   1, 2, 3+2, 1; ...   % Exceed threshold
          6, 4,   3, 1, 2,   1; ...   % Set due to window length
          7, 4,   2, 3, 1,   1; ...   % Set due to window margin
          8, 4,   1, 2, 3,   0; ...
          9, 4,   3, 1, 2,   0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5) + 100;
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, bdm, hkBiasCurrent, ...
          0, ...        % Time threshold
          2, ...        % Window length
          3, ...        % Threshold
          1.1e-9, ...   % windowMarginSec
          expIsSweeping)
      end

      if ALL_ENABLED
        % Complex test
        DATA = [ ...
          1, 0,   2, 3, 1,   0; ...
          2, 4,   1, 2, 3,   1; ...   % BDM=4
          3, 0,   3, 1, 2,   0; ...
          4, 0,   2, 3, 1-2, 0; ...   % Exceed threshold (BDM=0)
          5, 0,   1, 2, 3,   0; ...
          6, 0,   3, 1, 2,   0; ...
          7, 4,   2, 3, 1,   1; ...
          8, 4,   1, 2, 3+2, 1; ...   % BDM=4. Exceed threshold
          9, 4,   3, 1, 2,   1; ...
          10, 4,   2, 3, 1,   0; ...
          11, 0,   1, 2, 3+2, 0; ...   % BDM=0. Exceed threshold
          12, 4,   3, 1, 2,   0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5) + 100;
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, bdm, hkBiasCurrent, ...
          6, ...   % Time threshold
          2, ...   % Window length
          3, ...   % Threshold
          0, ...   % windowMarginSec
          expIsSweeping)
      end

    end    % function



  end    % methods(Test)



end
