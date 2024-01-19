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

      function test(...
          tt2000, bdm, hkBiasCurrentFpa, ...
          bdm4TrickEndTt2000, windowLengthPts, currentMinMaxDiffThresholdTm, windowMarginSec, ...
          expIsSweeping)

        assert(issorted(tt2000), 'ascend')
        assert(isa(expIsSweeping, 'double'))
        bdmFpa           = bicas.utils.FPArray.floatNan2int(bdm, 'uint8');
        hkBiasCurrentFpa = bicas.utils.FPArray.floatNan2int(hkBiasCurrentFpa, 'uint16');

        Bso = bicas.create_default_BSO();
        Bso.override_value('PROCESSING.L2.AUTODETECT_SWEEPS.END_MUX4_TRICK_UTC',              spdfbreakdowntt2000(bdm4TrickEndTt2000), 'test');
        Bso.override_value('PROCESSING.L2.AUTODETECT_SWEEPS.WINDOW_LENGTH_PTS',               windowLengthPts,                         'test');
        Bso.override_value('PROCESSING.L2.AUTODETECT_SWEEPS.WINDOW_MINMAX_DIFF_THRESHOLD_TM', currentMinMaxDiffThresholdTm,            'test');
        Bso.override_value('PROCESSING.L2.AUTODETECT_SWEEPS.WINDOW_MARGIN_SEC',               windowMarginSec,                         'test');
        Bso.make_read_only()

        % CALL TESTED FUNCTION
        actIsSweepingFpa = bicas.proc.L1L2.autodetect_sweeps(tt2000, bdmFpa, hkBiasCurrentFpa, Bso);

        actIsSweeping = actIsSweepingFpa.logical2doubleNan();
        %[actIsSweeping, expIsSweeping]
        testCase.verifyEqual(actIsSweeping, expIsSweeping)
      end

      %===================================================================

      function test2(bdm, bdm4TrickEndTt2000, currentMinMaxDiffThresholdTm, expIsSweeping)
        % NOTE: Only test window lengths equal or shorter than length of
        % data.
        % Constant BDM over time.
        %
        for windowLengthPts = [1, 3, 6]
          test(...
            [0:1000:5000]', ...       % tt2000
            [1,1,1,1,1,1]' * bdm, ....
            [                         % hkBiasCurrentFpa
            1, 2, 3; ...
            3, 1, 2; ...
            2, 3, 1; ...
            1, 2, 3; ...
            3, 1, 2; ...
            2, 3, 1; ...
            ], ...
            bdm4TrickEndTt2000, ...
            windowLengthPts, ...      % Iterated over
            currentMinMaxDiffThresholdTm, ...
            0, ...                    % windowMarginSec
            [1 1 1 1 1 1]' * expIsSweeping)
        end
      end

      ALL_ENABLED = true;

      if ALL_ENABLED
        test2(4,  10000, 3, 1)
        test2(4,  10000, 1, 1)

        test2(0,  10000, 3, 0)
        test2(0,  10000, 1, 0)

        test2(4, -10000, 3, 0)
        test2(4, -10000, 1, 1)

        test2(0, -10000, 3, 0)
        test2(0, -10000, 1, 0)
      end

      if ALL_ENABLED
        % Window longer than data. ==> No sweep.
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
