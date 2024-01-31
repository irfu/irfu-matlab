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
          sbdaEndTt2000, windowLengthPts, currentMmDiffMinTm, windowMarginSec, ...
          expIsSweeping)

        assert(issorted(tt2000), 'ascend')
        assert(isa(expIsSweeping, 'double'))
        bdmFpa           = bicas.utils.FPArray.floatNan2int(bdm, 'uint8');
        hkBiasCurrentFpa = bicas.utils.FPArray.floatNan2int(hkBiasCurrent, 'uint16');

        Bso = bicas.create_default_BSO();
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC',                       spdfbreakdowntt2000(sbdaEndTt2000), 'test');
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS',             windowLengthPts,                    'test');
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MINMAX_DIFF_MINIMUM_TM', currentMmDiffMinTm,                 'test');
        Bso.override_value('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MARGIN_SEC',             windowMarginSec,                    'test');
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
      function test2(bdm, sbdaEndTt2000, currentMmDiffMinTm, expIsSweeping)
        % NOTE: Only test window lengths equal or shorter than length of
        % data.
        %
        assert(isfloat(bdm          ) && isscalar(bdm          ))
        assert(isfloat(expIsSweeping) && isscalar(expIsSweeping))

        % NOTE: Test window length up until exact number of records.
        for windowLengthPts = [2, 5, 6]
          test(...
            [0:1000:5000]', ...       % tt2000
            [1,1,1,1,1,1]' * bdm, ....
            [                         % hkBiasCurrent. Always MM diff == 2, for all windows.
            1, 2, 3; ...
            1, 2, 5; ...
            1, 4, 5; ...
            3, 4, 5; ...
            3, 4, 3; ...
            3, 2, 3; ...
            ], ...
            sbdaEndTt2000, ...
            windowLengthPts, ...      % Iterated over
            currentMmDiffMinTm, ...
            0, ...                    % windowMarginSec
            [1 1 1 1 1 1]' * expIsSweeping)
        end
      end

      %===================================================================

      ALL_ENABLED = true;

      if ALL_ENABLED
        % Window short or equal to number of records.
        % function test2(bdm, sbdaEndTt2000, currentMmDiffMinTm, expIsSweeping)

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
        % Window longer than data. ==> No SCDA sweep.

        % function test(...
        %     tt2000, bdm, hkBiasCurrent, ...
        %     sbdaEndTt2000, windowLengthPts, currentMmDiffMinTm, windowMarginSec, ...
        %     expIsSweeping)
        test(...
          [0:1000:2000]', ...
          [4 4 4]', ....
          [ ...
          100, 200, 300; ...
          300, 100, 200; ...
          200, 300, 100; ...
          ], ...
          -10000, ...   % Ensure SCDA is used.
          4, ...
          100, ...      % currentMmDiffMinTm
          0, ...        % windowMarginSec
          [0, 0, 0]')
      end

      if ALL_ENABLED
        % All SCDA, BDM=4. Windows triggered for separate bias currents.
        % MM diffs between channels (should not be relevant) > MM diff minimum
        DATA = [ ...
          6,  4,   1, 4, 7,     0; ...
          7,  4,   1, 4, 7,     1; ...
          8,  4,   4, 4, 7,     1; ...
          9,  4,   4, 4, 7,     0; ...
          10, 4,   4, 4, 7,     1; ...
          11, 4,   4, 1, 7,     1; ...
          12, 4,   4, 1, 7,     0; ...
          13, 4,   4, 1, 7,     1; ...
          14, 4,   4, 1, 4,     1; ...
          15, 4,   4, 1, 4,     0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5);
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, ...
          bdm, ....
          hkBiasCurrent, ...
          0, ...   % Time threshold. Ensure SCDA.
          2, ...   % Window length
          3, ...   % MM diff minimum
          0, ...   % windowMarginSec
          expIsSweeping)
      end

      if ALL_ENABLED
        % SCDA. BDM=0 (non-sweep). Exceed threshold. ==> Still no sweep
        DATA = [ ...
          10, 4,   2, 3, 1,   0; ...
          11, 0,   2, 3, 5,   0; ...
          12, 4,   3, 7, 5,   0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5);
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, bdm, hkBiasCurrent, ...
          6, ...   % Time threshold
          2, ...   % Window length
          3, ...   % MM diff minimum
          0, ...   % windowMarginSec
          expIsSweeping)
      end

      if ALL_ENABLED
        % SCDA, BDM=4.
        % Test window margin
        DATA = [ ...
          1, 4,   1, 4, 7,   0; ...
          2, 4,   1, 4, 7,   1; ...   % Set due to window margin
          3, 4,   1, 4, 7,   1; ...   % Set due to window length
          4, 4,   1, 4, 7,   1; ...   % Exceed MM diff cf next
          5, 4,   1, 6, 7,   1; ...   % Exceed MM diff cf previous
          6, 4,   1, 6, 7,   1; ...   % Set due to window length
          7, 4,   1, 6, 7,   1; ...   % Set due to window margin
          8, 4,   1, 6, 7,   0; ...
          9, 4,   1, 6, 7,   0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5) + 100;
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, bdm, hkBiasCurrent, ...
          0, ...        % Time threshold
          3, ...        % Window length
          2, ...        % MM diff minimum
          1.1e-9, ...   % windowMarginSec
          expIsSweeping)
      end

      if ALL_ENABLED
        % Complex test
        DATA = [ ...
          1,  0,   1, 2, 3,     0; ...
          2,  4,   3, 0, 5,     1; ...   % BDM=4
          3,  0,   1, 2, 3,     0; ...
          4,  0,   1, 2, 6,     0; ...   % Exceed threshold (BDM=0)
          5,  0,   1, 2, 5,     0; ...
          6,  0,   1, 2, 4,     0; ...
          7,  4,   1, 2, 3,     1; ...
          8,  4,   1, 2, 6,     1; ...   % BDM=4. Exceed threshold
          9,  4,   1, 2, 5,     0; ...
          10, 4,   1, 2, 4,     0; ...
          11, 0,   1, 2, 1,     0; ...   % BDM=0. Exceed threshold
          12, 4,   1, 2, 4,     0; ...
          ];
        tt2000        = int64(DATA(:, 1));
        bdm           =       DATA(:, 2);
        hkBiasCurrent =       DATA(:, 3:5) + 100;
        expIsSweeping =       DATA(:, 6);
        test(...
          tt2000, bdm, hkBiasCurrent, ...
          6, ...   % Time threshold
          2, ...   % Window length
          3, ...   % MM diff minimum
          0, ...   % windowMarginSec
          expIsSweeping)
      end

    end    % function



  end    % methods(Test)



end
