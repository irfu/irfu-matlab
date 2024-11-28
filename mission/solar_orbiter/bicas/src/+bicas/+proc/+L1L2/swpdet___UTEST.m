%
% matlab.unittest automatic test code for bicas.proc.L1L2.swpdet.
%
% NOTE: This is effectively mostly a test on
% bicas.proc.L1L2.swpdet.SBDA_SCDA_with_margins() only, for historical reasons.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef swpdet___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Abbreviation for currentMmDiffMinTm.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_window_shorter_equal_than_data(testCase)

      % Test using hard-coded data except for specified arguments.
      % Hardcoded arguments: tt2000, hkBiasCurrent, windowLengthHkCdfRecords
      %
      % NOTE: Only test window lengths equal or shorter than length of data.
      %
      % ARGUMENTS
      % =========
      % bdm, expIsSweeping
      %       Scalar constant values. FPAs as float/NaN.
      %
      function test(bdm, sbdaEndTt2000, currentMmDiffMinTm, expIsSweeping)
        assert(isfloat(bdm          ) && isscalar(bdm          ))
        assert(isfloat(expIsSweeping) && isscalar(expIsSweeping))

        % NOTE: Test window length up until exact number of records.
        for windowLengthHkCdfRecords = [2, 5, 6]
          % NOTE: Always MM diff == 2, for all windows.
          DATA = [...
            0,      1, 2, 3; ...
            1000,   1, 2, 5; ...
            2000,   1, 4, 5; ...
            3000,   3, 4, 5; ...
            4000,   3, 4, 3; ...
            5000,   3, 2, 3; ...
            ];

          S.tt2000                   = DATA(:, 1);
          S.bdm                      = [1,1,1,1,1,1]' * bdm;
          S.hkBiasCurrent            = DATA(:, 2:4);
          S.sbdaEndTt2000            = sbdaEndTt2000;
          S.windowLengthHkCdfRecords = windowLengthHkCdfRecords;
          S.currentMmDiffMinTm       = currentMmDiffMinTm;
          S.windowMarginSec          = 0;
          S.expIsSweeping            = [1 1 1 1 1 1]' * expIsSweeping;

          % Assertion on test data.
          mmDiff = max(S.hkBiasCurrent, [], 1) - min(S.hkBiasCurrent, [], 1);
          assert(all(mmDiff <= 2))

          testCase.test(S)
        end    % for
      end    % function

      %===================================================================

      % Window short or equal to number of records.
      % function test2(bdm, sbdaEndTt2000, currentMmDiffMinTm, expIsSweeping)

      % SBDA
      test(4,  10000, 3, 1)
      test(4,  10000, 1, 1)

      test(0,  10000, 3, 0)
      test(0,  10000, 1, 0)

      % SCDA
      test(4, -10000, 3, 0)
      test(4, -10000, 1, 1)

      test(0, -10000, 3, 0)
      test(0, -10000, 1, 0)
    end



    function test_window_longer_than_data(testCase)
      % Window longer than data. ==> No SCDA sweep.
      DATA = [ ...
        0,    4,   100, 200, 300,   0;
        1000, 4,   300, 100, 200,   0;
        2000, 4,   200, 300, 100,   0;
        ];
      S.tt2000                   = int64(DATA(:, 1));
      S.bdm                      =       DATA(:, 2);
      S.hkBiasCurrent            =       DATA(:, 3:5);
      S.sbdaEndTt2000            = -10000;
      S.windowLengthHkCdfRecords = 4;
      S.currentMmDiffMinTm       = 100;
      S.windowMarginSec          = 0;
      S.expIsSweeping            =       DATA(:, 6);
      testCase.test(S)
    end



    % All SCDA, BDM=4. Windows triggered for separate bias currents.
    % MM diffs between channels (should not be relevant) > MM diff minimum
    function test_SCDA_BDM4_trigger_from_any_channel(testCase)
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

      S.tt2000                   = int64(DATA(:, 1));
      S.bdm                      =       DATA(:, 2);
      S.hkBiasCurrent            =       DATA(:, 3:5);
      S.sbdaEndTt2000            = 0;   % Time threshold. Ensure SCDA.
      S.windowLengthHkCdfRecords = 2;
      S.currentMmDiffMinTm       = 3;
      S.windowMarginSec          = 0;
      S.expIsSweeping            =       DATA(:, 6);
      testCase.test(S)
    end



    % SCDA. BDM=0 (non-sweep). Exceed threshold. ==> Still no sweep
    function test_SCDA_BDM0_exceed_MM_diff(testCase)
      DATA = [ ...
        10, 4,   2, 3, 1,   0; ...
        11, 0,   2, 3, 5,   0; ...
        12, 4,   3, 7, 5,   0; ...
        ];
      S.tt2000                   = int64(DATA(:, 1));
      S.bdm                      =       DATA(:, 2);
      S.hkBiasCurrent            =       DATA(:, 3:5);
      S.sbdaEndTt2000            = 6;
      S.windowLengthHkCdfRecords = 2;
      S.currentMmDiffMinTm       = 3;
      S.windowMarginSec          = 0;
      S.expIsSweeping            =       DATA(:, 6);
      testCase.test(S)
    end



    % SCDA, BDM=4.
    % Test window margin.
    function test_SCDA_BDM4_window_margin(testCase)
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
      S.tt2000                   = int64(DATA(:, 1));
      S.bdm                      =       DATA(:, 2);
      S.hkBiasCurrent            =       DATA(:, 3:5) + 100;
      S.sbdaEndTt2000            = 0;
      S.windowLengthHkCdfRecords = 3;
      S.currentMmDiffMinTm       = 2;
      S.windowMarginSec          = 1.1e-9;
      S.expIsSweeping            =       DATA(:, 6);
      testCase.test(S)
    end



    function test_complex(testCase)
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
      S.tt2000                   = int64(DATA(:, 1));
      S.bdm                      =       DATA(:, 2);
      S.hkBiasCurrent            =       DATA(:, 3:5) + 100;
      S.sbdaEndTt2000            = 6;    % In the middle of data.
      S.windowLengthHkCdfRecords = 2;
      S.currentMmDiffMinTm       = 3;
      S.windowMarginSec          = 0;
      S.expIsSweeping            =       DATA(:, 6);
      testCase.test(S)
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Generic test. All input and expected output is specified in arguments.
    %
    % IMPLEMENTATION NOTE: All arguments are passed in a single struct to
    % force the caller to explicitly name the separate arguments. ==> Saves
    % confusion and comments.
    %
    %
    % ARGUMENTS
    % =========
    % S
    %       Struct.
    %       .expIsSweeping
    %           logical FPA as double/NaN
    %       .bdm, S.hkBiasCurrent
    %           FPAs as float/NaN
    %       etc.
    %
    function test(testCase, S)
      assert(issorted(S.tt2000), 'ascend')
      assert(isa(S.expIsSweeping, 'double'))

      bdmFpa           = bicas.utils.FPArray.floatNan2int(S.bdm, 'uint8');
      hkBiasCurrentFpa = bicas.utils.FPArray.floatNan2int(S.hkBiasCurrent, 'uint16');

      Bso = bicas.create_default_BSO();
      Bso.override_value('PROCESSING.L2.SWEEP_DETECTION.SBDA_SCDA_BOUNDARY_UTC', spdfbreakdowntt2000(S.sbdaEndTt2000),           'test');
      Bso.override_value('PROCESSING.L2.SWEEP_DETECTION.SCDA.WINDOW_LENGTH_HK_CDF_RECORDS',          S.windowLengthHkCdfRecords, 'test');
      Bso.override_value('PROCESSING.L2.SWEEP_DETECTION.SCDA.WINDOW_MINMAX_DIFF_MINIMUM_TM',         S.currentMmDiffMinTm,       'test');
      Bso.override_value('PROCESSING.L2.SWEEP_DETECTION.SCDA.WINDOW_MARGIN_SEC',                     S.windowMarginSec,          'test');
      Bso.make_read_only()

      % CALL TESTED FUNCTION
      actIsSweepingFpa = bicas.proc.L1L2.swpdet.SBDA_SCDA_with_margins(S.tt2000, bdmFpa, hkBiasCurrentFpa, Bso);

      actIsSweeping = actIsSweepingFpa.logical2doubleNan();
      %[actIsSweeping, expIsSweeping]
      testCase.assertEqual(actIsSweeping, S.expIsSweeping)
    end



  end    % methods(Static, Access=private)



end
