%
% matlab.unittest automatic test code for bicas.proc.L1L2.Saturation.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Saturation___UTEST < matlab.unittest.TestCase
% PROPOSAL: Split up into multiple files.
%   CON: init_object() is used by tests for multiple methods.
%   CON: Might want to share instance field S = bicas.proc.L1L2.const.C.SSID_DICT;



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    A = bicas.proc.L1L2.const.C.ASID_DICT;
    S = bicas.proc.L1L2.const.C.SSID_DICT;
    BDM0_DLR0_ROW = ["DC_V1", "DC_V12", "DC_V23", "AC_V12", "AC_V23"];
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_VSTB(testCase)

      % NOTE: Changes order of arguments to get easier-to-read hardcoded
      %       calls.
      function test(satArgsCa, ssidStr, isAchg, samplesAVolt, expVstbAr)
        ssid      = testCase.S(ssidStr);
        expVstbAr  = logical(expVstbAr);
        isAchgFpa = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat = testCase.init_object(satArgsCa{:});

        actVstbAr = Sat.get_VSTB(samplesAVolt, ssid, isAchgFpa);

        testCase.assertEqual(actVstbAr, expVstbAr)
      end

      function main()
        for isAchg = [0, 1, NaN]
          % DC single, 1x6
          test(...
            {99, 0.6, 3, 99,99,99}, "DC_V1",  isAchg, ...
            [-5, -1, 0, 1, 5, NaN], ...
            [ 1,  0, 0, 0, 1,   0])

          % DC diff, 3x2
          test(...
            {99, 0.6, 99, 3,99,99}, "DC_V12", isAchg, ...
            [-5, -1, 0; 1, 5, NaN], ...
            [ 1,  0, 0; 0, 1,   0])
        end

        % AC, low gain, 6x1
        for achgThreshold = [1, 7]
          test(...
            {99, 0.6, 99, 99, 3, achgThreshold}, "AC_V12", 0, ...
            [-5, -1, 0, 1, 5, NaN]', ...
            [ 1,  0, 0, 0, 1,   0]')
        end

        % AC, high gain, 1x6
        for aclgThreshold = [1, 7]
          test(...
            {99, 0.6, 99, 99, aclgThreshold,  3}, "AC_V23", 1, ...
            [-5, -1, 0, 1, 5, NaN], ...
            [ 1,  0, 0, 0, 1,   0])
        end

        % AC, unknown gain, 1x6
        for aclgThreshold = [1, 7]
          test({99, 0.6, 99, 99, aclgThreshold,  3}, "AC_V23", NaN, ...
            [-5, -1, 0, 1, 5, NaN], ...
            [ 0,  0, 0, 0, 0,   0])
        end
      end



      main()
    end



    function test_get_snapshot_VSQB(testCase)

      function test(satArgsCa, samplesAVolt, ssidStr, isAchg, expVsqb)
        ssid      = testCase.S(ssidStr);
        isAchgFpa = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat       = testCase.init_object(satArgsCa{:});

        actVsqb = Sat.get_snapshot_VSQB(samplesAVolt, ssid, isAchgFpa);

        testCase.assertEqual(actVsqb, expVsqb)
      end

      function main()
        % IMPLEMENTATION NOTE: Maybe somewhat overkill since
        % test_get_VSTB() tests the different types of thresholds. The
        % extensive tests exist for historical reasons.


        for isAchg = [0, 1, NaN]
          % Empty snapshot.
          test({99, 0.6, 3, 99,99,99}, zeros(1,0), "DC_V1", isAchg, false)
          % Length-one snapshot.
          test({99, 0.6, 3, 99,99,99}, [5],        "DC_V1", isAchg, true)

          % DC single
          test({99, 0.6, 3, 99,99,99}, [ 0, 0], "DC_V1", isAchg, false)
          test({99, 0.6, 3, 99,99,99}, [ 0, 5], "DC_V1", isAchg, false)
          test({99, 0.6, 3, 99,99,99}, [-5, 5], "DC_V1", isAchg, true)

          % DC diff
          test({99, 0.6, 99, 3, 99,99}, [0,  0], "DC_V12", isAchg, false)
          test({99, 0.6, 99, 3, 99,99}, [5, -5], "DC_V13", isAchg, true)
        end

        for unusedGainThreshold = [1, 7]
          % AC, low gain
          test({99, 0.6, 99, 99, 3, unusedGainThreshold}, [0,  0], "AC_V12", 0, false)
          test({99, 0.6, 99, 99, 3, unusedGainThreshold}, [5, -5], "AC_V13", 0, true)

          % AC, high gain
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [0,  0], "AC_V23", 1, false)
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [5, -5], "AC_V23", 1, true)
        end

        % NaN samples
        test({99, 0.6, 3, 99,99,99}, [ NaN, NaN], "DC_V1", 0, false)
        test({99, 0.6, 3, 99,99,99}, [ NaN,   5], "DC_V1", 0, false)
        test({99, 0.6, 3, 99,99,99}, [-5,     5], "DC_V1", 0, true)
      end

      main()
    end



    function test_get_snapshot_VSQB_many(testCase)

      function test(...
          satArgsCa, ...
          zvNValidSamplesPerRecord, samplesAVolt, ...
          ssidStr, isAchg, expVsqbAr)

        ssid      = testCase.S(ssidStr);
        expVsqbAr = logical(expVsqbAr);
        isAchgFpa = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat       = testCase.init_object(satArgsCa{:});

        actVsqbAr = Sat.get_snapshot_VSQB_many(...
          zvNValidSamplesPerRecord, samplesAVolt, ssid, isAchgFpa);

        testCase.assertEqual(actVsqbAr, expVsqbAr)
      end

      function main()
        % IMPLEMENTATION NOTE: Maybe somewhat overkill since
        % test_get_VSTB() tests the different types of thresholds. The
        % extensive tests exist for historical reasons.

        for isAchg = [0, 1, NaN]
          % DC single
          % ---------
          % Zero snapshots
          test({99, 0.6, 99,99,99,99}, zeros(0,1), zeros(0,0),   "DC_V1", isAchg, false(0,1))

          % One snapshot of varying length.
          test({99, 0.6, 3, 99,99,99}, [0], [5, 0, 5], "DC_V1", isAchg, [0])
          test({99, 0.6, 3, 99,99,99}, [1], [5, 0, 5], "DC_V1", isAchg, [1])
          test({99, 0.6, 3, 99,99,99}, [2], [5, 0, 5], "DC_V1", isAchg, [0])
          test({99, 0.6, 3, 99,99,99}, [3], [5, 0, 5], "DC_V1", isAchg, [1])

          % Many snapshots
          test({99, 0.6, 3, 99,99,99}, [0;1;2;3], [5 0 5; 5 0 5; 5 0 5; 5 0 5], "DC_V1", isAchg, [0; 1; 0; 1])

          % DC diff
          test({99, 0.6, 99, 3, 99,99}, 2, [0,  0], "DC_V12", isAchg, [0])
          test({99, 0.6, 99, 3, 99,99}, 2, [5, -5], "DC_V13", isAchg, [1])
        end

        for unusedGainThreshold = [1, 7]
          % AC, low gain
          test({99, 0.6, 99, 99, 3, unusedGainThreshold}, [1; 1], [0  0; 5 -5], "AC_V12", 0, [0; 1])

          % AC, high gain
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [1],    [0,  0],      "AC_V23", 1, [0])
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [1],    [5, -5],      "AC_V23", 1, [1])
        end

        % NaN samples
        test({99, 0.6, 3, 99,99,99}, [2], [ NaN,   5], "DC_V1", 0, false)
        test({99, 0.6, 3, 99,99,99}, [2], [-5,     5], "DC_V1", 0, true)
      end

      main()
    end



    %##############
    % Zero records
    %##############
    function test_get_VSQB___zero_records(testCase)
      for hasSwfFormat = 0:1
        AsrSamplesAVoltSrm = bicas.utils.SameRowsMap("uint8", 0, 'EMPTY');
        AsrSamplesAVoltSrm.add(testCase.A("DC_V1"), zeros(0, 1))
        AsrSamplesAVoltSrm.add(testCase.A("DC_V2"), zeros(0, 1))
        AsrSamplesAVoltSrm.add(testCase.A("DC_V3"), zeros(0, 1))

        testCase.test_get_VSQB(testCase, ...
          cwfSlidingWindowLengthSec = 10, ...
          vstbFractionThreshold     = 0.5, ...
          thresholdAVoltDcSingle    = 1, ...
          thresholdAVoltDcDiff      = 2, ...
          thresholdAVoltAclg        = 3, ...
          thresholdAVoltAchg        = 4, ...
          ...
          tt2000Ar                  = zeros(0, 1), ...
          AsrSamplesAVoltSrm        = AsrSamplesAVoltSrm, ...
          zvNValidSamplesPerRecord  = zeros(0, 1), ...
          bltsSsidAr                = string.empty(0, bicas.const.N_BLTS), ...
          isAchg                    = zeros(0, 1), ...
          hasSwfFormat              = logical(hasSwfFormat), ...
          ...
          expVsqbAr                 = false(0,1));
      end
    end


    %#####
    % CWF
    %#####
    % One subsequence
    % CWF, BDM=0, LRX=1/DC diff.
    function test_get_VSQB___CWF_DC_diff(testCase)

      AsrSamplesAVoltSrm        = bicas.utils.SameRowsMap("uint8", 3, 'EMPTY');
      AsrSamplesAVoltSrm.add(testCase.A("DC_V1"),  [2 0 2]')
      AsrSamplesAVoltSrm.add(testCase.A("DC_V12"), [2 0 2]'+1)
      AsrSamplesAVoltSrm.add(testCase.A("DC_V23"), [2 0 2]'+1)

      testCase.test_get_VSQB(testCase, ...
        cwfSlidingWindowLengthSec = 2.1, ...
        vstbFractionThreshold     = 0.4, ...
        thresholdAVoltDcSingle    = 1, ...
        thresholdAVoltDcDiff      = 2, ...
        thresholdAVoltAclg        = 3, ...
        thresholdAVoltAchg        = 4, ...
        ...
        tt2000Ar                  = [10 11 12]' * 1e9, ...
        AsrSamplesAVoltSrm        = AsrSamplesAVoltSrm, ...
        zvNValidSamplesPerRecord  = [1 1 1]', ...
        bltsSsidAr                = repmat(testCase.BDM0_DLR0_ROW, 3, 1), ...
        isAchg                    = [0 0 0]', ...
        hasSwfFormat              = false, ...
        ...
        expVsqbAr                 = [1 1 1]');
    end



    % One subsequence
    % Separate episodes of saturation on different channels.
    % CWF, BDM=0, LRX=1/AC diff.
    function test_get_VSQB___CWF_AC_diff(testCase)
      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap("uint8", 10, 'EMPTY');
      AsrSamplesAVoltSrm.add(testCase.A("DC_V1"),  [2 2 0 0 0 0 0 0 0 0]')
      AsrSamplesAVoltSrm.add(testCase.A("AC_V12"), [0 0 0 2 2 0 0 0 0 0]'+2)
      AsrSamplesAVoltSrm.add(testCase.A("AC_V23"), [0 0 0 0 0 0 2 2 0 0]'+2)

      testCase.test_get_VSQB(testCase, ...
        cwfSlidingWindowLengthSec = 2.1, ...
        vstbFractionThreshold     = 0.6, ...
        thresholdAVoltDcSingle    = 1, ...
        thresholdAVoltDcDiff      = 2, ...
        thresholdAVoltAclg        = 3, ...
        thresholdAVoltAchg        = 4, ...
        ...
        tt2000Ar                  = [10:19]' * 1e9, ...
        AsrSamplesAVoltSrm        = AsrSamplesAVoltSrm, ...
        zvNValidSamplesPerRecord  = [1 1 1 1 1 1 1 1 1 1]', ...
        bltsSsidAr                = repmat(testCase.BDM0_DLR0_ROW, 10, 1), ...
        isAchg                    = [0 0 0 0 0 0 0 0 0 0]', ...
        hasSwfFormat              = false, ...
        ...
        expVsqbAr                 = [1 1 0 1 1 0 1 1 0 0]');
    end



    % Multiple subsequences
    % Separate episodes of saturation on different channels, with different
    % thresholds, which together combine saturation always set.
    % CWF, BDM=0, LRX=1/AC diff.
    %
    % NOTE: Tests correct function behaviour, but unrealistic input data.
    % Function should never simultaneously receive DC and AC diffs (non-NaN
    % values), but the function is also not aware of LRX and should not care.
    function test_get_VSQB___mult_subsequences(testCase)
      BDM4_DLR0_ROW = ["DC_V1", "DC_V2", "DC_V3", "AC_V12", "AC_V23"];

      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap("uint8", 8, 'EMPTY');
      AsrSamplesAVoltSrm.add(testCase.A("DC_V1"),  [0 0 2 2 0 0 0 0]')
      AsrSamplesAVoltSrm.add(testCase.A("DC_V12"), [0 0 0 0 0 0 3 3]')
      AsrSamplesAVoltSrm.add(testCase.A("DC_V23"), [0 0 0 0 0 0 0 0]')
      AsrSamplesAVoltSrm.add(testCase.A("DC_V2"),  [3 3 0 0 0 0 0 0]')
      AsrSamplesAVoltSrm.add(testCase.A("DC_V3"),  [0 0 0 0 2 2 0 0]')

      testCase.test_get_VSQB(testCase, ...
        cwfSlidingWindowLengthSec = 2.1, ...
        vstbFractionThreshold     = 0.6, ...
        thresholdAVoltDcSingle    = 1, ...
        thresholdAVoltDcDiff      = 2, ...
        thresholdAVoltAclg        = 3, ...
        thresholdAVoltAchg        = 4, ...
        ...
        tt2000Ar                  = [10:17]' * 1e9, ...
        AsrSamplesAVoltSrm        = AsrSamplesAVoltSrm, ...
        zvNValidSamplesPerRecord  = [1 1 1 1 1 1 1 1]', ...
        bltsSsidAr                = [
          repmat(testCase.BDM0_DLR0_ROW, 4, 1); ...
          repmat(         BDM4_DLR0_ROW, 4, 1) ...
        ], ...
        isAchg                    = [0 0 0 0 0 0 0 0]', ...
        hasSwfFormat              = false, ...
        ...
        expVsqbAr                 = [1 1 1 1 1 1 1 1]');
    end



    %#####
    % SWF
    %#####
    % One subsequence
    % SWF, BDM=0, LRX=0/AC diff
    % Saturated DC single
    function test_get_VSQB___SWF_1(testCase)
      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap("uint8", 2, 'EMPTY');
      AsrSamplesAVoltSrm.add(testCase.A("DC_V1"), ...
        [
        [2 2 2 0 0 0 0] + 0; ...   % Saturated
        [0 0 0 2 2 2 2] + 0;
        ]);
      AsrSamplesAVoltSrm.add(testCase.A("AC_V12"), ...
        [
        [0 0 0 0 0 0 0] + 4; ...
        [0 0 0 2 2 2 2] + 6;
        ]+1);
      AsrSamplesAVoltSrm.add(testCase.A("AC_V23"), ...
        [
        [0 0 0 0 0 0 0] + 4; ...
        [0 0 0 2 2 2 2] + 6;
        ]+1);

      testCase.test_get_VSQB(testCase, ...
        cwfSlidingWindowLengthSec = 2.1, ...   % Should be irrelevant since SWF.
        vstbFractionThreshold     = 0.5, ...
        thresholdAVoltDcSingle    = 1, ...
        thresholdAVoltDcDiff      = 3, ...
        thresholdAVoltAclg        = 5, ...
        thresholdAVoltAchg        = 7, ...
        ...
        tt2000Ar                  = [10 11]' * 1e9, ...
        AsrSamplesAVoltSrm        = AsrSamplesAVoltSrm, ...
        zvNValidSamplesPerRecord  = [5 5]', ...
        bltsSsidAr                = repmat(testCase.BDM0_DLR0_ROW, 2, 1), ...
        isAchg                    = [0 1]', ...
        hasSwfFormat              = true, ...
        ...
        expVsqbAr                 = [1 0]');
    end



    % One subsequence
    % SWF, BDM=0, LRX=0/AC diff, LG+HG.
    % Varying snapshot lengths.
    % Saturated AC diffs
    function test_get_VSQB___SWF_2(testCase)

      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap("uint8", 4, 'EMPTY');
      AsrSamplesAVoltSrm.add(testCase.A("DC_V1"), ...
        [
        [0 0 0 0 0 0 0] + 0; ...
        [0 2 2 2 2 2 2] + 0; ...
        [0 0 0 0 0 0 0] + 0; ...
        [0 0 0 2 2 2 2] + 0;
        ])
      AsrSamplesAVoltSrm.add(testCase.A("AC_V12"), ...
        [
        [2 0 0 0 0 0 0] + 4; ...   % Saturated.
        [0 2 2 2 2 2 2] + 4; ...
        [2 2 2 0 0 0 0] + 6; ...   % Saturated.
        [0 0 0 2 2 2 2] + 6;
        ]+1)
      AsrSamplesAVoltSrm.add(testCase.A("AC_V23"), ...
        [
        [0 0 0 0 0 0 0] + 4; ...
        [0 2 2 2 2 2 2] + 4; ...
        [0 0 0 0 0 0 0] + 6; ...
        [0 0 0 2 2 2 2] + 6;
        ]+1)

      testCase.test_get_VSQB(testCase, ...
        cwfSlidingWindowLengthSec = 2.1, ...   % Should be irrelevant since SWF.
        vstbFractionThreshold     = 0.5, ...
        thresholdAVoltDcSingle    = 1, ...
        thresholdAVoltDcDiff      = 3, ...
        thresholdAVoltAclg        = 5, ...
        thresholdAVoltAchg        = 7, ...
        ...
        tt2000Ar                  = [10:13]' * 1e9, ...
        AsrSamplesAVoltSrm        = AsrSamplesAVoltSrm, ...
        zvNValidSamplesPerRecord  = [1 1 3 3]', ...
        bltsSsidAr                = repmat(testCase.BDM0_DLR0_ROW, 4, 1), ...
        isAchg                    = [0 0 1 1]', ...
        hasSwfFormat              = true, ...
        ...
        expVsqbAr                 = [1 0 1 0]');
    end



    function test_get_ASR_CWF_channel_VSTB___DC(testCase)
      testCase.test_get_ASR_CWF_channel_VSTB(testCase, ...
        thresholdAVoltDcSingle = 1, ...
        thresholdAVoltDcDiff   = 3, ...
        thresholdAVoltAclg     = 5, ...
        thresholdAVoltAchg     = 7, ...
        ssidStr                = "DC_V1", ...
        hasSwfFormat           = false, ...
        isAchg                 = [0,   0]', ...
        samplesAVolt           = [0.9, 1.1, ]', ...
        expVstbAr              = [0,   1]')
    end



    function test_get_ASR_CWF_channel_VSTB___AC(testCase)
      testCase.test_get_ASR_CWF_channel_VSTB(testCase, ...
        thresholdAVoltDcSingle = 1, ...
        thresholdAVoltDcDiff   = 3, ...
        thresholdAVoltAclg     = 5, ...
        thresholdAVoltAchg     = 7, ...
        ssidStr                = "AC_V12", ...
        hasSwfFormat           = false, ...
        isAchg                 = [0,   0,   1,   1  ]', ...
        samplesAVolt           = [4.9, 5.1, 6.9, 7.1]', ...
        expVstbAr              = [0,   1,   0,   1]')
    end



    function test_get_ASR_SWF_channel_VSQB___AC(testCase)
      testCase.test_get_ASR_SWF_channel_VSQB(testCase, ...
        vstbFractionThreshold    = 0.4, ...
        thresholdAVoltDcSingle   = 1, ...
        thresholdAVoltDcDiff     = 3, ...
        thresholdAVoltAclg       = 5, ...
        thresholdAVoltAchg       = 7, ...
        ssidStr                  = "AC_V12", ...
        hasSwfFormat             = true, ...
        isAchg                   = [0; 0; 1; 1], ...
        samplesAVolt             = [
        4.9, 4.9, 9, 9;
        4.9, 5.1, 9, 9;
        6.9, 7.1, 0, 0;
        7.1, 7.1, 0, 0;
        ], ...
        zvNValidSamplesPerRecord = [2; 2; 2; 2], ...
        expVstbAr                = [0; 1; 1; 1])
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function test_get_VSQB(testCase, V)
      arguments
        testCase
        V.cwfSlidingWindowLengthSec
        V.vstbFractionThreshold
        V.thresholdAVoltDcSingle
        V.thresholdAVoltDcDiff
        V.thresholdAVoltAclg
        V.thresholdAVoltAchg
        %
        V.tt2000Ar
        V.AsrSamplesAVoltSrm
        V.zvNValidSamplesPerRecord
        V.bltsSsidAr
        V.isAchg
        V.hasSwfFormat
        %
        V.expVsqbAr
      end

      % Modify/normalize arguments.
      V.tt2000Ar     = int64(V.tt2000Ar);
      V.bltsSsidAr   = testCase.S(V.bltsSsidAr);
      V.isAchgFpa    = bicas.utils.FPArray(logical(V.isAchg), 'FILL_POSITIONS', isnan(V.isAchg));
      V.hasSwfFormat = logical(V.hasSwfFormat);

      V.expVsqbAr = logical(V.expVsqbAr);

      L = bicas.Logger('NO_STDOUT', false);
      Sat = testCase.init_object(...
        V.cwfSlidingWindowLengthSec, ...
        V.vstbFractionThreshold, ...
        V.thresholdAVoltDcSingle, ...
        V.thresholdAVoltDcDiff, ...
        V.thresholdAVoltAclg, ...
        V.thresholdAVoltAchg);

      % CALL FUNCTION
      actVsqbAr = Sat.get_VSQB(...
        V.tt2000Ar, V.AsrSamplesAVoltSrm, V.zvNValidSamplesPerRecord, ...
        V.bltsSsidAr, V.isAchgFpa, V.hasSwfFormat, L);

      testCase.assertEqual(actVsqbAr, V.expVsqbAr)
    end



    function test_get_ASR_CWF_channel_VSTB(testCase, V)
      arguments
        testCase
        V.cwfSlidingWindowLengthSec
        V.vstbFractionThreshold
        V.thresholdAVoltDcSingle
        V.thresholdAVoltDcDiff
        V.thresholdAVoltAclg
        V.thresholdAVoltAchg
        V.ssidStr
        V.hasSwfFormat
        V.isAchg
        V.samplesAVolt
        V.expVstbAr
      end

      % Should be irrelevant for this test.
      cwfSlidingWindowLengthSec = 99;
      vstbFractionThreshold     = 0.5;

      % Modify/normalize arguments.
      V.isAchgFpa = bicas.utils.FPArray(logical(V.isAchg), 'FILL_POSITIONS', isnan(V.isAchg));
      V.expVstbAr  = logical(V.expVstbAr);

      Sat = testCase.init_object(...
        cwfSlidingWindowLengthSec, ...
        vstbFractionThreshold, ...
        V.thresholdAVoltDcSingle, ...
        V.thresholdAVoltDcDiff, ...
        V.thresholdAVoltAclg, ...
        V.thresholdAVoltAchg);

      actVstbAr = Sat.get_ASR_CWF_channel_VSTB(...
        testCase.S(V.ssidStr), V.isAchgFpa, V.samplesAVolt);

      testCase.assertEqual(actVstbAr, V.expVstbAr)
    end



    function test_get_ASR_SWF_channel_VSQB(testCase, V)
      arguments
        testCase
        V.cwfSlidingWindowLengthSec
        V.vstbFractionThreshold
        V.thresholdAVoltDcSingle
        V.thresholdAVoltDcDiff
        V.thresholdAVoltAclg
        V.thresholdAVoltAchg
        V.ssidStr
        V.hasSwfFormat
        V.isAchg
        V.samplesAVolt
        V.zvNValidSamplesPerRecord
        V.expVstbAr
      end

      % Should be irrelevant for this test.
      cwfSlidingWindowLengthSec = 99;

      % Modify/normalize arguments.
      V.isAchgFpa = bicas.utils.FPArray(logical(V.isAchg), 'FILL_POSITIONS', isnan(V.isAchg));
      V.expVstbAr = logical(V.expVstbAr);

      Sat = testCase.init_object(...
        cwfSlidingWindowLengthSec, ...
        V.vstbFractionThreshold, ...
        V.thresholdAVoltDcSingle, ...
        V.thresholdAVoltDcDiff, ...
        V.thresholdAVoltAclg, ...
        V.thresholdAVoltAchg);

      actVstbAr = Sat.get_ASR_SWF_channel_VSQB(...
        bicas.proc.L1L2.const.C.SSID_DICT(V.ssidStr), V.isAchgFpa, ...
        V.samplesAVolt, V.zvNValidSamplesPerRecord);

      testCase.assertEqual(actVstbAr, V.expVstbAr)
    end



    function S = init_object(...
        cwfSlidingWindowLengthSec, ...
        vstbFractionThreshold, ...
        thresholdAVoltDcSingle, ...
        thresholdAVoltDcDiff, ...
        thresholdAVoltAclg, ...
        thresholdAVoltAchg)

      Bso = bicas.create_default_BSO();
      Bso.override_value('PROCESSING.SATURATION.CWF_SLIDING_WINDOW_LENGTH_SEC',            cwfSlidingWindowLengthSec, 'test');
      Bso.override_value('PROCESSING.SATURATION.VSTB_FRACTION_THRESHOLD',                  vstbFractionThreshold,     'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.SINGLE',         thresholdAVoltDcSingle,    'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.DIFF',           thresholdAVoltDcDiff,      'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.LOW_GAIN',  thresholdAVoltAclg,        'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.HIGH_GAIN', thresholdAVoltAchg,        'test');
      Bso.make_read_only();

      S = bicas.proc.L1L2.Saturation(Bso);
    end



  end    % methods(Static, Access=private)



end
