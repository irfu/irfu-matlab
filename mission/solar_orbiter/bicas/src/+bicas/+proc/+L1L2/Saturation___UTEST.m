%
% matlab.unittest automatic test code for bicas.proc.L1L2.Saturation.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Saturation___UTEST < matlab.unittest.TestCase
% PROPOSAL: Split up into multiple files.
%   CON: init_object() is used by tests for multiple methods.
%   CON: Might want to share isntance field S = bicas.sconst.C.S_SSID_DICT;



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    S = bicas.sconst.C.S_SSID_DICT;
    BDM0_DLR0_ROW = ["DC_V1", "DC_V12", "DC_V23", "AC_V12", "AC_V23"];
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_TSF(testCase)

      % NOTE: Changes order of arguments to get easier-to-read hardcoded
      %       calls.
      function test(satArgsCa, sSsid, isAchg, samplesAVolt, expTsfAr)
        Ssid      = testCase.S(sSsid);
        expTsfAr  = logical(expTsfAr);
        isAchgFpa = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat = testCase.init_object(satArgsCa{:});

        actTsfAr = Sat.get_TSF(samplesAVolt, Ssid, isAchgFpa);

        testCase.assertEqual(actTsfAr, expTsfAr)
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



    function test_get_snapshot_saturation(testCase)

      function test(satArgsCa, samplesAVolt, sSsid, isAchg, expIsSaturated)
        Ssid      = testCase.S(sSsid);
        isAchgFpa = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat       = testCase.init_object(satArgsCa{:});

        actIsSaturated = Sat.get_snapshot_saturation(samplesAVolt, Ssid, isAchgFpa);

        testCase.assertEqual(actIsSaturated, expIsSaturated)
      end

      function main()
        % IMPLEMENTATION NOTE: Maybe somewhat overkill since
        % test_get_TSF() tests the different types of thresholds. The
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



    function test_get_snapshot_saturation_many(testCase)

      function test(...
          satArgsCa, ...
          zvNValidSamplesPerRecord, samplesAVolt, ...
          sSsid, isAchg, expIsSaturatedAr)

        Ssid             = testCase.S(sSsid);
        expIsSaturatedAr = logical(expIsSaturatedAr);
        isAchgFpa        = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat              = testCase.init_object(satArgsCa{:});

        actIsSaturatedAr = Sat.get_snapshot_saturation_many(...
          zvNValidSamplesPerRecord, samplesAVolt, Ssid, isAchgFpa);

        testCase.assertEqual(actIsSaturatedAr, expIsSaturatedAr)
      end

      function main()
        % IMPLEMENTATION NOTE: Maybe somewhat overkill since
        % test_get_TSF() tests the different types of thresholds. The
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
    function test_get_voltage_saturation_quality_bit___zero_records(testCase)
      for hasSwfFormat = 0:1
        V = [];
        V.cwfSlidingWindowLengthSec = 10;
        V.tsfFractionThreshold      = 0.5;
        V.thresholdAVoltDcSingle    = 1;
        V.thresholdAVoltDcDiff      = 2;
        V.thresholdAVoltAclg        = 3;
        V.thresholdAVoltAchg        = 4;
        %
        V.tt2000Ar                  = zeros(0, 1);
        V.AsrSamplesAVoltSrm        = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", 0, 'EMPTY');
        V.AsrSamplesAVoltSrm.add(testCase.S("DC_V1").Asid, zeros(0, 1))
        V.AsrSamplesAVoltSrm.add(testCase.S("DC_V2").Asid, zeros(0, 1))
        V.AsrSamplesAVoltSrm.add(testCase.S("DC_V3").Asid, zeros(0, 1))
        V.zvNValidSamplesPerRecord  = zeros(0, 1);
        V.bltsKSsidAr               = string.empty(0, bicas.const.N_BLTS);
        V.isAchg                    = zeros(0, 1);
        V.hasSwfFormat              = logical(hasSwfFormat);
        %
        V.expIsSaturatedAr          = false(0,1);

        testCase.test_get_voltage_saturation_quality_bit(testCase, V);
      end
    end


    %#####
    % CWF
    %#####
    % One subsequence
    % CWF, BDM=0, LRX=1/DC diff.
    function test_get_voltage_saturation_quality_bit___CWF_DC_diff(testCase)
      V = [];
      V.cwfSlidingWindowLengthSec = 2.1;
      V.tsfFractionThreshold      = 0.4;
      V.thresholdAVoltDcSingle    = 1;
      V.thresholdAVoltDcDiff      = 2;
      V.thresholdAVoltAclg        = 3;
      V.thresholdAVoltAchg        = 4;
      %
      V.tt2000Ar                  = [10 11 12]' * 1e9;
      V.AsrSamplesAVoltSrm        = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", 3, 'EMPTY');
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V1").Asid,  [2 0 2]')
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V12").Asid, [2 0 2]'+1)
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V23").Asid, [2 0 2]'+1)
      V.zvNValidSamplesPerRecord  = [1 1 1]';
      V.bltsKSsidAr               = repmat(testCase.BDM0_DLR0_ROW, 3, 1);
      V.isAchg                    = [0 0 0]';
      V.hasSwfFormat              = false;
      %
      V.expIsSaturatedAr          = [1 1 1]';

      testCase.test_get_voltage_saturation_quality_bit(testCase, V);
    end



    % One subsequence
    % Separate episodes of saturation on different channels.
    % CWF, BDM=0, LRX=1/AC diff.
    function test_get_voltage_saturation_quality_bit___CWF_AC_diff(testCase)
      V = [];
      V.cwfSlidingWindowLengthSec = 2.1;
      V.tsfFractionThreshold      = 0.6;
      V.thresholdAVoltDcSingle    = 1;
      V.thresholdAVoltDcDiff      = 2;
      V.thresholdAVoltAclg        = 3;
      V.thresholdAVoltAchg        = 4;
      %
      V.tt2000Ar                  = [10:19]' * 1e9;
      V.AsrSamplesAVoltSrm        = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", 10, 'EMPTY');
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V1").Asid,  [2 2 0 0 0 0 0 0 0 0]')
      V.AsrSamplesAVoltSrm.add(testCase.S("AC_V12").Asid, [0 0 0 2 2 0 0 0 0 0]'+2)
      V.AsrSamplesAVoltSrm.add(testCase.S("AC_V23").Asid, [0 0 0 0 0 0 2 2 0 0]'+2)
      V.zvNValidSamplesPerRecord  = [1 1 1 1 1 1 1 1 1 1]';
      V.bltsKSsidAr               = repmat(testCase.BDM0_DLR0_ROW, 10, 1);
      V.isAchg                    = [0 0 0 0 0 0 0 0 0 0]';
      V.hasSwfFormat              = false;
      %
      V.expIsSaturatedAr          = [1 1 0 1 1 0 1 1 0 0]';

      testCase.test_get_voltage_saturation_quality_bit(testCase, V);
    end



    % Multiple subsequences
    % Separate episodes of saturation on different channels, with different
    % thresholds, which together combine saturation always set.
    % CWF, BDM=0, LRX=1/AC diff.
    %
    % NOTE: Tests correct function behaviour, but unrealistic input data.
    % Function should never simultaneously receive DC and AC diffs (non-NaN
    % values), but the function is also not aware of LRX and should not care.
    function test_get_voltage_saturation_quality_bit___mult_subsequences(testCase)
      BDM4_DLR0_ROW = ["DC_V1", "DC_V2", "DC_V3", "AC_V12", "AC_V23"];

      V = [];
      V.cwfSlidingWindowLengthSec = 2.1;
      V.tsfFractionThreshold      = 0.6;
      V.thresholdAVoltDcSingle    = 1;
      V.thresholdAVoltDcDiff      = 2;
      V.thresholdAVoltAclg        = 3;
      V.thresholdAVoltAchg        = 4;
      %
      V.tt2000Ar                  = [10:17]' * 1e9;
      V.AsrSamplesAVoltSrm        = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", 8, 'EMPTY');
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V1").Asid,  [0 0 2 2 0 0 0 0]')
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V12").Asid, [0 0 0 0 0 0 3 3]')
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V23").Asid, [0 0 0 0 0 0 0 0]')
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V2").Asid,  [3 3 0 0 0 0 0 0]')
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V3").Asid,  [0 0 0 0 2 2 0 0]')
      V.zvNValidSamplesPerRecord  = [1 1 1 1 1 1 1 1]';
      V.bltsKSsidAr               = [
        repmat(testCase.BDM0_DLR0_ROW, 4, 1);
        repmat(         BDM4_DLR0_ROW, 4, 1);
      ];
      V.isAchg                    = [0 0 0 0 0 0 0 0]';
      V.hasSwfFormat              = false;
      %
      V.expIsSaturatedAr          = [1 1 1 1 1 1 1 1]';

      testCase.test_get_voltage_saturation_quality_bit(testCase, V);
    end



    %#####
    % SWF
    %#####
    % One subsequence
    % SWF, BDM=0, LRX=0/AC diff
    % Saturated DC single
    function test_get_voltage_saturation_quality_bit___SWF_1(testCase)
      V = [];
      V.cwfSlidingWindowLengthSec = 2.1;   % Should be irrelevant since SWF.
      V.tsfFractionThreshold      = 0.5;
      V.thresholdAVoltDcSingle    = 1;
      V.thresholdAVoltDcDiff      = 3;
      V.thresholdAVoltAclg        = 5;
      V.thresholdAVoltAchg        = 7;
      %
      V.tt2000Ar                  = [10 11]' * 1e9;
      V.AsrSamplesAVoltSrm        = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", 2, 'EMPTY');
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V1").Asid, ...
        [
        [2 2 2 0 0 0 0] + 0; ...   % Saturated
        [0 0 0 2 2 2 2] + 0;
        ])
      V.AsrSamplesAVoltSrm.add(testCase.S("AC_V12").Asid, ...
        [
        [0 0 0 0 0 0 0] + 4; ...
        [0 0 0 2 2 2 2] + 6;
        ]+1)
      V.AsrSamplesAVoltSrm.add(testCase.S("AC_V23").Asid, ...
        [
        [0 0 0 0 0 0 0] + 4; ...
        [0 0 0 2 2 2 2] + 6;
        ]+1)
      V.zvNValidSamplesPerRecord  = [5 5]';
      V.bltsKSsidAr               = repmat(testCase.BDM0_DLR0_ROW, 2, 1);
      V.isAchg                    = [0 1]';
      V.hasSwfFormat              = true;
      %
      V.expIsSaturatedAr          = [1 0]';

      testCase.test_get_voltage_saturation_quality_bit(testCase, V);
    end



    % One subsequence
    % SWF, BDM=0, LRX=0/AC diff, LG+HG.
    % Varying snapshot lengths.
    % Saturated AC diffs
    function test_get_voltage_saturation_quality_bit___SWF_2(testCase)
      V = [];
      V.cwfSlidingWindowLengthSec = 2.1;   % Should be irrelevant since SWF.
      V.tsfFractionThreshold      = 0.5;
      V.thresholdAVoltDcSingle    = 1;
      V.thresholdAVoltDcDiff      = 3;
      V.thresholdAVoltAclg        = 5;
      V.thresholdAVoltAchg        = 7;
      %
      V.tt2000Ar                  = [10:13]' * 1e9;
      V.AsrSamplesAVoltSrm        = bicas.utils.SameRowsMap("bicas.proc.L1L2.AntennaSignalId", 4, 'EMPTY');
      V.AsrSamplesAVoltSrm.add(testCase.S("DC_V1").Asid, ...
        [
        [0 0 0 0 0 0 0] + 0; ...
        [0 2 2 2 2 2 2] + 0; ...
        [0 0 0 0 0 0 0] + 0; ...
        [0 0 0 2 2 2 2] + 0;
        ])
      V.AsrSamplesAVoltSrm.add(testCase.S("AC_V12").Asid, ...
        [
        [2 0 0 0 0 0 0] + 4; ...   % Saturated.
        [0 2 2 2 2 2 2] + 4; ...
        [2 2 2 0 0 0 0] + 6; ...   % Saturated.
        [0 0 0 2 2 2 2] + 6;
        ]+1)
      V.AsrSamplesAVoltSrm.add(testCase.S("AC_V23").Asid, ...
        [
        [0 0 0 0 0 0 0] + 4; ...
        [0 2 2 2 2 2 2] + 4; ...
        [0 0 0 0 0 0 0] + 6; ...
        [0 0 0 2 2 2 2] + 6;
        ]+1)
      V.zvNValidSamplesPerRecord  = [1 1 3 3]';
      V.bltsKSsidAr               = repmat(testCase.BDM0_DLR0_ROW, 4, 1);
      V.isAchg                    = [0 0 1 1]';
      V.hasSwfFormat              = true;
      %
      V.expIsSaturatedAr          = [1 0 1 0]';

      testCase.test_get_voltage_saturation_quality_bit(testCase, V);
    end



    function test_get_one_ASR_CWF_channel_TSF_bit_array___DC(testCase)
      testCase.test_get_one_ASR_CWF_channel_TSF_bit_array(testCase, ...
        thresholdAVoltDcSingle = 1, ...
        thresholdAVoltDcDiff   = 3, ...
        thresholdAVoltAclg     = 5, ...
        thresholdAVoltAchg     = 7, ...
        sSsid                  = "DC_V1", ...
        hasSwfFormat           = false, ...
        isAchg                 = [0,   0]', ...
        samplesAVolt           = [0.9, 1.1, ]', ...
        expTsfAr               = [0,   1]')
    end



    function test_get_one_ASR_CWF_channel_TSF_bit_array___AC(testCase)
      testCase.test_get_one_ASR_CWF_channel_TSF_bit_array(testCase, ...
        thresholdAVoltDcSingle = 1, ...
        thresholdAVoltDcDiff   = 3, ...
        thresholdAVoltAclg     = 5, ...
        thresholdAVoltAchg     = 7, ...
        sSsid                  = "AC_V12", ...
        hasSwfFormat           = false, ...
        isAchg                 = [0,   0,   1,   1  ]', ...
        samplesAVolt           = [4.9, 5.1, 6.9, 7.1]', ...
        expTsfAr               = [0,   1,   0,   1]')
    end



    function test_get_one_ASR_SWF_channel_saturation_bit_array___AC(testCase)
      testCase.test_get_one_ASR_SWF_channel_saturation_bit_array(testCase, ...
        tsfFractionThreshold     = 0.4, ...
        thresholdAVoltDcSingle   = 1, ...
        thresholdAVoltDcDiff     = 3, ...
        thresholdAVoltAclg       = 5, ...
        thresholdAVoltAchg       = 7, ...
        sSsid                    = "AC_V12", ...
        hasSwfFormat             = true, ...
        isAchg                   = [0; 0; 1; 1], ...
        samplesAVolt             = [
        4.9, 4.9, 9, 9;
        4.9, 5.1, 9, 9;
        6.9, 7.1, 0, 0;
        7.1, 7.1, 0, 0;
        ], ...
        zvNValidSamplesPerRecord = [2; 2; 2; 2], ...
        expTsfAr                 = [0; 1; 1; 1])
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function test_get_voltage_saturation_quality_bit(testCase, V)
      % Modify/normalize arguments.
      V.tt2000Ar     = int64(V.tt2000Ar);
      V.bltsKSsidAr  = bicas.sconst.C.SSID_S_K_DICT(V.bltsKSsidAr);
      V.isAchgFpa    = bicas.utils.FPArray(logical(V.isAchg), 'FILL_POSITIONS', isnan(V.isAchg));
      V.hasSwfFormat = logical(V.hasSwfFormat);

      V.expIsSaturatedAr = logical(V.expIsSaturatedAr);

      L = bicas.Logger('NO_STDOUT', false);
      Sat = testCase.init_object(...
        V.cwfSlidingWindowLengthSec, ...
        V.tsfFractionThreshold, ...
        V.thresholdAVoltDcSingle, ...
        V.thresholdAVoltDcDiff, ...
        V.thresholdAVoltAclg, ...
        V.thresholdAVoltAchg);

      % CALL FUNCTION
      actIsSaturatedAr = Sat.get_voltage_saturation_quality_bit(...
        V.tt2000Ar, V.AsrSamplesAVoltSrm, V.zvNValidSamplesPerRecord, ...
        V.bltsKSsidAr, V.isAchgFpa, V.hasSwfFormat, L);

      testCase.assertEqual(actIsSaturatedAr, V.expIsSaturatedAr)
    end



    function test_get_one_ASR_CWF_channel_TSF_bit_array(testCase, V)
      arguments
        testCase
        V.cwfSlidingWindowLengthSec
        V.tsfFractionThreshold
        V.thresholdAVoltDcSingle
        V.thresholdAVoltDcDiff
        V.thresholdAVoltAclg
        V.thresholdAVoltAchg
        V.sSsid
        V.hasSwfFormat
        V.isAchg
        V.samplesAVolt
        V.expTsfAr
      end

      % Should be irrelevant for this test.
      cwfSlidingWindowLengthSec = 99;
      tsfFractionThreshold      = 0.5;

      % Modify/normalize arguments.
      V.isAchgFpa = bicas.utils.FPArray(logical(V.isAchg), 'FILL_POSITIONS', isnan(V.isAchg));
      V.expTsfAr  = logical(V.expTsfAr);

      Sat = testCase.init_object(...
        cwfSlidingWindowLengthSec, ...
        tsfFractionThreshold, ...
        V.thresholdAVoltDcSingle, ...
        V.thresholdAVoltDcDiff, ...
        V.thresholdAVoltAclg, ...
        V.thresholdAVoltAchg);

      actTsfAr = Sat.get_one_ASR_CWF_channel_TSF_bit_array(...
        bicas.sconst.C.S_SSID_DICT(V.sSsid), V.isAchgFpa, V.samplesAVolt);

      testCase.assertEqual(actTsfAr, V.expTsfAr)
    end



    function test_get_one_ASR_SWF_channel_saturation_bit_array(testCase, V)
      arguments
        testCase
        V.cwfSlidingWindowLengthSec
        V.tsfFractionThreshold
        V.thresholdAVoltDcSingle
        V.thresholdAVoltDcDiff
        V.thresholdAVoltAclg
        V.thresholdAVoltAchg
        V.sSsid
        V.hasSwfFormat
        V.isAchg
        V.samplesAVolt
        V.zvNValidSamplesPerRecord
        V.expTsfAr
      end

      % Should be irrelevant for this test.
      cwfSlidingWindowLengthSec = 99;

      % Modify/normalize arguments.
      V.isAchgFpa = bicas.utils.FPArray(logical(V.isAchg), 'FILL_POSITIONS', isnan(V.isAchg));
      V.expTsfAr  = logical(V.expTsfAr);

      Sat = testCase.init_object(...
        cwfSlidingWindowLengthSec, ...
        V.tsfFractionThreshold, ...
        V.thresholdAVoltDcSingle, ...
        V.thresholdAVoltDcDiff, ...
        V.thresholdAVoltAclg, ...
        V.thresholdAVoltAchg);

      actTsfAr = Sat.get_one_ASR_SWF_channel_saturation_bit_array(...
        bicas.sconst.C.S_SSID_DICT(V.sSsid), V.isAchgFpa, ...
        V.samplesAVolt, V.zvNValidSamplesPerRecord);

      testCase.assertEqual(actTsfAr, V.expTsfAr)
    end



    function S = init_object(...
        cwfSlidingWindowLengthSec, ...
        tsfFractionThreshold, ...
        thresholdAVoltDcSingle, ...
        thresholdAVoltDcDiff, ...
        thresholdAVoltAclg, ...
        thresholdAVoltAchg)

      Bso = bicas.create_default_BSO();
      Bso.override_value('PROCESSING.SATURATION.CWF_SLIDING_WINDOW_LENGTH_SEC',            cwfSlidingWindowLengthSec, 'test');
      Bso.override_value('PROCESSING.SATURATION.TSF_FRACTION_THRESHOLD',                   tsfFractionThreshold,      'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.SINGLE',         thresholdAVoltDcSingle,    'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.DIFF',           thresholdAVoltDcDiff,      'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.LOW_GAIN',  thresholdAVoltAclg,        'test');
      Bso.override_value('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.HIGH_GAIN', thresholdAVoltAchg,        'test');
      Bso.make_read_only();

      S = bicas.proc.L1L2.Saturation(Bso);
    end



  end    % methods(Static, Access=private)



end
