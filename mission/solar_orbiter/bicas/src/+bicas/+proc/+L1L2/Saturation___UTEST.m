%
% matlab.unittest automatic test code for bicas.proc.L1L2.Saturation.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Saturation___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_TSF(testCase)

      % NOTE: Changes order of arguments to get easier-to-read hardcoded
      %       calls.
      function test(satArgsCa, Ssid, isAchg, samplesAVolt, expTsfAr)
        expTsfAr = logical(expTsfAr);
        isAchgFpa = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat = bicas.proc.L1L2.Saturation___UTEST.init_object(satArgsCa{:});

        actTsfAr = Sat.get_TSF(samplesAVolt, Ssid, isAchgFpa);

        testCase.assertEqual(actTsfAr, expTsfAr)
      end

      function main()
        S = bicas.proc.L1L2.SignalSourceId.C;

        for isAchg = [0, 1, NaN]
          % DC single, 1x6
          test(...
            {99, 0.6, 3, 99,99,99}, S.DC_V1,  isAchg, ...
            [-5, -1, 0, 1, 5, NaN], ...
            [ 1,  0, 0, 0, 1,   0])

          % DC diff, 3x2
          test(...
            {99, 0.6, 99, 3,99,99}, S.DC_V12, isAchg, ...
            [-5, -1, 0; 1, 5, NaN], ...
            [ 1,  0, 0; 0, 1,   0])
        end

        % AC, low gain, 6x1
        for achgThreshold = [1, 7]
          test(...
            {99, 0.6, 99, 99, 3, achgThreshold}, S.AC_V12, 0, ...
            [-5, -1, 0, 1, 5, NaN]', ...
            [ 1,  0, 0, 0, 1,   0]')
        end

        % AC, high gain, 1x6
        for aclgThreshold = [1, 7]
          test(...
            {99, 0.6, 99, 99, aclgThreshold,  3}, S.AC_V23, 1, ...
            [-5, -1, 0, 1, 5, NaN], ...
            [ 1,  0, 0, 0, 1,   0])
        end

        % AC, unknown gain, 1x6
        for aclgThreshold = [1, 7]
          test({99, 0.6, 99, 99, aclgThreshold,  3}, S.AC_V23, NaN, ...
            [-5, -1, 0, 1, 5, NaN], ...
            [ 0,  0, 0, 0, 0,   0])
        end
      end

      main()
    end



    function test_get_snapshot_saturation(testCase)

      function test(satArgsCa, samplesAVolt, Ssid, isAchg, expIsSaturated)
        isAchgFpa = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat = bicas.proc.L1L2.Saturation___UTEST.init_object(satArgsCa{:});

        actIsSaturated = Sat.get_snapshot_saturation(samplesAVolt, Ssid, isAchgFpa);

        testCase.assertEqual(actIsSaturated, expIsSaturated)
      end

      function main()
        S = bicas.proc.L1L2.SignalSourceId.C;

        % IMPLEMENTATION NOTE: Maybe somewhat overkill since
        % test_get_TSF() tests the different types of thresholds. The
        % extensive tests exist for historical reasons.


        for isAchg = [0, 1, NaN]
          % Empty snapshot.
          test({99, 0.6, 3, 99,99,99}, zeros(1,0), S.DC_V1, isAchg, false)
          % Length-one snapshot.
          test({99, 0.6, 3, 99,99,99}, [5],        S.DC_V1, isAchg, true)

          % DC single
          test({99, 0.6, 3, 99,99,99}, [ 0, 0], S.DC_V1, isAchg, false)
          test({99, 0.6, 3, 99,99,99}, [ 0, 5], S.DC_V1, isAchg, false)
          test({99, 0.6, 3, 99,99,99}, [-5, 5], S.DC_V1, isAchg, true)

          % DC diff
          test({99, 0.6, 99, 3, 99,99}, [0,  0], S.DC_V12, isAchg, false)
          test({99, 0.6, 99, 3, 99,99}, [5, -5], S.DC_V13, isAchg, true)
        end

        for unusedGainThreshold = [1, 7]
          % AC, low gain
          test({99, 0.6, 99, 99, 3, unusedGainThreshold}, [0,  0], S.AC_V12, 0, false)
          test({99, 0.6, 99, 99, 3, unusedGainThreshold}, [5, -5], S.AC_V13, 0, true)

          % AC, high gain
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [0,  0], S.AC_V23, 1, false)
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [5, -5], S.AC_V23, 1, true)
        end

        % NaN samples
        test({99, 0.6, 3, 99,99,99}, [ NaN, NaN], S.DC_V1, 0, false)
        test({99, 0.6, 3, 99,99,99}, [ NaN,   5], S.DC_V1, 0, false)
        test({99, 0.6, 3, 99,99,99}, [-5,     5], S.DC_V1, 0, true)
      end

      main()
    end



    function test_get_snapshot_saturation_many(testCase)

      function test(...
          satArgsCa, ...
          zvNValidSamplesPerRecord, samplesAVolt, ...
          Ssid, isAchg, expIsSaturatedAr)

        expIsSaturatedAr = logical(expIsSaturatedAr);
        isAchgFpa  = bicas.utils.FPArray.floatNan2logical(isAchg);
        Sat = bicas.proc.L1L2.Saturation___UTEST.init_object(satArgsCa{:});

        actIsSaturatedAr = Sat.get_snapshot_saturation_many(...
          zvNValidSamplesPerRecord, samplesAVolt, Ssid, isAchgFpa);

        testCase.assertEqual(actIsSaturatedAr, expIsSaturatedAr)
      end

      function main()
        S = bicas.proc.L1L2.SignalSourceId.C;

        % IMPLEMENTATION NOTE: Maybe somewhat overkill since
        % test_get_TSF() tests the different types of thresholds. The
        % extensive tests exist for historical reasons.

        for isAchg = [0, 1, NaN]
          % DC single
          % ---------
          % Zero snapshots
          test({99, 0.6, 99,99,99,99}, zeros(0,1), zeros(0,0),   S.DC_V1, isAchg, false(0,1))

          % One snapshot of varying length.
          test({99, 0.6, 3, 99,99,99}, [0],        [5, 0, 5],    S.DC_V1, isAchg, [0])
          test({99, 0.6, 3, 99,99,99}, [1],        [5, 0, 5],    S.DC_V1, isAchg, [1])
          test({99, 0.6, 3, 99,99,99}, [2],        [5, 0, 5],    S.DC_V1, isAchg, [0])
          test({99, 0.6, 3, 99,99,99}, [3],        [5, 0, 5],    S.DC_V1, isAchg, [1])

          % Many snapshots
          test({99, 0.6, 3, 99,99,99}, [0;1;2;3], [5 0 5; 5 0 5; 5 0 5; 5 0 5],    S.DC_V1, isAchg, [0; 1; 0; 1])

          % DC diff
          test({99, 0.6, 99, 3, 99,99}, 2, [0,  0], S.DC_V12, isAchg, [0])
          test({99, 0.6, 99, 3, 99,99}, 2, [5, -5], S.DC_V13, isAchg, [1])
        end

        for unusedGainThreshold = [1, 7]
          % AC, low gain
          test({99, 0.6, 99, 99, 3, unusedGainThreshold}, [1; 1], [0  0; 5 -5], S.AC_V12, 0, [0; 1])

          % AC, high gain
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [1],    [0,  0],      S.AC_V23, 1, [0])
          test({99, 0.6, 99, 99, unusedGainThreshold, 3}, [1],    [5, -5],      S.AC_V23, 1, [1])
        end

        % NaN samples
        test({99, 0.6, 3, 99,99,99}, [2], [ NaN,   5], S.DC_V1, 0, false)
        test({99, 0.6, 3, 99,99,99}, [2], [-5,     5], S.DC_V1, 0, true)
      end

      main()
    end



    function test_get_voltage_saturation_quality_bit(testCase)

      function test(V, expIsSaturatedAr)
        % Modify/normalize arguments.
        V.tt2000Ar     = int64(V.tt2000Ar);
        V.dlrFpa       = bicas.utils.FPArray(logical(V.dlr),    'FILL_POSITIONS', isnan(V.dlr));
        V.bdmFpa       = bicas.utils.FPArray(uint8(V.bdm),      'FILL_POSITIONS', isnan(V.bdm));
        V.isAchgFpa    = bicas.utils.FPArray(logical(V.isAchg), 'FILL_POSITIONS', isnan(V.isAchg));
        V.hasSwfFormat = logical(V.hasSwfFormat);

        expIsSaturatedAr = logical(expIsSaturatedAr);

        L = bicas.Logger('NO_STDOUT', false);
        S = bicas.proc.L1L2.Saturation___UTEST.init_object(...
          V.cwfSlidingWindowLengthSec, ...
          V.tsfFractionThreshold, ...
          V.thresholdAVoltDcSingle, ...
          V.thresholdAVoltDcDiff, ...
          V.thresholdAVoltAclg, ...
          V.thresholdAVoltAchg);

        % CALL FUNCTION
        actIsSaturatedAr = S.get_voltage_saturation_quality_bit(...
          V.tt2000Ar, V.AsrSamplesAVoltSrm, V.zvNValidSamplesPerRecord, ...
          V.bdmFpa, V.dlrFpa, V.lrx, V.isAchgFpa, V.hasSwfFormat, L);

        testCase.assertEqual(actIsSaturatedAr, expIsSaturatedAr)
      end

      C = bicas.proc.L1L2.SignalSourceId.C;
      ALL_ENABLED = true;
      %ALL_ENABLED = false;

      if ALL_ENABLED
        % Zero records.
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
          V.AsrSamplesAVoltSrm.add(C.DC_V1.Asid, zeros(0, 1))
          V.AsrSamplesAVoltSrm.add(C.DC_V2.Asid, zeros(0, 1))
          V.AsrSamplesAVoltSrm.add(C.DC_V3.Asid, zeros(0, 1))
          V.zvNValidSamplesPerRecord  = zeros(0, 1);
          V.bdm                       = zeros(0, 1);
          V.dlr                       = zeros(0, 1);
          V.lrx                       = zeros(0, 1);
          V.isAchg                    = zeros(0, 1);
          V.hasSwfFormat              = logical(hasSwfFormat);

          test(V, false(0,1));
        end
      end

      %#####
      % CWF
      %#####
      % One subsequence
      % CWF, BDM=0, LRX=1/DC diff.
      if ALL_ENABLED
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
        V.AsrSamplesAVoltSrm.add(C.DC_V1.Asid,  [2 0 2]')
        V.AsrSamplesAVoltSrm.add(C.DC_V12.Asid, [2 0 2]'+1)
        V.AsrSamplesAVoltSrm.add(C.DC_V23.Asid, [2 0 2]'+1)
        V.zvNValidSamplesPerRecord  = [1 1 1]';
        V.bdm                       = [0 0 0]';
        V.dlr                       = [0 0 0]';
        V.lrx                       = [1 1 1]';
        V.isAchg                    = [0 0 0]';
        V.hasSwfFormat              = false;

        test(V, [1 1 1]');
      end

      % One subsequence
      % Separate episodes of saturation on different channels.
      % CWF, BDM=0, LRX=1/AC diff.
      if ALL_ENABLED
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
        V.AsrSamplesAVoltSrm.add(C.DC_V1.Asid,  [2 2 0 0 0 0 0 0 0 0]')
        V.AsrSamplesAVoltSrm.add(C.AC_V12.Asid, [0 0 0 2 2 0 0 0 0 0]'+2)
        V.AsrSamplesAVoltSrm.add(C.AC_V23.Asid, [0 0 0 0 0 0 2 2 0 0]'+2)
        V.zvNValidSamplesPerRecord  = [1 1 1 1 1 1 1 1 1 1]';
        V.bdm                       = [0 0 0 0 0 0 0 0 0 0]';
        V.dlr                       = [0 0 0 0 0 0 0 0 0 0]';
        V.lrx                       = [0 0 0 0 0 0 0 0 0 0]';
        V.isAchg                    = [0 0 0 0 0 0 0 0 0 0]';
        V.hasSwfFormat              = false;

        test(V, [1 1 0 1 1 0 1 1 0 0]');
      end

      % Multiple subsequences
      % Separate episodes of saturation on different channels.
      % CWF, BDM=0, LRX=1/AC diff.
      if ALL_ENABLED
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
        V.AsrSamplesAVoltSrm.add(C.DC_V1.Asid,  [0 0 2 2 0 0 0 0]')
        V.AsrSamplesAVoltSrm.add(C.DC_V12.Asid, [0 0 0 0 0 0 2 2]'+2)
        V.AsrSamplesAVoltSrm.add(C.DC_V23.Asid, [0 0 0 0 0 0 0 0]'+2)
        V.AsrSamplesAVoltSrm.add(C.DC_V2.Asid,  [2 2 0 0 0 0 0 0]')
        V.AsrSamplesAVoltSrm.add(C.DC_V3.Asid,  [0 0 0 0 2 2 0 0]')
        V.zvNValidSamplesPerRecord  = [1 1 1 1 1 1 1 1]';
        V.bdm                       = [0 0 0 0 4 4 4 4]';
        V.dlr                       = [0 0 0 0 0 0 0 0]';
        V.lrx                       = [1 1 1 1 1 1 1 1]';
        V.isAchg                    = [0 0 0 0 0 0 0 0]';
        V.hasSwfFormat              = false;

        test(V, [0 0 1 1 1 1 0 0]');
      end

      %#####
      % SWF
      %#####
      % One subsequence
      % SWF, BDM=0, LRX=0/AC diff
      % Saturated DC single
      if ALL_ENABLED
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
        V.AsrSamplesAVoltSrm.add(C.DC_V1.Asid, ...
          [
          [2 2 2 0 0 0 0] + 0; ...   % Saturated
          [0 0 0 2 2 2 2] + 0;
          ])
        V.AsrSamplesAVoltSrm.add(C.AC_V12.Asid, ...
          [
          [0 0 0 0 0 0 0] + 4; ...
          [0 0 0 2 2 2 2] + 6;
          ]+1)
        V.AsrSamplesAVoltSrm.add(C.AC_V23.Asid, ...
          [
          [0 0 0 0 0 0 0] + 4; ...
          [0 0 0 2 2 2 2] + 6;
          ]+1)
        V.zvNValidSamplesPerRecord  = [5 5]';
        V.bdm                       = [0 0]';
        V.dlr                       = [0 0]';
        V.lrx                       = [0 0]';
        V.isAchg                    = [0 1]';
        V.hasSwfFormat              = true;

        test(V, [1 0]');
      end

      % One subsequence
      % SWF, BDM=0, LRX=0/AC diff, LG+HG.
      % Varying snapshot lengths.
      % Saturated AC diffs
      if ALL_ENABLED
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
        V.AsrSamplesAVoltSrm.add(C.DC_V1.Asid, ...
          [
          [0 0 0 0 0 0 0] + 0; ...
          [0 2 2 2 2 2 2] + 0; ...
          [0 0 0 0 0 0 0] + 0; ...
          [0 0 0 2 2 2 2] + 0;
          ])
        V.AsrSamplesAVoltSrm.add(C.AC_V12.Asid, ...
          [
          [2 0 0 0 0 0 0] + 4; ...   % Saturated.
          [0 2 2 2 2 2 2] + 4; ...
          [2 2 2 0 0 0 0] + 6; ...   % Saturated.
          [0 0 0 2 2 2 2] + 6;
          ]+1)
        V.AsrSamplesAVoltSrm.add(C.AC_V23.Asid, ...
          [
          [0 0 0 0 0 0 0] + 4; ...
          [0 2 2 2 2 2 2] + 4; ...
          [0 0 0 0 0 0 0] + 6; ...
          [0 0 0 2 2 2 2] + 6;
          ]+1)
        V.zvNValidSamplesPerRecord  = [1 1 3 3]';
        V.bdm                       = [0 0 0 0]';
        V.dlr                       = [0 0 0 0]';
        V.lrx                       = [0 0 0 0]';
        V.isAchg                    = [0 0 1 1]';
        V.hasSwfFormat              = true;

        test(V, [1 0 1 0]');
      end

    end    % test_get_voltage_saturation_quality_bit



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



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
