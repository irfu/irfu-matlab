%
% UNFINISHED
%
% matlab.unittest automatic test code for bicas.proc.L1L2.Saturation.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Saturation___UTEST < matlab.unittest.TestCase



    properties(Constant)
    end



    properties(TestParameter)
    end



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_get_TSF(testCase)
            function test(satArgsCa, samplesAVolt, Ssid, isAchg, expTsfAr)
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
                    test({99, 0.6, 3, 99,99,99}, [-5, -1, 0, 1, 5, NaN], S.DC_V1,  isAchg, [1, 0, 0, 0, 1, 0])

                    % DC diff, 3x2
                    test({99, 0.6, 99, 3,99,99}, [-5, -1, 0; 1, 5, NaN], S.DC_V12, isAchg, [1, 0, 0; 0, 1, 0])
                end

                % AC, low gain, 6x1
                for achgThreshold = [1, 7]
                    test({99, 0.6, 99, 99, 3, achgThreshold}, [-5, -1, 0, 1, 5, NaN]', S.AC_V12, 0,    [1, 0, 0, 0, 1, 0]')
                end

                % AC, high gain, 1x6
                for aclgThreshold = [1, 7]
                    test({99, 0.6, 99, 99, aclgThreshold,  3}, [-5, -1, 0, 1, 5, NaN],  S.AC_V23, 1,    [1, 0, 0, 0, 1, 0] )
                end

                % AC, unknown gain, 1x6
                for aclgThreshold = [1, 7]
                    test({99, 0.6, 99, 99, aclgThreshold,  3}, [-5, -1, 0, 1, 5, NaN],  S.AC_V23, NaN,  [0, 0, 0, 0, 0, 0] )
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
                    test({99, 0.6, 99,99,99,99}, ones(0,1), ones(0,0),    S.DC_V1, isAchg, false(0,1))

                    % One snapshot of varying length.
                    test({99, 0.6, 3, 99,99,99}, [0],       [5, 0, 5],    S.DC_V1, isAchg, [0])
                    test({99, 0.6, 3, 99,99,99}, [1],       [5, 0, 5],    S.DC_V1, isAchg, [1])
                    test({99, 0.6, 3, 99,99,99}, [2],       [5, 0, 5],    S.DC_V1, isAchg, [0])
                    test({99, 0.6, 3, 99,99,99}, [3],       [5, 0, 5],    S.DC_V1, isAchg, [1])

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



%         function test_get_voltage_saturation_quality_bit(testCase)
%             % TODO
%         end



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
