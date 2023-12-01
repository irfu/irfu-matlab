%
% matlab.unittest automatic test code for bicas.proc.L1L2.qual.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_modify_quality_filter(testCase)
            % IMPLEMENTATION NOTE: Test exist partly for historical reasons.
            % The "real tests" test the two functions called by
            % bicas.proc.L1L2.qual.modify_quality_filter().
            %
            % PROPOSAL: QUALITY_FLAG FPs.

            % One output variable.
            function test(ZvIn, isLfr, NsoTable, S, expZvOut)
                % NOTE: Does not need to handle PROCESSING.ZV_QUALITY_FLAG_MAX.
                % That is handled by bicas.write_dataset_CDF().

                Bso = bicas.create_default_BSO();
                Bso.override_value('PROCESSING.L2.REMOVE_DATA.MUX_MODES',             S.bdmRemoveArray,  'test')
                Bso.override_value('PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S', S.lfrBdmMarginSec, 'test')
                Bso.override_value('PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S', S.tdsBdmMarginSec, 'test')
                Bso.make_read_only()
                L = bicas.Logger('human-readable', false);

                % CALL FUNCTION
                [actZvUfv, actZv_QUALITY_FLAG, actZv_L2_QUALITY_BITMASK] = ...
                    bicas.proc.L1L2.qual.modify_quality_filter(ZvIn, isLfr, NsoTable, Bso, L);

                testCase.assertEqual(actZvUfv,                 expZvOut.ufv)
                testCase.assertEqual(actZv_QUALITY_FLAG,       expZvOut.QUALITY_FLAG_Fpa)
                testCase.assertEqual(actZv_L2_QUALITY_BITMASK, expZvOut.L2_QUALITY_BITMASK)
            end

            %===================================================================

            ENA           = zeros(0, 1);
            ECA           = cell(0, 1);
            ENA_UINT8_FPA = bicas.utils.FPArray(uint8(ENA));
            EmptyNsoTable    = bicas.NsoTable(int64(ENA), int64(ENA), ECA);
            NonemptyNsoTable = bicas.NsoTable(...
                int64([1, 2]'*1e9), ...
                int64([2, 3]'*1e9), ...
                {bicas.const.NSOID.PARTIAL_SATURATION, ...
                 bicas.const.NSOID.FULL_SATURATION}');

            for isLfr = [false, true]
                for NsoTable = [EmptyNsoTable, NonemptyNsoTable]
                    if 1
                        %===============
                        % "Simple test"
                        %===============
                        % Empty data
                        % LFR/TDS, empty/non-empty NSO table
                        Settings = struct(...
                            'bdmRemoveArray',  [1, 2], ...
                            'lfrBdmMarginSec', 1.5, ...
                            'tdsBdmMarginSec', 2.5);
                        ZvIn = struct(...
                            'Epoch',               int64(ENA), ...
                            'bdmFpa',              ENA_UINT8_FPA, ...
                            'QUALITY_FLAG_Fpa',    bicas.utils.FPArray(uint8(ENA)), ...
                            'ufv',                 logical(ENA));
                        expZvOut = struct(...
                            'QUALITY_FLAG_Fpa',    bicas.utils.FPArray(uint8(ENA)), ...
                            'L2_QUALITY_BITMASK',  uint16(ENA), ...
                            'ufv',                 false(0, 1));
                        test(ZvIn, isLfr, NsoTable, Settings, expZvOut);
                    end    % if
                end   % for
            end    % for

        end
        
        
        
        function test_get_quality_by_NSOs(testCase)
            
            function test(NsoidSettingsMap, NsoTable, Epoch, exp_QUALITY_FLAG_doubleNaN, exp_L2_QUALITY_BITMASK)
                Epoch = int64(Epoch(:));
                exp_QUALITY_FLAG_Fpa   = bicas.utils.FPArray(exp_QUALITY_FLAG_doubleNaN(:), 'FILL_VALUE', NaN).cast('uint8');
                exp_L2_QUALITY_BITMASK = uint16(exp_L2_QUALITY_BITMASK(:));

                L = bicas.Logger('human-readable', false);

                [act_QUALITY_FLAG_Fpa, act_L2_QUALITY_BITMASK] = bicas.proc.L1L2.qual.get_quality_by_NSOs(...
                    NsoidSettingsMap, NsoTable, Epoch, L);
            
                testCase.assertEqual(act_QUALITY_FLAG_Fpa,   exp_QUALITY_FLAG_Fpa)
                testCase.assertEqual(act_L2_QUALITY_BITMASK, exp_L2_QUALITY_BITMASK)
            end

            ENA = zeros(0, 1);
            ECA = cell(0, 1);

            EMPTY_NSO_TABLE = bicas.NsoTable(int64(ENA), int64(ENA), ECA);
            EMPTY_NSO_SETTINGS_MAP = containers.Map;

            %====================================
            % Empty NSO table. Various Epoch ZVs
            %====================================
            EPOCH_DOUBLE_CA = {zeros(0,1), [10], [10;20;30]};
            EPOCH_DOUBLE_CA = EPOCH_DOUBLE_CA(2);
            for i = 1:numel(EPOCH_DOUBLE_CA)
                Epoch_double = EPOCH_DOUBLE_CA{i};
                
                test(...
                    EMPTY_NSO_SETTINGS_MAP, ...
                    EMPTY_NSO_TABLE, ...
                    Epoch_double, ...
                    3*ones(size(Epoch_double)), ...
                    zeros(size(Epoch_double)) ...
                )
            end

            %=====================================
            % Non-overlapping NSOs, one at a time
            %=====================================
            % Nontrivial NSO settings.
            NSOID_SETTINGS = containers.Map;
            NSOID_SETTINGS('N' )    = bicas.proc.L1L2.NsoidSetting(uint8(3), uint16(0));     % "Neutral" = Do nothing
            NSOID_SETTINGS('QF')    = bicas.proc.L1L2.NsoidSetting(uint8(1), uint16(0));     % Cap QUALITY_FLAG
            NSOID_SETTINGS('BM')    = bicas.proc.L1L2.NsoidSetting(uint8(3), uint16(1+2));   % Set L2_QUALITY_BITMASK bits
            NSOID_SETTINGS('QF_BM') = bicas.proc.L1L2.NsoidSetting(uint8(2), uint16(2+4));   % Cap QUALITY_FLAG & set L2_QUALITY_BITMASK bits
            NSO_TABLE = bicas.NsoTable(...
                int64([1, 4, 7, 10]'*1e9), ...
                int64([2, 5, 8, 11]'*1e9), ...
                {'N', 'QF', 'BM', 'QF_BM'}');
            test(...
                NSOID_SETTINGS, NSO_TABLE, ...
                [0:3]*1e9, ...
                [3, 3, 3, 3], ...
                [0, 0, 0, 0] ...
            );
            test(...
                NSOID_SETTINGS, NSO_TABLE, ...
                [3:6]*1e9, ...
                [3, 1, 1, 3], ...
                [0, 0, 0, 0] ...
            );
            test(...
                NSOID_SETTINGS, ...
                NSO_TABLE, ...
                [6:9]'*1e9, ...
                [3, 3, 3, 3], ...
                [0, 3, 3, 0] ...
            );
            test(...
                NSOID_SETTINGS, ...
                NSO_TABLE, ...
                [9:12]'*1e9, ...
                [3, 2, 2, 3], ...
                [0, 6, 6, 0] ...
            );
            % NSO overlaps Epoch boundaries.
            test(...
                NSOID_SETTINGS, ...
                NSO_TABLE, ...
                [10.1, 10.9]'*1e9, ...
                [2, 2], ...
                [6, 6] ...
            );
        
            %==================
            % Overlapping NSOs
            %==================
            NSOID_SETTINGS = containers.Map;
            NSOID_SETTINGS('A') = bicas.proc.L1L2.NsoidSetting(uint8(2), uint16(1+2));
            NSOID_SETTINGS('B') = bicas.proc.L1L2.NsoidSetting(uint8(1), uint16(2+4));
            NSO_TABLE = bicas.NsoTable(...
                int64([1, 2]'*1e9), ...
                int64([2, 3]'*1e9), ...
                {'A', 'B'}');
            test(...
                NSOID_SETTINGS, NSO_TABLE, ...
                [0:4]*1e9, ...
                [3, 2, 1, 1, 3], ...
                [0, 3, 7, 6, 0] ...
            );
            % Epoch does not overlap with any NSOs.
            test(...
                NSOID_SETTINGS, NSO_TABLE, ...
                [-1, 4]*1e9, ...
                [3, 3], ...
                [0, 0] ...
            );
        end



        function test_set_voltage_current_FV(testCase)

            function test(...
                    zv_Epoch, zvUfv, zvAsrSamplesAVoltSrm, zvCurrentAAmpere, ...
                    expZvAsrSamplesAVoltSrm, expZvCurrentAAmpere)

                nRows = irf.assert.sizes(...
                    zv_Epoch,         [-1], ...
                    zvCurrentAAmpere, [-1, 3]);
                assert(zvAsrSamplesAVoltSrm.nRows == nRows)
                L = bicas.Logger('none', false);

                % NOTE: Modifies argument zvAsrSamplesAVoltSrm (handle object).
                actZvCurrentAAmpere = bicas.proc.L1L2.qual.set_voltage_current_FV(...
                    zv_Epoch, zvAsrSamplesAVoltSrm, zvCurrentAAmpere, zvUfv, L);

                actZvAsrSamplesAVoltSrm = zvAsrSamplesAVoltSrm;
                testCase.verifyEqual(actZvAsrSamplesAVoltSrm,  expZvAsrSamplesAVoltSrm)
                testCase.verifyEqual(actZvCurrentAAmpere, expZvCurrentAAmpere)
            end

            %===================================================================

            % Empty data.
            test( ...
                int64(zeros(0, 1)), ...
                false(0, 1), ...
                bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(0, 1)), ...
                zeros(0, 3),  ...
                bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(0, 1)), ...
                zeros(0, 3)  ...
            )

            % Non-empty input data that is not altered.
            test( ...
                int64([10, 11, 12, 13, 14]'), ...
                false(5, 1), ...
                bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(5, 1)), ...
                zeros(5, 3),  ...
                bicas.proc.L1L2.qual___UTEST.AsSrm(zeros(5, 1)), ...
                zeros(5, 3)  ...
            )

            % Non-empty input data that is altered: 1 unaltered + 1 altered
            test( ...
                int64([10, 11]'), ...
                logical([0, 1]'), ...
                bicas.proc.L1L2.qual___UTEST.AsSrm([1, 2]'), ...
                [1:3; 11:13],  ...
                bicas.proc.L1L2.qual___UTEST.AsSrm([1, NaN]'), ...
                [1:3; NaN(1, 3)] ...
            )
            % Non-empty input data that is altered: 2 unaltered + (1+2) altered
            test( ...
                int64([10, 11, 12, 13, 14]'), ...
                logical([0, 1, 0, 1, 1]'), ...
                bicas.proc.L1L2.qual___UTEST.AsSrm([1, 2, 3, 4, 5]'), ...
                [1:3; 11:13; 21:23; 31:33; 41:43],  ...
                bicas.proc.L1L2.qual___UTEST.AsSrm([1, NaN, 3, NaN, NaN]'), ...
                [1:3; NaN(1,3); 21:23; NaN(2,3)] ...
            )

        end



        function test_get_UFV_records_from_settings(testCase)

            function actZvUfv = call_func(zv_Epoch, zvBdmDoubleNan, isLfr, bdmRemoveArray, lfrBdmMarginSec, tdsBdmMarginSec)
                assert(isfloat(zvBdmDoubleNan))
                zvBdmFpa = bicas.utils.FPArray(zvBdmDoubleNan, 'FILL_VALUE', NaN).cast('int8');

                Bso = bicas.create_default_BSO();
                Bso.override_value('PROCESSING.L2.REMOVE_DATA.MUX_MODES',             bdmRemoveArray,  'test')
                Bso.override_value('PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S', lfrBdmMarginSec, 'test')
                Bso.override_value('PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S', tdsBdmMarginSec, 'test')
                Bso.make_read_only()

                L = bicas.Logger('none', false);
                
                actZvUfv = bicas.proc.L1L2.qual.get_UFV_records_from_settings(...
                    zv_Epoch, zvBdmFpa, isLfr, Bso, L);
            end



            function test(zv_Epoch, zvBdm, isLfr, bdmRemoveArray, lfrBdmMarginSec, tdsBdmMarginSec, expZvUfv)
                irf.assert.sizes(...
                    zv_Epoch, [-1, 1], ...
                    zvBdm,    [-1, 1], ...
                    expZvUfv, [-1, 1] ...
                )

                actZvUfv = call_func(zv_Epoch, zvBdm, isLfr, bdmRemoveArray, lfrBdmMarginSec, tdsBdmMarginSec);

                testCase.verifyEqual(actZvUfv, expZvUfv)
            end
            %===================================================================



            if 1
                for isLfr = [false, true]
                    % NOTE: Test PROCESSING.L2.REMOVE_DATA.MUX_MODES = [] (0x0)
                    % as is likely to be set to.
                    % Iterate over bdmRemoveArray values without BDM=0.
                    for bdmRemoveArrayCa = {[], zeros(0,1), [2]', [1,2,3]'}
                        bdmRemoveArray = bdmRemoveArrayCa{1};

                        % Zero records
                        test(...
                            int64(zeros(0, 1)), ...
                            zeros(0, 1), ...
                            isLfr, bdmRemoveArray, 10, 20, false(0, 1))

                        % No UFV
                        test(...
                            int64([10, 11, 12]'), ...
                            [0, 0, 0]', ...
                            isLfr, bdmRemoveArray, 10, 20, [false, false, false]' ...
                        )
                        % All UFV
                        test(...
                            int64([10, 11, 12]'), ...
                            [0, 0, 0]', ...
                            isLfr, bdmRemoveArray, 10, 20, [false, false, false]' ...
                        )
                    end
                end
            end

            % LFR
            % Test BDM=unknown/FV
            test(...
                int64([10:20]') * 1e9, ...
                [2, 2, 0, 0, 0, NaN, 0, 0, 0, 3, 3]', ...
                true, [2, 3], 1.5, 2.5, ...
                logical([1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1])' ...
            )

            % TDS
            % Test BDM=unknown/FV
            test(...
                int64([10:20]') * 1e9, ...
                [2, 2, 0, 0, 0, NaN, 0, 0, 0, 3, 3]', ...
                false, [2, 3], 1.5, 2.5, ...
                logical([1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1])' ...
            )
        end



        function test_sliding_window_over_fraction(testCase)

            % Whether to also test calling the tested function with the same
            % arguments modified to represent the opposite direction.
            % Useful for disabling while debugging failed tests.
            TEST_REVERSE_ORDER = true;                
            
            % "Raw test", without any additional functionality. Not intended to
            % be called with hardcoded arguments.
            function test_raw(...
                    tt2000Ar, bFlag1Ar, expBFlagAr, flagFractionThreshold, intervalLengthSec)

                actBFlagAr = bicas.proc.L1L2.qual.sliding_window_over_fraction(...
                    tt2000Ar, bFlag1Ar, flagFractionThreshold, intervalLengthSec);

                testCase.verifyEqual(actBFlagAr, expBFlagAr)
            end

            % Test function intended to be called with hardcoded arguments.
            function test(...
                    timeSecAr, bFlag1Ar, expBFlagAr, ...
                    flagFractionThreshold, intervalLengthSec)

                assert(all(ismember(bFlag1Ar,   [0,1])))
                assert(all(ismember(expBFlagAr, [0,1])))

                % Convert arguments to simplify calls.
                % NOTE: Normalize 0x0, transpose, change types
                timeSecAr  = int64(  timeSecAr (:)) * 1e9;
                bFlag1Ar   = logical(bFlag1Ar  (:));
                expBFlagAr = logical(expBFlagAr(:));

                % Run test with arguments as they are.
                test_raw(...
                    timeSecAr, ...
                    bFlag1Ar, expBFlagAr, ...
                    flagFractionThreshold, intervalLengthSec)

                % Run test with arguments modified to correspond to a reversed
                % order of samples.
                if TEST_REVERSE_ORDER
                    test_raw(...
                        flipud(-timeSecAr), ...   % NOTE: Negation
                        flipud(bFlag1Ar), flipud(expBFlagAr), ...
                        flagFractionThreshold, intervalLengthSec)
                end
            end

            % TODO: Identical timestamps.       -- DONE
            % TODO: Data gaps.                  -- DONE
            % TODO: Varying sampling rate.      -- DONE
            % TODO: Total time < window length  -- DONE
            % TODO: Test equality in interval length, fraction.
            
            ENABLE_ALL = 1;
            if ENABLE_ALL
                % ============
                % Zero samples
                % ============
                test( [], [], [], 0, 1 )
                test( [], [], [], 1, 1 )

                % ==========
                % One sample
                % ==========
                test(...
                    [10], ...
                    [ 0], ...
                    [ 0], ...
                    1, 2)
                test(...
                    [10], ...
                    [ 1], ...
                    [ 1], ...
                    1, 2)
            
                % ======================
                % Multiple samples
                % Constant sampling rate
                % ======================
                
                % Change length of window
                % -----------------------
                % Window that always covers 2 samples.
                test(...
                    [10, 11, 12], ...
                    [ 0,  1,  0], ...
                    [ 0,  0,  0], ...
                    1, 2)
                test(...
                    [10, 11, 12], ...
                    [ 1,  0,  0], ...
                    [ 0,  0,  0], ...
                    1, 2)
                % Window that always covers 1 sample.
                test(...
                    [10, 11, 12], ...
                    [ 0,  1,  0], ...
                    [ 0,  1,  0], ...
                    1, 1)
                test(...
                    [10, 11, 12], ...
                    [ 1,  0,  0], ...
                    [ 1,  0,  0], ...
                    1, 1)
                
                % Change fraction required
                % ------------------------
                % Window that always covers 2 samples. Threshold not reached.
                test(...
                    [10, 11, 12], ...
                    [ 0,  1,  0], ...
                    [ 0,  0,  0], ...
                    0.55, 2.01)
                test(...
                    [10, 11, 12], ...
                    [ 1,  0,  0], ...
                    [ 0,  0,  0], ...
                    0.55, 2)
                % Window that always covers 2 samples. Threshold reached.
                test(...
                    [10, 11, 12], ...
                    [ 0,  1,  0], ...
                    [ 1,  1,  1], ...
                    0.45, 2)
                test(...
                    [10, 11, 12], ...
                    [ 1,  0,  0], ...
                    [ 1,  1,  0], ...
                    0.45, 2)

                % Different thresholds +-0.5
                % --------------------------
                test(...
                    [10:22], ...
                    [0 0 1 0 1 1 1 0 0 0 1 0 0], ...
                    [0 1 1 1 1 1 1 1 0 1 1 1 0], ...
                    0.49, 2)
                test(...
                    [10:22], ...
                    [0 0 1 0 1 1 1 0 0 0 1 0 0], ...
                    [0 0 0 0 1 1 1 0 0 0 0 0 0], ...
                    0.51, 2)
            end
            
            if ENABLE_ALL
                % Data gap
                % --------
                % 2/3 samples in window satisfies threshold.
                test(...
                    [10 11 12 13   15 16 17 18], ...
                    [ 1  1  1  1    1  0  0  0], ...
                    [ 1  1  1  1    1  0  0  0], ...
                    0.66, 3)
                test(...
                    [10 11 12 13   16 17 18], ...
                    [ 1  1  1  1    0  0  0], ...
                    [ 1  1  1  1    0  0  0], ...
                    0.66, 3)
                % 1/3 samples in window satisfies threshold.
                test(...
                    [10 11 12 13   15 16 17 18], ...
                    [ 1  1  1  1    0  0  0  0], ...
                    [ 1  1  1  1    1  0  0  0], ...
                    0.33, 3)
            end
            
            if ENABLE_ALL
                % =====================
                % Varying sampling rate
                % =====================
                test(...
                    [10 11 12 13  15 17 19 21], ...
                    [ 1  0  1  0   1  0  1  0], ...
                    [ 1  1  1  1   1  0  1  0], ...
                    0.5, 2)
                test(...
                    [10 11 12 13  15 17 19 21], ...
                    [ 1  0  1  0   1  0  0  0], ...
                    [ 1  1  1  1   1  1  0  0], ...
                    0.5, 4)
                test(...
                    [10 11 12 13  15 17 19 21], ...
                    [ 0  0  0  0   1  0  0  0], ...
                    [ 0  0  0  1   1  1  0  0], ...
                    0.5, 4)

                % ====================
                % Identical timestamps
                % ====================
                test(...
                    [10 11 12 13 13 14 15 16], ...
                    [ 0  0  1  0  0  1  0  0], ...
                    [ 0  0  1  1  1  1  0  0], ...
                    0.66, 3)
                test(...
                    [10 11 12 13 13 14 15 16], ...
                    [ 0  0  0  1  1  0  0  0], ...
                    [ 0  0  0  0  0  0  0  0], ...
                    0.01, 3)
            end
            
            % ================================
            % Window same length as total time
            % ================================
            test(...
                [10 11 12], ...
                [ 0  1  0], ...
                [ 1  1  1], ...
                0.33, 3)
            test(...
                [10 11 12], ...
                [ 0  1  0], ...
                [ 0  0  0], ...
                0.34, 3)

            % ====================================
            % Window length longer than total time
            % ====================================
            test(...
                [10 11 12], ...
                [ 0  1  0], ...
                [ 0  0  0], ...
                0.30, 4)
            test(...
                [10 11 12], ...
                [ 1  1  0], ...
                [ 1  1  1], ...
                0.30, 4)
        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Utility function for creating a simplified AsrSamplesAVoltSrm data
        % structure with identical values on every channel.
        function AsrSamplesAVoltSrm = AsSrm(v)
            assert(iscolumn(v))
            AsrSamplesAVoltSrm = bicas.utils.SameRowsMap(...
                'char', size(v, 1), 'CONSTANT', v, ...
                bicas.proc.L1L2.AntennaSignalId.C.ALL_ASID_NAMES_CA(:));
        end



    end    % methods(Static, Access=private)



end
