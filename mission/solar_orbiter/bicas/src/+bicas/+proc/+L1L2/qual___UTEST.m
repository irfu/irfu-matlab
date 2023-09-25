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
            % PROPOSAL: Test for specific sets of NSOs (one or multiple) in one timestamp/sample.
            %   PRO: Can simplify a lot using helper function.
            %   PRO: More useful if having many NSOIDs.

            % One output variable.
            function test(ZvIn, isLfr, NsoTable, S, expZvOut)
                % NOTE: Does not need to handle PROCESSING.ZV_QUALITY_FLAG_MAX.
                % That is handled by bicas.write_dataset_CDF().

                SETTINGS = bicas.create_default_SETTINGS();
                SETTINGS.override_value('PROCESSING.L2.REMOVE_DATA.MUX_MODES',             S.rmMuxModesArray,     'test')
                SETTINGS.override_value('PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S', S.lfrMuxModeMarginSec, 'test')
                SETTINGS.override_value('PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S', S.tdsMuxModeMarginSec, 'test')
                SETTINGS.make_read_only()

                L = bicas.Logger('human-readable', false);
                [actZvUfv, actZv_QUALITY_FLAG, actZv_L2_QUALITY_BITMASK] = ...
                    bicas.proc.L1L2.qual.modify_quality_filter(ZvIn, isLfr, NsoTable, SETTINGS, L);
                testCase.verifyEqual(actZvUfv,                 expZvOut.ufv)
                testCase.verifyEqual(actZv_QUALITY_FLAG,       expZvOut.QUALITY_FLAG)
                testCase.verifyEqual(actZv_L2_QUALITY_BITMASK, expZvOut.L2_QUALITY_BITMASK)
            end

            %===================================================================

            ENA = zeros(0, 1);
            ECA = cell(0, 1);
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
                        % LFR/TDS, empty/nont empty NSO table
                        Settings = struct(...
                            'rmMuxModesArray',     [1, 2], ...
                            'lfrMuxModeMarginSec', 1.5, ...
                            'tdsMuxModeMarginSec', 2.5);
                        ZvIn = struct(...
                            'Epoch',               int64(ENA), ...
                            'MUX_SET',             ENA, ...
                            'QUALITY_FLAG',        ENA, ...
                            'ufv',                 logical(ENA));
                        expZvOut = struct(...
                            'QUALITY_FLAG',        ENA, ...
                            'L2_QUALITY_BITMASK',  uint16(ENA), ...
                            'ufv',                 false(0, 1));
                        test(ZvIn, isLfr, NsoTable, Settings, expZvOut);
                    end
                end
            end

            if 1
                %================
                % "Complex test"
                %================
                % (1) UFV set by mux modes (LFR)
                % (2) thruster firings.
                NsoTable = bicas.NsoTable(int64(7e9), int64(8e9), {bicas.const.NSOID.THRUSTER_FIRING});
                Settings = struct(...
                    'rmMuxModesArray',     [1, 2], ...
                    'lfrMuxModeMarginSec', 1.5, ...
                    'tdsMuxModeMarginSec', 2.5);
                ZvIn = struct(...
                    'Epoch',               int64([0:9]'*1e9), ...
                    'MUX_SET',             [0, 0, 1, 2, 0, 0, 0, 0, 0, 0]', ...
                    'QUALITY_FLAG',        2*ones(10, 1), ...
                    'ufv',                 logical([0, 0, 0, 0, 0, 0, 0, 0, 0, 0]'));
                expZvOut = struct(...
                    'QUALITY_FLAG',        [2, 2, 2, 2, 2, 2, 2, 1, 1, 2]', ...
                    'L2_QUALITY_BITMASK',  uint16( [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]'), ...
                    'ufv',                 logical([0, 1, 1, 1, 1, 0, 0, 0, 0, 0]'));
                test(ZvIn, true, NsoTable, Settings, expZvOut);
            end


            if 1
                %================
                % "Complex test"
                %================
                % (1) PARTIAL_SATURATION
                % (2) FULL_SATURATION with
                % (3) time overlap between NSO events.
                % (LFR)
                % Bit constants.
                PS = bicas.const.L2QBM_PARTIAL_SATURATION;
                FS = bicas.const.L2QBM_FULL_SATURATION;

                NsoTable = bicas.NsoTable(...
                    int64([1, 2]'*1e9), ...
                    int64([2, 3]'*1e9), ...
                    {bicas.const.NSOID.PARTIAL_SATURATION, ...
                     bicas.const.NSOID.FULL_SATURATION}');
                Settings = struct(...
                    'rmMuxModesArray',     [1, 2], ...
                    'lfrMuxModeMarginSec', 1.5, ...
                    'tdsMuxModeMarginSec', 2.5);
                ZvIn = struct(...
                    'Epoch',               int64([0:4]'*1e9), ...
                    'MUX_SET',             [0, 0,  0,     0,     0]', ...
                    'QUALITY_FLAG',        [3, 3,  3,     3,     3]', ...
                    'ufv',                 logical([0, 0,  0,     0,     0]'));
                expZvOut = struct(...
                    'QUALITY_FLAG',        [3, 1,  0,     0,     3]', ...
                    'L2_QUALITY_BITMASK',  uint16( [0, PS, PS+FS, PS+FS, 0]'), ...
                    'ufv',                 logical([0, 0,  0,     0,     0]'));
                test(ZvIn, true, NsoTable, Settings, expZvOut);
            end

        end



        function test_set_voltage_current_fill_value(testCase)

            function test(...
                    zv_Epoch, zvUfv, zvDemuxerOutput, zvCurrentAAmpere, ...
                    expZvDemuxerOutput, expZvCurrentAAmpere)

                irf.assert.sizes(...
                    zv_Epoch,         [-1], ...
                    zvCurrentAAmpere, [-1, 3])
                L = bicas.Logger('none', false);

                [actZvDemuxerOutput, actZvCurrentAAmpere] = bicas.proc.L1L2.qual.set_voltage_current_fill_value(...
                    zv_Epoch, zvDemuxerOutput, zvCurrentAAmpere, zvUfv, L);

                testCase.verifyEqual(actZvDemuxerOutput,  expZvDemuxerOutput)
                testCase.verifyEqual(actZvCurrentAAmpere, expZvCurrentAAmpere)
            end

            %===================================================================

            % Empty data.
            test( ...
                int64(zeros(0, 1)), ...
                false(0, 1), ...
                bicas.proc.L1L2.qual___UTEST.DemuxerOutput(zeros(0, 1)), ...
                zeros(0, 3),  ...
                bicas.proc.L1L2.qual___UTEST.DemuxerOutput(zeros(0, 1)), ...
                zeros(0, 3)  ...
            )

            % Non-empty input data that is not altered.
            test( ...
                int64([10, 11, 12, 13, 14]'), ...
                false(5, 1), ...
                bicas.proc.L1L2.qual___UTEST.DemuxerOutput(zeros(5, 1)), ...
                zeros(5, 3),  ...
                bicas.proc.L1L2.qual___UTEST.DemuxerOutput(zeros(5, 1)), ...
                zeros(5, 3)  ...
            )

            % Non-empty input data that is altered.
            test( ...
                int64([10, 11]'), ...
                logical([0, 1]'), ...
                bicas.proc.L1L2.qual___UTEST.DemuxerOutput([1, 2]'), ...
                [1:3; 11:13],  ...
                bicas.proc.L1L2.qual___UTEST.DemuxerOutput([1, NaN]'), ...
                [1:3; NaN(1, 3)] ...
            )

        end



        function test_get_UFV_records_from_settings(testCase)

            function actZvUfv = call_func(zv_Epoch, zv_MUX_SET, isLfr, rmMuxModesArray, lfrMuxModeMarginSec, tdsMuxModeMarginSec)
                SETTINGS = bicas.create_default_SETTINGS();
                SETTINGS.override_value('PROCESSING.L2.REMOVE_DATA.MUX_MODES',             rmMuxModesArray,     'test')
                SETTINGS.override_value('PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S', lfrMuxModeMarginSec, 'test')
                SETTINGS.override_value('PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S', tdsMuxModeMarginSec, 'test')
                SETTINGS.make_read_only()

                L = bicas.Logger('none', false);

                actZvUfv = bicas.proc.L1L2.qual.get_UFV_records_from_settings(...
                    zv_Epoch, zv_MUX_SET, isLfr, SETTINGS, L);
            end

            function test(zv_Epoch, zv_MUX_SET, isLfr, rmMuxModesArray, lfrMuxModeMarginSec, tdsMuxModeMarginSec, expZvUfv)
                irf.assert.sizes(...
                    zv_Epoch,   [-1, 1], ...
                    zv_MUX_SET, [-1, 1], ...
                    expZvUfv,   [-1, 1] ...
                )

                actZvUfv = call_func(zv_Epoch, zv_MUX_SET, isLfr, rmMuxModesArray, lfrMuxModeMarginSec, tdsMuxModeMarginSec);

                testCase.verifyEqual(actZvUfv, expZvUfv)
            end
            %===================================================================

            if 1
                for isLfr = [false, true]
                    % NOTE: Test PROCESSING.L2.REMOVE_DATA.MUX_MODES = [] (0x0)
                    % as is likely to be set to.
                    for rmMuxModesArrayCa = {[], zeros(0,1), [2]', [1,2,3]'}
                        rmMuxModesArray = rmMuxModesArrayCa{1};

                        test(int64(zeros(0, 1)), zeros(0, 1), isLfr, rmMuxModesArray, 10, 20, false(0, 1))

                        % No UFV
                        test(...
                            int64([10, 11, 12]'), ...
                            [0, 0, 0]', ...
                            isLfr, rmMuxModesArray, 10, 20, [false, false, false]' ...
                        )
                        % All UFV
                        test(...
                            int64([10, 11, 12]'), ...
                            [0, 0, 0]', ...
                            isLfr, rmMuxModesArray, 10, 20, [false, false, false]' ...
                        )
                    end
                end
            end

            % LFR
            test(...
                int64([10:20]') * 1e9, ...
                [2, 2, 0, 0, 0, 0, 0, 0, 0, 3, 3]', ...
                true, [2, 3], 1.5, 2.5, ...
                logical([1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1])' ...
            )

            % TDS
            test(...
                int64([10:20]') * 1e9, ...
                [2, 2, 0, 0, 0, 0, 0, 0, 0, 3, 3]', ...
                false, [2, 3], 1.5, 2.5, ...
                logical([1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1])' ...
            )
        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Utility function for creating a simplified DemuxerOutput data
        % structure with identical values on every channel.
        function S = DemuxerOutput(v)
            assert(iscolumn(v))
            S = struct();
            for asidCa = bicas.proc.L1L2.AntennaSignalId.C.ALL_ASID_NAMES_CA(:)'
                S.(asidCa{1}) = v;
            end
        end



    end    % methods(Static, Access=private)



end
