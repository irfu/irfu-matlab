%
% matlab.unittest automatic test code for bicas.proc.qual.
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



        function test_NSO_table_to_QRC_flag_arrays(testCase)

            function test(allQrcidCa, NsoTable, Epoch, ExpQrcFlagsMap)
                % Normalize/modify arguments
                Epoch = int64(Epoch(:));

                L = bicas.Logger('human-readable', false);

                % CALL TESTED FUNCTION
                ActQrcFlagsMap = bicas.proc.qual.NSO_table_to_QRC_flag_arrays(...
                    allQrcidCa, NsoTable, Epoch, L);

                % ASSERT EXPECTED RESULT
                % ----------------------
                % IMPLEMENTATION NOTE: testCase.assertEqual() (and isequaln())
                % can handle containers.Map, but that is not very helpful for
                % debugging by understanding any found difference between the
                % two maps. Therefore explicitly comparing the map subcomponents.
                testCase.assertEqual(...
                    sort(ActQrcFlagsMap.keys), ...
                    sort(ExpQrcFlagsMap.keys))

                qrcidCa = ActQrcFlagsMap.keys;
                for i = 1:numel(qrcidCa)
                    qrcid = qrcidCa{i};
                    testCase.assertEqual(...
                        ActQrcFlagsMap(qrcid), ...
                        ExpQrcFlagsMap(qrcid))
                end
            end

            ALL_ENABLED = true;
            %ALL_ENABLED = false;

            %====================================
            % Empty NSO table. Various Epoch ZVs
            %====================================
            if ALL_ENABLED
                EMPTY_NSO_TABLE = bicas.NsoTable(...
                    int64(zeros(0, 1)), ...
                    int64(zeros(0, 1)), ...
                    cell(0, 1));
                EPOCH_DOUBLE_CA = {zeros(0,1), [10], [10;20;30]};
                for i = 1:numel(EPOCH_DOUBLE_CA)
                    Epoch_double = EPOCH_DOUBLE_CA{i};
                    ExpQrcFlagsMap = containers.Map();
                    test(...
                        {}, EMPTY_NSO_TABLE, ...
                        Epoch_double, ...
                        ExpQrcFlagsMap ...
                    )
                end
            end

            %=========================================
            % Two non-overlapping NSOs, one at a time
            %=========================================
            % Nontrivial NSO settings.
            ALL_QRCID_CA = {'QRCID1', 'QRCID2'};
            NSO_TABLE = bicas.NsoTable(...
                int64([1, 4]'*1e9), ...
                int64([2, 5]'*1e9), ...
                {'QRCID1', 'QRCID2'}');

            if ALL_ENABLED
                % Time interval is superset of NSO 1/2.
                ExpQrcFlagsMap = containers.Map();
                ExpQrcFlagsMap('QRCID1') = logical([0 1 1 0]');
                ExpQrcFlagsMap('QRCID2') = false(4,1);
                test(...
                    ALL_QRCID_CA, NSO_TABLE, ...
                    [0:3]*1e9, ...
                    ExpQrcFlagsMap ...
                );
            end

            if ALL_ENABLED
                % Time interval is superset of NSO 2/2.
                ExpQrcFlagsMap = containers.Map();
                ExpQrcFlagsMap('QRCID1') = false(4,1);
                ExpQrcFlagsMap('QRCID2') = logical([0 1 1 0]');
                test(...
                    ALL_QRCID_CA, NSO_TABLE, ...
                    [3:6]*1e9, ...
                    ExpQrcFlagsMap ...
                );
            end

            if ALL_ENABLED
                % Time interval from middle of NSO 1 to middle of NSO 2.
                ExpQrcFlagsMap = containers.Map();
                ExpQrcFlagsMap('QRCID1') = logical([1 0 0]');
                ExpQrcFlagsMap('QRCID2') = logical([0 0 1]');
                test(...
                    ALL_QRCID_CA, NSO_TABLE, ...
                    [2:4]'*1e9, ...
                    ExpQrcFlagsMap ...
                );
            end

            %========================================
            % Two overlapping NSOs, one unused QRCID
            %========================================
            ALL_QRCID_CA = {'QRCID1', 'QRCID2', 'QRCID3'};
            NSO_TABLE = bicas.NsoTable(...
                int64([1, 2]'*1e9), ...
                int64([2, 3]'*1e9), ...
                {'QRCID1', 'QRCID2'}');

            if ALL_ENABLED
                % Time interval covers all NSOs.
                ExpQrcFlagsMap = containers.Map();
                ExpQrcFlagsMap('QRCID1') = logical([0 1 1 0 0]');
                ExpQrcFlagsMap('QRCID2') = logical([0 0 1 1 0]');
                ExpQrcFlagsMap('QRCID3') = logical([0 0 0 0 0]');
                test(...
                    ALL_QRCID_CA, NSO_TABLE, ...
                    [0:4]*1e9, ...
                    ExpQrcFlagsMap ...
                );
            end

            if ALL_ENABLED
                % Epoch does not overlap with any NSOs (though time interval does).
                ExpQrcFlagsMap = containers.Map();
                ExpQrcFlagsMap('QRCID1') = logical([0 0]');
                ExpQrcFlagsMap('QRCID2') = logical([0 0]');
                ExpQrcFlagsMap('QRCID3') = logical([0 0]');
                test(...
                    ALL_QRCID_CA, NSO_TABLE, ...
                    [-1, 4]*1e9, ...
                    ExpQrcFlagsMap ...
                );
            end
        end



        function test_QRC_flag_arrays_to_quality_ZVs(testCase)

            function test(nRec, QrcFlagsMap, QrcSettingsMap, ...
                    exp_QUALITY_FLAG, exp_L2_QUALITY_BITMASK)
                exp_QUALITY_FLAG       = uint8( exp_QUALITY_FLAG(:));
                exp_L2_QUALITY_BITMASK = uint16(exp_L2_QUALITY_BITMASK(:));

                % CALL TESTED FUNCTION
                [act_QUALITY_FLAG, act_L2_QUALITY_BITMASK] = ...
                    bicas.proc.qual.QRC_flag_arrays_to_quality_ZVs(...
                        nRec, QrcFlagsMap, QrcSettingsMap);

                testCase.assertEqual(act_QUALITY_FLAG,       exp_QUALITY_FLAG)
                testCase.assertEqual(act_L2_QUALITY_BITMASK, exp_L2_QUALITY_BITMASK)
            end

            % =======================
            % Zero QRCIDs are defined
            % =======================
            QrcSettingsMap = containers.Map();
            QrcFlagsMap    = containers.Map();

            % Zero records
            test(0, QrcFlagsMap, QrcSettingsMap, ...
                [], [] ...
            )

            % Non-zero records
            test(3, QrcFlagsMap, QrcSettingsMap, ...
                4*ones(3,1), zeros(3,1) ...
            )

            % ==========================
            % Several QRCIDs are defined
            % ==========================
            QrcSettingsMap = containers.Map();
            QrcSettingsMap('QRCID1') = bicas.proc.QrcSetting(uint8(2), uint16(2));
            QrcSettingsMap('QRCID2') = bicas.proc.QrcSetting(uint8(3), uint16(4));

            % Zero records
            QrcFlagsMap = containers.Map();
            QrcFlagsMap('QRCID1') = false(0, 1);
            QrcFlagsMap('QRCID2') = false(0, 1);
            test(0, QrcFlagsMap, QrcSettingsMap, ...
                [], [] ...
            )

            % Non-zero records
            QrcFlagsMap = containers.Map();
            QrcFlagsMap('QRCID1') = logical([0 0 1 1]');
            QrcFlagsMap('QRCID2') = logical([0 1 0 1]');
            test(4, QrcFlagsMap, QrcSettingsMap, ...
                [4 3 2 2], [0 4 2 4+2] ...
            )
        end



    end    % methods(Test)



end
