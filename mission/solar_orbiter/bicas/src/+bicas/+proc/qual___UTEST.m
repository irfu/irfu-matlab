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



    function test_NSO_table_to_QRCB_arrays___empty_NSO_table(testCase)
      % Empty NSO table. Various Epoch ZVs, QRCIDs.

      EMPTY_NSO_TABLE = bicas.NsoTable(...
        int64(zeros(0, 1)), ...
        int64(zeros(0, 1)), ...
        strings(0, 1));
      REQUESTED_QRCID_AR_CA = {
        string.empty(0, 1),    % Zero QRCIDs
        ["QRCID1"],
        ["QRCID1"; "QRCID2"]
        };
      EPOCH_AR_CA = {
        zeros(0,1),   % No timestamps
        [10],         % One timestamp
        [10;20;30]    % Multiple timestamps
        };

      for i = 1:numel(EPOCH_AR_CA)
        for j = 1:numel(REQUESTED_QRCID_AR_CA)
          tt2000Ar         = EPOCH_AR_CA{i};
          requestedQrcidAr = REQUESTED_QRCID_AR_CA{j};

          ExpQrcbMap = containers.Map();
          for qrcid = requestedQrcidAr'
            ExpQrcbMap(qrcid) = false(size(tt2000Ar));
          end

          testCase.test_NSO_table_to_QRCB_arrays(...
            requestedQrcidAr, EMPTY_NSO_TABLE, ...
            tt2000Ar, ...
            ExpQrcbMap ...
            )
        end
      end
    end



    function test_NSO_table_to_QRCB_arrays___nonoverlapping_events(testCase)
      % Two non-overlapping NSO events.

      % Nontrivial NSO table.
      ALL_QRCID_AR = ["QRCID1", "QRCID2"];
      NSO_TABLE = bicas.NsoTable(...
        int64([1, 4]'*1e9), ...
        int64([2, 5]'*1e9), ...
        ["QRCID1", "QRCID2"]');

      % Time interval is superset of NSO event 1/2.
      ExpQrcbMap = containers.Map();
      ExpQrcbMap("QRCID1") = logical([0 1 1 0]');
      ExpQrcbMap("QRCID2") = logical([0 0 0 0]');
      testCase.test_NSO_table_to_QRCB_arrays(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [0:3]*1e9, ...
        ExpQrcbMap ...
        );

      % Time interval is superset of NSO event 2/2.
      ExpQrcbMap = containers.Map();
      ExpQrcbMap("QRCID1") = logical([0 0 0 0]');
      ExpQrcbMap("QRCID2") = logical([0 1 1 0]');
      testCase.test_NSO_table_to_QRCB_arrays(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [3:6]*1e9, ...
        ExpQrcbMap ...
        );

      % Time interval from middle of NSO event 1 to middle of NSO event 2.
      ExpQrcbMap = containers.Map();
      ExpQrcbMap("QRCID1") = logical([1 0 0]');
      ExpQrcbMap("QRCID2") = logical([0 0 1]');
      testCase.test_NSO_table_to_QRCB_arrays(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [2:4]'*1e9, ...
        ExpQrcbMap ...
        );
    end



    function test_NSO_table_to_QRCB_arrays___overlapping_events_nonexistent(testCase)
      % Two overlapping NSOs, one requested non-existing QRCID.

      ALL_QRCID_AR = ["QRCID1", "QRCID2", "QRCID3"];
      NSO_TABLE = bicas.NsoTable(...
        int64([1, 2]'*1e9), ...
        int64([2, 3]'*1e9), ...
        ["QRCID1", "QRCID2"]');

      % Time interval covers all NSO events.
      ExpQrcbMap = containers.Map();
      ExpQrcbMap("QRCID1") = logical([0 1 1 0 0]');
      ExpQrcbMap("QRCID2") = logical([0 0 1 1 0]');
      ExpQrcbMap("QRCID3") = logical([0 0 0 0 0]');
      testCase.test_NSO_table_to_QRCB_arrays(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [0:4]*1e9, ...
        ExpQrcbMap ...
        );

      % Epoch does not overlap with any NSO events (though time interval does).
      ExpQrcbMap = containers.Map();
      ExpQrcbMap("QRCID1") = logical([0 0]');
      ExpQrcbMap("QRCID2") = logical([0 0]');
      ExpQrcbMap("QRCID3") = logical([0 0]');
      testCase.test_NSO_table_to_QRCB_arrays(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [-1, 4]*1e9, ...
        ExpQrcbMap ...
        );
    end



    function test_QRCB_arrays_to_quality_ZVs(testCase)

      function test(nRec, QrcbMap, QrcsMap, ...
          exp_QUALITY_FLAG, exp_L2_QUALITY_BITMASK)
        exp_QUALITY_FLAG       = uint8( exp_QUALITY_FLAG(:));
        exp_L2_QUALITY_BITMASK = uint16(exp_L2_QUALITY_BITMASK(:));

        % CALL TESTED FUNCTION
        [act_QUALITY_FLAG, act_L2_QUALITY_BITMASK] = ...
          bicas.proc.qual.QRCB_arrays_to_quality_ZVs(...
          nRec, QrcbMap, QrcsMap);

        testCase.assertEqual(act_QUALITY_FLAG,       exp_QUALITY_FLAG)
        testCase.assertEqual(act_L2_QUALITY_BITMASK, exp_L2_QUALITY_BITMASK)
      end

      % =======================
      % Zero QRCIDs are defined
      % =======================
      QrcsMap     = containers.Map();
      QrcbMap = containers.Map();

      % Zero records
      test(0, QrcbMap, QrcsMap, ...
        [], [] ...
        )

      % Non-zero records
      test(3, QrcbMap, QrcsMap, ...
        4*ones(3,1), zeros(3,1) ...
        )

      % ==========================
      % Several QRCIDs are defined
      % ==========================
      QrcsMap = containers.Map();
      QrcsMap('QRCID1') = bicas.proc.QrcSetting(uint8(2), uint16(2));
      QrcsMap('QRCID2') = bicas.proc.QrcSetting(uint8(3), uint16(4));

      % Zero records
      QrcbMap = containers.Map();
      QrcbMap('QRCID1') = false(0, 1);
      QrcbMap('QRCID2') = false(0, 1);
      test(0, QrcbMap, QrcsMap, ...
        [], [] ...
        )

      % Non-zero records
      QrcbMap = containers.Map();
      QrcbMap('QRCID1') = logical([0 0 1 1]');
      QrcbMap('QRCID2') = logical([0 1 0 1]');
      test(4, QrcbMap, QrcsMap, ...
        [4 3 2 2], [0 4 2 4+2] ...
        )
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



      function test_NSO_table_to_QRCB_arrays(...
          testCase, requestedQrcidAr, NsoTable, tt2000Ar, ExpQrcbMap)

        % Normalize/modify arguments
        requestedQrcidAr = requestedQrcidAr(:);
        tt2000Ar         = int64(tt2000Ar(:));

        L = bicas.Logger('HUMAN_READABLE', false);

        % CALL TESTED FUNCTION
        ActQrcbMap = bicas.proc.qual.NSO_table_to_QRCB_arrays(...
          requestedQrcidAr, NsoTable, tt2000Ar, L);

        % ASSERT EXPECTED RESULT
        % ----------------------
        % IMPLEMENTATION NOTE: testCase.assertEqual() (and isequaln()) can
        % handle containers.Map, but that is not very helpful for debugging by
        % understanding any found difference between the two maps. Therefore
        % explicitly comparing the map subcomponents.
        testCase.assertEqual(...
          sort(ActQrcbMap.keys), ...
          sort(ExpQrcbMap.keys))

        qrcidCa = ActQrcbMap.keys;
        for i = 1:numel(qrcidCa)
          qrcid = qrcidCa{i};
          testCase.assertEqual(...
            ActQrcbMap(qrcid), ...
            ExpQrcbMap(qrcid))
        end
      end
    end    % methods(Access=private)



end
