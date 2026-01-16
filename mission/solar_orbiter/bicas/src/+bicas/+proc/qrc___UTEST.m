%
% matlab.unittest automatic test code for bicas.proc.qrc.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qrc___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_NSO_table_to_QRCBM___empty_NSO_table(T)
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

          ExpQrcbm = bicas.proc.QrcbMap(numel(tt2000Ar));
          for qrcid = requestedQrcidAr'
            ExpQrcbm.add(qrcid, false(size(tt2000Ar)));
          end

          T.test_NSO_table_to_QRCBM(...
            requestedQrcidAr, EMPTY_NSO_TABLE, ...
            tt2000Ar, ...
            ExpQrcbm ...
            )
        end
      end
    end



    function test_NSO_table_to_QRCBM___nonoverlapping_events(T)
      % Two non-overlapping NSO events.

      % Nontrivial NSO table.
      ALL_QRCID_AR = ["QRCID1", "QRCID2"];
      NSO_TABLE = bicas.NsoTable(...
        int64([1, 4]'*1e9), ...
        int64([2, 5]'*1e9), ...
        ["QRCID1", "QRCID2"]');

      % Time interval is superset of NSO event 1/2.
      ExpQrcbm = bicas.proc.QrcbMap(4);
      ExpQrcbm.add("QRCID1", logical([0 1 1 0]'));
      ExpQrcbm.add("QRCID2", logical([0 0 0 0]'));
      T.test_NSO_table_to_QRCBM(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [0:3]*1e9, ...
        ExpQrcbm ...
        );

      % Time interval is superset of NSO event 2/2.
      ExpQrcbm = bicas.proc.QrcbMap(4);
      ExpQrcbm.add("QRCID1", logical([0 0 0 0]'));
      ExpQrcbm.add("QRCID2", logical([0 1 1 0]'));
      T.test_NSO_table_to_QRCBM(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [3:6]*1e9, ...
        ExpQrcbm ...
        );

      % Time interval from middle of NSO event 1 to middle of NSO event 2.
      ExpQrcbm = bicas.proc.QrcbMap(3);
      ExpQrcbm.add("QRCID1", logical([1 0 0]'));
      ExpQrcbm.add("QRCID2", logical([0 0 1]'));
      T.test_NSO_table_to_QRCBM(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [2:4]'*1e9, ...
        ExpQrcbm ...
        );
    end



    function test_NSO_table_to_QRCBM___overlapping_events_nonexistent(T)
      % Two overlapping NSOs, one requested non-existing QRCID.

      ALL_QRCID_AR = ["QRCID1", "QRCID2", "QRCID3"];
      NSO_TABLE = bicas.NsoTable(...
        int64([1, 2]'*1e9), ...
        int64([2, 3]'*1e9), ...
        ["QRCID1", "QRCID2"]');

      % Time interval covers all NSO events.
      ExpQrcbm = bicas.proc.QrcbMap(5);
      ExpQrcbm.add("QRCID1", logical([0 1 1 0 0]'));
      ExpQrcbm.add("QRCID2", logical([0 0 1 1 0]'));
      ExpQrcbm.add("QRCID3", logical([0 0 0 0 0]'));
      T.test_NSO_table_to_QRCBM(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [0:4]*1e9, ...
        ExpQrcbm ...
        );

      % Epoch does not overlap with any NSO events (though time interval does).
      ExpQrcbm = bicas.proc.QrcbMap(2);
      ExpQrcbm.add("QRCID1", logical([0 0]'));
      ExpQrcbm.add("QRCID2", logical([0 0]'));
      ExpQrcbm.add("QRCID3", logical([0 0]'));
      T.test_NSO_table_to_QRCBM(...
        ALL_QRCID_AR, NSO_TABLE, ...
        [-1, 4]*1e9, ...
        ExpQrcbm ...
        );
    end


    function test_NSO_table_to_QRCBM___QRCID_req_omitted_req_nonexist(T)
      % (1) Request QRCID which does not exist in NSO table.
      % (2) Omit to request QRCID in NSO table.

      NSO_TABLE = bicas.NsoTable(...
        int64([1, 2]'*1e9), ...
        int64([2, 3]'*1e9), ...
        ["QRCID1", "QRCID2"]');

      ExpQrcbm = bicas.proc.QrcbMap(6);
      ExpQrcbm.add("QRCID1", logical([0 1 1 0 0 0]'));
      ExpQrcbm.add("QRCID3", logical([0 0 0 0 0 0]'));
      T.test_NSO_table_to_QRCBM(...
        string(ExpQrcbm.qrcidAr), NSO_TABLE, ...
        [0:5]*1e9, ...
        ExpQrcbm ...
        );
    end



    function test_QRCB_arrays_to_quality_ZVs(T)

      % Test function shared between tests.
      function test(Qrcbm, Qrcsm, lxqbmName, expQfl, expLxqbm)
        expQfl   = uint8( expQfl(:));
        expLxqbm = uint16(expLxqbm(:));

        % CALL TESTED FUNCTION
        [actQfl, actLxqbm] = bicas.proc.qrc.QRCB_arrays_to_quality_ZVs( ...
          Qrcbm, Qrcsm, lxqbmName);

        T.assertEqual(actQfl,   expQfl)
        T.assertEqual(actLxqbm, expLxqbm)
      end



      % Zero QRCIDs defined
      function test_zero_QRCIDs()
        Qrcbm = bicas.proc.QrcbMap(0);
        Qrcsm = bicas.proc.QrcSettingsMap();

        % Zero records
        test(Qrcbm, Qrcsm, "L2_QUALITY_BITMASK", [], [])

        % Non-zero records
        Qrcbm = bicas.proc.QrcbMap(3);
        test(Qrcbm, Qrcsm, "L2_QUALITY_BITMASK", ...
          4*ones(3,1), zeros(3,1))
      end

      % Several QRCIDs are defined
      function test_nonzero_QRCIDs()
        Qrcsm = bicas.proc.QrcSettingsMap();

        Qrcs = bicas.proc.QrcSettingL2(qfl=uint8(2), l2qbm=uint16(2));
        Qrcsm.add("QRCID_1", Qrcs);

        Qrcs = bicas.proc.QrcSettingL2(qfl=uint8(3), l2qbm=uint16(4));
        Qrcsm.add("QRCID_2", Qrcs);



        % Zero records
        Qrcbm = bicas.proc.QrcbMap(0);
        Qrcbm.add("QRCID_1", false(0, 1));
        Qrcbm.add("QRCID_2", false(0, 1));
        test(Qrcbm, Qrcsm, "L2_QUALITY_BITMASK", ...
          [], [])

        % Non-zero records
        Qrcbm = bicas.proc.QrcbMap(4);
        Qrcbm.add("QRCID_1", logical([0 0 1 1]'));
        Qrcbm.add("QRCID_2", logical([0 1 0 1]'));
        test(Qrcbm, Qrcsm, "L2_QUALITY_BITMASK", ...
          [4 3 2 2]', [0 4 2 4+2]')
      end



      test_zero_QRCIDs()
      test_nonzero_QRCIDs()
    end



    function test_QRCB_arrays_to_GA_CAVEATS___zero_QRCs_zero_records(T)
      Qrcbm = bicas.proc.QrcbMap(0);
      Qrcsm = bicas.proc.QrcSettingsMap();
      actGaCaveats = bicas.proc.qrc.QRCB_arrays_to_GA_CAVEATS(Qrcbm, Qrcsm);
      T.assertEqual(actGaCaveats, string.empty(0, 1))
    end



    function test_QRCB_arrays_to_GA_CAVEATS___zero_QRCs_nonzero_records(T)
      Qrcbm = bicas.proc.QrcbMap(3);
      Qrcsm = bicas.proc.QrcSettingsMap();
      actGaCaveats = bicas.proc.qrc.QRCB_arrays_to_GA_CAVEATS(Qrcbm, Qrcsm);
      T.assertEqual(actGaCaveats, string.empty(0, 1))
    end



    function test_QRCB_arrays_to_GA_CAVEATS___empty_QRCS_CAVEATS(T)
      Qrcbm = bicas.proc.QrcbMap(3);
      Qrcbm.add("QRCID_1", logical([0; 1; 0]))
      Qrcbm.add("QRCID_2", logical([0; 0; 0]))
      Qrcbm.add("QRCID_3", logical([1; 0; 0]))

      Qrcsm = bicas.proc.QrcSettingsMap();
      Qrcsm.add("QRCID_1", bicas.proc.QrcSettingL2())
      Qrcsm.add("QRCID_2", bicas.proc.QrcSettingL2())
      Qrcsm.add("QRCID_3", bicas.proc.QrcSettingL2())

      actGaCaveats = bicas.proc.qrc.QRCB_arrays_to_GA_CAVEATS(Qrcbm, Qrcsm);
      T.assertEqual(actGaCaveats, string.empty(0, 1))
    end



    % Complex test
    function test_QRCB_arrays_to_GA_CAVEATS___complex(T)
      Qrcbm = bicas.proc.QrcbMap(3);
      Qrcbm.add("QRCID_1", logical([0; 1; 0]))
      Qrcbm.add("QRCID_2", logical([0; 0; 1]))
      Qrcbm.add("QRCID_3", logical([1; 0; 0]))
      Qrcbm.add("QRCID_4", logical([0; 0; 0]))

      Qrcsm = bicas.proc.QrcSettingsMap();
      % NOTE: Test (1) using >=2 CAVEATS strings per QRCS, (2) sorting of
      % combined CAVEATS strings, (3) empty CAVEATS (for matching QRC).
      gaCaveats1 = ["d CAVEATS 1a"; "b CAVEATS 1b"];
      gaCaveats2 = ["c CAVEATS 2a"; "a CAVEATS 2b"];
      Qrcsm.add("QRCID_1", bicas.proc.QrcSettingL2(gaCaveats=gaCaveats1))
      Qrcsm.add("QRCID_2", bicas.proc.QrcSettingL2(gaCaveats=gaCaveats2))
      Qrcsm.add("QRCID_3", bicas.proc.QrcSettingL2())
      Qrcsm.add("QRCID_4", bicas.proc.QrcSettingL3(gaCaveats=["Should not be used."]))

      actGaCaveats = bicas.proc.qrc.QRCB_arrays_to_GA_CAVEATS(Qrcbm, Qrcsm);
      T.assertEqual(actGaCaveats, sort([gaCaveats1; gaCaveats2]))
    end



    function test_LxQBM_to_bit_positions(T)
      function test(lxqbm, expBitPosAr)
        lxqbm       = uint16(lxqbm);

        actBitPosAr = bicas.proc.qrc.LxQBM_to_bit_positions(lxqbm);
        T.assertEqual(actBitPosAr, expBitPosAr)
      end

      % Zero set bits.
      test(0, zeros(0, 1))

      % One set bit.
      test(1, [0])
      test(2, [1])
      test(32768, [15])

      % Multiple set bits
      test(1+4+32768, [0 2 15]')
    end



    function test_filter_saturation_QRCBs___zero_QRCBs(T)
      Qrcbm = bicas.proc.QrcbMap(3);

      for saturationQualitySchemeId = ["GLOBAL_SATURATION", "CHANNEL_SATURATION"]
        ActQrcbm = bicas.proc.qrc.filter_saturation_QRCBs( ...
          Qrcbm, saturationQualitySchemeId);

        T.assertEqual(ActQrcbm, Qrcbm)
        T.assertFalse(ActQrcbm == Qrcbm)   % Test reference inequality.
      end
    end



    function test_filter_saturation_QRCBs___some_QRCBs_global_saturation(T)
      Qrcbm = bicas.proc.QrcbMap(3);
      Qrcbm.add("FULL_SATURATION",    true(3, 1));
      Qrcbm.add("PARTIAL_SATURATION", true(3, 1));
      Qrcbm.add("SATURATION_ZV_V3",   true(3, 1));
      Qrcbm.add("SWEEP",              true(3, 1));

      ExpQrcbm = bicas.proc.QrcbMap(3);
      ExpQrcbm.add("FULL_SATURATION",    true( 3, 1));
      ExpQrcbm.add("PARTIAL_SATURATION", true( 3, 1));
      ExpQrcbm.add("SATURATION_ZV_V3",   false(3, 1));
      ExpQrcbm.add("SWEEP",              true( 3, 1));

      ActQrcbm = bicas.proc.qrc.filter_saturation_QRCBs(Qrcbm, "GLOBAL_SATURATION");

      T.assertEqual(ActQrcbm, ExpQrcbm)
      T.assertFalse(ActQrcbm == Qrcbm)   % Test reference inequality.
    end



    function test_filter_saturation_QRCBs___some_QRCBs_channel_saturation(T)
      Qrcbm = bicas.proc.QrcbMap(3);
      Qrcbm.add("FULL_SATURATION",    true(3, 1));
      Qrcbm.add("PARTIAL_SATURATION", true(3, 1));
      Qrcbm.add("SATURATION_ZV_V3",   true(3, 1));
      Qrcbm.add("SWEEP",              true(3, 1));

      ExpQrcbm = bicas.proc.QrcbMap(3);
      ExpQrcbm.add("FULL_SATURATION",    false(3, 1));
      ExpQrcbm.add("PARTIAL_SATURATION", false(3, 1));
      ExpQrcbm.add("SATURATION_ZV_V3",   true( 3, 1));
      ExpQrcbm.add("SWEEP",              true( 3, 1));

      ActQrcbm = bicas.proc.qrc.filter_saturation_QRCBs(Qrcbm, "CHANNEL_SATURATION");

      T.assertEqual(ActQrcbm, ExpQrcbm)
      T.assertFalse(ActQrcbm == Qrcbm)   % Test reference inequality.
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function test_NSO_table_to_QRCBM(...
        T, requestedQrcidAr, NsoTable, tt2000Ar, ExpQrcbm)

      assert(isa(ExpQrcbm, "bicas.proc.QrcbMap"))

      % Normalize/modify arguments
      requestedQrcidAr = requestedQrcidAr(:);
      tt2000Ar         = int64(tt2000Ar(:));

      L = bicas.Logger('HUMAN_READABLE', false);

      % CALL TESTED FUNCTION
      ActQrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        requestedQrcidAr, NsoTable, tt2000Ar, L);

      % ASSERT EXPECTED RESULT
      % ----------------------
      % IMPLEMENTATION NOTE: T.assertEqual() (and isequaln()) can
      % handle containers.Map, but that is not very helpful for debugging by
      % understanding any found difference between the two maps. Therefore
      % explicitly comparing the map subcomponents.
      T.assertEqual(...
        sort(ActQrcbm.qrcidAr), ...
        sort(ExpQrcbm.qrcidAr))

      for qrcid = ActQrcbm.qrcidAr'
        T.assertEqual(...
          ActQrcbm.get(qrcid), ...
          ExpQrcbm.get(qrcid))
      end
    end
  end    % methods(Access=private)



end
