%
% matlab.unittest automatic test code for bicas.proc.L2L3.qrc.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qrc___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Separate test class for tests of get_synthetic_L2_QFL().
  %   NOTE: File total 371 rows.



  %############
  %############
  % PROPERTIES
  %############
  %############
  % Additional properties of testCase objects. Needed for setup and teardown
  % methods which store/read their own data from the testCase object.
  properties
    L
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(T)
      T.L = bicas.Logger('HUMAN_READABLE', false);
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_L2QBM_to_QRCBs(testCase)

      function test(l2QbmAr, Qrcsm, ExpQrcbm)
        ActQrcbm = bicas.proc.L2L3.qrc.L2QBM_to_QRCBs(l2QbmAr, Qrcsm);

        testCase.assertEqual(ActQrcbm, ExpQrcbm);
      end

      % Zero QRCSs, zero-length L2QBM
      test( ...
        uint16.empty(0, 1), ...
        bicas.proc.QrcSettingsMap(), ...
        bicas.proc.QrcbMap(0))

      % Zero QRCSs, non-zero length L2QBM
      test( ...
        uint16([0:15]'), ...
        bicas.proc.QrcSettingsMap(), ...
        bicas.proc.QrcbMap(16))

      % Non-empty QRCSM, non-zero length L2QBM
      Qrcsm = bicas.proc.QrcSettingsMap();

      Qrcs = bicas.proc.QrcSettingL2(l2qbm=uint16(1));
      Qrcsm.add("QRCID_1", Qrcs);
      Qrcs = bicas.proc.QrcSettingL2(l2qbm=uint16(4));
      Qrcsm.add("QRCID_2", Qrcs);

      ExpQrcbm = bicas.proc.QrcbMap(8);
      ExpQrcbm.add("QRCID_1", logical([0 1 0 1 0 1 0 1]'))
      ExpQrcbm.add("QRCID_2", logical([0 0 0 0 1 1 1 1]'))

      test( ...
        uint16(0:7)', ...
        Qrcsm, ...
        ExpQrcbm)
    end



    %======================
    %======================
    % get_synthetic_L2_QFL
    %======================
    %======================
    % ANT3_FAILING:     QFL<=1
    % SATURATION_ZV_V*: QFL =0
    % THRUSTER_FIRING:  QFL<=1
    % ANT3_FAILING:     QFL<=1
    % --
    % NOTE: Not testing QFL fill positions.



    function test_get_synthetic_L2_QFL___nonexception_QRC(T)
      tt2000Ar = int64([11 12 13 14 15   31 32 33 34 35]');
      nRec     = numel(tt2000Ar);
      L2QflFpa = bicas.utils.FPArray( ...
        uint8([4 3 2 1 0   4 3 2 1 0]'), 'FILL_VALUE', uint8(255));

      ChannelSaturationQrcbm = T.get_channel_saturation_QRCBM_false(nRec);

      NsoTable = bicas.NsoTable( ...
        int64([10; 30]), ...
        int64([20; 40]), ...
        ["THRUSTER_FIRING" "SWEEP"]');

      % CALL TESTED FUNCTION
      SyntheticL2QflFpa = bicas.proc.L2L3.qrc.get_synthetic_L2_QFL( ...
        tt2000Ar, NsoTable, ChannelSaturationQrcbm, L2QflFpa, T.L);

      T.assertEqual(SyntheticL2QflFpa.array(uint8(0)), uint8(  [4 3 2 1 0   4 3 2 1 0]'))
      T.assertEqual(SyntheticL2QflFpa.fpAr,            logical([0 0 0 0 0   0 0 0 0 0]'))
    end



    % NOTE: L2 QFL=1 and QFL=0 during ANT3_FAILING (QFL<=1) test whether
    % QFL=0 is recognized as timestamp at which to not use synthetic value.
    function test_get_synthetic_L2_QFL___exception_QRCs(T)
      tt2000Ar = int64([11 12 13 14 15   31 32 33 34 35]');
      nRec     = numel(tt2000Ar);
      L2QflFpa = bicas.utils.FPArray( ...
        uint8([4 3 2 1 0   4 3 2 1 0]'), 'FILL_VALUE', uint8(255));

      ChannelSaturationQrcbm = T.get_channel_saturation_QRCBM_false(nRec);

      NsoTable = bicas.NsoTable( ...
        int64([10; 30]), ...
        int64([20; 40]), ...
        ["ANT3_FAILING", "SATURATION_ZV_V1"]');

      % CALL TESTED FUNCTION
      SyntheticL2QflFpa = bicas.proc.L2L3.qrc.get_synthetic_L2_QFL( ...
        tt2000Ar, NsoTable, ChannelSaturationQrcbm, L2QflFpa, T.L);

      T.assertEqual(SyntheticL2QflFpa.array(uint8(0)), uint8(  [4 4 4 4 0   4 4 4 4 4]'))
      T.assertEqual(SyntheticL2QflFpa.fpAr,            logical([0 0 0 0 0   0 0 0 0 0]'))
    end



    % Test absence of QRCs.
    % QFL=FP
    function test_get_synthetic_L2_QFL___no_QRC_events(T)
      tt2000Ar = int64([11 12 13 14 15   20]');
      nRec     = numel(tt2000Ar);
      L2QflFpa = bicas.utils.FPArray( ...
        uint8([4 3 2 1 0   255]'), 'FILL_VALUE', uint8(255));

      ChannelSaturationQrcbm = T.get_channel_saturation_QRCBM_false(nRec);

      NsoTable = bicas.NsoTable( ...
        int64.empty( 0, 1), ...
        int64.empty( 0, 1), ...
        string.empty(0, 1));

      % CALL TESTED FUNCTION
      SyntheticL2QflFpa = bicas.proc.L2L3.qrc.get_synthetic_L2_QFL( ...
        tt2000Ar, NsoTable, ChannelSaturationQrcbm, L2QflFpa, T.L);

      T.assertEqual(SyntheticL2QflFpa.array(uint8(0)), uint8(  [4 3 2 1 0   0]'))
      T.assertEqual(SyntheticL2QflFpa.fpAr,            logical([0 0 0 0 0   1]'))
    end



    % Test every case (?) of one QRC at a time.
    function test_get_synthetic_L2_QFL___every_QRCID_one_simultaneous(T)
      tt2000Ar = int64([11 12   31 32   51 52   100 101]');
      nRec     = numel(tt2000Ar);
      L2QflFpa = bicas.utils.FPArray( ...
        uint8([1 0   1 0   1 0   1 0]'), 'FILL_VALUE', uint8(255));
      % NOTE: L2 QFL=1 and QFL=0 during ANT3_FAILING (QFL<=1) to test whether
      % QFL=0 is recognized to not use synthetic values for.

      ChannelSaturationQrcbm = T.get_channel_saturation_QRCBM_false(nRec);
      ChannelSaturationQrcbm.set(...
        "SATURATION_ZV_V1", logical([0 0   0 0   0 0   1 1])')

      NsoTable = bicas.NsoTable( ...
        int64([10; 30; 50]), ...
        int64([20; 40; 60]), ...
        ["THRUSTER_FIRING", "ANT3_FAILING", "SATURATION_ZV_V1"]');

      % CALL TESTED FUNCTION
      SyntheticL2QflFpa = bicas.proc.L2L3.qrc.get_synthetic_L2_QFL( ...
        tt2000Ar, NsoTable, ChannelSaturationQrcbm, L2QflFpa, T.L);

      T.assertEqual(SyntheticL2QflFpa.array(uint8(0)), uint8(  [1 0   4 0   4 4   4 4]'))
      T.assertEqual(SyntheticL2QflFpa.fpAr,            logical([0 0   0 0   0 0   0 0]'))
    end



    % 1x non-exception + 2x exception QRC simultaneously.
    function test_get_synthetic_L2_QFL___mult_simultaneous_QRC(T)
      tt2000Ar = int64([11 12 13 14 15]');
      nRec     = numel(tt2000Ar);
      L2QflFpa = bicas.utils.FPArray( ...
        uint8([4 3 2 1 0]'), 'FILL_VALUE', uint8(255));

      ChannelSaturationQrcbm = T.get_channel_saturation_QRCBM_false(nRec);

      NsoTable = bicas.NsoTable( ...
        int64([10; 10; 10]), ...
        int64([20; 20; 20]), ...
        ["THRUSTER_FIRING", "ANT3_FAILING", "SATURATION_ZV_V1"]');

      % CALL TESTED FUNCTION
      SyntheticL2QflFpa = bicas.proc.L2L3.qrc.get_synthetic_L2_QFL( ...
        tt2000Ar, NsoTable, ChannelSaturationQrcbm, L2QflFpa, T.L);

      % THRUSTER_FIRING overrides.
      T.assertEqual(SyntheticL2QflFpa.array(uint8(0)), uint8(  [1 1 1 1 1]'))
      T.assertEqual(SyntheticL2QflFpa.fpAr,            logical([0 0 0 0 0]'))
    end



    %=======
    %=======
    % Other
    %=======
    %=======



    function test_L2QBM_to_channel_saturation_QRCBs(T)

      function test(l2qbmAr, saturationQualitySchemeId, ExpQrcbm)
        ActQrcbm = bicas.proc.L2L3.qrc.L2QBM_to_channel_saturation_QRCBs(...
          l2qbmAr, saturationQualitySchemeId);

        if 0
          % Separately compare parts of QRCBM to make debugging easier.
          T.assertEqual(ActQrcbm.qrcidAr, ExpQrcbm.qrcidAr)
          for qrcid = ActQrcbm.qrcidAr'
            qrcid
            T.assertEqual(ActQrcbm.get(qrcid), ExpQrcbm.get(qrcid))
          end
        end

        T.assertEqual(ActQrcbm, ExpQrcbm)
      end



      function test_GLOBAL_SATURATION()
        l2qbmAr = uint16(2.^[0:15]');
        ExcQrcbm = bicas.proc.QrcbMap(numel(l2qbmAr));
        ExcQrcbm.add_false(bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR)
        test(l2qbmAr, "GLOBAL_SATURATION", ExcQrcbm)
      end

      function test_CHANNEL_SATURATION()
        % NOTE: Tests assumes locations of quality bits. Could be fixed, but
        % then it still becomes hard to test non-saturation bits.
        l2qbmAr = uint16(2.^[4:6]');
        ExcQrcbm = bicas.proc.QrcbMap(numel(l2qbmAr));
        ExcQrcbm.add_false(bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR)
        ExcQrcbm.set("SATURATION_ZV_V13", logical([1 0 0]'))
        ExcQrcbm.set("SATURATION_ZV_V23", logical([0 1 0]'))
        test(l2qbmAr, "CHANNEL_SATURATION", ExcQrcbm)
      end

      test_GLOBAL_SATURATION()
      test_CHANNEL_SATURATION()
    end



    function test_set_FPA_samples_FP(T)
      % PROPOSAL: Split up in separate test functions.

      function test(v, Qrcbm, Qrcsm, qrcsFieldName, expV)
        FV     = cast(-1, class(v));
        Fpa    = bicas.utils.FPArray(v,    'FILL_VALUE', FV);
        ExpFpa = bicas.utils.FPArray(expV, 'FILL_VALUE', FV);

        ActFpa = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
          Fpa, Qrcbm, Qrcsm, qrcsFieldName);

        T.assertEqual(ActFpa, ExpFpa)
      end

      function test_zero_rows()
        Qrcbm = bicas.proc.QrcbMap(0);
        Qrcbm.add("QRCID_1", logical.empty(0, 1))
        Qrcbm.add("QRCID_2", logical.empty(0, 1))
        v    = zeros(0, 3);
        expV = zeros(0, 3);

        Qrcsm = bicas.proc.QrcSettingsMap();
        Qrcsm.add("QRCID_1", bicas.proc.QrcSettingL3(vdcFvIndexAr=[1 2]'))
        Qrcsm.add("QRCID_2", bicas.proc.QrcSettingL3(vdcFvIndexAr=[2 3]'))

        test(v, Qrcbm, Qrcsm, "vdcFvIndexAr", expV)
      end

      function test_zero_QRCBs_QRCSs()
        DATA = [
          -1  2  3; ...
          2 -1  4; ...
          3  4  5; ...
          4  5  6; ...
          ];
        Qrcbm = bicas.proc.QrcbMap(size(DATA, 1));
        v    = DATA(:,  1:3);
        expV = v;

        Qrcsm = bicas.proc.QrcSettingsMap();

        test(v, Qrcbm, Qrcsm, "vdcFvIndexAr", expV)
      end

      function test_complex(mc)
        DATA = [
          0 1   -1  2  3   -1  2 -1; ...
          1 0    2 -1  4   -1 -1  4; ...
          1 0    3  4  5   -1 -1  5; ...
          1 1    4  5  6   -1 -1 -1; ...
          ];
        Qrcbm = bicas.proc.QrcbMap(size(DATA, 1));
        Qrcbm.add("QRCID_1", logical(DATA(:, 1)))
        Qrcbm.add("QRCID_2", logical(DATA(:, 2)))
        v    = cast(DATA(:, 3:5), mc);
        expV = cast(DATA(:, 6:8), mc);

        Qrcsm = bicas.proc.QrcSettingsMap();
        Qrcsm.add("QRCID_1", bicas.proc.QrcSettingL3(vdcFvIndexAr=[1 2]'))
        Qrcsm.add("QRCID_2", bicas.proc.QrcSettingL3(vdcFvIndexAr=[  3]'))

        test(v, Qrcbm, Qrcsm, "vdcFvIndexAr", expV)
      end

      test_zero_rows()
      test_zero_QRCBs_QRCSs()
      test_complex("double")
      test_complex("single")
      test_complex("int8")
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function ChannelSaturationQrcbm = get_channel_saturation_QRCBM_false(nRec)
      qrcbFalseAr = false(nRec, 1);

      ChannelSaturationQrcbm = bicas.proc.QrcbMap(nRec);
      ChannelSaturationQrcbm.add("SATURATION_ZV_V1",  qrcbFalseAr)
      ChannelSaturationQrcbm.add("SATURATION_ZV_V2",  qrcbFalseAr)
      ChannelSaturationQrcbm.add("SATURATION_ZV_V3",  qrcbFalseAr)
      ChannelSaturationQrcbm.add("SATURATION_ZV_V12", qrcbFalseAr)
      ChannelSaturationQrcbm.add("SATURATION_ZV_V13", qrcbFalseAr)
      ChannelSaturationQrcbm.add("SATURATION_ZV_V23", qrcbFalseAr)
    end



  end    % methods(Static, Access=private)



end
