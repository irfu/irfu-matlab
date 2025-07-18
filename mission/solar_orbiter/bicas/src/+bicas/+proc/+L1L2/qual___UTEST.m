%
% matlab.unittest automatic test code for bicas.proc.L1L2.qual, except for
% method bicas.proc.L1L2.qual.sliding_window_over_fraction().
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



    function test_get_saturation_QRCBs(testCase)
      % NOTE: Does not (truly) test the call to
      % bicas.proc.L1L2.qual.sliding_window_over_fraction().

      DATA_AR = [...
        0 0 0  0 0 0  0 0 0 ; ...
        1 0 0  0 0 0  0 0 0 ; ...
        0 1 0  0 0 0  0 0 0 ; ...
        0 0 1  0 0 0  0 0 0 ; ...
        0 0 0  0 0 0  0 0 0 ; ...
        0 0 0  1 0 0  0 0 0 ; ...
        0 0 0  0 1 0  0 0 0 ; ...
        0 0 0  0 0 1  0 0 0 ; ...
        0 0 0  0 0 0  0 0 0 ; ...
        0 0 0  0 0 0  1 0 0 ; ...
        0 0 0  0 0 0  0 1 0 ; ...
        0 0 0  0 0 0  0 0 1 ; ...
        0 0 0  0 0 0  0 0 0 ; ...
        1 1 1  1 1 1  0 0 0 ; ...
        % Below: Both DC_V12 and AC_V12 saturation, but still only one bit set.
        0 0 0  1 0 0  1 0 0 ; ...
        0 0 0  0 1 0  0 1 0 ; ...
        0 0 0  0 0 1  0 0 1 ; ...
      ];
      N_ROWS    = size(DATA_AR, 1);
      % NOTE: (Almost) only the size is important.
      TT2000_AR = int64((1:N_ROWS) * 1e9)';
      % (iRec, iSdid)
      VSIB_AR_ALL = logical(DATA_AR(:, 1:9));
      CHANNEL_SATURATION_AR = [...
        VSIB_AR_ALL(:, 1:3), ...
        VSIB_AR_ALL(:, 4:6) | VSIB_AR_ALL(:, 7:9)];
      FULL_SATURATION_AR    = any(VSIB_AR_ALL, 2);



      VsibZvm = bicas.utils.ZvMap(N_ROWS);
      function add_channel_VSIB(ssidStr, iVsib)
        ssid = bicas.proc.L1L2.const.C.SDID_DICT(ssidStr);
        VsibZvm.add(ssid, VSIB_AR_ALL(:, iVsib));
      end
      add_channel_VSIB("DC_V1",  1)
      add_channel_VSIB("DC_V2",  2)
      add_channel_VSIB("DC_V3",  3)
      add_channel_VSIB("DC_V12", 4)
      add_channel_VSIB("DC_V13", 5)
      add_channel_VSIB("DC_V23", 6)
      add_channel_VSIB("AC_V12", 7)
      add_channel_VSIB("AC_V13", 8)
      add_channel_VSIB("AC_V23", 9)

      ExpGlobalSaturationQrcbMap  = bicas.proc.QrcbMap(N_ROWS);

      ExpGlobalSaturationQrcbMap.add("FULL_SATURATION",    FULL_SATURATION_AR)
      ExpGlobalSaturationQrcbMap.add("PARTIAL_SATURATION", false(N_ROWS, 1))

      ExpChannelSaturationQrcbMap = bicas.proc.QrcbMap(N_ROWS);

      ExpChannelSaturationQrcbMap.add("FULL_SATURATION",    false(N_ROWS, 1))
      ExpChannelSaturationQrcbMap.add("PARTIAL_SATURATION", false(N_ROWS, 1))

      function add_channel_saturation_QRCB(qrcid, iCol)
        ExpGlobalSaturationQrcbMap.add( qrcid, false(N_ROWS, 1))

        qrcbAr = CHANNEL_SATURATION_AR(:, iCol);
        ExpChannelSaturationQrcbMap.add(qrcid, qrcbAr)

      end
      add_channel_saturation_QRCB("SATURATION_ZV_V1",  1)
      add_channel_saturation_QRCB("SATURATION_ZV_V2",  2)
      add_channel_saturation_QRCB("SATURATION_ZV_V3",  3)
      add_channel_saturation_QRCB("SATURATION_ZV_V12", 4)
      add_channel_saturation_QRCB("SATURATION_ZV_V13", 5)
      add_channel_saturation_QRCB("SATURATION_ZV_V23", 6)



      isSwf                     = false;
      vstbFractionThreshold     = 0.9;
      cwfSlidingWindowLengthSec = 1.01;



      % CALL TESTED FUNCTION
      ActQrcbMap = bicas.proc.L1L2.qual.get_saturation_QRCBs( ...
        TT2000_AR, "GLOBAL_SATURATION", VsibZvm, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec);

      testCase.assertEqual(ActQrcbMap, ExpGlobalSaturationQrcbMap)



      % CALL TESTED FUNCTION
      ActQrcbMap = bicas.proc.L1L2.qual.get_saturation_QRCBs( ...
      TT2000_AR, "CHANNEL_SATURATION", VsibZvm, isSwf, ...
      vstbFractionThreshold, cwfSlidingWindowLengthSec);

      testCase.assertEqual(ActQrcbMap, ExpChannelSaturationQrcbMap)
    end



    function test_set_5xBLTS_voltage_samples_FV(testCase)

      function test(samplesAvoltAr, ssidAr, QrcbMap, QrcsMap, bExpNan)
        expSamplesAvoltAr          = samplesAvoltAr;
        expSamplesAvoltAr(bExpNan) = NaN;

        actSamplesAvoltAr = bicas.proc.L1L2.qual.set_5xBLTS_voltage_samples_FV(...
          samplesAvoltAr, ssidAr, QrcbMap, QrcsMap);

        testCase.assertEqual(actSamplesAvoltAr, expSamplesAvoltAr)
      end

      S = bicas.proc.L1L2.const.C.SSID_DICT;



      function test_empty()
        QrcbMap = bicas.proc.QrcbMap(0);
        QrcsMap = containers.Map();
        test(...
          double.empty(0, 1), ...
          uint8.empty( 0, 1), ...
          QrcbMap, QrcsMap, ...
          logical.empty(0, 1))
      end

      function test_nonempty_samples_empty_maps()
        QrcbMap = bicas.proc.QrcbMap(2);
        QrcsMap = containers.Map();

        samplesAvoltAr = reshape(1:2*3, [2, 3]);
        ssidAr         = S(["DC_V1" "DC_V12" "DC_V23"; "DC_V1" "DC_V2" "DC_V3"]);
        test(...
          samplesAvoltAr, ...
          ssidAr, ...
          QrcbMap, QrcsMap, ...
          false(2, 3))
      end

      % Deliberately 3D arrays.
      function test_complex_3D()
        QrcbMap = bicas.proc.QrcbMap(4);
        QrcbMap.add("QRCID_1", logical([1 1 0 0]'))
        QrcbMap.add("QRCID_2", logical([0 1 1 0]'))

        QrcsMap = containers.Map();
        QrcsMap("QRCID_1") = bicas.proc.QrcSettingL1rL2(...
          bicas.const.QUALITY_FLAG_MAX, ...
          bicas.const.LxQBM_NONE, ...
          S(["DC_V1"; "DC_V12"]), ...
          bicas.const.QRCS_CURRENT_FV_NONE);
        QrcsMap("QRCID_2") = bicas.proc.QrcSettingL1rL2(...
          bicas.const.QUALITY_FLAG_MAX, ...
          bicas.const.LxQBM_NONE, ...
          S(["DC_V1"; "DC_V2" ]), ...
          bicas.const.QRCS_CURRENT_FV_NONE);

        samplesAvoltAr = reshape(1:4*3*2, [4, 3, 2]);
        ssidAr(:, :, 1) = S([...
          "DC_V1"  "DC_V2"  "DC_V3"; ...
          "DC_V12" "DC_V13" "DC_V23" ; ...
          "DC_V2"  "DC_V3"  "DC_V1"; ...
          "AC_V12" "AC_V13" "AC_V23"]);
        ssidAr(:, :, 2) = ssidAr(end:-1:1, :, 1);   % Reverse rows.

        bExpNan(:, :, 1) = logical([1 0 0; 1 0 0; 1 0 1; 0 0 0]);
        bExpNan(:, :, 2) = logical([0 0 0; 1 0 1; 0 0 0; 0 0 0]);

        test(...
          samplesAvoltAr, ...
          ssidAr, ...
          QrcbMap, QrcsMap, ...
          bExpNan)
      end

      test_empty()
      test_nonempty_samples_empty_maps()
      test_complex_3D()
    end



    function test_set_current_samples_FV(testCase)

      function test(QrcbMap, QrcsMap, currentAr, expCurrentAr)
        irf.assert.sizes(currentAr, [QrcbMap.nRecords, 3]);

        actCurrentAr = bicas.proc.L1L2.qual.set_current_samples_FV(...
          currentAr, QrcbMap, QrcsMap);

        testCase.verifyEqual(actCurrentAr, expCurrentAr)
      end

      %===================================================================

      % Empty data.
      function test_empty()
        QrcbMap = bicas.proc.QrcbMap(0);
        QrcsMap = containers.Map();
        test( ...
          QrcbMap, QrcsMap, ...
          zeros(0, 3),  ...
          zeros(0, 3)  ...
          )
      end

      % Non-empty input data that is not altered.
      function test_data_unaltered()
        % NC = No change
        function test_NC()
          test( ...
            QrcbMap, QrcsMap, ...
            zeros(5, 3),  ...
            zeros(5, 3)  ...
            )
        end

        % Empty maps.
        QrcbMap = bicas.proc.QrcbMap(5);
        QrcsMap = containers.Map();
        test_NC()

        % QRCB=false, QRCS=all antennas FV
        QrcbMap.add("QRCID_1", false(5, 1))
        QrcsMap("QRCID_1") = bicas.proc.QrcSettingL1rL2(...
          bicas.const.QUALITY_FLAG_MAX, ...
          uint16(0), ...
          bicas.const.QRCS_VOLTAGE_FV_NONE, ...
          [1:3]');
        test_NC()

        % QRCB=true, QRCS=no antennas FV
        QrcbMap.set("QRCID_1", true(5, 1))
        QrcsMap("QRCID_1") = bicas.proc.QrcSettingL1rL2(...
          bicas.const.QUALITY_FLAG_MAX, ...
          uint16(0), ...
          bicas.const.QRCS_VOLTAGE_FV_NONE, ...
          bicas.const.QRCS_VOLTAGE_FV_NONE);
        test_NC()
      end

      function test_data_altered_unaltered()
        % Non-empty input data that is altered: 2 unaltered + (1+2) altered
        QrcbMap = bicas.proc.QrcbMap(5);
        QrcsMap = containers.Map();

        QrcbMap.add("QRCID_1",  logical([0 1 0 0 0]'))
        QrcbMap.add("QRCID_23", logical([0 0 0 1 1]'))


        QrcsMap("QRCID_1") = bicas.proc.QrcSettingL1rL2(...
          bicas.const.QUALITY_FLAG_MAX, ...
          uint16(0), ...
          bicas.const.QRCS_VOLTAGE_FV_NONE, ...
          [1]');
        QrcsMap("QRCID_23") = bicas.proc.QrcSettingL1rL2(...
          bicas.const.QUALITY_FLAG_MAX, ...
          uint16(0), ...
          bicas.const.QRCS_VOLTAGE_FV_NONE, ...
          [2 3]');

        N = NaN;
        test( ...
          QrcbMap, QrcsMap, ...
          [...
          1  2  3; ...
          11 12 13; ...
          21 22 23; ...
          31 32 33; ...
          41 42 43],  ...
          [ ...
          1  2  3; ...
          N  12 13; ...
          21 22 23; ...
          31 N  N; ...
          41 N  N])
      end

      test_empty()
      test_data_unaltered()
      test_data_altered_unaltered()
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Utility function for creating a simplified SamplesZvm with identical
    % values on every channel.
    function SamplesZvm = init_SamplesZvm(samplesAr)
      assert(iscolumn(samplesAr))

      SamplesZvm = bicas.utils.ZvMap(numel(samplesAr));
      for asrSdid = bicas.proc.L1L2.const.C.SDID_ASR_AR'
        SamplesZvm.add(asrSdid, samplesAr);
      end
    end



  end    % methods(Static, Access=private)



end
