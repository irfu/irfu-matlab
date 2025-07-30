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



    function test_get_saturation_QRCBs(testCase)
      % NOTE: Does not (truly) test the call to
      % bicas.utils.sliding_window_over_fraction().

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

      % Expected value for schema GLOBAL_SATURATION.
      ExpGlobalSaturationQrcbm  = bicas.proc.QrcbMap(N_ROWS);
      ExpGlobalSaturationQrcbm.add(      "FULL_SATURATION", FULL_SATURATION_AR)
      ExpGlobalSaturationQrcbm.add_false("PARTIAL_SATURATION")

      % Expected value for schema CHANNEL_SATURATION.
      ExpChannelSaturationQrcbm = bicas.proc.QrcbMap(N_ROWS);
      ExpChannelSaturationQrcbm.add_false("FULL_SATURATION")
      ExpChannelSaturationQrcbm.add_false("PARTIAL_SATURATION")

      function add_channel_saturation_QRCB(qrcid, iCol)
        ExpGlobalSaturationQrcbm.add_false(qrcid)

        qrcbAr = CHANNEL_SATURATION_AR(:, iCol);
        ExpChannelSaturationQrcbm.add(qrcid, qrcbAr)
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
      ActQrcbm = bicas.proc.L1L2.qual.get_saturation_QRCBs( ...
        TT2000_AR, "GLOBAL_SATURATION", VsibZvm, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec);

      testCase.assertEqual(ActQrcbm, ExpGlobalSaturationQrcbm)



      % CALL TESTED FUNCTION
      ActQrcbm = bicas.proc.L1L2.qual.get_saturation_QRCBs( ...
        TT2000_AR, "CHANNEL_SATURATION", VsibZvm, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec);

      testCase.assertEqual(ActQrcbm, ExpChannelSaturationQrcbm)
    end



    function test_set_5xBLTS_voltage_samples_FV(testCase)

      function test(samplesAvoltAr, ssidAr, Qrcbm, Qrcsm, bExpNan)
        expSamplesAvoltAr          = samplesAvoltAr;
        expSamplesAvoltAr(bExpNan) = NaN;

        actSamplesAvoltAr = bicas.proc.L1L2.qual.set_5xBLTS_voltage_samples_FV(...
          samplesAvoltAr, ssidAr, Qrcbm, Qrcsm);

        testCase.assertEqual(actSamplesAvoltAr, expSamplesAvoltAr)
      end

      S = bicas.proc.L1L2.const.C.SSID_DICT;



      function test_empty()
        Qrcbm = bicas.proc.QrcbMap(0);
        Qrcsm   = bicas.proc.QrcSettingsMap();
        test(...
          double.empty(0, 1), ...
          uint8.empty( 0, 1), ...
          Qrcbm, Qrcsm, ...
          logical.empty(0, 1))
      end

      function test_nonempty_samples_empty_maps()
        Qrcbm = bicas.proc.QrcbMap(2);
        Qrcsm   = bicas.proc.QrcSettingsMap();

        samplesAvoltAr = reshape(1:2*3, [2, 3]);
        ssidAr         = S(["DC_V1" "DC_V12" "DC_V23"; "DC_V1" "DC_V2" "DC_V3"]);
        test(...
          samplesAvoltAr, ...
          ssidAr, ...
          Qrcbm, Qrcsm, ...
          false(2, 3))
      end

      % Deliberately 3D arrays.
      function test_complex_3D()
        Qrcbm = bicas.proc.QrcbMap(4);
        Qrcbm.add("QRCID_1", logical([1 1 0 0]'))
        Qrcbm.add("QRCID_2", logical([0 1 1 0]'))

        Qrcsm  = bicas.proc.QrcSettingsMap();
        Qrcs = bicas.proc.QrcSettingL2(voltageFvSsidAr=S(["DC_V1"; "DC_V12"]));
        Qrcsm.add("QRCID_1", Qrcs);
        Qrcs = bicas.proc.QrcSettingL2(voltageFvSsidAr=S(["DC_V1"; "DC_V2" ]));
        Qrcsm.add("QRCID_2", Qrcs);

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
          Qrcbm, Qrcsm, ...
          bExpNan)
      end

      test_empty()
      test_nonempty_samples_empty_maps()
      test_complex_3D()
    end



    function test_set_current_samples_FV(testCase)

      function test(Qrcbm, Qrcsm, currentAr, expCurrentAr)
        irf.assert.sizes(currentAr, [Qrcbm.nRecords, 3]);

        actCurrentAr = bicas.proc.L1L2.qual.set_current_samples_FV(...
          currentAr, Qrcbm, Qrcsm);

        testCase.verifyEqual(actCurrentAr, expCurrentAr)
      end

      %===================================================================

      % Empty data.
      function test_empty_empty()
        Qrcbm = bicas.proc.QrcbMap(0);
        Qrcsm   = bicas.proc.QrcSettingsMap();

        test( ...
          Qrcbm, Qrcsm, ...
          zeros(0, 3),  ...
          zeros(0, 3))
      end

      % Non-empty input data that is not altered.
      function test_data_unaltered()

        % NED = Non-empty data
        % NC  = No change (in output)
        function test_NED_NC()
          test( ...
            Qrcbm, Qrcsm, ...
            zeros(5, 3),  ...
            zeros(5, 3))
        end

        % Zero QRCB arrays, zero QRCSs.
        Qrcbm = bicas.proc.QrcbMap(5);
        Qrcsm   = bicas.proc.QrcSettingsMap();
        test_NED_NC()

        % QRCB=false, QRCS=all antennas FV
        Qrcbm.add("QRCID_1", false(5, 1))
        Qrcs = bicas.proc.QrcSettingL2(currentFvIantAr=[1:3]');
        Qrcsm.add("QRCID_1", Qrcs);
        test_NED_NC()

        % QRCB=true, QRCS=no antennas FV
        Qrcbm.set("QRCID_1", true(5, 1))
        Qrcsm = bicas.proc.QrcSettingsMap();
        Qrcs  = bicas.proc.QrcSettingL2();    % Default ==> Does nothing.
        Qrcsm.add("QRCID_1", Qrcs);
        test_NED_NC()
      end

      function test_data_altered_unaltered()
        % Non-empty input data that is altered: 2 unaltered + (1+2) altered
        Qrcbm = bicas.proc.QrcbMap(5);
        Qrcsm   = bicas.proc.QrcSettingsMap();

        Qrcbm.add("QRCID_1",  logical([0 1 0 0 0]'))
        Qrcbm.add("QRCID_23", logical([0 0 0 1 1]'))

        Qrcs = bicas.proc.QrcSettingL2(currentFvIantAr=[1]');
        Qrcsm.add("QRCID_1",  Qrcs);
        Qrcs = bicas.proc.QrcSettingL2(currentFvIantAr=[2 3]');
        Qrcsm.add("QRCID_23", Qrcs);

        N = NaN;
        test( ...
          Qrcbm, Qrcsm, ...
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

      test_empty_empty()
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
