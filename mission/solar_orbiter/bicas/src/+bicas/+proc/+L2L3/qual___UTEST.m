%
% matlab.unittest automatic test code for bicas.proc.L2L3.qual.
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



    function test_L2QBM_to_QRCBs(testCase)

      function test(l2QbmAr, Qrcsm, ExpQrcbMap)
        ActQrcbMap = bicas.proc.L2L3.qual.L2QBM_to_QRCBs(l2QbmAr, Qrcsm);

        testCase.assertEqual(ActQrcbMap, ExpQrcbMap);
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
      % QrcsMap = containers.Map();

      Qrcs = bicas.proc.QrcSettingL2(L2_QUALITY_BITMASK=uint16(1));
      Qrcsm.add("QRCID_1", Qrcs);
      % QrcsMap("QRCID_1") = Qrcs;
      Qrcs = bicas.proc.QrcSettingL2(L2_QUALITY_BITMASK=uint16(4));
      Qrcsm.add("QRCID_2", Qrcs);
      % QrcsMap("QRCID_2") = Qrcs;

      ExpQrcbMap = bicas.proc.QrcbMap(8);
      ExpQrcbMap.add("QRCID_1", logical([0 1 0 1 0 1 0 1]'))
      ExpQrcbMap.add("QRCID_2", logical([0 0 0 0 1 1 1 1]'))

      test( ...
        uint16(0:7)', ...
        Qrcsm, ...
        ExpQrcbMap)
    end



    function test_L2QBM_to_saturation_QRCBs(T)

      function test(l2qbmAr, saturationQualitySchemeId, ExpQrcbMap)
        ActQrcbMap = bicas.proc.L2L3.qual.L2QBM_to_saturation_QRCBs(...
          l2qbmAr, saturationQualitySchemeId);

        % for qrcid = ActQrcbMap.qrcidAr'
        %   qrcid
        %   T.assertEqual(ActQrcbMap.get(qrcid), ExpQrcbMap.get(qrcid))
        % end
        T.assertEqual(ActQrcbMap, ExpQrcbMap)
      end

      function test_GLOBAL_SATURATION()
        l2qbmAr = uint16(2.^[0:15]');
        ExcQrcbMap = bicas.proc.QrcbMap(numel(l2qbmAr));
        ExcQrcbMap.add_false(bicas.const.Q.SATURATION_QRCID_AR)
        test(l2qbmAr, "GLOBAL_SATURATION", ExcQrcbMap)
      end

      function test_CHANNEL_SATURATION()
        % NOTE: Tests assumes locations of quality bits. Could be fixed, but
        % then it still becomes hard to test non-saturation bits.
        l2qbmAr = uint16(2.^[4:6]');
        ExcQrcbMap = bicas.proc.QrcbMap(numel(l2qbmAr));
        ExcQrcbMap.add_false(bicas.const.Q.SATURATION_QRCID_AR)
        ExcQrcbMap.set("SATURATION_ZV_V13", logical([1 0 0]'))
        ExcQrcbMap.set("SATURATION_ZV_V23", logical([0 1 0]'))
        test(l2qbmAr, "CHANNEL_SATURATION", ExcQrcbMap)
      end

      test_GLOBAL_SATURATION()
      test_CHANNEL_SATURATION()
    end



    function test_set_VDC_EDC_samples_FV(T)

      function test(VDC, EDC, QrcbMap, QrcsMap, Exp_VDC, Exp_EDC)
        FV = single(-1);
        VDC_Fpa     = bicas.utils.FPArray(single(VDC),     'FILL_VALUE', FV);
        EDC_Fpa     = bicas.utils.FPArray(single(EDC),     'FILL_VALUE', FV);
        Exp_VDC_Fpa = bicas.utils.FPArray(single(Exp_VDC), 'FILL_VALUE', FV);
        Exp_EDC_Fpa = bicas.utils.FPArray(single(Exp_EDC), 'FILL_VALUE', FV);

        [Act_VDC_Fpa, Act_EDC_Fpa] = bicas.proc.L2L3.qual.set_VDC_EDC_samples_FV(...
          VDC_Fpa, EDC_Fpa, QrcbMap, QrcsMap);

        T.assertEqual(Act_VDC_Fpa, Exp_VDC_Fpa)
        T.assertEqual(Act_EDC_Fpa, Exp_EDC_Fpa)
      end

      function test_zero_rows()
        QrcbMap = bicas.proc.QrcbMap(0);
        QrcbMap.add("QRCID_1", logical.empty(0, 1))
        QrcbMap.add("QRCID_2", logical.empty(0, 1))
        VDC     = zeros(0, 3);
        Exp_VDC = zeros(0, 3);
        EDC     = zeros(0, 3);
        Exp_EDC = zeros(0, 3);

        Qrcsm = bicas.proc.QrcSettingsMap();
        Qrcsm.add("QRCID_1", bicas.proc.QrcSettingL3(vdcFvIndexAr=[1 2]'))
        Qrcsm.add("QRCID_2", bicas.proc.QrcSettingL3(edcFvIndexAr=[2 3]'))

        test(VDC, EDC, QrcbMap, Qrcsm, Exp_VDC, Exp_EDC)
      end

      function test_zero_QRCBs_QRCSs()
        DATA = [
          -1  2  3    -1 3 4; ...
           2 -1  4     3 4 5; ...
           3  4  5     4 5 6; ...
           4  5  6     5 6 7; ...
        ];
        QrcbMap = bicas.proc.QrcbMap(size(DATA, 1));
        VDC     = DATA(:,  1:3);
        Exp_VDC = VDC;
        EDC     = DATA(:,  4:6);
        Exp_EDC = EDC;

        Qrcsm = bicas.proc.QrcSettingsMap();

        test(VDC, EDC, QrcbMap, Qrcsm, Exp_VDC, Exp_EDC)
      end

      function test_complex()
        DATA = [
          0 1   -1  2  3   -1  2  3   -1 3 4  -1 -1 -1; ...
          1 1    2 -1  4   -1 -1  4    3 4 5   3 -1 -1; ...
          1 1    3  4  5   -1 -1  5    4 5 6   4 -1 -1; ...
          1 0    4  5  6   -1 -1  6    5 6 7   5  6  7; ...
        ];
        QrcbMap = bicas.proc.QrcbMap(size(DATA, 1));
        QrcbMap.add("QRCID_1", logical(DATA(:, 1)))
        QrcbMap.add("QRCID_2", logical(DATA(:, 2)))
        VDC     = DATA(:,  3:5);
        Exp_VDC = DATA(:,  6:8);
        EDC     = DATA(:,  9:11);
        Exp_EDC = DATA(:, 12:14);

        Qrcsm = bicas.proc.QrcSettingsMap();
        Qrcsm.add("QRCID_1", bicas.proc.QrcSettingL3(vdcFvIndexAr=[1 2]'))
        Qrcsm.add("QRCID_2", bicas.proc.QrcSettingL3(edcFvIndexAr=[2 3]'))

        test(VDC, EDC, QrcbMap, Qrcsm, Exp_VDC, Exp_EDC)
      end

      test_zero_rows()
      test_zero_QRCBs_QRCSs()
      test_complex()
    end



  end    % methods(Test)



end
