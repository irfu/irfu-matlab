%
% matlab.unittest automatic test code for bicas.proc.QrcSettingsMap.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef QrcSettingsMap___UTEST < matlab.unittest.TestCase



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    QRCS_1 = bicas.proc.QrcSettingL2(QUALITY_FLAG=uint8(1));
    QRCS_2 = bicas.proc.QrcSettingL2(QUALITY_FLAG=uint8(2));
    QRCS_3 = bicas.proc.QrcSettingL2(QUALITY_FLAG=uint8(3));
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_add_get(T)
      Qrcsm = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"]);

      Qrcsm.add("QRCID_1", "PTID_1", T.QRCS_1);
      Qrcsm.add("QRCID_1", "PTID_2", T.QRCS_2);
      Qrcsm.add("QRCID_2", "PTID_1", T.QRCS_3);
      Qrcsm.add("QRCID_2", "PTID_2", T.QRCS_2);

      ActQrcs1 = Qrcsm.get("QRCID_1", "PTID_1");
      T.assertEqual(ActQrcs1, T.QRCS_1)
      ActQrcs2 = Qrcsm.get("QRCID_1", "PTID_2");
      T.assertEqual(ActQrcs2, T.QRCS_2)
      ActQrcs3 = Qrcsm.get("QRCID_2", "PTID_1");
      T.assertEqual(ActQrcs3, T.QRCS_3)

      % Implicitly test equality of QRCS class since these tests rely on that.
      T.assertNotEqual(T.QRCS_1, T.QRCS_2)
    end



    function test_get_QRCIDs(T)
      Qrcsm = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"]);

      T.assertEqual(Qrcsm.get_QRCIDs(), string.empty(0, 1))

      Qrcsm.add("QRCID_1", "PTID_1", T.QRCS_1);
      Qrcsm.add("QRCID_1", "PTID_2", T.QRCS_2);
      Qrcsm.add("QRCID_2", "PTID_1", T.QRCS_3);

      T.assertEqual(Qrcsm.get_QRCIDs(), ["QRCID_1"; "QRCID_2"])
    end



    function test_merge___empty(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap(string.empty(0, 1));
      Qrcsm2 = bicas.proc.QrcSettingsMap(string.empty(0, 1));

      Qrcsm1.merge(Qrcsm2)
      T.assertEqual(Qrcsm1, Qrcsm2)
    end



    function test_merge___nonempty(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"]);
      Qrcsm1.add("QRCID_1", "PTID_1", T.QRCS_1);
      Qrcsm1.add("QRCID_1", "PTID_2", T.QRCS_2);

      Qrcsm2 = bicas.proc.QrcSettingsMap(["PTID_2"; "PTID_3"]);
      Qrcsm2.add("QRCID_2", "PTID_3", T.QRCS_3);

      Qrcsm1.merge(Qrcsm2)

      ExpQrcsm = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"; "PTID_3"]);
      ExpQrcsm.add("QRCID_1", "PTID_1", T.QRCS_1);
      ExpQrcsm.add("QRCID_1", "PTID_2", T.QRCS_2);
      ExpQrcsm.add("QRCID_2", "PTID_3", T.QRCS_3);

      T.assertEqual(Qrcsm1, ExpQrcsm)

      % Implicitly check that can using the union om PTIDs.
      ExpQrcsm.add("QRCID_3", "PTID_1", T.QRCS_1);
      ExpQrcsm.add("QRCID_3", "PTID_2", T.QRCS_2);
      ExpQrcsm.add("QRCID_3", "PTID_3", T.QRCS_3);
    end



    function test_merge___collision(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"]);
      Qrcsm1.add("QRCID_1", "PTID_1", T.QRCS_1);
      Qrcsm1.add("QRCID_1", "PTID_2", T.QRCS_2);    % Collision

      Qrcsm2 = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"]);
      Qrcsm2.add("QRCID_1", "PTID_2", T.QRCS_2);    % Collision
      Qrcsm2.add("QRCID_2", "PTID_1", T.QRCS_3);

      T.assertError(...
          @() Qrcsm1.merge(Qrcsm2), ...
          ?MException)
    end



    function test_equality(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"]);
      Qrcsm2 = bicas.proc.QrcSettingsMap(["PTID_1"; "PTID_2"]);

      T.assertEqual(Qrcsm1, Qrcsm2)

      Qrcsm1.add("QRCID_1", "PTID_1", T.QRCS_1);
      Qrcsm1.add("QRCID_1", "PTID_2", T.QRCS_2);
      Qrcsm1.add("QRCID_2", "PTID_1", T.QRCS_3);

      T.assertNotEqual(Qrcsm1, Qrcsm2)

      % Add same QRCSs, but in reverse order.
      Qrcsm2.add("QRCID_2", "PTID_1", T.QRCS_3);
      Qrcsm2.add("QRCID_1", "PTID_2", T.QRCS_2);
      Qrcsm2.add("QRCID_1", "PTID_1", T.QRCS_1);

      T.assertEqual(Qrcsm1, Qrcsm2)
    end



  end    % methods(Test)



end
