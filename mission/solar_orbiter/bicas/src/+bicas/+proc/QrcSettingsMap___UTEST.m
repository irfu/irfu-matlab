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



    function test_qrcidAr(T)
      Qrcsm = bicas.proc.QrcSettingsMap();
      T.assertEqual(Qrcsm.qrcidAr, string.empty(0, 1))

      % IMPLEMENTATION NOTE: Deliberately add QRCIDs in non-alphanumeric order.
      Qrcsm.add("QRCID_2", T.QRCS_2);
      T.assertEqual(Qrcsm.qrcidAr, "QRCID_2")

      % Deliberately assert alphanumeric order for return value.
      Qrcsm.add("QRCID_1", T.QRCS_1);
      T.assertEqual(Qrcsm.qrcidAr, ["QRCID_1"; "QRCID_2"])
    end



    function test_add_get(T)
      Qrcsm = bicas.proc.QrcSettingsMap();

      Qrcsm.add("QRCID_1", T.QRCS_1);
      Qrcsm.add("QRCID_2", T.QRCS_2);

      ActQrcs1 = Qrcsm.get("QRCID_1");
      T.assertEqual(ActQrcs1, T.QRCS_1)

      ActQrcs2 = Qrcsm.get("QRCID_2");
      T.assertEqual(ActQrcs2, T.QRCS_2)

      % Implicitly test equality of QRCS class since these tests rely on that.
      T.assertNotEqual(T.QRCS_1, T.QRCS_2)
    end



    function test_add_QRCSM___empty(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap();
      Qrcsm2 = bicas.proc.QrcSettingsMap();

      Qrcsm1.add_QRCSM(Qrcsm2)
      T.assertEqual(Qrcsm1, Qrcsm2)
    end



    function test_add_QRCSM___nonempty(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap();
      Qrcsm1.add("QRCID_1", T.QRCS_1);
      Qrcsm1.add("QRCID_2", T.QRCS_2);

      Qrcsm2 = bicas.proc.QrcSettingsMap();
      Qrcsm2.add("QRCID_3", T.QRCS_3);

      Qrcsm1.add_QRCSM(Qrcsm2)

      ExpQrcsm = bicas.proc.QrcSettingsMap();
      ExpQrcsm.add("QRCID_1", T.QRCS_1);
      ExpQrcsm.add("QRCID_2", T.QRCS_2);
      ExpQrcsm.add("QRCID_3", T.QRCS_3);

      T.assertEqual(Qrcsm1, ExpQrcsm)
    end



    function test_add_QRCSM___collision(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap();
      Qrcsm1.add("QRCID_1", T.QRCS_1);
      Qrcsm1.add("QRCID_2", T.QRCS_2);    % Key collision

      Qrcsm2 = bicas.proc.QrcSettingsMap();
      Qrcsm2.add("QRCID_2", T.QRCS_3);    % key collision. Other QRCS.
      Qrcsm2.add("QRCID_3", T.QRCS_3);

      T.assertError(...
          @() Qrcsm1.add_QRCSM(Qrcsm2), ...
          ?MException)
    end



    function test_equality(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap();
      Qrcsm2 = bicas.proc.QrcSettingsMap();

      T.assertEqual(Qrcsm1, Qrcsm2)

      Qrcsm1.add("QRCID_1", T.QRCS_1);      T.assertNotEqual(Qrcsm1, Qrcsm2)
      Qrcsm1.add("QRCID_2", T.QRCS_2);      T.assertNotEqual(Qrcsm1, Qrcsm2)
      Qrcsm1.add("QRCID_3", T.QRCS_3);

      T.assertNotEqual(Qrcsm1, Qrcsm2)

      % Add same QRCSs, but in reverse order.
      Qrcsm2.add("QRCID_3", T.QRCS_3);      T.assertNotEqual(Qrcsm1, Qrcsm2)
      Qrcsm2.add("QRCID_2", T.QRCS_2);      T.assertNotEqual(Qrcsm1, Qrcsm2)
      Qrcsm2.add("QRCID_1", T.QRCS_1);

      T.assertEqual(Qrcsm1, Qrcsm2)
    end



    function test_remove_many(T)
      Qrcsm = bicas.proc.QrcSettingsMap();

      Qrcsm.add("QRCID_1", T.QRCS_1);
      Qrcsm.add("QRCID_2", T.QRCS_2);
      ExpQrcsm = copy(Qrcsm);
      Qrcsm.remove_many(string.empty(0, 1))
      T.assertEqual(Qrcsm, ExpQrcsm);

      Qrcsm.remove_many(["QRCID_2"; "QRCID_1"])
      ExpQrcsm = bicas.proc.QrcSettingsMap();
      T.assertEqual(Qrcsm, ExpQrcsm);
    end



    function test_copy(T)
      Qrcsm1 = bicas.proc.QrcSettingsMap();
      Qrcsm1.add("QRCID_1", T.QRCS_1);

      Qrcsm2 = copy(Qrcsm1);
      T.assertEqual(Qrcsm1, Qrcsm2);

      Qrcsm2.add("QRCID_2", T.QRCS_2);
      T.assertNotEqual(Qrcsm1, Qrcsm2);

      Qrcsm1.add("QRCID_2", T.QRCS_2);
      T.assertEqual(Qrcsm1, Qrcsm2);
    end



  end    % methods(Test)



end
