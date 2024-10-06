%
% matlab.unittest automatic test code for bicas.proc.L1L2.const.
%
% NOTE: Does not cover all functions, but probably the most important cases.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef const___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_is_ASID(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertTrue(bicas.proc.L1L2.const.is_ASID(C.ASID_DICT("DC_V1")))
      testCase.assertTrue(bicas.proc.L1L2.const.is_ASID(C.ASID_DICT("AC_V23")))

      testCase.assertTrue(bicas.proc.L1L2.const.is_ASID(C.ASID_DICT(...
        ["AC_V23", "DC_V12"; "DC_V1", "DC_V13"])))

      testCase.assertFalse(bicas.proc.L1L2.const.is_ASID(C.SSID_DICT("DC_V1")))
      testCase.assertFalse(bicas.proc.L1L2.const.is_ASID(C.SDID_DICT("AC_V23")))
    end



    function test_get_ASID_category(testCase)
      testCase.assertEqual(...
        bicas.proc.L1L2.const.get_ASID_category(bicas.proc.L1L2.const.C.ASID_DICT("DC_V3")), ...
        "DC_SINGLE")
      testCase.assertEqual(...
        bicas.proc.L1L2.const.get_ASID_category(bicas.proc.L1L2.const.C.ASID_DICT("AC_V13")), ...
        "AC_DIFF")
    end



    function test_ASID_to_SSID___SSID_ASR_to_ASID___SDID_ASR_to_ASID(testCase)
      C = bicas.proc.L1L2.const.C;

      for s = ["DC_V13", "AC_V23", "DC_V2"]
        % testCase.assertEqual() does check for MATLAB class.

        asid = C.ASID_DICT(s);
        ssid = C.SSID_DICT(s);
        sdid = C.SDID_DICT(s);

        testCase.assertEqual(bicas.proc.L1L2.const.ASID_to_SSID(    asid), ssid)
        testCase.assertEqual(bicas.proc.L1L2.const.SSID_ASR_to_ASID(ssid), asid)
        testCase.assertEqual(bicas.proc.L1L2.const.SDID_ASR_to_ASID(sdid), asid)
      end

      ssidUnknown = C.SSID_DICT("UNKNOWN");
      sdidUnknown = C.SDID_DICT("NOWHERE");
      testCase.assertError(...
          @() bicas.proc.L1L2.const.SSID_ASR_to_ASID(ssidUnknown), ...
          ?MException)
      testCase.assertError(...
          @() bicas.proc.L1L2.const.SDID_ASR_to_ASID(sdidUnknown), ...
          ?MException)
    end



    function test_is_SSID(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT("DC_V1")))
      testCase.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT("AC_V23")))
      testCase.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT("REF25V")))

      testCase.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT(...
        ["AC_V23", "DC_V12"; "DC_V1", "DC_V13"])))

      testCase.assertFalse(bicas.proc.L1L2.const.is_SSID(C.ASID_DICT("DC_V1")))
      testCase.assertFalse(bicas.proc.L1L2.const.is_SSID(C.SDID_DICT("AC_V23")))
    end



    function test_SSID_is_ASR(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertTrue( bicas.proc.L1L2.const.SSID_is_ASR(C.SSID_DICT("DC_V13")));
      testCase.assertFalse(bicas.proc.L1L2.const.SSID_is_ASR(C.SSID_DICT("GND")));
    end



    function test_SSID_is_AC(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertFalse(bicas.proc.L1L2.const.SSID_is_AC(C.SSID_DICT("DC_V3")));
      testCase.assertFalse(bicas.proc.L1L2.const.SSID_is_AC(C.SSID_DICT("DC_V13")));
      testCase.assertTrue( bicas.proc.L1L2.const.SSID_is_AC(C.SSID_DICT("AC_V13")));
    end



    function test_SSID_is_diff(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertFalse(bicas.proc.L1L2.const.SSID_is_diff(C.SSID_DICT("DC_V3")));
      testCase.assertTrue( bicas.proc.L1L2.const.SSID_is_diff(C.SSID_DICT("DC_V13")));
      testCase.assertTrue( bicas.proc.L1L2.const.SSID_is_diff(C.SSID_DICT("AC_V13")));
    end



    function test_is_SDID(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT("DC_V1")))
      testCase.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT("AC_V23")))
      testCase.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT("NOWHERE")))

      testCase.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT(...
        ["AC_V23", "DC_V12"; "DC_V1", "DC_V13"])))

      testCase.assertFalse(bicas.proc.L1L2.const.is_SDID(C.ASID_DICT("DC_V1")))
      testCase.assertFalse(bicas.proc.L1L2.const.is_SDID(C.SSID_DICT("AC_V23")))
    end



    function test_SDID_is_ASR(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertTrue( bicas.proc.L1L2.const.SDID_is_ASR(C.SDID_DICT("DC_V13")));
      testCase.assertFalse(bicas.proc.L1L2.const.SDID_is_ASR(C.SDID_DICT("NOWHERE")));
    end



    function test_SDID_is_nowhere(testCase)
      C = bicas.proc.L1L2.const.C;

      testCase.assertFalse(bicas.proc.L1L2.const.SDID_is_nowhere(C.SDID_DICT("DC_V13")));
      testCase.assertTrue( bicas.proc.L1L2.const.SDID_is_nowhere(C.SDID_DICT("NOWHERE")));
    end



  end    % methods(Test)



end
