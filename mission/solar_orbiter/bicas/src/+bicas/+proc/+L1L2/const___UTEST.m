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



    % Test multiple functions for converting between ASID, SSID, SDID.
    function test_ASID_to_SSID___SSID_ASR_to_ASID___SDID_ASR_to_ASID(T)
      C = bicas.proc.L1L2.const.C;

      for s = ["DC_V13", "AC_V23", "DC_V2"]
        % T.assertEqual() does check for MATLAB class.

        asid = C.ASID_DICT(s);
        ssid = C.SSID_DICT(s);
        sdid = C.SDID_DICT(s);

        T.assertEqual(bicas.proc.L1L2.const.ASID_to_SSID(    asid), ssid)
        T.assertEqual(bicas.proc.L1L2.const.SSID_ASR_to_ASID(ssid), asid)
        T.assertEqual(bicas.proc.L1L2.const.SDID_ASR_to_ASID(sdid), asid)
      end

      ssidUnknown = C.SSID_DICT("UNKNOWN");
      sdidNowhere = C.SDID_DICT("NOWHERE");
      T.assertError(...
        @() bicas.proc.L1L2.const.SSID_ASR_to_ASID(ssidUnknown), ...
        ?MException)
      T.assertError(...
        @() bicas.proc.L1L2.const.SDID_ASR_to_ASID(sdidNowhere), ...
        ?MException)
    end



    %======
    % ASID
    %======

    function test_is_ASID(T)
      C = bicas.proc.L1L2.const.C;

      T.assertTrue(bicas.proc.L1L2.const.is_ASID(C.ASID_DICT("DC_V1")))
      T.assertTrue(bicas.proc.L1L2.const.is_ASID(C.ASID_DICT("AC_V23")))

      % Test non-scalar
      T.assertTrue(bicas.proc.L1L2.const.is_ASID(C.ASID_DICT(...
        ["AC_V23", "DC_V12"; "DC_V1", "DC_V13"])))

      % Test non-ASID
      T.assertFalse(bicas.proc.L1L2.const.is_ASID(C.SSID_DICT("DC_V1")))
      T.assertFalse(bicas.proc.L1L2.const.is_ASID(C.SDID_DICT("AC_V23")))
    end



    function test_get_ASID_category(T)
      T.assertEqual(...
        bicas.proc.L1L2.const.get_ASID_category(bicas.proc.L1L2.const.C.ASID_DICT("DC_V3")), ...
        "DC_SINGLE")
      T.assertEqual(...
        bicas.proc.L1L2.const.get_ASID_category(bicas.proc.L1L2.const.C.ASID_DICT("AC_V13")), ...
        "AC_DIFF")
    end



    function test_ASID_is_AC(T)
      function test(asidStrAr, expB)
        T.test_array_to_logical( ...
          @(asidAr) bicas.proc.L1L2.const.ASID_is_AC(asidAr), ...
          bicas.proc.L1L2.const.C.ASID_DICT, asidStrAr, expB)
      end

      test("DC_V3", 0);
      test("DC_V13", 0);
      test("AC_V13", 1);
      test(...
        ["AC_V13", "DC_V13"; "DC_V1", "AC_V12"], ...
        [1 0; 0 1]);
    end



    function test_ASID_is_diff(T)
      function test(asidStrAr, expB)
        T.test_array_to_logical( ...
          @(asidAr) bicas.proc.L1L2.const.ASID_is_diff(asidAr), ...
          bicas.proc.L1L2.const.C.ASID_DICT, asidStrAr, expB)
      end

      test("DC_V3", 0);
      test("DC_V13", 1);
      test("AC_V13", 1);
      test(...
        ["AC_V13", "DC_V13"; "DC_V1", "AC_V12"], ...
        [1 1; 0 1]);
    end



    %======
    % SSID
    %======

    function test_is_SSID(T)
      C = bicas.proc.L1L2.const.C;

      T.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT("DC_V1")))
      T.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT("AC_V23")))
      T.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT("REF25V")))
      T.assertTrue(bicas.proc.L1L2.const.is_SSID(C.SSID_DICT(...
        ["AC_V23", "DC_V12"; "DC_V1", "DC_V13"; "UNKNOWN", "GND"])))

      % Not SSIDs.
      T.assertFalse(bicas.proc.L1L2.const.is_SSID(C.ASID_DICT("DC_V1")))
      T.assertFalse(bicas.proc.L1L2.const.is_SSID(C.SDID_DICT("AC_V23")))
    end



    function test_SSID_is_ASR(T)
      function test(ssidStrAr, expB)
        T.test_array_to_logical( ...
          @(ssidAr) bicas.proc.L1L2.const.SSID_is_ASR(ssidAr), ...
          bicas.proc.L1L2.const.C.SSID_DICT, ssidStrAr, expB)
      end

      test("DC_V13", true);
      test(...
        ["GND", "DC_V13"; "REF25V", "DC_V2"; "UNKNOWN", "AC_V23"], ...
        [0 1; 0 1; 0 1]);
    end



    function test_SSID_is_AC(T)
      function test(ssidStrAr, expB)
        T.test_array_to_logical( ...
          @(ssidAr) bicas.proc.L1L2.const.SSID_is_AC(ssidAr), ...
          bicas.proc.L1L2.const.C.SSID_DICT, ssidStrAr, expB)
      end

      test("DC_V3", 0);
      test("DC_V13", 0);
      test("AC_V13", 1);
      test(...
        ["AC_V13", "DC_V13"; "DC_V1", "AC_V12"; "REF25V", "UNKNOWN"], ...
        [1 0; 0 1; 0 0]);
    end



    function test_SSID_is_diff(T)
      function test(ssidStrAr, expB)
        T.test_array_to_logical( ...
          @(ssidAr) bicas.proc.L1L2.const.SSID_is_diff(ssidAr), ...
          bicas.proc.L1L2.const.C.SSID_DICT, ssidStrAr, expB)
      end

      test(["DC_V3"], false)
      test(["DC_V13"], true)
      test(["AC_V13"], true)

      test(...
        ["AC_V13", "DC_V12"; "DC_V1", "DC_V3"; "REF25V", "UNKNOWN"], ...
        [1 1; 0 0; 0 0])
    end



    %======
    % SDID
    %======

    function test_is_SDID(T)
      C = bicas.proc.L1L2.const.C;

      T.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT("DC_V1")))
      T.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT("AC_V23")))
      T.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT("NOWHERE")))

      T.assertTrue(bicas.proc.L1L2.const.is_SDID(C.SDID_DICT(...
        ["AC_V23", "DC_V12"; "DC_V1", "DC_V13"])))

      T.assertFalse(bicas.proc.L1L2.const.is_SDID(C.ASID_DICT("DC_V1")))
      T.assertFalse(bicas.proc.L1L2.const.is_SDID(C.SSID_DICT("AC_V23")))
    end



    function test_SDID_is_ASR(T)
      C = bicas.proc.L1L2.const.C;

      T.assertTrue( bicas.proc.L1L2.const.SDID_is_ASR(C.SDID_DICT("DC_V13")));
      T.assertFalse(bicas.proc.L1L2.const.SDID_is_ASR(C.SDID_DICT("NOWHERE")));
    end



    function test_SDID_is_nowhere(T)
      C = bicas.proc.L1L2.const.C;

      T.assertFalse(bicas.proc.L1L2.const.SDID_is_nowhere(C.SDID_DICT("DC_V13")));
      T.assertTrue( bicas.proc.L1L2.const.SDID_is_nowhere(C.SDID_DICT("NOWHERE")));
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################



  methods(Access=private)



    % Utility function for simplifying some tests.
    % Only somewhat useful.
    function test_array_to_logical(T, fh, ArgDict, arg, expRv)
      arg = ArgDict(arg);
      T.assertEqual(fh(arg), logical(expRv));
    end



  end    % methods(Access=private)



end
