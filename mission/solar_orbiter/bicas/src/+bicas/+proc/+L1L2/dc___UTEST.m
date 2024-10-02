%
% matlab.unittest automatic test code for bicas.proc.L1L2.dc().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef dc___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      S = bicas.sconst.C.SSID_S_K_DICT;
      D = bicas.sconst.C.SDID_S_K_DICT;

      SSID_BDM0_DLR0_ROW = S(["DC_V1", "DC_V12", "DC_V23", "AC_V12", "AC_V23"]);
      SDID_BDM0_DLR0_ROW = D(["DC_V1", "DC_V12", "DC_V23", "AC_V12", "AC_V23"]);

      SSID_BDM0_DLR1_ROW = S(["DC_V1", "DC_V13", "DC_V23", "AC_V13", "AC_V23"]);
      SDID_BDM0_DLR1_ROW = D(["DC_V1", "DC_V13", "DC_V23", "AC_V13", "AC_V23"]);

      SSID_BDM4_DLR1_ROW = S(["DC_V1", "DC_V2", "DC_V3", "AC_V13", "AC_V23"]);
      SDID_BDM4_DLR1_ROW = D(["DC_V1", "DC_V2", "DC_V3", "AC_V13", "AC_V23"]);

      testCase.test_get_KSSID_KSDID_arrays(...
        bicas.utils.FPArray(zeros(0,1,'uint8')), ...
        bicas.utils.FPArray(false(0,1,'logical')), ...
        zeros(0, 5), ...
        zeros(0, 5));

      testCase.test_get_KSSID_KSDID_arrays(...
        bicas.utils.FPArray(uint8(0)), ...
        bicas.utils.FPArray(false), ...
        SSID_BDM0_DLR0_ROW, ...
        SDID_BDM0_DLR0_ROW);

      testCase.test_get_KSSID_KSDID_arrays(...
        bicas.utils.FPArray(uint8(4)), ...
        bicas.utils.FPArray(true), ...
        SSID_BDM4_DLR1_ROW, ...
        SDID_BDM4_DLR1_ROW);

      % "Complex test"
      testCase.test_get_KSSID_KSDID_arrays(...
        bicas.utils.FPArray(uint8([0; 0; 0; 4])), ...
        bicas.utils.FPArray([false; false; true; true]), ...
        [ ...
        SSID_BDM0_DLR0_ROW;
        SSID_BDM0_DLR0_ROW;
        SSID_BDM0_DLR1_ROW;
        SSID_BDM4_DLR1_ROW;
        ], ...
        [ ...
        SDID_BDM0_DLR0_ROW;
        SDID_BDM0_DLR0_ROW;
        SDID_BDM0_DLR1_ROW;
        SDID_BDM4_DLR1_ROW;
        ] ...
        );
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function test_get_KSSID_KSDID_arrays(...
        testCase, bdmFpa, dlrFpa, expBltsKSsidArray, expBltsKSdidArray)

      [actBltsKSsidArray, actBltsKSdidArray] = ...
        bicas.proc.L1L2.dc.get_KSSID_KSDID_arrays(bdmFpa, dlrFpa);

      testCase.assertEqual(size(actBltsKSsidArray), size(actBltsKSdidArray))
      testCase.assertEqual(class(actBltsKSsidArray), 'uint8')
      testCase.assertEqual(class(actBltsKSdidArray), 'uint8')

      testCase.assertEqual(actBltsKSsidArray, uint8(expBltsKSsidArray))
      testCase.assertEqual(actBltsKSdidArray, uint8(expBltsKSdidArray))
    end



  end    % methods(Access=private)



end