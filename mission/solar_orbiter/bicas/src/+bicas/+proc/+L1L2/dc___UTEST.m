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



    function test_get_SSID_SDID_arrays(testCase)
      S = bicas.proc.L1L2.const.C.SSID_DICT;
      D = bicas.proc.L1L2.const.C.SDID_DICT;

      %=============================================================
      % Outputs for corresponding inputs (one CDF record at a time)
      %=============================================================
      SSID_BDM0_DLR0_ROW = S(["DC_V1", "DC_V12", "DC_V23", "AC_V12", "AC_V23"]);
      SDID_BDM0_DLR0_ROW = D(["DC_V1", "DC_V12", "DC_V23", "AC_V12", "AC_V23"]);

      SSID_BDM0_DLR1_ROW = S(["DC_V1", "DC_V13", "DC_V23", "AC_V13", "AC_V23"]);
      SDID_BDM0_DLR1_ROW = D(["DC_V1", "DC_V13", "DC_V23", "AC_V13", "AC_V23"]);

      SSID_BDM4_DLR1_ROW = S(["DC_V1", "DC_V2",  "DC_V3",  "AC_V13", "AC_V23"]);
      SDID_BDM4_DLR1_ROW = D(["DC_V1", "DC_V2",  "DC_V3",  "AC_V13", "AC_V23"]);

      % BDMX = Unknown BDM
      SSID_BDMX_DLR1_ROW = S(["UNKNOWN", "UNKNOWN",  "UNKNOWN",  "AC_V13", "AC_V23"]);
      SDID_BDMX_DLR1_ROW = D(["NOWHERE", "NOWHERE",  "NOWHERE",  "AC_V13", "AC_V23"]);

      %=======
      % Tests
      %=======
      testCase.test_get_SSID_SDID_arrays_helper(...
        bicas.utils.FPArray(zeros(0,1,'uint8')), ...
        bicas.utils.FPArray(false(0,1,'logical')), ...
        zeros(0, 5), ...
        zeros(0, 5));

      testCase.test_get_SSID_SDID_arrays_helper(...
        bicas.utils.FPArray(uint8(0)), ...
        bicas.utils.FPArray(false), ...
        SSID_BDM0_DLR0_ROW, ...
        SDID_BDM0_DLR0_ROW);

      testCase.test_get_SSID_SDID_arrays_helper(...
        bicas.utils.FPArray(uint8(4)), ...
        bicas.utils.FPArray(true), ...
        SSID_BDM4_DLR1_ROW, ...
        SDID_BDM4_DLR1_ROW);

      testCase.test_get_SSID_SDID_arrays_helper(...
        bicas.utils.FPArray(uint8(0), 'ONLY_FILL_POSITIONS'), ...
        bicas.utils.FPArray(true), ...
        SSID_BDMX_DLR1_ROW, ...
        SDID_BDMX_DLR1_ROW);

      % "Complex test"
      testCase.test_get_SSID_SDID_arrays_helper(...
        bicas.utils.FPArray(uint8([0; 0; 0; 4; 255]), 'FILL_VALUE', uint8(255)), ...
        bicas.utils.FPArray([false; false; true; true; true]), ...
        [ ...
        SSID_BDM0_DLR0_ROW;
        SSID_BDM0_DLR0_ROW;
        SSID_BDM0_DLR1_ROW;
        SSID_BDM4_DLR1_ROW;
        SSID_BDMX_DLR1_ROW;
        ], ...
        [ ...
        SDID_BDM0_DLR0_ROW;
        SDID_BDM0_DLR0_ROW;
        SDID_BDM0_DLR1_ROW;
        SDID_BDM4_DLR1_ROW;
        SDID_BDMX_DLR1_ROW;
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



    % Helper function for testing bicas.proc.L1L2.dc.get_SSID_SDID_arrays().
    function test_get_SSID_SDID_arrays_helper(...
        testCase, bdmFpa, dlrFpa, expBltsSsidArray, expBltsSdidArray)
      % PROPOSAL: Create FPAs inside this function instead of by caller.

      [actBltsSsidArray, actBltsSdidArray] = ...
        bicas.proc.L1L2.dc.get_SSID_SDID_arrays(bdmFpa, dlrFpa);

      testCase.assertEqual(size(actBltsSsidArray), size(actBltsSdidArray))
      testCase.assertEqual(size(actBltsSsidArray, 2), bicas.const.N_BLTS)
      testCase.assertEqual(class(actBltsSsidArray), 'uint8')
      testCase.assertEqual(class(actBltsSdidArray), 'uint8')

      testCase.assertEqual(actBltsSsidArray, uint8(expBltsSsidArray))
      testCase.assertEqual(actBltsSdidArray, uint8(expBltsSdidArray))
    end



  end    % methods(Access=private)



end
