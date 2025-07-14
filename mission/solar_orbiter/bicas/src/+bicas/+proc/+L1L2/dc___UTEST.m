%
% matlab.unittest automatic test code for bicas.proc.L1L2.dc().
%
% NOTE: Most functions are not tested due to being hard to test.
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



    function test_get_VSIB_5xBLTS_NEW___CWF(testCase)
      L = bicas.Logger('NO_STDOUT', false);

      SatSettings = struct();
      %SatSettings.cwfSlidingWindowLengthSec   = NaN;   % Not used by test.
      %SatSettings.vstbFractionThreshold       = NaN;   % Not used by test.
      SatSettings.upperThresholdAVoltDcSingle = 2;
      SatSettings.upperThresholdAVoltDcDiff   = 5;
      SatSettings.upperThresholdAVoltAclg     = 8;
      SatSettings.upperThresholdAVoltAchg     = 1;

      % (5 BLTS) x (3 records)
      bltsSamplesAVoltAr = single([...
        1 3 5;
        1 3 6;
        6 3 1;
        7 9 7;
        9 7 9]);
      bltsSamplesAVoltAr = permute(bltsSamplesAVoltAr, [2 3 1]);

      bltsSsidAr         = permute(bicas.proc.L1L2.const.C.SSID_DICT(...
        ["DC_V1" "DC_V12" "DC_V23" "AC_V12" "AC_V23"]), [3 2 1]);
      bltsSsidAr         = repmat(bltsSsidAr, [3 1 1]);

      isAchgFpa          = repmat(bicas.utils.FPArray(false, 'NO_FILL_POSITIONS'), [3, 1]);
      uspr               = repmat(1, [3, 1]);

      expBltsVsibAr      = permute(logical([ ...
        0 1 1;
        0 0 1;
        1 0 0;
        0 1 0;
        1 0 1;
        ]), [2 1]);

      % CALL TESTED CODE
      actBltsVsibAr = bicas.proc.L1L2.dc.get_VSIB_5xBLTS_NEW(...
        bltsSamplesAVoltAr, false, uspr, ...
        bltsSsidAr, isAchgFpa, SatSettings, L);

      testCase.assertEqual(actBltsVsibAr, expBltsVsibAr)
    end



    function test_get_VSIB_5xBLTS_NEW___SWF(testCase)
      % PROPOSAL: USPR test data that varies over records.
      %   CON: Must have multiple records. ==> More test data.
      L = bicas.Logger('NO_STDOUT', false);

      SatSettings = struct();
      %SatSettings.cwfSlidingWindowLengthSec   = NaN;
      SatSettings.vstbFractionThreshold       = 0.49;
      SatSettings.upperThresholdAVoltDcSingle = NaN;   % Not used by test.
      SatSettings.upperThresholdAVoltDcDiff   = 2;
      SatSettings.upperThresholdAVoltAclg     = NaN;   % Not used by test.
      SatSettings.upperThresholdAVoltAchg     = NaN;   % Not used by test.



      % (5 BLTS) x (1 record) x (4 spr)
      bltsSamplesAVoltAr = single([...
        1 1 1 1 0 0;
        1 1 3 1 0 0;
        1 3 3 1 0 0;
        1 3 3 3 0 0;
        3 3 3 3 0 0]);
      bltsSamplesAVoltAr = permute(bltsSamplesAVoltAr, [3 2 1]);
      bltsSsidAr         = repmat(bicas.proc.L1L2.const.C.SSID_DICT("DC_V12"),     [1, 5]);
      isAchgFpa          = repmat(bicas.utils.FPArray(false, 'NO_FILL_POSITIONS'), [1, 1]);   % Value should be irrelevant
      expBltsVsibAr      = permute(logical([0 0 1 1 1]), [3, 2, 1]);
      uspr               = repmat(4, [1, 1]);

      % CALL TESTED CODE
      actBltsVsibAr = bicas.proc.L1L2.dc.get_VSIB_5xBLTS_NEW(...
        bltsSamplesAVoltAr, true, uspr, ...
        bltsSsidAr, isAchgFpa, SatSettings, L);

      testCase.assertEqual(actBltsVsibAr, expBltsVsibAr)
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
