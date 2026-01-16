%
% matlab.unittest automatic test code for bicas.proc.L1L2.sat.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef sat___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Split up into multiple files.
  %   CON: init_object() is used by tests for multiple methods.
  %   CON: Might want to share instance field S = bicas.proc.L1L2.const.C.SSID_DICT;



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_get_VSTB(T)
      SatSettings.upperThresholdAvoltDcSingle = 2;
      SatSettings.upperThresholdAvoltDcDiff   = 3;
      SatSettings.upperThresholdAvoltAclg     = 12;
      SatSettings.upperThresholdAvoltAchg     = 13;

      function test(samplesAvolt, ssidStrAr, isAchgFpa, expVstbAr)
        actVstbAr = bicas.proc.L1L2.sat.get_VSTB(...
          SatSettings, samplesAvolt, ...
          bicas.proc.L1L2.const.C.SSID_DICT(ssidStrAr), isAchgFpa);

        T.assertEqual(actVstbAr, logical(expVstbAr))
      end

      % 0x0
      test([], [], bicas.utils.FPArray(logical([])), [])

      % Sample values inside and outside of all thresholds + NaN.
      SAMPLES_AVOLT_AR = [2, 3; 12, 13; 99, 99; NaN, NaN];
      FPA              = bicas.utils.FPArray(logical([0 1; 0 1; 0 1; 0 1]));
      SSID_STR_AR      = ["DC_V1" "DC_V12"; "AC_V12", "AC_V23"; "GND", "REF25V"; "DC_V1" "AC_V12"];

      test(-1.01 * SAMPLES_AVOLT_AR, SSID_STR_AR, FPA, [1, 1; 1, 1; 0, 0; 0, 0])
      test(-0.99 * SAMPLES_AVOLT_AR, SSID_STR_AR, FPA, [0, 0; 0, 0; 0, 0; 0, 0])
      test( 0.99 * SAMPLES_AVOLT_AR, SSID_STR_AR, FPA, [0, 0; 0, 0; 0, 0; 0, 0])
      test( 1.01 * SAMPLES_AVOLT_AR, SSID_STR_AR, FPA, [1, 1; 1, 1; 0, 0; 0, 0])
    end



    function test_get_upper_thresholds(T)
      SatSettings.upperThresholdAvoltDcSingle = 2;
      SatSettings.upperThresholdAvoltDcDiff   = 3;
      SatSettings.upperThresholdAvoltAclg     = 12;
      SatSettings.upperThresholdAvoltAchg     = 13;

      function test(ssidStrAr, isAchgFpa, expUpperThresholdAvoltAr)
        actUpperThresholdAvoltAr = bicas.proc.L1L2.sat.get_upper_thresholds(...
          SatSettings, bicas.proc.L1L2.const.C.SSID_DICT(ssidStrAr), isAchgFpa);

        T.assertEqual(actUpperThresholdAvoltAr, expUpperThresholdAvoltAr)
      end

      % Empty arrays
      test([], bicas.utils.FPArray(logical([])), [])

      % 2D arrays, all kinds of thresholds incl. undefined.
      test( ...
        ["DC_V3", "DC_V13"; ...
        "AC_V12", "AC_V23"; ...
        "REF25V", "UNKNOWN"; ...
        "GND",    "GND"], ...
        bicas.utils.FPArray(logical([0 1; 0 1; 0 1; 0 1])), ...
        [2 3; 12 13; NaN NaN; NaN NaN])

      % Iterate of over all possible values for isAchgFpa. Should not affect the
      % specified SSIDs (all except AC diffs).
      IS_ACHG_FPA_CA = {...
        bicas.utils.FPArray(false(3,2), 'NO_FILL_POSITIONS'), ...
        bicas.utils.FPArray(true( 3,2), 'NO_FILL_POSITIONS'), ...
        bicas.utils.FPArray(true( 3,2), 'ONLY_FILL_POSITIONS') ...
        };
      for i = 1:numel(IS_ACHG_FPA_CA)
        test( ...
          ["DC_V3", "DC_V13"; ...
          "REF25V", "UNKNOWN"; ...
          "GND",    "GND"], ...
          IS_ACHG_FPA_CA{i}, ...
          [2 3; NaN NaN; NaN NaN])
      end
    end



  end    % methods(Test)



end
