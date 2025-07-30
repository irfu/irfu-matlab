%
% matlab.unittest automatic test code for bicas.proc.L1L2.dc().
%
% NOTE: Most functions are not tested due to being hard to test.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef dc___UTEST < matlab.unittest.TestCase



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    L = bicas.Logger('HUMAN_READABLE', false);
    S = bicas.proc.L1L2.const.C.SSID_DICT;
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Test CWF NaN data.
    %
    % NOTE: bicas.tf.apply_TF() itself splits by NaN (by default).
    %
    function test_calibrate_voltages_1xBLTS___CWF_NaN(testCase)
      Bso = testCase.get_simple_TF_BSO();

      Vcds = bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierTest( ...
        itfBiasAvpiv=@(omegaRps) (ones(size(omegaRps))*2), ...
        offsetAvolt=1, ...
        itfLfrAvpiv =@(omegaRps) (ones(size(omegaRps))*3));
      Vcal = bicas.proc.L1L2.cal.VoltageCalibration(Vcds, true, Bso);

      SAMPLES_TM = [1 2 3 NaN 5 6 NaN 8 NaN]';   % NOTE: Pre-existing NaN!
      I_LSF      = 1;
      FREQ_HZ    = solo.hwzv.const.LSF_HZ(I_LSF);

      actSamplesAVolt = bicas.proc.L1L2.dc.calibrate_voltages_1xBLTS( ...
        ... % Variables which DO NOT VARY over CDF records at all.
        Vcal         = Vcal, ...
        L            = testCase.L, ...
        iBlts        = 1, ...
        isLfr        = true, ...
        isTdsCwf     = false, ...
        hasSwfFormat = false, ...
        ... % Variables which DO VARY over CDF records.
        tt2000       = uint64( [1:9]' * 1e9 / FREQ_HZ), ...
        ssid         = testCase.S(repmat(["DC_V2"], 9, 1)), ...
        samplesTm    = single(SAMPLES_TM), ...
        uspr         = repmat(1,       9, 1), ...
        freqHz       = repmat(FREQ_HZ, 9, 1), ...
        iLsf         = repmat(I_LSF,   9, 1), ...
        isAchgFpa    = bicas.utils.FPArray(repmat(logical(0), 9, 1), 'NO_FILL_POSITIONS'), ...
        NbriFpa      = bicas.utils.FPArray(repmat(uint8(  1), 9, 1), 'NO_FILL_POSITIONS'), ...
        NbciFpa      = bicas.utils.FPArray(repmat(uint8(  0), 9, 1), 'NO_FILL_POSITIONS'));

      % Emulated calibration
      expSamplesAvolt    = SAMPLES_TM*2*3 + 1;
      % NOTE: Setting one extra NaN due to
      % PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES.
      expSamplesAvolt(8) = NaN;
      testCase.assertEqual(actSamplesAVolt, expSamplesAvolt, AbsTol=1e-14)
    end



    % Test SWF NaN data.
    function test_calibrate_voltages_1xBLTS___SWF_NaN(testCase)
      Bso = testCase.get_simple_TF_BSO();

      Vcds = bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierTest( ...
        itfBiasAvpiv=@(omegaRps) (ones(size(omegaRps))*2), ...
        offsetAvolt=1, ...
        itfLfrAvpiv =@(omegaRps) (ones(size(omegaRps))*3));
      Vcal = bicas.proc.L1L2.cal.VoltageCalibration(Vcds, true, Bso);

      SAMPLES_TM = [...
        1:9;  ...
        2:10; ...
        3 NaN 5 6 NaN 8 9 10 11;   % NaN (not excluded by USPR).
        NaN(1, 9)];
      USPR = [5 6 7 8]';
      I_LSF      = 2;
      FREQ_HZ    = solo.hwzv.const.LSF_HZ(I_LSF);

      % NOTE: Split data due to NaN. bicas.tf.apply_TF() splits by NaN (by
      % default).
      actSamplesAVolt = bicas.proc.L1L2.dc.calibrate_voltages_1xBLTS( ...
        ... % Variables which DO NOT VARY over CDF records at all.
        Vcal         = Vcal, ...
        L            = testCase.L, ...
        iBlts        = 1, ...
        isLfr        = true, ...
        isTdsCwf     = false, ...
        hasSwfFormat = true, ...
        ... % Variables which DO VARY over CDF records.
        tt2000       = uint64( [1:4]' * 1e9 / FREQ_HZ), ...
        ssid         = testCase.S(repmat(["DC_V2"], 4, 1)), ...
        samplesTm    = single(SAMPLES_TM), ...
        uspr         = USPR, ...
        freqHz       = repmat(FREQ_HZ, 4, 1), ...
        iLsf         = repmat(I_LSF,   4, 1), ...
        isAchgFpa    = bicas.utils.FPArray(repmat(logical(0), 4, 1), 'NO_FILL_POSITIONS'), ...
        NbriFpa      = bicas.utils.FPArray(repmat(uint8(  1), 4, 1), 'NO_FILL_POSITIONS'), ...
        NbciFpa      = bicas.utils.FPArray(repmat(uint8(  0), 4, 1), 'NO_FILL_POSITIONS'));

      % Emulated calibration
      expSamplesAvolt                 = SAMPLES_TM*2*3 + 1;
      expSamplesAvolt(1, USPR(1)+1:9) = NaN;
      expSamplesAvolt(2, USPR(2)+1:9) = NaN;
      expSamplesAvolt(3, USPR(3)+1:9) = NaN;
      % NOTE: Setting one extra NaN due to
      % PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES.
      expSamplesAvolt(3, 1)           = NaN;
      testCase.assertEqual(actSamplesAVolt, expSamplesAvolt, AbsTol=1e-14)
    end



    function test_get_SSID_SDID_arrays(testCase)
      S = testCase.S;
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



    function test_get_VSIB_5xBLTS___CWF(testCase)
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
      actBltsVsibAr = bicas.proc.L1L2.dc.get_VSIB_5xBLTS(...
        bltsSamplesAVoltAr, false, uspr, ...
        bltsSsidAr, isAchgFpa, SatSettings, L);

      testCase.assertEqual(actBltsVsibAr, expBltsVsibAr)
    end



    function test_get_VSIB_5xBLTS___SWF(testCase)
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
      actBltsVsibAr = bicas.proc.L1L2.dc.get_VSIB_5xBLTS(...
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



    function Bso = get_simple_TF_BSO(testCase)
      % IMPLEMENTATION NOTE: The default configuration makes it hard to predict
      % the behaviour for bicas.tf.apply_TF() even for simple TFs and data.
      % Therefore deactivating multiple features.
      Bso = bicas.create_default_BSO();

      Bso.override_value('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE', -1,    'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED',    false, 'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE', -1,    'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.FV_SPLITTING.ENABLED',      true,  'test');
      % NOTE: FV_SPLITTING.MIN_SAMPLES is independent of FV_SPLITTING.ENABLED.
      Bso.override_value('PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES',  2,     'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.HIGH_FREQ_LIMIT_FRACTION',  Inf,   'test');

      Bso.make_read_only()
    end



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
