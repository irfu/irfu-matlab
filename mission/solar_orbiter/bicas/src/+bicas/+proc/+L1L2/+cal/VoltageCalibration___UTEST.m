%
% matlab.unittest automatic test code for
% bicas.proc.L1L2.cal.VoltageCalibration.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef VoltageCalibration___UTEST < matlab.unittest.TestCase



  %#####################
  %#####################
  % CONSTANT PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    SSID = bicas.proc.L1L2.const.C.SSID_DICT;
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_calibrate_voltage_TM_to_avolt___nominal(testCase)
      SAMPLES_TM_CA = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (2*3*x+1), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      testCase.test_basic(...
        bltsSamplesTmCa  =SAMPLES_TM_CA, ...
        expSamplesAvoltCa=expSamplesAvoltCa, ...
        ufv              =false)
    end



    function test_calibrate_voltage_TM_to_avolt___scalar_BIAS(testCase)
      SAMPLES_TM_CA     = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (4*3*x+1), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      testCase.test_basic(...
        bltsSamplesTmCa  =SAMPLES_TM_CA, ...
        expSamplesAvoltCa=expSamplesAvoltCa, ...
        useBiasTfScalar  =true)
    end



    function test_calibrate_voltage_TM_to_avolt___BIAS_offsets_disabled(testCase)
      SAMPLES_TM_CA     = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (2*3*x), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      testCase.test_basic(...
        bltsSamplesTmCa    =SAMPLES_TM_CA, ...
        expSamplesAvoltCa  =expSamplesAvoltCa, ...
        biasOffsetsDisabled=true)
    end



    % (LFR) AC
    % Also tests that special case for AC does not crash.
    function test_calibrate_voltage_TM_to_avolt___AC(testCase)
      SAMPLES_TM_CA     = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (2*3*x+1), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      testCase.test_basic(ufv=false, ...
        bltsSamplesTmCa  =SAMPLES_TM_CA, ...
        expSamplesAvoltCa=expSamplesAvoltCa, ...
        ssid             =testCase.SSID("AC_V12"))
    end



    function test_calibrate_voltage_TM_to_avolt___ASR_UFV(testCase)
      SAMPLES_TM_CA     = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (NaN*x), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      testCase.test_basic(...
        bltsSamplesTmCa  =SAMPLES_TM_CA, ...
        expSamplesAvoltCa=expSamplesAvoltCa, ...
        ufv              =true)
    end



  % Non-ASR SSIDs, with/without UFV.
    function test_calibrate_voltage_TM_to_avolt___Non_ASR_w_wo_UFV(testCase)
      SAMPLES_TM_CA = {[1 2 3]'; [4 5 6 7]'};

      for ssid = testCase.SSID(["GND", "REF25V", "UNKNOWN"])
        expSamplesAvoltCa = cellfun(@(x) (1*x), SAMPLES_TM_CA, ...
          'UniformOutput', false);
        testCase.test_basic(...
          bltsSamplesTmCa  =SAMPLES_TM_CA, ...
          expSamplesAvoltCa=expSamplesAvoltCa, ...
          ssid             =ssid)

        % +UFV
        expSamplesAvoltCa = cellfun(@(x) (NaN*x), SAMPLES_TM_CA, ...
          'UniformOutput', false);
        testCase.test_basic(...
          bltsSamplesTmCa  =SAMPLES_TM_CA, ...
          expSamplesAvoltCa=expSamplesAvoltCa, ...
          ssid             =ssid, ...
          ufv              =true)
      end
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % "Basic test" which can be altered to test changing some specific
    % variables not relating to switching between LFR, TDS-CWF, TDS-RSWF ITFs.
    %
    % Concrete tests assume that the caller knows the hardcoded values.
    %
    function test_basic(testCase, A)
      arguments
        testCase
        A.bltsSamplesTmCa
        A.expSamplesAvoltCa
        A.ufv                 = false
        A.ssid                = testCase.SSID("DC_V1");
        A.biasOffsetsDisabled = false
        A.useBiasTfScalar     = false
      end

      Bso = testCase.get_simple_TF_BSO(...
        biasOffsetsDisabled=A.biasOffsetsDisabled, ...
        useBiasTfScalar    =A.useBiasTfScalar);

      dtSec       = 1/solo.hwzv.const.LFR_F2_HZ * ones(size(A.bltsSamplesTmCa));

      isLfr       = true;
      isTdsCwf    = false;
      zvcti       = [1, 1];

      iBlts       = 1;
      isAchg      = true;
      iCalibTimeL = 1;
      iCalibTimeH = 1;
      iLsf        = 2;
      CalSettings = bicas.proc.L1L2.CalibrationSettings(...
        iBlts, A.ssid, isAchg, iCalibTimeL, iCalibTimeH, iLsf);

      Vcds = bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierTest( ...
        itfBiasAvpiv =@(omegaRps) (ones(size(omegaRps)) * 2), ...
        offsetAvolt  =1, ...
        kItfBiasAvpiv=4, ...   % NOTE: Different factor from in itfBiasAvpiv.
        itfLfrAvpiv  =@(omegaRps) (ones(size(omegaRps)) * 3));
      Vcal = bicas.proc.L1L2.cal.VoltageCalibration(Vcds, true, Bso);

      % CALL TESTED CODE
      actSamplesAvoltCa = Vcal.calibrate_voltage_TM_to_avolt( ...
        dtSec, A.bltsSamplesTmCa, isLfr, isTdsCwf, CalSettings, zvcti, A.ufv);

      testCase.assertEqual(actSamplesAvoltCa, A.expSamplesAvoltCa)
    end



  end    % methods(Access=private)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function Bso = get_simple_TF_BSO(A)
      arguments
        A.biasOffsetsDisabled = false
        A.useBiasTfScalar     = false
      end
      Bso = bicas.create_default_BSO();

      if A.useBiasTfScalar
        biasTfType = 'SCALAR';
      else
        biasTfType = 'FULL';
      end

      % IMPLEMENTATION NOTE: The default configuration makes it hard to predict
      % the behaviour for bicas.tf.apply_TF() even for simple TFs and data.
      % Therefore deactivating multiple features.
      Bso.override_value('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE', -1,    'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED',    false, 'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE', -1,    'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.FV_SPLITTING.ENABLED',      false, 'test');
      % NOTE: FV_SPLITTING.MIN_SAMPLES is independent of FV_SPLITTING.ENABLED.
      Bso.override_value('PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES',  1,     'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.HIGH_FREQ_LIMIT_FRACTION',  Inf,   'test');

      Bso.override_value('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED', A.biasOffsetsDisabled, 'test');
      Bso.override_value('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF',               biasTfType,            'test');

      Bso.make_read_only()
    end



  end    % methods(Static, Access=private)



end
