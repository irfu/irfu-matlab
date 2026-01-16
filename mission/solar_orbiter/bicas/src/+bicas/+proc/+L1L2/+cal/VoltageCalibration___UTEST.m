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



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  % Technically, additional properties of testCase objects with cell array
  % default values. Test methods with arguments with the same name will be
  % called once for every element in the cell arrays.
  properties(TestParameter)
    % All legal values for setting "PROCESSING.CALIBRATION.TF.METHOD".
    TF_METHOD = {'FFT'; 'KERNEL'}
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_calibrate_voltage_TM_to_avolt___nominal(T, TF_METHOD)
      SAMPLES_TM_CA = {[1 2 3]'; [4 5 6 7]'; zeros(0, 1)};
      expSamplesAvoltCa = cellfun(@(x) (2*3*x+1), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      T.test_basic(...
        bltsSamplesTmCa   = SAMPLES_TM_CA, ...
        expSamplesAvoltCa = expSamplesAvoltCa, ...
        tfMethod          = TF_METHOD)
    end



    % Checks behaviour w.r.t. sample=NaN.
    % NOTE: Basic check on the splitting of 1D arrays of samples which should
    %       be used under the hood.
    function test_calibrate_voltage_TM_to_avolt___NaN(T, TF_METHOD)
      N = NaN;

      SAMPLES_TM_CA = {[N N N]'; [4 N 6 7 N N 10 11 12]'};
      expSamplesAvoltCa = cellfun(@(x) (2*3*x+1), SAMPLES_TM_CA, ...
        'UniformOutput', false);
      % Sample which PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES will
      % set to NaN.
      expSamplesAvoltCa{2}(1) = NaN;    % NOTE: WORKS despite syntax!

      T.test_basic(...
        bltsSamplesTmCa   = SAMPLES_TM_CA, ...
        expSamplesAvoltCa = expSamplesAvoltCa, ...
        tfMethod          = TF_METHOD)
    end



    function test_calibrate_voltage_TM_to_avolt___scalar_BIAS(T, TF_METHOD)
      SAMPLES_TM_CA     = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (4*3*x+1), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      T.test_basic(...
        bltsSamplesTmCa   = SAMPLES_TM_CA, ...
        expSamplesAvoltCa = expSamplesAvoltCa, ...
        useBiasTfScalar   = true, ...
        tfMethod          = TF_METHOD)
    end



    function test_calibrate_voltage_TM_to_avolt___BIAS_offsets_disabled(T, TF_METHOD)
      SAMPLES_TM_CA     = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (2*3*x), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      T.test_basic(...
        bltsSamplesTmCa     = SAMPLES_TM_CA, ...
        expSamplesAvoltCa   = expSamplesAvoltCa, ...
        biasOffsetsDisabled = true, ...
        tfMethod            = TF_METHOD)
    end



    % (LFR) AC
    % Also tests that special case for AC does not crash.
    function test_calibrate_voltage_TM_to_avolt___AC(T, TF_METHOD)
      SAMPLES_TM_CA     = {[1 2 3]'; [4 5 6 7]'};
      expSamplesAvoltCa = cellfun(@(x) (2*3*x+1), SAMPLES_TM_CA, ...
        'UniformOutput', false);

      T.test_basic(...
        bltsSamplesTmCa   = SAMPLES_TM_CA, ...
        expSamplesAvoltCa = expSamplesAvoltCa, ...
        ssid              = T.SSID("AC_V12"), ...
        tfMethod          = TF_METHOD)
    end



    % Non-ASR SSIDs.
    function test_calibrate_voltage_TM_to_avolt___Non_ASR(T, TF_METHOD)
      SAMPLES_TM_CA = {[1 2 3]'; [4 5 6 7]'};

      for ssid = T.SSID(["GND", "REF25V", "UNKNOWN"])
        expSamplesAvoltCa = cellfun(@(x) (1*x), SAMPLES_TM_CA, ...
          'UniformOutput', false);
        T.test_basic(...
          bltsSamplesTmCa   = SAMPLES_TM_CA, ...
          expSamplesAvoltCa = expSamplesAvoltCa, ...
          ssid              = ssid, ...
          tfMethod            = TF_METHOD)
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
    function test_basic(T, A)
      arguments
        T
        A.bltsSamplesTmCa
        A.expSamplesAvoltCa
        A.ssid                = T.SSID("DC_V1");
        A.biasOffsetsDisabled = false
        A.useBiasTfScalar     = false
        A.tfMethod
      end

      Bso = T.get_simple_TF_BSO(...
        biasOffsetsDisabled = A.biasOffsetsDisabled, ...
        useBiasTfScalar     = A.useBiasTfScalar, ...
        tfMethod            = A.tfMethod);

      dtSec       = 1/solo.hwzv.const.LFR_F2_HZ;

      isLfr       = true;
      isTdsCwf    = false;
      NbriFpa     = bicas.utils.FPArray(1, 'NO_FILL_POSITIONS');
      NbciFpa     = bicas.utils.FPArray(1, 'NO_FILL_POSITIONS');

      iBlts       = 1;
      isAchg      = true;
      iCalibTimeL = 1;
      iCalibTimeH = 1;
      iLsf        = 2;
      CalSettings = bicas.proc.L1L2.cal.CalibrationSettings(...
        iBlts, A.ssid, isAchg, iCalibTimeL, iCalibTimeH, iLsf);

      Vcds = bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierTest( ...
        itfBiasAvpiv  = @(omegaRps) (ones(size(omegaRps)) * 2), ...
        offsetAvolt   = 1, ...
        kItfBiasAvpiv = 4, ...   % NOTE: Different factor from in itfBiasAvpiv.
        itfLfrAvpiv   = @(omegaRps) (ones(size(omegaRps)) * 3));
      Vcal = bicas.proc.L1L2.cal.VoltageCalibration(Vcds, true, Bso);

      % CALL TESTED CODE
      actSamplesAvoltCa = Vcal.calibrate_voltage_TM_to_avolt( ...
        dtSec, A.bltsSamplesTmCa, isLfr, isTdsCwf, CalSettings, NbriFpa, NbciFpa);

      T.assertEqual(actSamplesAvoltCa, A.expSamplesAvoltCa)
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
        A.tfMethod
      end
      Bso = bicas.create_default_BSO();

      if A.useBiasTfScalar
        biasTfType = 'SCALAR';
      else
        biasTfType = 'FULL';
      end

      Bso.override_value('PROCESSING.CALIBRATION.TF.METHOD',                     A.tfMethod, 'test');

      % IMPLEMENTATION NOTE: The default configuration makes it hard to predict
      % the behaviour for bicas.tf.apply_TF() even for simple TFs and data.
      % Therefore deactivating multiple features.
      Bso.override_value('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE',  -1,    'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED',     false, 'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE',  -1,    'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.FV_SPLITTING.ENABLED',       true,  'test');
      % NOTE: FV_SPLITTING.MIN_SAMPLES is independent of FV_SPLITTING.ENABLED.
      Bso.override_value('PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES',   2,     'test');
      Bso.override_value('PROCESSING.CALIBRATION.TF.HIGH_FREQ_LIMIT_FRACTION',   Inf,   'test');

      % NOTE: Using a Hann window (which only applies to KERNEL method)
      % effectively scales down also simple TFs (constants) for very short
      % sequences of samples which affects the tests. Therefore disabling
      % functionality to simplify tests.
      Bso.override_value('PROCESSING.CALIBRATION.TF.KERNEL.HANN_WINDOW_ENABLED', false, 'test');

      Bso.override_value('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED', A.biasOffsetsDisabled, 'test');
      Bso.override_value('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF',               biasTfType,            'test');

      Bso.make_read_only()
    end



  end    % methods(Static, Access=private)



end
