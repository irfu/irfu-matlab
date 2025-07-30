%
% Class that collects functionality related to *DETECTING* saturation.
%
% NOTE: Excludes functionality for setting saturation quality bits and
% QUALITY_FLAG via NSO table.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef sat
  % PROPOSAL: Reorder methods to make easier to follow calls recursively within
  %           class.
  %
  % PROPOSAL: Better terminology than "upperTreshold".
  %   NOTE: Applies to variables and settings keys.
  %   ~abs, absolute, magnitude
  %   ~positive
  %   ~threshold
  %   PROPOSAL: Drop any other prefix/suffix.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Extract relevant settings to struct which can be used internally.
    %
    % NOTE: Useful for tests.
    function SatSettings = from_BSO_extract_saturation_settings(Bso)
      S = struct();

      % How long the sliding window should be when using CDF data.
      S.cwfSlidingWindowLengthSec   = Bso.get_fv('PROCESSING.SATURATION.CWF_SLIDING_WINDOW_LENGTH_SEC');
      % Threshold for the sample-length weighted fraction of VSTB-labelled
      % samples within either (1) a sliding window (CWF), or (2) snapshot. If
      % fraction of VSTB-labelled samples excedes this fraction, then the entire
      % sliding window or snapshot is labelled as saturated.
      S.vstbFractionThreshold       = Bso.get_fv('PROCESSING.SATURATION.SAMPLE_FRACTION_THRESHOLD');

      % Higher thresholds for saturation. Sample values above these values, or
      % below the negated value, count as threshold-saturated (VSIB).
      S.upperThresholdAVoltDcSingle = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.SINGLE');
      S.upperThresholdAVoltDcDiff   = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.DIFF');
      S.upperThresholdAVoltAclg     = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.LOW_GAIN');
      S.upperThresholdAVoltAchg     = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.HIGH_GAIN');

      % ==========
      % ASSERTIONS
      % ==========
      function assert_positive_float(x)
        % NOTE: Positive, not non-negative.
        assert(isfinite(x) && isscalar(x) && isfloat(x) && (x > 0))
      end

      assert_positive_float(S.cwfSlidingWindowLengthSec)
      assert(...
        isfinite(S.vstbFractionThreshold) && ...
        isscalar(S.vstbFractionThreshold) && ...
        isfloat( S.vstbFractionThreshold) && ...
        (0 <= S.vstbFractionThreshold) && (S.vstbFractionThreshold <= 1))

      assert_positive_float(S.upperThresholdAVoltDcSingle)
      assert_positive_float(S.upperThresholdAVoltDcDiff)
      assert_positive_float(S.upperThresholdAVoltAclg)
      assert_positive_float(S.upperThresholdAVoltAchg)

      % Rename return value.
      SatSettings = S;
    end



    % Derive VSTBs. Vectorized.
    function vstbAr = get_VSTB(...
        SatSettings, samplesAVoltAr, ssidAr, isAchgFpa)

      upperThresholdAVolt = bicas.proc.L1L2.sat.get_upper_thresholds(...
        SatSettings, ssidAr, isAchgFpa);

      vstbAr = abs(samplesAVoltAr) > upperThresholdAVolt;
    end



    % Get saturation thresholds. Vectorized.
    %
    % ARGUMENTS
    % =========
    % ssidAr
    % isAchgFpa
    %       Must have the same size as ssidAr.
    %
    % RETURN VALUE
    % ============
    % upperThresholdAVoltAr
    %       Same size as ssidAr. Non-negative threshold values.
    %       Non-ASR SSIDs have threshold NaN.
    %
    function upperThresholdAVoltAr = get_upper_thresholds(...
        SatSettings, ssidAr, isAchgFpa)
      % Tmk = bicas.utils.Timekeeper('get_upper_thresholds', L);

      assert(bicas.proc.L1L2.const.is_SSID(ssidAr))
      assert(isa(isAchgFpa, 'bicas.utils.FPArray') & strcmp(isAchgFpa.mc, "logical"))
      assert(isequal(size(ssidAr), size(isAchgFpa)))

      bIsAsr    = bicas.proc.L1L2.const.SSID_is_ASR( ssidAr);
      bIsDiff   = bicas.proc.L1L2.const.SSID_is_diff(ssidAr);
      bIsAc     = bicas.proc.L1L2.const.SSID_is_AC(  ssidAr);

      % NOTE: Must use bIsAsr to ensure that one excludes non-ASRs.
      bDcSingle = bIsAsr & ~bIsDiff;
      bDcDiff   = bIsAsr &  bIsDiff & ~bIsAc;
      bAcDiff   = bIsAsr &  bIsDiff &  bIsAc;

      % NOTE: Excluding fill positions in the FPA (i.e. bAclg and bAchg do not
      %       cover all elements).
      bAclg = bAcDiff & ~isAchgFpa.array(true);
      bAchg = bAcDiff &  isAchgFpa.array(false);

      % NOTE: Threshold set to NaN for elements without a set threshold (e.g.
      %       SSID=GND).
      %   NOTE: Inequality with NaN always gives false(!)
      upperThresholdAVoltAr            = NaN(size(ssidAr));
      upperThresholdAVoltAr(bDcSingle) = SatSettings.upperThresholdAVoltDcSingle;
      upperThresholdAVoltAr(bDcDiff)   = SatSettings.upperThresholdAVoltDcDiff;
      upperThresholdAVoltAr(bAclg)     = SatSettings.upperThresholdAVoltAclg;
      upperThresholdAVoltAr(bAchg)     = SatSettings.upperThresholdAVoltAchg;

      % Tmk.stop_log(numel(ssidAr), 'element')
    end



  end    % methods(Static)



end
