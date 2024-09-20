%
% Class that collects functionality related to *DETECTING* saturation.
%
% NOTE: Excludes functionality for setting saturation quality bits and
% QUALITY_FLAG via NSO table.
%
%
% IMPLEMENTATION NOTE
% ===================
% Class is designed as an instantiable class in order to:
% (1) reduce the number of arguments (eliminates arguments that configure
%     the saturation criteria),
% (2) only extract saturation criteria from BSO once, to possibly increase
%     performance.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Saturation
  % PROPOSAL: Merge higherThresholdAVolt* to a struct somehow.
  %
  % PROBLEM: Can not deduce saturation limits for channels reconstructed from
  %          other channels.
  %   Ex: BDM=4 ==> Have DC_V1/V2/V3 ==> Derives DC_V12/V13/V23 ==> Saturation
  %       on DC diffs can only be deduced from the samples from which the
  %       samples originate.
  % PROBLEM: get_voltage_saturation_quality_bit() can not correctly handle
  %          AsrSamplesAVoltSrm data which derives from non-ASR channels.
  %   Ex: BDM=5-7 ==> 2.5V Ref/GND stored in AsrSamplesAVoltSrm, but are
  %       represented by ASIDs.
  %
  % PROPOSAL: Only detect saturation in BLTSs (which are true antenna signals).
  %   Have saturation bits propagate to all signals in (not-yet-implemented)
  %   _Sdid_SamplesAVoltSrm in the same way as signals do.
  %   PROPOSAL: Separate SDID SRM for quality bits.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(GetAccess=private, SetAccess=immutable)
    % How long the sliding window should be when using CDF data.
    cwfSlidingWindowLengthSec

    % Threshold for the sample-length weighted fraction of TSF-labelled samples
    % within either (1) a sliding window (CWF), or (2) snapshot. If fraction of
    % TSF-labelled samples excedes this fraction, then the entire sliding window
    % or snapshot is labelled as saturated.
    tsfFractionThreshold

    % Higher thresholds for saturation. Sample values above these values, or
    % below the negated value, count as threshold-saturated (TSF).
    higherThresholdAVoltDcSingle
    higherThresholdAVoltDcDiff
    higherThresholdAVoltAclg
    higherThresholdAVoltAchg
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = Saturation(Bso)
      obj.cwfSlidingWindowLengthSec    = Bso.get_fv('PROCESSING.SATURATION.CWF_SLIDING_WINDOW_LENGTH_SEC');
      obj.tsfFractionThreshold         = Bso.get_fv('PROCESSING.SATURATION.TSF_FRACTION_THRESHOLD');

      obj.higherThresholdAVoltDcSingle = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.SINGLE');
      obj.higherThresholdAVoltDcDiff   = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.DIFF');
      obj.higherThresholdAVoltAclg     = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.LOW_GAIN');
      obj.higherThresholdAVoltAchg     = Bso.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.HIGH_GAIN');



      % ==========
      % ASSERTIONS
      % ==========
      function assert_positive_float(x)
        % NOTE: Positive, not non-negative.
        assert(isfinite(x) && isscalar(x) && isfloat(x) && (x > 0))
      end

      assert_positive_float(obj.cwfSlidingWindowLengthSec)
      assert(...
        isfinite(obj.tsfFractionThreshold) && ...
        isscalar(obj.tsfFractionThreshold) && ...
        isfloat( obj.tsfFractionThreshold) && ...
        (0 <= obj.tsfFractionThreshold) && (obj.tsfFractionThreshold <= 1))

      assert_positive_float(obj.higherThresholdAVoltDcSingle)
      assert_positive_float(obj.higherThresholdAVoltDcDiff)
      assert_positive_float(obj.higherThresholdAVoltAclg)
      assert_positive_float(obj.higherThresholdAVoltAchg)
    end



    % Given an arbitrary-size ARRAY of samples, get TSF bits for every
    % sample.
    %
    % NOTE: The data may refer to both CWF data and SWF data, but the
    % function itself makes no distinction between the two. The caller has
    % to make distinctions between those two if needed. For example, this
    % function returns TSFs for each sample in a snapshot, but the caller
    % might one to condense this to one saturation bit per snapshot
    % according to some algorithm that has no analogue for CWF data.
    %
    %
    % ARGUMENTS
    % =========
    % samplesAVolt
    %       Arbitrary-size array. May contain NaN.
    %
    %
    % RETURN VALUE
    % ============
    % tsfAr
    %       Float. Same size as samplesAVolt. Whether corresponding elements
    %       in samplesAVolt are deemed to be outside the relevant
    %       thresholds. False is returned for all input elements if there
    %       are no thresholds for this kind of data (e.g. for non-ASR
    %       sources). False is returned for NaN input elements.
    %
    function tsfAr = get_TSF(obj, samplesAVolt, Ssid, isAchgFpa)
      % PROPOSAL: Better name.
      %   ~sample-to-TSF
      %   ~threshold_saturation

      assert(isfloat(samplesAVolt))
      assert(isa(Ssid, 'bicas.proc.L1L2.SignalSourceId'))
      assert(isa(isAchgFpa, 'bicas.utils.FPArray') && isscalar(isAchgFpa))

      % Default value that used if there are no thresholds.
      tsfAr = false(size(samplesAVolt));

      if ~Ssid.is_ASR()
        return
      end

      % CASE: ASR (i.e. no non-plasma/unknown signal, no special case)

      % ====================
      % Determine thresholds
      % ====================
      if Ssid.Asid.is_diff()
        % CASE: DC/AC diff
        % ----------------

        isAchg = isAchgFpa.logical2doubleNan();
        if Ssid.Asid.is_AC()
          % CASE: AC diff
          % -------------
          if isAchg == 0
            highThresholdAVolt = obj.higherThresholdAVoltAclg;
          elseif isAchg == 1
            highThresholdAVolt = obj.higherThresholdAVoltAchg;
          else
            return
          end
        else
          % CASE: DC diff
          % -------------
          highThresholdAVolt = obj.higherThresholdAVoltDcDiff;
        end
      else
        % CASE: DC single
        % ---------------
        % NOTE: Not using terms "min" and "max" since they are
        % ambiguous (?).
        highThresholdAVolt = obj.higherThresholdAVoltDcSingle;
      end
      lowerThresholdAVolt = -highThresholdAVolt;

      % =========================================
      % Use thresholds on array to determine TSFs
      % =========================================
      % NOTE: Has to be able ignore NaN.
      tsfAr = (samplesAVolt < lowerThresholdAVolt) | (highThresholdAVolt < samplesAVolt);
    end



    % Determine whether ONE snapshot should be labelled as saturated.
    %
    % ARGUMENTS
    % =========
    % samplesAVolt
    %   Snapshot samples. (1, iSampleInSnapshot) = row vector.
    %   NOTE: Should only contain the length of the snapshot. No padding at
    %         the end of array.
    %
    % RETURN VALUE
    % ============
    % isSaturated
    %       Logical. Scalar.
    %
    function isSaturated = get_snapshot_saturation(obj, samplesAVolt, Ssid, isAchg)
      irf.assert.sizes(samplesAVolt, [1, NaN, 1])     % Row vector.

      tsfAr = obj.get_TSF(samplesAVolt, Ssid, isAchg);

      isSaturated = (sum(tsfAr, 'all') / numel(samplesAVolt)) > obj.tsfFractionThreshold;
    end



    % Determine whether multiple snapshots (with same settings) are
    % saturated. Uses ZV-like variables.
    %
    % ARGUMENTS
    % =========
    % zvNValidSamplesPerRecord
    %       ZV-like array. (iCdfRecord). Length of separate snapshots.
    % zvSamplesAVolt
    %       ZV-like array. (iCdfRecord, iSampleInSnapshot)
    function isSaturatedAr = get_snapshot_saturation_many(obj, ...
        zvNValidSamplesPerRecord, zvSamplesAVolt, Ssid, isAchgFpa)

      nRecs = irf.assert.sizes(...
        zvNValidSamplesPerRecord, [-1],  ...
        zvSamplesAVolt,           [-1, NaN, 1]);

      isSaturatedAr = false(nRecs, 1);
      for iRec = 1:nRecs
        isSaturatedAr(iRec) = obj.get_snapshot_saturation(...
          zvSamplesAVolt(iRec, 1:zvNValidSamplesPerRecord(iRec)), ...
          Ssid, isAchgFpa);
      end
    end



    % Given ZV-like variables, get saturation bits for quality bitmask.
    %
    % NOTE: Applies to both CWF and SWF data.
    %
    % PROBLEM: Function is conceptually bad (buggy) for edge cases (non-antenna
    % signals, reconstructed signals). See unofficial class comments.
    %
    %
    % RETURN VALUE
    % ============
    % isSaturatedAr
    %       (iCdfRecords). Logical. Quality bit for saturation.
    %
    function isSaturatedAr = get_voltage_saturation_quality_bit(...
        obj, tt2000Ar, AsrSamplesAVoltSrm, zvNValidSamplesPerRecord, ...
        bltsKSsidAr, isAchgFpa, hasSwfFormat, L)
      % PROPOSAL: Vectorize. Obtain vectors of thresholds for each channel. Then
      %           look for saturation.
      %   NOTE: Only ACHG influences the calibration thresholds for each channel
      %         (SDID/ASR). Could otherwise have scalar values per channel.
      %   PRO: Easier to keep track of what thresholds are a function of.

      % ASSERTIONS
      bicas.utils.assert_ZV_Epoch(tt2000Ar)
      assert(islogical(hasSwfFormat) && isscalar(hasSwfFormat))
      assert(isa(bltsKSsidAr, 'uint8'))
      nRows = irf.assert.sizes(...
        tt2000Ar,                 [-1], ...
        zvNValidSamplesPerRecord, [-1], ...
        bltsKSsidAr,              [-1, bicas.const.N_BLTS]);
      assert(isa(AsrSamplesAVoltSrm, "bicas.utils.SameRowsMap"))
      assert(AsrSamplesAVoltSrm.nRows == nRows)



      L.logf('info', ...
        ['Detecting threshold saturation (voltages) -', ...
        ' One sequence of records with identical settings at a time.'])
      Tmk = bicas.utils.Timekeeper('get_voltage_saturation_quality_bit', L);

      % IMPLEMENTATION NOTE: Below code for cases CWF and SWF do ~duplicate
      % code, but it is difficult to use the same implementation for both
      % without (1) making the impleementation harder to understand and (2)
      % having one particular variable with different meanings in the two cases.
      if ~hasSwfFormat
        %===========
        % CASE: CWF
        %===========
        tsfAr = false(nRows, 1);
        for Asid = AsrSamplesAVoltSrm.keys'
          asidTsfAr = obj.get_one_ASR_CWF_channel_TSF_bit_array(...
            bicas.proc.L1L2.SignalSourceId(Asid), isAchgFpa, ...
            AsrSamplesAVoltSrm(Asid));

          % Merge (OR) bits over ASIDs.
          tsfAr = any([tsfAr, asidTsfAr], 2);
        end

        isSaturatedAr = bicas.proc.L1L2.qual.sliding_window_over_fraction(...
          tt2000Ar, tsfAr, ...
          obj.tsfFractionThreshold, obj.cwfSlidingWindowLengthSec);
      else
        %===========
        % CASE: SWF
        %===========
        isSaturatedAr = false(nRows, 1);
        for Asid = AsrSamplesAVoltSrm.keys'
          asidIsSaturatedAr = obj.get_one_ASR_SWF_channel_saturation_bit_array(...
            bicas.proc.L1L2.SignalSourceId(Asid), isAchgFpa, ...
            AsrSamplesAVoltSrm(Asid), zvNValidSamplesPerRecord);

          % Merge (OR) bits over ASIDs.
          isSaturatedAr = any([isSaturatedAr, asidIsSaturatedAr], 2);
        end
      end



      if hasSwfFormat
        Tmk.stop_log(nRows, 'CDF record')
      else
        % Log some saturation statistics which may help tell whether how
        % much the saturation varies over time, which may
        % influence/explain if the above processing is slow. Should only
        % be relevant for CWF.
        % NOTE: Only reflects the behaviour of the final saturation bit,
        % not the TSF.
        nSaturationChanges = numel(find(isSaturatedAr(1:end-1) ~= isSaturatedAr(2:end)));
        Tmk.stop_log(nRows, 'CDF record', nSaturationChanges, 'sat. flag change')
        L.logf('debug', 'SPEED -- %g [CDF rows/sat. flag change]', nRows/nSaturationChanges)
      end

    end    % function



    % Return TSF for CWF data.
    function tsfAr = get_one_ASR_CWF_channel_TSF_bit_array(obj, Ssid, isAchgFpa, samplesAVolt)
      nRows = irf.assert.sizes( ...
        Ssid,         [1], ...
        isAchgFpa,    [-1], ...
        samplesAVolt, [-1]);

      % NOTE: Splits into subsequences also when ACHG does not matter (DC).
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        isAchgFpa.logical2doubleNan());

      tsfAr = false(nRows, 1);

      for iSs = 1:nSs
        iRec1   = iRec1Ar(iSs);
        iRec2   = iRec2Ar(iSs);

        tsfAr(iRec1:iRec2) = obj.get_TSF(...
          samplesAVolt(iRec1:iRec2), Ssid, isAchgFpa(iRec1));
      end
    end



    % Return final saturation bit for SWF data.
    function saturationBitAr = get_one_ASR_SWF_channel_saturation_bit_array(...
        obj, Ssid, isAchgFpa, samplesAVolt, zvNValidSamplesPerRecord)
      [nRows, ~] = irf.assert.sizes( ...
        Ssid,                     [1], ...
        isAchgFpa,                [-1], ...
        samplesAVolt,             [-1, -2], ...
        zvNValidSamplesPerRecord, [-1]);

      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        isAchgFpa.logical2doubleNan());

      saturationBitAr = false(nRows, 1);

      for iSs = 1:nSs
        iRec1   = iRec1Ar(iSs);
        iRec2   = iRec2Ar(iSs);

        saturationBitAr(iRec1:iRec2) = obj.get_snapshot_saturation_many(...
          zvNValidSamplesPerRecord(iRec1:iRec2), ...
          samplesAVolt(            iRec1:iRec2, :), ...
          Ssid, isAchgFpa(iRec1));
      end

    end



  end    % methods(Access=public)



end
