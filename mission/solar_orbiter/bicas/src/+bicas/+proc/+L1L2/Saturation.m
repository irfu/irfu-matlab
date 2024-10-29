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
  % PROBLEM: get_VSQB() can not correctly handle
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

    % Threshold for the sample-length weighted fraction of VSTB-labelled samples
    % within either (1) a sliding window (CWF), or (2) snapshot. If fraction of
    % VSTB-labelled samples excedes this fraction, then the entire sliding window
    % or snapshot is labelled as saturated.
    vstbFractionThreshold

    % Higher thresholds for saturation. Sample values above these values, or
    % below the negated value, count as threshold-saturated (VSTB).
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
      obj.vstbFractionThreshold        = Bso.get_fv('PROCESSING.SATURATION.VSTB_FRACTION_THRESHOLD');

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
        isfinite(obj.vstbFractionThreshold) && ...
        isscalar(obj.vstbFractionThreshold) && ...
        isfloat( obj.vstbFractionThreshold) && ...
        (0 <= obj.vstbFractionThreshold) && (obj.vstbFractionThreshold <= 1))

      assert_positive_float(obj.higherThresholdAVoltDcSingle)
      assert_positive_float(obj.higherThresholdAVoltDcDiff)
      assert_positive_float(obj.higherThresholdAVoltAclg)
      assert_positive_float(obj.higherThresholdAVoltAchg)
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
    % vsqbAr
    %       (iCdfRecords). Logical.
    %
    function vsqbAr = get_VSQB(...
        obj, tt2000Ar, AsrSamplesAVoltSrm, zvNValidSamplesPerRecord, ...
        bltsSsidAr, isAchgFpa, hasSwfFormat, L)
      % PROPOSAL: Vectorize. Obtain vectors of thresholds for each channel. Then
      %           look for saturation.
      %   NOTE: Only ACHG influences the calibration thresholds for each channel
      %         (SDID/ASR). Could otherwise have scalar values per channel.
      %   PRO: Easier to keep track of what thresholds are a function of.

      % ASSERTIONS
      bicas.utils.assert_ZV_Epoch(tt2000Ar)
      assert(islogical(hasSwfFormat) && isscalar(hasSwfFormat))
      assert(bicas.proc.L1L2.const.is_SSID(bltsSsidAr))
      nRows = irf.assert.sizes(...
        tt2000Ar,                 [-1], ...
        zvNValidSamplesPerRecord, [-1], ...
        bltsSsidAr,               [-1, bicas.const.N_BLTS]);
      assert(isa(AsrSamplesAVoltSrm, "bicas.utils.SameRowsMap"))
      assert(AsrSamplesAVoltSrm.nRows == nRows)



      L.logf('info', ...
        ['Detecting threshold saturation (voltages) -', ...
        ' One sequence of records with identical settings at a time.'])
      Tmk = bicas.utils.Timekeeper('get_VSQB', L);

      % IMPLEMENTATION NOTE: Below code for cases CWF and SWF do ~duplicate
      % code, but it is difficult to use the same implementation for both
      % without (1) making the implementation harder to understand and (2)
      % having one particular variable with different meanings in the two cases.
      if ~hasSwfFormat
        %===========
        % CASE: CWF
        %===========
        vstbAr = false(nRows, 1);
        for asid = AsrSamplesAVoltSrm.keys'
          asidVstbAr = obj.get_ASR_CWF_channel_VSTB(...
            bicas.proc.L1L2.const.ASID_to_SSID(asid), isAchgFpa, ...
            AsrSamplesAVoltSrm(asid));

          % Merge (OR) bits over ASIDs.
          vstbAr = any([vstbAr, asidVstbAr], 2);
        end

        vsqbAr = bicas.proc.L1L2.qual.sliding_window_over_fraction(...
          tt2000Ar, vstbAr, ...
          obj.vstbFractionThreshold, obj.cwfSlidingWindowLengthSec);
      else
        %===========
        % CASE: SWF
        %===========
        vsqbAr = false(nRows, 1);
        for asid = AsrSamplesAVoltSrm.keys'
          asidVsqbAr = obj.get_ASR_SWF_channel_VSQB(...
            bicas.proc.L1L2.const.ASID_to_SSID(asid), isAchgFpa, ...
            AsrSamplesAVoltSrm(asid), zvNValidSamplesPerRecord);

          % Merge (OR) bits over ASIDs.
          vsqbAr = any([vsqbAr, asidVsqbAr], 2);
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
        % not the VSTB.
        nSaturationChanges = numel(find(vsqbAr(1:end-1) ~= vsqbAr(2:end)));
        Tmk.stop_log(nRows, 'CDF record', nSaturationChanges, 'sat. flag change')
        L.logf('debug', 'SPEED -- %g [CDF rows/sat. flag change]', nRows/nSaturationChanges)
      end

    end    % function



    % Return VSTB (not VSQB) for one channel of CWF data.
    function vstbAr = get_ASR_CWF_channel_VSTB(obj, ssid, isAchgFpa, samplesAVolt)
      nRows = irf.assert.sizes( ...
        ssid,         [ 1], ...
        isAchgFpa,    [-1], ...
        samplesAVolt, [-1]);

      % NOTE: Splits into subsequences also when ACHG does not matter (DC).
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        isAchgFpa.logical2doubleNan());

      vstbAr = false(nRows, 1);

      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        vstbAr(iRec1:iRec2) = obj.get_VSTB(...
          samplesAVolt(iRec1:iRec2), ssid, isAchgFpa(iRec1));
      end
    end



    % Return VSQB (not VSTB) for one (ASR) channel of SWF data.
    function vsqbAr = get_ASR_SWF_channel_VSQB(...
        obj, ssid, isAchgFpa, samplesAVolt, zvNValidSamplesPerRecord)
      [nRows, ~] = irf.assert.sizes( ...
        ssid,                     [ 1], ...
        isAchgFpa,                [-1], ...
        samplesAVolt,             [-1, -2], ...
        zvNValidSamplesPerRecord, [-1]);

      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        isAchgFpa.logical2doubleNan());

      vsqbAr = false(nRows, 1);

      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        vsqbAr(iRec1:iRec2) = obj.get_snapshot_VSQB_many(...
          zvNValidSamplesPerRecord(iRec1:iRec2), ...
          samplesAVolt(            iRec1:iRec2, :), ...
          ssid, isAchgFpa(iRec1));
      end

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
    %
    function vsqbAr = get_snapshot_VSQB_many(obj, ...
        zvNValidSamplesPerRecord, zvSamplesAVolt, ssid, isAchgFpa)

      nRecs = irf.assert.sizes(...
        zvNValidSamplesPerRecord, [-1],  ...
        zvSamplesAVolt,           [-1, NaN, 1]);

      vsqbAr = false(nRecs, 1);
      for iRec = 1:nRecs
        vsqbAr(iRec) = obj.get_snapshot_VSQB(...
          zvSamplesAVolt(iRec, 1:zvNValidSamplesPerRecord(iRec)), ...
          ssid, isAchgFpa);
      end
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
    % vsqb
    %       Logical. Scalar.
    %
    function vsqb = get_snapshot_VSQB(obj, samplesAVolt, ssid, isAchg)
      assert(isrow(samplesAVolt))     % Row vector(!).

      vstbAr = obj.get_VSTB(samplesAVolt, ssid, isAchg);

      vsqb = (sum(vstbAr, 'all') / numel(samplesAVolt)) > obj.vstbFractionThreshold;
    end



    % Given an arbitrary-size ARRAY of samples, get VSTB bits for every
    % sample.
    %
    % NOTE: The data may refer to both CWF data and SWF data, but the
    % function itself makes no distinction between the two. The caller has
    % to make distinctions between those two if needed. For example, this
    % function returns VSTBs for each sample in a snapshot, but the caller
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
    % vstbAr
    %       Float. Same size as samplesAVolt. Whether corresponding elements
    %       in samplesAVolt are deemed to be outside the relevant
    %       thresholds. False is returned for all input elements if there
    %       are no thresholds for this kind of data (e.g. for non-ASR
    %       sources). False is returned for NaN input elements.
    %
    function vstbAr = get_VSTB(obj, samplesAVolt, ssid, isAchgFpa)
      % PROPOSAL: Better name.
      %   ~sample-to-VSTB
      %   ~threshold_saturation

      assert(isfloat(samplesAVolt))
      assert(bicas.proc.L1L2.const.is_SSID(ssid) & isscalar(ssid))
      assert(isa(isAchgFpa, 'bicas.utils.FPArray') && isscalar(isAchgFpa))

      % Default value that is used if there are no thresholds.
      vstbAr = false(size(samplesAVolt));

      if ~bicas.proc.L1L2.const.SSID_is_ASR(ssid)
        return
      end

      % CASE: ASR (i.e. no non-plasma/unknown signal, no special case)

      % ====================
      % Determine thresholds
      % ====================
      if bicas.proc.L1L2.const.SSID_is_diff(ssid)
        % CASE: DC/AC diff
        % ----------------

        isAchg = isAchgFpa.logical2doubleNan();
        if bicas.proc.L1L2.const.SSID_is_AC(ssid)
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

      % ==========================================
      % Use thresholds on array to determine VSTBs
      % ==========================================
      % NOTE: Has to be able ignore NaN.
      vstbAr = (samplesAVolt < lowerThresholdAVolt) | (highThresholdAVolt < samplesAVolt);
    end



  end    % methods(Access=public)



end
