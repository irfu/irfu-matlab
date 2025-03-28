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
  % PROPOSAL: Reorg. to non-instantiatble class. -- IMPLEMENTED
  %   PRO: Makes no sense for the caller to instantiate.
  %     PRO: Only instantiated once in code, followed by one method call.
  %     PRO: Constructor only takes one argument.
  %   CON: Constructor argument BSO is split into six relevant values stored in
  %        instance variables.
  %     PROPOSAL: Derive struct/class which is passed around internally between
  %               methods.
  %
  % PROPOSAL: Reorder methods to make easier to follow calls recursively within
  %           class.
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



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Extract relevant settings to struct which can be used internally.
    %
    % NOTE: Useful for tests.
    function SatSettings = convert_BSO_to_struct(Bso)
      S = struct();

      % How long the sliding window should be when using CDF data.
      S.cwfSlidingWindowLengthSec   = Bso.get_fv('PROCESSING.SATURATION.CWF_SLIDING_WINDOW_LENGTH_SEC');
      % Threshold for the sample-length weighted fraction of VSIB-labelled
      % samples within either (1) a sliding window (CWF), or (2) snapshot. If
      % fraction of VSIB-labelled samples excedes this fraction, then the entire
      % sliding window or snapshot is labelled as saturated.
      S.vsibFractionThreshold       = Bso.get_fv('PROCESSING.SATURATION.SAMPLE_FRACTION_THRESHOLD');

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
        isfinite(S.vsibFractionThreshold) && ...
        isscalar(S.vsibFractionThreshold) && ...
        isfloat( S.vsibFractionThreshold) && ...
        (0 <= S.vsibFractionThreshold) && (S.vsibFractionThreshold <= 1))

      assert_positive_float(S.upperThresholdAVoltDcSingle)
      assert_positive_float(S.upperThresholdAVoltDcDiff)
      assert_positive_float(S.upperThresholdAVoltAclg)
      assert_positive_float(S.upperThresholdAVoltAchg)

      % Rename return value.
      SatSettings = S;
    end



    % Given an arbitrary-size ARRAY of samples, with constant SSID and isAchg,
    % get VSIBs for every sample.
    %
    % NOTE: The data may refer to both CWF data and SWF data, but the
    % function itself makes no distinction between the two. The caller has
    % to make distinctions between those two if needed. For example, this
    % function returns VSIBs for each sample in a snapshot, but the caller
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
    % vsibAr
    %       Float. Same size as samplesAVolt. Whether corresponding elements
    %       in samplesAVolt are deemed to be outside the relevant
    %       thresholds. False is returned for all input elements if there
    %       are no thresholds for this kind of data (e.g. for non-ASR
    %       sources). False is returned for NaN input elements.
    %
    function vsibAr = get_VSIB(SatSettings, samplesAVoltAr, ssid, isAchgFpa)
      % PROPOSAL: Better name.
      %   ~sample-to-VSIB
      %   ~threshold_saturation
      %
      % PROPOSAL/TODO: Replace with get_VSIB_NEW().

      assert(isfloat(samplesAVoltAr))
      assert(bicas.proc.L1L2.const.is_SSID(ssid) & isscalar(ssid))
      assert(isa(isAchgFpa, 'bicas.utils.FPArray') && isscalar(isAchgFpa))

      ssidAr    = repmat(ssid,      size(samplesAVoltAr));
      isAchgFpa = repmat(isAchgFpa, size(samplesAVoltAr));
      vsibAr    = bicas.proc.L1L2.sat.get_VSIB_NEW(...
        SatSettings, samplesAVoltAr, ssidAr, isAchgFpa);
    end



    % Vectorized.
    function vsibAr = get_VSIB_NEW(SatSettings, samplesAVoltAr, ssidAr, isAchgFpa)
      upperThresholdAVolt = bicas.proc.L1L2.sat.get_upper_thresholds(...
        SatSettings, ssidAr, isAchgFpa);

      vsibAr = abs(samplesAVoltAr) > upperThresholdAVolt;
    end



    % Get saturation thresholds. Vectorized.
    %
    % ARGUMENTS
    % =========
    % ssidAr
    % isAchgFpa
    %       Must be the same size as ssidAr.
    %
    % RETURN VALUE
    % ============
    % upperThresholdAVoltAr
    %       Same size as ssidAr. Non-negative threshold values.
    %
    function upperThresholdAVoltAr = get_upper_thresholds(...
        SatSettings, ssidAr, isAchgFpa)

      assert(bicas.proc.L1L2.const.is_SSID(ssidAr))
      assert(isa(isAchgFpa, 'bicas.utils.FPArray') & strcmp(isAchgFpa.mc, "logical"))
      assert(isequal(size(ssidAr), size(isAchgFpa)))

      bIsAsr    = bicas.proc.L1L2.const.SSID_is_ASR( ssidAr);
      bIsDiff   = bicas.proc.L1L2.const.SSID_is_diff(ssidAr);
      bIsAc     = bicas.proc.L1L2.const.SSID_is_AC(  ssidAr);

      bDcSingle = bIsAsr & ~bIsDiff;
      bDcDiff   = bIsAsr &  bIsDiff & ~bIsAc;
      bAcDiff   = bIsAsr &  bIsDiff &  bIsAc;

      % NOTE: Excluding fill positions in the FPA (i.e. bAclg and bAchg do not
      %       cover all elements).
      bAclg = bAcDiff & ~isAchgFpa.array(true);
      bAchg = bAcDiff &  isAchgFpa.array(false);

      % NOTE: Threshold set to NaN for element without known threshold (e.g.
      %       SSID=GND).
      %   NOTE: Inequality with NaN always gives false(!)
      upperThresholdAVoltAr = NaN(size(ssidAr));
      upperThresholdAVoltAr(bDcSingle) = SatSettings.upperThresholdAVoltDcSingle;
      upperThresholdAVoltAr(bDcDiff)   = SatSettings.upperThresholdAVoltDcDiff;
      upperThresholdAVoltAr(bAclg)     = SatSettings.upperThresholdAVoltAclg;
      upperThresholdAVoltAr(bAchg)     = SatSettings.upperThresholdAVoltAchg;
    end



    % Determine whether ONE snapshot should be labelled as saturated.
    %
    % ARGUMENTS
    % =========
    % samplesAVolt
    %   Snapshot samples. (1, iSampleInSnapshot) = row vector.
    %   NOTE: Should only contain the actual snapshot data, i.e. no padding at
    %         the end of array.
    %
    % RETURN VALUE
    % ============
    % vsqb
    %       Logical. Scalar.
    %
    function vsqb = get_snapshot_VSQB(...
        SatSettings, samplesAVolt, ssid, isAchg)

      assert(isrow(samplesAVolt))     % Row vector(!).

      vsibAr = bicas.proc.L1L2.sat.get_VSIB(...
        SatSettings, samplesAVolt, ssid, isAchg);

      vsqb = (sum(vsibAr, 'all') / numel(samplesAVolt)) > SatSettings.vsibFractionThreshold;
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
    function vsqbAr = get_snapshot_VSQB_many(...
        SatSettings, zvNValidSamplesPerRecord, zvSamplesAVolt, ssid, isAchgFpa)

      nRecs = irf.assert.sizes(...
        zvNValidSamplesPerRecord, [-1],  ...
        zvSamplesAVolt,           [-1, NaN, 1]);

      vsqbAr = false(nRecs, 1);
      for iRec = 1:nRecs
        vsqbAr(iRec) = bicas.proc.L1L2.sat.get_snapshot_VSQB(...
          SatSettings, ...
          zvSamplesAVolt(iRec, 1:zvNValidSamplesPerRecord(iRec)), ...
          ssid, isAchgFpa);
      end
    end



    % Given ZV-like variables (incl. all ASR ASID channels, CWF/SWF), get
    % saturation bits (Nx1) for quality bitmask.
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
        Bso, tt2000Ar, AsrSamplesAVoltSrm, zvNValidSamplesPerRecord, ...
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

      SatSettings = bicas.proc.L1L2.sat.convert_BSO_to_struct(Bso);



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
        vsibAr = false(nRows, 1);
        for asid = AsrSamplesAVoltSrm.keys'
          asidSamplesAr = bicas.proc.L1L2.sat.set_reconstructed_samples_to_NaN(...
            AsrSamplesAVoltSrm(asid), asid, bltsSsidAr);

          asidVsibAr = bicas.proc.L1L2.sat.get_ASR_CWF_channel_VSIB(...
            SatSettings, bicas.proc.L1L2.const.ASID_to_SSID(asid), ...
            isAchgFpa, asidSamplesAr);

          % L.logf('debug', 'unique(asidVsibAr)         = [%s]', strjoin(string(unique(asidVsibAr)), ','))

          % Merge (OR) bits, record-by-record, over ASIDs.
          vsibAr = any([vsibAr, asidVsibAr], 2);
        end
        % L.logf('debug', 'unique(vsibAr) = [%s]', strjoin(string(unique(vsibAr)), ','))

        vsqbAr = bicas.proc.L1L2.qual.sliding_window_over_fraction(...
          tt2000Ar, vsibAr, ...
          SatSettings.vsibFractionThreshold, ...
          SatSettings.cwfSlidingWindowLengthSec);
      else
        %===========
        % CASE: SWF
        %===========
        vsqbAr = false(nRows, 1);
        for asid = AsrSamplesAVoltSrm.keys'
          asidSamplesAr = bicas.proc.L1L2.sat.set_reconstructed_samples_to_NaN(...
            AsrSamplesAVoltSrm(asid), asid, bltsSsidAr);

          asidIsSaturatedAr = bicas.proc.L1L2.sat.get_ASR_SWF_channel_VSQB(...
            SatSettings, ...
            bicas.proc.L1L2.const.ASID_to_SSID(asid), isAchgFpa, ...
            asidSamplesAr, zvNValidSamplesPerRecord);

          % L.logf('debug', 'unique(asidIsSaturatedAr) = [%s]', strjoin(string(unique(asidIsSaturatedAr)), ','))

          % Merge (OR) bits, record-by-record, over ASIDs.
          vsqbAr = any([vsqbAr, asidIsSaturatedAr], 2);
        end

        % L.logf('debug', 'unique(isSaturatedAr) = [%s]', strjoin(string(unique(isSaturatedAr)), ','))
      end



      if hasSwfFormat
        Tmk.stop_log(nRows, 'CDF record')
      else
        % Log some saturation statistics which may help tell whether how
        % much the saturation varies over time, which may
        % influence/explain if the above processing is slow. Should only
        % be relevant for CWF.
        % NOTE: Only reflects the behaviour of the final VSQB, not the VSIB.
        nSaturationChanges = numel(find(vsqbAr(1:end-1) ~= vsqbAr(2:end)));
        Tmk.stop_log(nRows, 'CDF record', nSaturationChanges, 'sat. flag change')
        L.logf('debug', 'SPEED -- %g [CDF rows/sat. flag change]', nRows/nSaturationChanges)
      end

    end    % function



    % Return VSIB (not VSQB) for one channel of CWF data over which isAchgFpa
    % may vary.
    function vsibAr = get_ASR_CWF_channel_VSIB(...
        SatSettings, ssid, isAchgFpa, samplesAVolt)

      nRows = irf.assert.sizes( ...
        ssid,         [ 1], ...
        isAchgFpa,    [-1], ...
        samplesAVolt, [-1]);

      % NOTE: Splits into subsequences also when ACHG does not matter (DC).
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        isAchgFpa.logical2doubleNan());

      vsibAr = false(nRows, 1);

      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        vsibAr(iRec1:iRec2) = bicas.proc.L1L2.sat.get_VSIB(...
          SatSettings, samplesAVolt(iRec1:iRec2), ssid, isAchgFpa(iRec1));
      end
    end



    % Return VSQB (not VSIB) for one (ASR) channel of SWF data over which
    % isAchgFpa may vary.
    function vsqbAr = get_ASR_SWF_channel_VSQB(...
        SatSettings, ssid, isAchgFpa, samplesAVolt, zvNValidSamplesPerRecord)

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

        vsqbAr(iRec1:iRec2) = bicas.proc.L1L2.sat.get_snapshot_VSQB_many(...
          SatSettings, ...
          zvNValidSamplesPerRecord(iRec1:iRec2), ...
          samplesAVolt(            iRec1:iRec2, :), ...
          ssid, isAchgFpa(iRec1));
      end

    end



    % Function to simplify other function get_VSQB().
    %
    % For a given ASID, determine which rows contain samples which originate
    % from L1R (i.e. which were not reconstructed). Set all other samples to
    % NaN.
    function asidSamplesAr = set_reconstructed_samples_to_NaN(...
        asidSamplesAr, asid, bltsSsidAr)

      ssid = bicas.proc.L1L2.const.ASID_to_SSID(asid);
      %bUse = any(repmat(ssid, [nRows, 1]) == bltsSsidAr, 2);
      bUse = any(ssid == bltsSsidAr, 2);
      asidSamplesAr(~bUse, :) = NaN;
    end



  end    % methods(Static)



end
