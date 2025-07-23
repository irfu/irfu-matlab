%
% Collection of code relating to quality ZVs for L1/L1R to L2 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Determine quality ZVs.
    % Overwrite selected data in selected CDF records with fill values/NaN.
    %
    % NOTE: Does not handle PROCESSING.ZV_QUALITY_FLAG_MAX. That is handled
    %       by bicas.write_dataset_CDF().
    % NOTE: Does not seem able to ever set zv_L2_QUALITY_BITMASK to fill
    %       value.
    %
    % RETURN VALUES
    % =============
    % zvUfv
    %       Logical array. UFV records.
    % QUALITY_FLAG
    %       Cap for output dataset ZV QUALITY_FLAG (cap for input dataset
    %       QUALITY_FLAG).
    % L2_QUALITY_BITMASK
    %       Quality ZV to use.
    %
    % function [zvUfv, QUALITY_FLAG, L2_QUALITY_BITMASK] = ...
    %     get_UFV_quality_ZVs(InZv, isLfr, NsoTable, Bso, L)
    %   % PROPOSAL: Replace InZv-->Separate arguments.
    %
    %   % ASSERTIONS
    %   irf.assert.struct(InZv, {'Epoch', 'bdmFpa', 'autodetectedVsqb'}, {})
    %   irf.assert.sizes( ...
    %     InZv.Epoch,            [-1], ...
    %     InZv.bdmFpa,           [-1], ...
    %     InZv.autodetectedVsqb, [-1]);
    %   assert(isscalar(isLfr) && islogical(isLfr))
    %
    %   Epoch            = InZv.Epoch;
    %   BdmFpa           = InZv.bdmFpa;
    %   autodetectedVsqb = InZv.autodetectedVsqb;
    %   clear InZv
    %
    %   %============================================
    %   % Find CDF records to remove due to settings
    %   %============================================
    %   zvUfv = bicas.proc.L1L2.qual.get_UFV_from_removing_BDMs(...
    %     Epoch, BdmFpa, isLfr, Bso, L);    % CALL TO DELETED FUNCTION
    %
    %   %============================================
    %   % Create quality ZVs based on
    %   % (1) NSO events table, and
    %   % (2) processing-generated QRCs (saturation)
    %   %============================================
    %   [QUALITY_FLAG, L2_QUALITY_BITMASK] = bicas.proc.L1L2.qual.get_quality_ZVs(...
    %     bicas.const.QRCS_L2_MAP, NsoTable, Epoch, autodetectedVsqb, L);
    % end



    % Derive QUALITY_FLAG *cap* and L2_QUALITY_FLAG *bits*. Return values
    % are then supposed to be used for creating global versions of the
    % actual ZVs.
    %
    % More custom arguments for more custom QRCs should be added as needed.
    %
    % ARGUMENTS
    % =========
    % fullSaturationQrcbAr
    %       Autodetected full saturation array (as opposed to derived from NSO
    %       table).
    %
    %
    % RETURN VALUES
    % =============
    % QUALITY_FLAG
    %       Highest allowed value for ZV QUALITY_FLAG.
    %       NOTE: Will never have FPs.
    % L2_QUALITY_BITMASK
    %       Array. L2_QUALITY_BITMASK bits set based on NSOs only. Should be
    %       merged (OR:ed) with pre-existing global L2_QUALITY_BITMASK.
    %
    % function [QUALITY_FLAG, L2_QUALITY_BITMASK] = get_quality_ZVs(...
    %     QrcsL2Map, NsoTable, Epoch, fullSaturationQrcbAr, L)
    %
    %   assert(islogical(fullSaturationQrcbAr))
    %
    %   Qrcbm = bicas.proc.qual.NSO_table_to_QRCB_arrays(...
    %     bicas.const.Q.ALL_QRCID_AR, NsoTable, Epoch, L);
    %
    %   % Remove QRCIDs which this function can not handle (and should not
    %   % need to) since they are not intended for L2_QUALITY_BITMASK.
    %   Qrcbm.remove("BAD_DENSITY");
    %
    %   % Add autodetected saturation.
    %   qrcbAr = Qrcbm("FULL_SATURATION");
    %   qrcbAr = qrcbAr | fullSaturationQrcbAr;
    %   Qrcbm("FULL_SATURATION") = qrcbAr;
    %
    %   % Call generic function for setting QUALITY_FLAG and *_QUALITY_BITMASK.
    %   [QUALITY_FLAG, L2_QUALITY_BITMASK] = ...
    %     bicas.proc.qual.QRCB_arrays_to_quality_ZVs(...
    %     size(Epoch, 1), Qrcbm, QrcsL2Map);
    % end



    % NOTE: Always sets the GLOBAL_SATURATION and CHANNEL_SATURATION QRCIDs but
    % sets the QRCB arrays differently depending on "saturationQualitySchemeId".
    %
    % ARGUMENTS
    % =========
    % VsibZvm
    %       ZVM with VSIBs (not VSTBs, not VSQBs, not QRCBs) for the respective
    %       ASR channels.
    %
    function SaturationQrcbm = get_saturation_QRCBs(...
        tt2000Ar, saturationQualitySchemeId, ...
        VsibZvm, isSwf, vstbFractionThreshold, cwfSlidingWindowLengthSec)

      % PROPOSAL: Change order of arguments.

      SaturationQrcbm = bicas.proc.QrcbMap(numel(tt2000Ar));
      SaturationQrcbm.add_false(bicas.const.Q.SATURATION_QRCID_AR)

      %---------------------------------
      % Obtain CHANNEL_SATURATION QRCBs
      %---------------------------------
      ChannelSaturationQrcbm = ...
        bicas.proc.L1L2.qual.get_QRCBs_channel_saturation(...
        VsibZvm, tt2000Ar, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec);

      switch(saturationQualitySchemeId)
        case "GLOBAL_SATURATION"
        %------------------------------------------------------
        % CHANNEL_SATURATION QRCBs --> GLOBAL_SATURATION QRCBs
        %------------------------------------------------------
        GlobalSaturationQrcbm = ...
          bicas.proc.L1L2.qual.channel_saturation_to_global_saturation_QRCBs(...
          ChannelSaturationQrcbm, numel(tt2000Ar));
          SaturationQrcbm.union(GlobalSaturationQrcbm)

        case "CHANNEL_SATURATION"
          SaturationQrcbm.union(ChannelSaturationQrcbm)

        otherwise
          error("BICAS:ConfigurationBug", ...
            "Illegal argument saturationQualitySchemeId=""%s"".", ...
            saturationQualitySchemeId)
      end

      assert(isequal( ...
        sort(SaturationQrcbm.qrcidAr), ...
        sort(bicas.const.Q.SATURATION_QRCID_AR)))
    end



    % Given VSIBs for 9x ASRs (incl. reconstruction), derive channel saturation
    % QRCBs.
    %
    % NOTE: One can only obtain channel saturation on the 9x ASRs *AFTER*
    %       reconstruction, using the already set VSIBs contained within
    %       the SCHDs.
    %       Ex: V1 and V2 are saturated by being high ==> V12 = 0 (approx.)
    %           ==> V12 does not appear saturated, but the reconstructed V12
    %               value is still affected by the saturation on V1 and V2.
    %
    %
    % IMPLEMENTATION NOTE: Using moving window algo. (CWF) on ASRs instead of
    %                      on BLTSs
    % =======================================================================
    % This function applies the moving window algorithm
    % (bicas.proc.L1L2.qual.sliding_window_over_fraction()) on the ASR (CWF)
    % channels rather than on the BLTSs, since one can then avoid possible bad
    % VSQB behaviour.
    %
    % Ex 1: Consider a reconstructed channel when a VSQB (after moving window)
    % ends on one underlying source channel and begins on another, at roughly
    % the same time: If there is a non-saturated hole (VSQB=0) between the end
    % and beginning, then the reconstructed channel will have a non-saturated
    % hole (VSQB=0) that might never have existed if moving window was applied
    % on the reconstructed channel's threshold saturation bits instead of the
    % two source channels.
    %
    % Ex 2: If the two source channels for a reconstructed channel separately
    % contain too few threshold saturated samples for the moving window algo. to
    % set VSQB=1 (e.g. 30% for a window fraction 50%), but enough when combined
    % (e.g. 30%+30% > 50%), then the reconstructed channel's saturation bits
    % would be zero, despite being very much affected by saturation.
    %
    function Qrcbm = get_QRCBs_channel_saturation(...
        VsibZvm, tt2000Ar, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec)

      assert(isa(VsibZvm, "bicas.utils.ZvMap"))
      assert(isscalar(isSwf) & islogical(isSwf))

      % IMPLEMENTATION NOTE: containers.Map does not support string-valued keys
      % (sic!)
      Qrcbm = bicas.proc.QrcbMap(numel(tt2000Ar));



      % Update Qrcbm wrt. the corresponding arguments.
      function handle_channel(sdidStr, channelSaturationQrcid)
        sdid       = bicas.proc.L1L2.const.C.SDID_DICT(sdidStr);
        vsibSdidAr = VsibZvm.get(sdid);

        if ~isSwf
          vsibSdidAr = bicas.proc.L1L2.qual.sliding_window_over_fraction(...
            tt2000Ar, vsibSdidAr, ...
            vstbFractionThreshold, cwfSlidingWindowLengthSec);
        end

        % IMPLEMENTATION NOTE: The (nested) function is called for the same
        % QRCID up to two times.
        if Qrcbm.has_QRCID(channelSaturationQrcid)
          Qrcbm.set(channelSaturationQrcid, ...
            Qrcbm.get(channelSaturationQrcid) | vsibSdidAr)
        else
          Qrcbm.add(channelSaturationQrcid, vsibSdidAr);
        end
      end



      handle_channel("DC_V1",  "SATURATION_ZV_V1")
      handle_channel("DC_V2",  "SATURATION_ZV_V2")
      handle_channel("DC_V3",  "SATURATION_ZV_V3")
      handle_channel("DC_V12", "SATURATION_ZV_V12")
      handle_channel("DC_V13", "SATURATION_ZV_V13")
      handle_channel("DC_V23", "SATURATION_ZV_V23")
      handle_channel("AC_V12", "SATURATION_ZV_V12")
      handle_channel("AC_V13", "SATURATION_ZV_V13")
      handle_channel("AC_V23", "SATURATION_ZV_V23")
    end



    % Convert channel saturation QRCBs to global saturation QRCBs.
    %
    % NOTE: Always sets QRCID="PARTIAL_SATURATION" QRCB=false since it can not
    %       be autodetected.
    % NOTE: Produces QRCB values using both the old scheme ("GLOBAL_SATURATION")
    %       and new scheme ("CHANNEL_SATURATION") for backward-compatibility
    %       during development and testing. "GLOBAL_SATURATION" should be phased
    %       out eventually. /Erik P G Johansson, 2025-07-08
    %
    % ARGUMENTS
    % =========
    % ChannelSaturationQrcbm
    %       Must contain at least all the CHANNEL_SATURATION QRCIDs, but may
    %       contain more which are then ignored.
    %
    function GlobalSaturationQrcbm = ...
        channel_saturation_to_global_saturation_QRCBs( ...
        ChannelSaturationQrcbm, nRecords)

      fullSaturationQrcbAr = false(nRecords, 1);
      for qrcid = bicas.const.Q.CHANNEL_SATURATION_QRCID_AR'
        fullSaturationQrcbAr = ...
          fullSaturationQrcbAr | ChannelSaturationQrcbm.get(qrcid);
      end

      GlobalSaturationQrcbm = bicas.proc.QrcbMap(nRecords);

      GlobalSaturationQrcbm.add("FULL_SATURATION",    fullSaturationQrcbAr);
      % NOTE: Always false since partial saturation (QRC) can not be
      % autodetected.
      GlobalSaturationQrcbm.add("PARTIAL_SATURATION", false(nRecords, 1));
    end



    % Set samples based on NSO table and SSIDs.
    %
    % ARGUMENTS
    % =========
    % voltageAr
    %       Float. 1D or more dimensions. First dimension is CDF records.
    %       NOTE: Does not have to have any particular unit.
    % ssidAr
    %       Same size as samplesAr. Same number of rows as Qrcbm.
    %
    function voltageAr = set_5xBLTS_voltage_samples_FV(...
        voltageAr, ssidAr, Qrcbm, Qrcsm)

      % IMPLEMENTATION NOTE: Input arrays samplesAr & ssidAr must have same
      % arbitrary size (not arbitrary for first dimension).
      %   PRO: It simplifies testing.
      %   PRO: It simplifies the implementation.
      %   PRO: It makes the function more generic.
      % IMPLEMENTATION NOTE: Argument Qrcsm is there instead of global
      % constant to make testing simpler & more robust.

      % ASSERTIONS
      assert(isfloat(voltageAr))
      assert(isequal(size(voltageAr), size(ssidAr)))
      assert(isa(Qrcsm, "bicas.proc.QrcSettingsMap"))
      assert(Qrcbm.nRecords == size(voltageAr, 1))    % Nbr. of records.

      sizeAr = size(voltageAr);
      bFv    = false(sizeAr);
      for qrcid = Qrcbm.qrcidAr'
        qrcbAr = Qrcbm.get(qrcid);    % (nRecords, 1)
        Qrcs   = Qrcsm.get(qrcid);

        % if isempty(Qrcs)
        %   continue
        % end

        % Arrays of the same size as voltageAr.
        bQrcbAr    = repmat(qrcbAr, [1, sizeAr(2:end)]);
        bSsidMatch = ismember(ssidAr, Qrcs.voltageFvSsidAr);

        bFv        = bFv | (bQrcbAr & bSsidMatch);
      end
      voltageAr(bFv) = NaN;
    end



    % Overwrite records of current with FVs as specified in QRCBs.
    % Cf. bicas.proc.L1L2.qual.set_5xBLTS_voltage_samples_FV().
    %
    % ARGUMENTS
    % =========
    % currentAr
    %       Float. Size (nRecords, 3).
    %
    function currentAr = set_current_samples_FV(currentAr, Qrcbm, Qrcsm)
      % IMPLEMENTATION NOTE: Argument Qrcsm is there (instead of using a
      % global constant) to make test code simpler & more robust.

      assert(isfloat(currentAr))
      assert(isa(Qrcbm, 'bicas.proc.QrcbMap'))
      assert(isa(Qrcsm, 'bicas.proc.QrcSettingsMap'))
      irf.assert.sizes(currentAr, [Qrcbm.nRecords, 3])

      % PROPOSAL: Create and use bAntennas = logical size (1, 3) + repmat.
      %   PRO: Faster?

      iAntennaAr = repmat([1:3], [Qrcbm.nRecords, 1]);
      bFv        = false(size(currentAr));
      for qrcid = Qrcbm.qrcidAr'
        qrcbAr = Qrcbm.get(qrcid);    % (nRecords, 1)
        Qrcs  = Qrcsm.get(qrcid);

        if isempty(Qrcs)
          continue
        end

        % Arrays of the same size as currentAr.
        bQrcbAr       = repmat(qrcbAr, [1, 3]);
        bAntennaMatch = ismember(iAntennaAr, Qrcs.currentFvIantAr);

        bFv           = bFv | (bQrcbAr & bAntennaMatch);
      end
      currentAr(bFv) = NaN;
    end



    % Given a (timestamped) 1D array of flagged samples (bool), label all
    % samples positions which are part of a sliding window (of a specified
    % length) with a fraction of flagged samples which is above a
    % specified threshold.
    %
    % IMPORTANT NOTE: The function's performance (speed) is not great for
    % large datasets (several millions of samples) with mostly saturation
    % (VSIB=1) and length greater than the window. See
    % bicas.proc.L1L2.qual___sliding_window_over_fraction_speedTest.
    % This has been observed to make up almost half of the execution time.
    % Ex: solo_L2_rpw-lfr-surv-cwf-e_20231116_V01.cdf.2024-01-05T22.26.49.log
    % Based on understanding as of 2024-01-08:
    % T ~ (nSamples-windowLengthInSamples) * windowLengthInSamples for
    % nSamples > windowLengthInSamples
    % It is therefore not a crippling problem, but it is a problem.
    %
    % TODO: Above problem should probably be fixed!
    %
    %
    % DETAILS
    % =======
    % * Data gaps effectively count as being filled with samples which are
    %   not flagged.
    % * The algorithm should be fast for data without or few flagged samples.
    % * Individual samples are weighted by their estimated inverse sampling
    %   rate which is the distance to the nearest sample (if there are
    %   >=2 samples). This should make the algorithm handle varying sampling
    %   rate sensibly.
    % * Due to how the algorithm estimates the weight for each sample, if
    %   a window contains samples with a slightly varying sampling rate, the
    %   found fraction of flagged samples within the window will be slightly
    %   lower than one might expect. For that reason, a flagged fraction
    %   threshold of 1 might not trigger (flag) windows within which all
    %   samples are flagged.
    % * Due to how the algorithm estimates the weight for each sample,
    %   samples with identical timestamps count as having weight zero.
    % * Timestamps must increase (but not strictly increase).
    %
    %
    % ARGUMENTS
    % =========
    % tt2000Ar
    %       Column array of TT2000 values.
    % bFlag1Ar
    %       Column array of logical. Samples which are "flagged", e.g. for
    %       saturation.
    % minFlaggedFraction
    %       Minimum fraction of (weighted) flagged samples within a window,
    %       for all samples within the window to be flagged.
    % windowLengthSec
    %       Length of sliding window. Individual samples count as having a
    %       length equal to the inverse sampling frequency (simplification).
    %
    %
    % RETURN VALUE
    % ============
    % bFlag2Ar
    %       Column array of logical. Modified version of bFlag1Ar such that,
    %       every sliding window (of length intervalLengthSec) contains at
    %       least a fraction minFlaggedFraction of weighted flagged samples.
    %
    function bFlag2Ar = sliding_window_over_fraction(...
        tt2000Ar, bFlag1Ar, minFlaggedFraction, windowLengthSec)
      % PROPOSAL: Better name
      %   "sliding_window_exceding_fraction"
      %   "sliding_window"
      %   moving window
      %   interval over fraction
      %   density (of set bits/flags)
      %   bit, flag
      %
      % PROPOSAL: Move to bicas.utils.
      %   PRO: More generic than quality variables.
      %   PRO: Independent of L1/L1R-L2 processing in principle.

      %============
      % ASSERTIONS
      %============
      % Sizes:
      irf.assert.sizes(...
        tt2000Ar, [-1], ...
        bFlag1Ar, [-1] ...
        );
      assert(isscalar(minFlaggedFraction))
      assert(isscalar(windowLengthSec))
      % Types/classes:
      assert(isa(tt2000Ar, 'int64'))
      assert(islogical(bFlag1Ar))
      assert(isfloat(minFlaggedFraction))
      assert(isfloat(windowLengthSec))
      % NOTE: Algorithm requires that timestamps increase (nut not
      %       strictly increase).
      assert(issorted(tt2000Ar, 'ascend'))
      assert((0 <= minFlaggedFraction) && (minFlaggedFraction <= 1), ...
        'flagFractionThreshold = %d is not a legal value.', minFlaggedFraction)
      assert(windowLengthSec >= 0)

      %===========================
      % ALGORITHM / SPECIAL CASES
      %===========================
      if all(~bFlag1Ar)
        % CASE: (1) All samples are false, or
        %       (2) there are zero samples.
        bFlag2Ar = false(size(bFlag1Ar));

      elseif isscalar(bFlag1Ar)
        % CASE: There is exactly one sample.

        % NOTE: Algorithm can not handle this case since STL becomes
        % infinite. Therefore special case.
        bFlag2Ar = bFlag1Ar;

      else
        % CASE: (1) There is at least one flagged sample, and
        %       (2) There are at least two samples.

        timeSecAr = double(tt2000Ar) / 1e9;

        bFlag2ForwardAr = bicas.proc.L1L2.qual.sliding_window_over_fraction_forward_pass(...
          timeSecAr, bFlag1Ar, ...
          minFlaggedFraction, windowLengthSec);

        % NOTE: Same call as above, except that (1) reversing the order of
        % timestamps and samples, and (2) negating the timestamps (so
        % that they increment despite their order being reversed).
        bFlag2BackwardAr = bicas.proc.L1L2.qual.sliding_window_over_fraction_forward_pass(...
          -timeSecAr(end:-1:1), bFlag1Ar(end:-1:1), ...
          minFlaggedFraction, windowLengthSec);

        bFlag2BackwardAr = bFlag2BackwardAr(end:-1:1);

        bFlag2Ar = bFlag2ForwardAr | bFlag2BackwardAr;
      end
    end    % function



    % Effectively internal function to
    % bicas.proc.L1L2.qual.sliding_window_over_fraction() to simplify its
    % implementation. Runs one "pass" in the forward direction.
    %
    function bFlag2Ar = sliding_window_over_fraction_forward_pass(...
        timeSecAr, bFlag1Ar, minFlaggedFraction, windowLengthSec)
      % PROPOSAL: Better name
      %   algorithm
      %   pass
      %
      % PROPOSAL: Use smallest window length that is equal to or greater than
      %           the specified one (instead of the largest window length
      %           that is equal to or less than the specified one).
      %   CON: If there is a data gap, then the difference can be very
      %        large, and the window could become too much large.
      % PROPOSAL: Always use argument for window length when calculating
      %           fraction.
      %   PRO: Prevents window from becoming too small before a data gap
      %        that is longer than the argument window length.

      % Naming conventions
      % ==================
      % STL  = Sample Time Length. Length of time assigned to each sample.
      %        Equal to twice the longest distance to the nearest sample.
      %        Intended for (1) weighing samples with different sampling
      %        rate, and (2) for including half in the window length for
      %        samples at the beginning and end of window.
      % STLW = STL-Weighted

      DEBUG_ENABLED = 0;

      % DEBUG
      if DEBUG_ENABLED
        fprintf('--------sliding_window_over_fraction_forward_pass\n')
      end

      n = numel(bFlag1Ar);
      assert(n >= 2)

      % Pre-allocate
      bFlag2Ar  = false(size(bFlag1Ar));

      diffSecAr = [Inf; diff(timeSecAr); Inf];
      % NOTE: Returns Inf for array length == 1 which must therefore be
      %       avoided.
      stlSecAr  = min([diffSecAr(1:end-1), diffSecAr(2:end)], [], 2);

      % Modified cumulative sum so that a difference between indices i and
      % i+1 represents the STL of sample i.
      cumulStlwFlagAr = [0; cumsum(bFlag1Ar .* stlSecAr)];

      % ==================================================================
      % Iterate over flagged samples and use these to find "windows"
      % (approximately fixed-length time intervals) which begin at those
      % samples
      % ==================================================================
      iFlagSet1Ar        = find(bFlag1Ar);
      i1                 = iFlagSet1Ar(1) - 1;
      prevSetWindowFlags = false;
      % Iterate over starting indices: i0
      for i0 = iFlagSet1Ar'
        % CASE: i0 = Index to a flagged sample.

        % =============================================================
        % Obtain window that begins with i0 (already set) and ends with
        % i1 (to be determined)
        % =============================================================
        while true
          % If no more sample can be added to the window, then keep
          % the window size as it is.
          if i1+1 > n
            break
          end
          % CASE: i1+1 <= n (i.e. one can safely use i1+1 as an index)

          % If a one sample larger window is too large, then keep the
          % current window size.
          % PROPOSAL: Derive arrays of time of beginnings and end of
          %           every sample and use that instead.
          edgesStlSec              = stlSecAr(i0)/2 + stlSecAr(i1+1)/2;
          candidateWindowLengthSec = timeSecAr(i1+1) - timeSecAr(i0) + edgesStlSec;
          if candidateWindowLengthSec > windowLengthSec
            break
          end
          % CASE: A one sample larger window is not too large.

          i1 = i1 + 1;
        end
        % CASE: i1 is the highest value for which
        %       (1) i0 <= i1 <= n, AND
        %       (2) cumulTimeSecAr(i1) < cumulTimeSecAr+intervalLengthSec.
        %       i0:i1 = Range of indices which define the window.



        % edgesStlSec       = stlSecAr(i0)/2 + stlSecAr(i1)/2;
        windowStlwFlagSec = cumulStlwFlagAr(i1+1) - cumulStlwFlagAr(i0);
        % IMPLEMENTATION NOTE: Using the argument window length rather
        % than window length calculated from the sample/index range
        % prevents the window from becoming too small (1) before a data
        % gap that is longer than the argument window length, and (2)
        % before the end of samples.
        fractionStlwFlag  = windowStlwFlagSec / windowLengthSec;

        % IMPLEMENTATION NOTE: Threshold should count as lower value
        % (equality) so that one can require all elements to be
        % flagged by setting minFlaggedFraction=1.
        setWindowFlags = (fractionStlwFlag >= minFlaggedFraction);
        if setWindowFlags
          % IMPLEMENTATION NOTE: Setting bFlag2Ar(...) = true slows
          % down the entire algorithm. Therefore sets as few indices
          % as possible by avoiding indices set in the previous window
          % (if set). ==> Speeds up algorithm.
          if prevSetWindowFlags
            bFlag2Ar(max(i0, prevI1+1):i1) = true;
          else
            bFlag2Ar(i0:i1)                = true;
          end
        end
        prevSetWindowFlags = setWindowFlags;
        %prevSetWindowFlags = false;    % TEST. Disable speedup.
        prevI1             = i1;

        if DEBUG_ENABLED
          fprintf('Found window i0:i1 = %i:%i\n', i0, i1)
          fprintf('    prevSetWindowFlags  = %g\n', prevSetWindowFlags)
          fprintf('    prevI1              = %g\n', prevI1)
          fprintf('    timeSecAr([i0, i1]) = %g - %g\n', timeSecAr(i0), timeSecAr(i1))
          % fprintf('    edgesStlSec         = %g\n', edgesStlSec)
          fprintf('    windowLengthSec     = %g\n', windowLengthSec)
          fprintf('    windowStlwFlagSec   = %g\n', windowStlwFlagSec)
          fprintf('    fractionStlwFlag    = %g\n', fractionStlwFlag)
          fprintf('    ==> setWindowFlags = %d\n', setWindowFlags)
        end

        % If future windows can not be larger due to lack of
        % samples, then exit function.
        % IMPLEMENTATION NOTE: This prevents the algorithm from
        % evaluating (calculating fractions for) unnecessarily small
        % windows at the high timestamps end.
        if i1+1 > n
          return
        end
      end    % for

    end    % function



  end    % methods(Static)



end
