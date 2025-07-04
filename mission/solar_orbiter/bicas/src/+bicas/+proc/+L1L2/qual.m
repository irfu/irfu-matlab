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
    %     Epoch, BdmFpa, isLfr, Bso, L);
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
    % isFullSaturation
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
    %     QrcsL2Map, NsoTable, Epoch, isFullSaturation, L)
    %
    %   assert(islogical(isFullSaturation))
    %
    %   QrcFlagsMap = bicas.proc.qual.NSO_table_to_QRC_flag_arrays(...
    %     bicas.const.ALL_QRCID_AR, NsoTable, Epoch, L);
    %
    %   % Remove QRCIDs which this function can not handle (and should not
    %   % need to) since they are not intended for L2_QUALITY_BITMASK.
    %   QrcFlagsMap.remove("BAD_DENSITY");
    %
    %   % Add autodetected saturation.
    %   b = QrcFlagsMap("FULL_SATURATION");
    %   b = b | isFullSaturation;
    %   QrcFlagsMap("FULL_SATURATION") = b;
    %
    %   % Call generic function for setting QUALITY_FLAG and *_QUALITY_BITMASK.
    %   [QUALITY_FLAG, L2_QUALITY_BITMASK] = ...
    %     bicas.proc.qual.QRC_flag_arrays_to_quality_ZVs(...
    %     size(Epoch, 1), QrcFlagsMap, QrcsL2Map);
    % end



    % Given 9x ASR VSIBs (channel-specific VSIBs), derive quality variables
    % wrt. channel saturation.
    %
    % NOTE: Can only be applied to 9x ASRs *AFTER* reconstruction, using
    %       the already set VSIBs contained within the SCHDs.
    % Ex: V1 and V2 are saturated by being high ==> V12 = 0 (approx.)
    %     ==> V12 does not appear saturated.
    %
    % NOTE: Can produce quality variable values using either old scheme
    % ("GLOBAL_SATURATION") and new scheme ("CHANNEL_SATURATION") for
    % backward-compatibility during development and testing. "GLOBAL_SATURATION"
    % should be phased out eventually. /Erik P G Johansson, 2025-07-2
    %
    %
    % ARGUMENTS
    % =========
    % saturationQualitySchemeId
    %       String constant. Which scheme to use for setting quality zVariables
    %       (which QRCIDs are actually used).
    %
    %
    % IMPLEMENTATION NOTE: Channels to use for QRCS
    % =============================================
    % Uses ASR channel (SDID) VSIBs to determine QRCSs for channels in this
    % function, despite that there are no SSIDs and that the ASR channels really
    % correspond to ZVs, which may be non-ASR samples in rare cases, e.g. GND
    % (BSM=6/7). This is necessary since reconstructed samples never have SSIDs.
    % Therefore, this operation can not be performed on the BLTSs.
    %
    %
    % IMPLEMENTATION NOTE: Using moving window algo. on ASRs instead of BLTSs
    % =======================================================================
    % This function applies the moving window algorithm
    % (bicas.proc.L1L2.qual.sliding_window_over_fraction()) on the ASR (CWF)
    % channels rather than on the BLTSs, since one can then avoid possible bad
    % VSQB behaviour.
    % --
    % Ex 1: Consider a reconstructed channel when a VSQB (after moving window)
    % ends on one underlying source channel and begins on another, at roughly
    % the same time: If there is a non-saturated hole between the end and
    % beginning, then the reconstructed channel may have a non-saturated hole
    % that would never exist if moving window was applied on the reconstructed
    % channel's threshold saturation bits instead of the two source channels.
    % --
    % Ex 2: If the two source channels for a reconstructed channel separately
    % contain too few threshold saturated samples for the moving window algo.
    % (e.g. 30% for a window fraction 50%), but enough when combined (e.g.
    % 30%+30%>50%), then the reconstructed channel's saturation bits would be
    % zero, despite being very much affected by saturation.
    %
    function [QUALITY_FLAG, L2_QUALITY_BITMASK] = get_quality_ZVs_channel_saturation(...
        Cdac, tt2000Ar, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec, saturationQualitySchemeId)
      % PROPOSAL: Separate function for one ASR/channel at a time.
      % PROPOSAL: Somehow only submit VSIBs, not SCHDs with channel samples.
      %   PRO: Easier to have separate functions for applying moving window
      %        algo. to VSIBs (for all channels) in a separate function (not
      %        called from here).
      % PROPOSAL: Reorg. to
      %   (1) function for converting samples to channel saturation QRC flags,
      %   (2) function for converting channel saturation QRC flags to
      %       global saturation QRC flags
      %   (3) function for selecting channel/global saturation QRC flags
      %       depending on saturation scheme,
      %   (4) using generic function for QRC flags-->quality variables,
      %       bicas.proc.qual.QRC_flag_arrays_to_quality_ZVs().

      assert(isa(Cdac, "bicas.proc.L1L2.ChannelDataAsrCollection"))
      assert(isscalar(isSwf) & islogical(isSwf))
      assert(isscalar(saturationQualitySchemeId) & isstring(saturationQualitySchemeId))

      nRec = Cdac.nRecords;

      % Default, "empty", quality ZVs which will be modified for different
      % combinations of SDID and QRCID.
      QUALITY_FLAG       = repmat(bicas.const.QUALITY_FLAG_MAX, nRec, 1);
      L2_QUALITY_BITMASK = zeros(nRec, 1, 'uint16');

      function handle_channel(sdidStr, channelSaturationQrcid)
        switch(saturationQualitySchemeId)
          case "GLOBAL_SATURATION"
            qrcid = 'FULL_SATURATION';
          case "CHANNEL_SATURATION"
            qrcid = channelSaturationQrcid;
          otherwise
            error("BICAS:ConfigurationBug", ...
              "Illegal argument saturationQualitySchemeId=""%s"".", ...
              saturationQualitySchemeId)
        end

        Qrcs = bicas.const.QRCS_L2_MAP(qrcid);
        sdid = bicas.proc.L1L2.const.C.SDID_DICT(sdidStr);
        Schd = Cdac.get_channel(sdid);

        if isSwf
          vsqbSdidAr = Schd.vsibAr;
        else
          vsqbSdidAr = bicas.proc.L1L2.qual.sliding_window_over_fraction(...
            tt2000Ar, Schd.vsibAr, ...
            vstbFractionThreshold, cwfSlidingWindowLengthSec);
        end

        QUALITY_FLAG_sdid = ...
            uint8( vsqbSdidAr) * Qrcs.QUALITY_FLAG ...
          + uint8(~vsqbSdidAr) * bicas.const.QUALITY_FLAG_MAX;
        L2_QUALITY_BITMASK_sdid = uint16(vsqbSdidAr) * Qrcs.Lx_QUALITY_BITMASK;

        QUALITY_FLAG       = min(QUALITY_FLAG, QUALITY_FLAG_sdid);
        L2_QUALITY_BITMASK = bitor(L2_QUALITY_BITMASK, L2_QUALITY_BITMASK_sdid);
      end

      handle_channel('DC_V1',  'SATURATION_ZV_V1')
      handle_channel('DC_V2',  'SATURATION_ZV_V2')
      handle_channel('DC_V3',  'SATURATION_ZV_V3')
      handle_channel('DC_V12', 'SATURATION_ZV_V12')
      handle_channel('DC_V13', 'SATURATION_ZV_V13')
      handle_channel('DC_V23', 'SATURATION_ZV_V23')
      handle_channel('AC_V12', 'SATURATION_ZV_V12')
      handle_channel('AC_V13', 'SATURATION_ZV_V13')
      handle_channel('AC_V23', 'SATURATION_ZV_V23')
    end



    % Overwrite records of voltage & current with FVs as specified in arbitrary
    % array.
    %
    % ARGUMENTS
    % =========
    % zv_Epoch
    %       NOTE: Only needed for logging.
    % zvAsrSamplesAVoltSrm
    %       ASR samples.
    %       NOTE: Handle object which is MODIFIED in-place.
    %
    function zvCurrentAAmpere = set_voltage_current_FV(...
        zv_Epoch, zvAsrSamplesAVoltSrm, zvCurrentAAmpere, zvUfv, L)
      % PROPOSAL: Separate functions for ASR samples and bias currents.

      assert(islogical(zvUfv))
      assert(isa(zvAsrSamplesAVoltSrm, 'bicas.utils.SameRowsMap'))
      assert(isa(zvAsrSamplesAVoltSrm, 'handle'))

      % Log
      logHeaderStr = sprintf(...
        ['Interval(s) of CDF records for which data should be set', ...
        ' to fill values (i.e. removed), regardless of reason.\n']);
      bicas.proc.L1L2.qual.log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)

      % ====================================
      % Set CURRENT values to fill value/NaN
      % ====================================
      zvCurrentAAmpere(zvUfv, :) = NaN;

      % ====================================
      % Set VOLTAGE values to fill value/NaN
      % ====================================
      % NOTE: Should really use future bicas.utils.SameSizeTypeMap here
      %       (if implemented) which contains size on other dimensions.
      keyArray = zvAsrSamplesAVoltSrm.keys;
      nSpr     = size(zvAsrSamplesAVoltSrm(keyArray(1)), 2);

      % IMPLEMENTATION NOTE: bicas.utils.SameRowsMap.set_rows() can not
      % handle logical indexing.
      iUfv = find(zvUfv);
      nanArray = NaN(size(iUfv, 1), nSpr);
      tempSrm = bicas.utils.SameRowsMap(...
        "uint8", size(nanArray, 1), ...
        'CONSTANT', nanArray, zvAsrSamplesAVoltSrm.keys);
      zvAsrSamplesAVoltSrm.set_rows(tempSrm, iUfv);
    end



    % Find CDF records to remove (set to fill value) based on settings (not
    % data itself, almost, since BDM is data).
    %
    % Ex: Sweeps
    %
    % NOTE: It is not obvious that data should be set to FV instead of having
    % quality bitmask/flag modified. Nonetheless, I think setting data to fill
    % value was requested by YK many years ago.
    % /Erik P G Johansson 2023-11-28
    %
    % NOTE: This function is a historical remnant from old functionality for
    % removing sweeps using BDM=4 when BDM=0 was the nominal BDM. This
    % functionality has not been used for a long time and could potentially be
    % removed. Setting PROCESSING.L2.REMOVE_DATA.MUX_MODES would be removed with
    % it.
    % /Erik P G Johansson 2025-01-16
    %
    %
    % ARGUMENTS
    % =========
    % zvBdmFpa
    %       Demultiplexer data, from BIAS HK or LFR.
    %       Fill positions are not matched agains BDMs stored in settings.
    %
    function zvUfv = get_UFV_from_removing_BDMs(...
        zv_Epoch, zvBdmFpa, isLfr, Bso, L)
      % PROPOSAL: Separate function for logging which records that should be removed.
      % PROPOSAL: Arguments for settings.
      %   CON: Logs the settings keys.
      %   CON: Settings used depends on argument isLfr.
      % PROPOSAL: Create QRCID and QRC flags and use generic functionality for
      %   setting data to FV (functionality which does not yet exist).

      bicas.utils.assert_ZV_Epoch(zv_Epoch)
      assert(islogical(isLfr));
      assert(isa(zvBdmFpa, 'bicas.utils.FPArray'))

      %===============
      % Read settings
      %===============
      [bdmRemoveArray, settingBdmRemoveKey] = Bso.get_fv(...
        'PROCESSING.L2.REMOVE_DATA.MUX_MODES');
      bdmRemoveArray = bdmRemoveArray(:);
      if     isLfr,   settingMarginKey = 'PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_SEC';    % LFR
      else,           settingMarginKey = 'PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_SEC';    % TDS
      end
      [removeMarginSec, settingMarginKey] = Bso.get_fv(settingMarginKey);

      %==========================================
      % Find exact indices/CDF records to remove
      %==========================================
      % NOTE: ismember(NaN, nan) == false
      zvUfv = irf.utils.true_with_margin(...
        zv_Epoch, ...
        ismember(zvBdmFpa.int2doubleNan(), bdmRemoveArray), ...
        removeMarginSec * 1e9, removeMarginSec * 1e9);

      %=====
      % Log
      %=====
      logHeaderStr = sprintf(...
        ['Found interval(s) of CDF records for which data should be set to', ...
        ' fill values (i.e. removed) based on settings.\n', ...
        '    NOTE: This may not be all CDF records which will be removed.\n', ...
        '    Setting %s = [%s]\n', ...
        '    Setting %s = %f\n'], ...
        settingBdmRemoveKey, ...
        strjoin(irf.str.sprintf_many('%g', bdmRemoveArray), ', '), ...
        settingMarginKey, ...
        removeMarginSec);
      bicas.proc.L1L2.qual.log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)
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



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Log UFV records
    %
    % NOTE: Only logs (including header) if there are records to remove.
    %
    function log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)
      % PROPOSAL: Redefine, rework to function that can be used for
      % logging separate UFVs obtained in different ways.
      %   Ex: UFVs due to excluding BDMs.
      %   Ex: UFVs due to automatically detected sweeps.
      %   Ex: UFVs due to detected sweeps via QUALITY_BITMASK (future).

      LL = 'info';    % LL = Log Level

      [i1Array, i2Array] = irf.utils.split_by_false(zvUfv);
      nUfvIntervals = numel(i1Array);
      if nUfvIntervals > 0

        %==============
        % Log settings
        %==============
        L.logf(LL, logHeaderStr)

        %===============
        % Log intervals
        %===============
        for iRi = 1:nUfvIntervals
          iCdfRecord1 = i1Array(iRi);
          iCdfRecord2 = i2Array(iRi);
          utcStr1 = bicas.utils.TT2000_to_UTC_str(zv_Epoch(iCdfRecord1), 9);
          utcStr2 = bicas.utils.TT2000_to_UTC_str(zv_Epoch(iCdfRecord2), 9);
          L.logf(LL, '    Records %8i-%8i, %s -- %s', ...
            iCdfRecord1, iCdfRecord2, utcStr1, utcStr2);
        end
      end

    end



  end    % methods(Static, Access=private)



end
