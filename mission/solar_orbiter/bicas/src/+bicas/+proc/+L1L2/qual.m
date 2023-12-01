%
% Collection of code relating to quality variables for L1/L1R to L2 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual
    % PROPOSAL: Redefine as applying to quality variables for all processing
    %           (between any archiving levels).



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Set quality zVariables.
        % Overwrite selected data in selected CDF records with fill values/NaN.
        %
        % NOTE: Does not handle PROCESSING.ZV_QUALITY_FLAG_MAX. That is handled
        %       by bicas.write_dataset_CDF().
        % NOTE: Does not seem able to ever set zv_L2_QUALITY_BITMASK to fill
        %       value.
        function [zvUfv, QUALITY_FLAG_Fpa, L2_QUALITY_BITMASK] = ...
                modify_quality_filter(InZv, isLfr, NsoTable, SETTINGS, L)

            irf.assert.struct(InZv, {'Epoch', 'ufv', 'bdmFpa', 'QUALITY_FLAG_Fpa'}, {})
            Epoch            = InZv.Epoch;
            zvUfv            = InZv.ufv;
            zvBdmFpa         = InZv.bdmFpa;
            QUALITY_FLAG_Fpa = InZv.QUALITY_FLAG_Fpa;
            clear ZvIn

            % ASSERTIONS
            assert(isscalar(isLfr) && islogical(isLfr))
            assert(isa(QUALITY_FLAG_Fpa, 'bicas.utils.FPArray') && strcmp(QUALITY_FLAG_Fpa.mc, 'uint8'))
            nRecords = irf.assert.sizes( ...
                Epoch,            [-1], ...
                zvBdmFpa,         [-1], ...
                QUALITY_FLAG_Fpa, [-1], ...
                zvUfv,            [-1]);

            % Pre-allocate
            L2_QUALITY_BITMASK = zeros(nRecords, 1, 'uint16');


            
            %============================================================
            % Find CDF records to remove due to settings
            %============================================================
            zvUfvSettings = bicas.proc.L1L2.qual.get_UFV_records_from_settings(...
                Epoch, zvBdmFpa, isLfr, SETTINGS, L);

            zvUfv = zvUfv | zvUfvSettings;

            %==============================================
            % Modify quality ZVs based on NSO events table
            %==============================================
            [QUALITY_FLAG_NsoFpa, L2_QUALITY_BITMASK_nso] = bicas.proc.L1L2.qual.get_quality_by_NSOs(...
                bicas.const.NSOID_SETTINGS, NsoTable, Epoch, L);
            QUALITY_FLAG_Fpa   = QUALITY_FLAG_Fpa.min(QUALITY_FLAG_NsoFpa);
            L2_QUALITY_BITMASK = bitor(L2_QUALITY_BITMASK, L2_QUALITY_BITMASK_nso);



            assert(isa(L2_QUALITY_BITMASK, 'uint16'))
        end    % modify_quality_filter
        
        
        
        % Derive QUALITY_FLAG *cap* and L2_QUALITY_FLAG *bits* set due to NSO
        % table by itself. Return values are then supposed to be used for
        % creating global versions of the actual ZVs.
        %
        %
        % RETURN VALUES
        % =============
        % QUALITY_FLAG_CapFpa
        %       FPA. Cap (highest allowed value) for ZV QUALITY_FLAG.
        %       NOTE: Will never have FPs.
        % L2_QUALITY_BITMASK
        %       Array. L2_QUALITY_BITMASK bits set based on NSOs only. Should be
        %       merged (OR:ed) with global L2_QUALITY_BITMASK.
        %
        function [QUALITY_FLAG_Fpa, L2_QUALITY_BITMASK] = ...
                get_quality_by_NSOs(NsoidSettingsMap, NsoTable, Epoch, L)

            % Variable naming conventions:
            % ----------------------------
            % GE = Global Event = NSO event in global NSO event table.
            % CE = CDF Event    = NSO event that overlaps with CDF records.
            % NA                = Numeric Array

            % NOTE: iCdfEventNa = CDF events as indices to global events.
            [bCeRecordsCa, ceNsoidCa, iCeNa] = NsoTable.get_NSO_timestamps(Epoch);
            nCe = numel(ceNsoidCa);
            nGe = numel(NsoTable.evtNsoidCa);
            L.logf('info', ...
                ['Searched non-standard operations (NSO) table.', ...
                ' Found %i relevant NSO events out of a total of %i NSO events.'], ...
                nCe, nGe);

            % Pre-allocate
            % NOTE: QUALITY_FLAG is set to max value.
            QUALITY_FLAG_Fpa = bicas.utils.FPArray(...
                bicas.const.QUALITY_FLAG_MAX * ones(size(Epoch), 'uint8'));
            L2_QUALITY_BITMASK = zeros(size(Epoch), 'uint16');
            
            

            % Iterate over index into LOCAL/CDF NSO events table.
            for kCe = 1:nCe

                % Index into GLOBAL NSO events table.
                iGe = iCeNa(kCe);
                eventNsoid = ceNsoidCa{kCe};
                % Indices into ZVs.
                bCeRecords = bCeRecordsCa{kCe};

                %===========================================================
                % Log the relevant NSO event in the GLOBAL NSO events table
                %===========================================================
                L.logf('info', '    %s -- %s %s', ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGe)), ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array( iGe)), ...
                    eventNsoid);

                %================================
                % Take action depending on NSOID
                %================================
                if NsoidSettingsMap.isKey(eventNsoid)
                    Setting = NsoidSettingsMap(eventNsoid);
                    
                    % NOTE: Variables contain SCALAR values (i.e. they are not
                    %       arrays).
                    QUALITY_FLAG_nsoid       = Setting.QUALITY_FLAG;
                    L2_QUALITY_BITMASK_nsoid = Setting.L2_QUALITY_BITMASK;
                else
                    % ASSERTION
                    % NOTE: Not perfect assertion on legal NSOIDs since code
                    % only checks those relevant for the data (time interval)
                    % currently processed. (Therefore also checks all NSOIDs
                    % when reads NSO table.)
                    error('Can not interpret RCS NSOID "%s".', eventNsoid)
                end

                QUALITY_FLAG_CeFpa              = QUALITY_FLAG_Fpa(bCeRecords, 1);
                QUALITY_FLAG_Fpa(bCeRecords, 1) = QUALITY_FLAG_CeFpa.min(...
                    bicas.utils.FPArray(QUALITY_FLAG_nsoid));
            
                L2_QUALITY_BITMASK( bCeRecords, 1) = bitor(...
                    L2_QUALITY_BITMASK(bCeRecords, 1), ...
                    L2_QUALITY_BITMASK_nsoid);
                
            end    % for

        end



        % Overwrite selected records of voltage & current with FVs.
        %
        % ARGUMENTS
        % =========
        % zv_Epoch
        %       NOTE: Only needed for logging.
        % zvAsrSamplesAVoltSrm
        %       ASR samples.
        %       NOTE: Handle object which is MODIFIED.
        function zvCurrentAAmpere = set_voltage_current_FV(...
                zv_Epoch, zvAsrSamplesAVoltSrm, zvCurrentAAmpere, zvUfv, L)
            % PROPOSAL: Separate functions for ASR samples and bias currents.
            
            assert(islogical(zvUfv))
            assert(isa(zvAsrSamplesAVoltSrm, 'bicas.utils.SameRowsMap'))

            % Log
            logHeaderStr = sprintf(...
                ['Interval(s) of CDF records for which data should be set', ...
                ' to fill values (i.e. removed), regardless of reason.\n']);
            bicas.proc.L1L2.qual.log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)

            % Set current values to fill value/NaN.
            zvCurrentAAmpere(zvUfv, :) = NaN;

            % ====================================
            % Set voltage values to fill value/NaN
            % ====================================
            % NOTE: Should really use future bicas.utils.SameSizeTypeMap here
            %       which contains size on other dimensions.
            keysCa = zvAsrSamplesAVoltSrm.keys;
            nSpr   = size(zvAsrSamplesAVoltSrm(keysCa{1}), 2);

            % IMPLEMENTATION NOTE: bicas.utils.SameRowsMap.setRows() can not
            % handle logical indexing.
            iUfv = find(zvUfv);
            nanArray = NaN(size(iUfv, 1), nSpr);
            tempSrm = bicas.utils.SameRowsMap(...
                'char', size(nanArray, 1), ...
                'CONSTANT', nanArray, zvAsrSamplesAVoltSrm.keys);
            zvAsrSamplesAVoltSrm.setRows(tempSrm, iUfv);
        end



        % Find CDF records to remove (set to fill value) based on settings (not
        % data itself, almost, since MUX mode is data).
        %
        % Ex: Sweeps
        %
        % NOTE: It is not obvious that data should be set to FV instead of
        % having quality bitmask/flag modified. Nonetheless, I think setting
        % data to fill value was requested by YK many years ago. /Erik P G
        % Johansson 2023-11-28
        %
        %
        % ARGUMENTS
        % =========
        % zvBdmFpa
        %       Demultiplexer data, from BIAS HK or LFR.
        %       Fill positions are not matched agains BDMs stored in settings.
        %
        function zvUfv = get_UFV_records_from_settings(...
                zv_Epoch, zvBdmFpa, isLfr, SETTINGS, L)
            % PROPOSAL: Separate function for logging which records that should be removed.
            % PROPOSAL: Arguments for settings.
            %   CON: Logs the settings keys.
            %   CON: Settings used depends on argument isLfr.

            bicas.utils.assert_ZV_Epoch(zv_Epoch)
            assert(islogical(isLfr));
            assert(isa(zvBdmFpa, 'bicas.utils.FPArray'))

            %===============
            % Read settings
            %===============
            [bdmRemoveArray, settingBdmRemoveKey] = SETTINGS.get_fv(...
                'PROCESSING.L2.REMOVE_DATA.MUX_MODES');
            bdmRemoveArray = bdmRemoveArray(:);
            if     isLfr   settingMarginKey = 'PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S';    % LFR
            else           settingMarginKey = 'PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S';    % TDS
            end
            [removeMarginSec, settingMarginKey] = SETTINGS.get_fv(settingMarginKey);

            %==========================================
            % Find exact indices/CDF records to remove
            %==========================================
            % NOTE: ismember(NaN, nan) == false
            zvUfv = irf.utils.true_with_margin(...
                zv_Epoch, ...
                ismember(zvBdmFpa.int2doubleNan(), bdmRemoveArray), ...
                removeMarginSec * 1e9);

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



        % NOTE: 2023-11-30: CURRENTLY UNUSED, BUT IS PLANNED TO BE USED.
        %
        function tsfAr = determine_TSF(...
                samplesAVolt, Ssid, hasSnapshotFormat, biasHighGain)
            % PROPOSAL: Better name.
            %   ~sample-to-TSF
            %       PRO: Can use same scheme for TSF-to-SWSF function.
            % PROPOSAL: Retrieve thresholds from settings.
            % TODO-DEC: How handle NaN/missing samples?
            %   Ex: Channels set to NaN when there is no BLTS due to LRX setting.
            %
            % TODO-DEC: How handle snapshots? If sliding window is not applied
            %           here, then when should the decision how to set SWSF for
            %           snapshots be done?
            %   PROPOSAL: Return TSFs for every snapshot sample for sequence of
            %             samples, regardless of whether snapshots or not. Other
            %             function handles snapshots specifically.
            %   PROPOSAL: Separate function for snapshots and which returns
            %             value with other meaning (not apply sliding window
            %             to). That function uses this function in its
            %             implementation though.
            
            assert(isfloat(samplesAVolt))
            assert(isa(Ssid, 'bicas.proc.L1L2.SignalSourceId'))
            
            if ~Ssid.is_ASR()
                tsfAr = false(size(samplesAVolt));
                return
            end
            
            % CASE: ASR (no non-plasma signal, no special case)
            
            % ====================
            % Determine thresholds
            % ====================
            if Ssid.Asid.is_diff()
                if Ssid.Asid.is_AC()
                    % CASE: AC diff
                    % -------------
                    assert(isfinite(biasHighGain))   % Not NaN.
                    % TODO-DEC: Too high thresholds?
                    if biasHighGain
                        higherThresholdAVolt = 0.3/20;
                    else
                        higherThresholdAVolt = 0.3;
                    end
                else
                    % CASE: DC diff
                    % -------------
                    higherThresholdAVolt = 2;
                end
            else
                % CASE: Single
                % ------------
                % NOTE: Not using terms "min" and "max" since they are
                % ambiguous (?).
                higherThresholdAVolt =  40;
            end
            lowerThreshold = -higherThresholdAVolt;
            
            % NOTE: Has to be able handle NaN.
            if ~hasSnapshotFormat
                tsfAr = (samplesAVolt < lowerThresholdAVolt) || (higherThresholdAVolt < samplesAVolt);
            else
                % size(ssBltsSamplesTm) = [1        2048           5];
                error
            end
        end



        % NOTE: 2023-11-30: CURRENTLY UNUSED, BUT IS PLANNED TO BE USED.
        %
        % Detect saturation in an arbitrary (one-channel/scalar) signal.
        %
        % Primarily intended for detecting saturation in CWF diffs, but could
        % potentially be used for singles too. Only intended for detecting "full
        % saturation" (not "partial saturation") to be used for setting in
        % quality bits directly.
        %
        % ARGUMENTS
        % =========
        % tt2000Ar
        % y
        %       Column vector of samples over time.
        %       Missing data is represented with NaN (not fill value).
        %       NOTE: Does not specify unit. Can be both physical units or TM,
        %       as long as it is consistent with the threshold arguments.
        % lowerThreshold, higherThreshold
        %       Scalar thresholds. If y exceeds any of these, and does so
        %       "sufficiently often", then the corresponding values are labelled
        %       as saturated.
        %
        %
        % RETURN VALUE
        % ============
        % bSwsf
        %       Array of logical where elements are set for positions where
        %       argument y is saturated (SWSF).
        %
        function bSwsf = determine_SWSF(...
                tt2000Ar, y, ...
                lowerThreshold, higherThreshold, ...
                saturationFractionThreshold, intervalLengthSec)
            % PROPOSAL: Remove assertions and rely on sub-function's assertion.
            % PROPOSAL: Separate threshold saturation from sliding window saturation (the output).
            %   PRO: Can detect threshold saturation separately for subsequences,
            %        using different thresholds, and then sliding window
            %        saturation.
            %        later for entire sequence.
            % TODO-DEC: How handle NaN/missing samples?
            %   Ex: Channels set to NaN when there is no BLTS due to LRX setting.

            %============
            % ASSERTIONS
            %============
            irf.assert.sizes(...
                tt2000Ar, [-1], ...
                y,        [-1] ...
            )
            assert(isa(tt2000Ar, 'int64'))
            assert(isfloat(y))
            assert(isscalar(lowerThreshold))
            assert(isscalar(higherThreshold))
            % NOTE: No assertions on saturationFractionThreshold and
            % intervalLengthSec which are checked by
            % bicas.proc.qual.sliding_window_over_fraction() further below.

            bTsf = (y < lowerThreshold) || (higherThreshold < y);

            bSwsf = bicas.proc.qual.sliding_window_over_fraction(...
                tt2000Ar, bTsf, saturationFractionThreshold, intervalLengthSec);
        end



        % NOTE: 2023-11-30: CURRENTLY UNUSED, BUT IS PLANNED TO BE USED.
        %
        % Given a (timestamped) 1D array of flagged samples (bool), label all
        % samples positions which are part of a sliding window (of a specified
        % length) with a fraction of flagged samples which is above a
        % specified threshold.
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
        % * Timestamps must increase (but not strictly increase). This is
        % important as there is late be data which
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
        %       Length of sliding window.
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
            %   smooth
            %   density (of set bits/flags)
            %   bit, flag
            %
            % PROPOSAL: Move to bicas.utils.
            %   PRO: More generic that quality variables.
            %   PRO: Independent of L1/L1R-L2 proessing in principle.
            %
            % TODO-NI: Distinguishing name for set bits before & after?
            %   PROPOSAL: Suffix 1 & 2
            %   PROPOSAL: Suffix before & after
            %   PROPOSAL: rawFlag vs slidingWindowFlag
            %   PROPOSAL: TSF=Threshold Saturation Flag,
            %             SWSF=Sliding Window Saturation Flag.
            %       CON: Function is generic. Should not make reference to
            %            saturation.
            %
            % TODO-DEC: Exact algorithm to use? How implement?
            %   NOTE: Most data is not saturated.
            %       PROPOSAL: Faster to iterate over saturated samples than
            %                 non-satured.
            %   PROPOSAL: More than x percent saturation within moving/rolling time
            %             period of length t. ==> Label entire period.
            %       TODO-DEC: How handle data that is shorter than window time
            %                 interval?
            %           PROPOSAL: Ignore time interval. Apply fraction to all
            %                     data.
            %       PROPOSAL: Iterate over intervals which are entirely
            %                 threshold saturated or not.
            %       PROPOSAL: Iterate over every length-t interval.
            %           CON: Slow?
            %           NOTE: Must still iterate in both directions.

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
                    timeSecAr, bFlag1Ar, minFlaggedFraction, windowLengthSec);

                % NOTE: Same call as above, except that (1) reversing the order of
                % timestamps and samples, and (2) negating the timestamps (so
                % that they increment despite their order being reversed).
                bFlag2BackwardAr = bicas.proc.L1L2.qual.sliding_window_over_fraction_forward_pass(...
                    -timeSecAr(end:-1:1), bFlag1Ar(end:-1:1), minFlaggedFraction, windowLengthSec);

                bFlag2BackwardAr = bFlag2BackwardAr(end:-1:1);
                
                bFlag2Ar = bFlag2ForwardAr | bFlag2BackwardAr;

%                 bFlag2Ar = bicas.proc.L1L2.qual.sliding_window_over_fraction_algorithm(...
%                     timeSecAr, bFlag1Ar, minFlagFraction, maxWindowLengthSec);
            end
        end    % function
        
        
        
        % NOTE: 2023-11-30: CURRENTLY UNUSED, BUT IS PLANNED TO BE USED.
        %
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
            %        Intended for (1) weighing sections with different sampling
            %        rate, and (2) for including half in the window length.
            % STLW = STL-Weighted
            
            DEBUG_ENABLED = 0;

            % DEBUG
            if DEBUG_ENABLED
                fprintf('--------sliding_window_over_fraction_forward_pass\n')
            end

            n = numel(bFlag1Ar);
            assert(n >= 2)

            % Pre-allocate
            bFlag2Ar = false(size(bFlag1Ar));

            diffSecAr = [Inf; diff(timeSecAr); Inf];
            % NOTE: Returns Inf for array length == 1 which must therefore be
            %       avoided.
            stlSecAr  = min([diffSecAr(1:end-1), diffSecAr(2:end)], [], 2);

            % Modified cumulative sum so that a difference between indices i and
            % i+1 represents the STL of sample i.
            cumulStlwFlagAr = [0; cumsum(bFlag1Ar .* stlSecAr)];

            % =======================================
            % Iterate over time intervals ("windows")
            % =======================================
            iFlagSet1Ar = find(bFlag1Ar);
            i1 = iFlagSet1Ar(1);
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
                    bFlag2Ar(i0:i1) = true;
                end

                if DEBUG_ENABLED
                    fprintf('Found interval i0:i1 = %i:%i\n', i0, i1)
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
        function log_UFV_records(zv_Epoch, zvUfv, logHeaderStr, L)
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
                    utc1  = irf.cdf.TT2000_to_UTC_str(zv_Epoch(iCdfRecord1));
                    utc2  = irf.cdf.TT2000_to_UTC_str(zv_Epoch(iCdfRecord2));
                    L.logf(LL, '    Records %7i-%7i, %s -- %s', ...
                        iCdfRecord1, iCdfRecord2, utc1, utc2);
                end
            end

        end



    end    % methods(Static, Access=private)



end
