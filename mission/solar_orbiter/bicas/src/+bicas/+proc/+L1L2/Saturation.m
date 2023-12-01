%
% Class that collects functionality related to *DETECTING* saturation.
%
% NOTE: Excludes functionality for setting saturation quality bits and
% QUALITY_FLAG via NSO table.
%
% NOTE: 2023-12-01: Class is partly complete. Class is currently unused, but is
% planned to be used.
%
%
% IMPLEMENTATION NOTE
% ===================
% Class is designed as an instantiable class in order to:
% (1) reduce the number of arguments (eliminates arguments that configure
%     the saturation criteria),
% (2) only extract saturation criteria from SETTINGS once, to possibly increase
%     performance.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Saturation
    % PROPOSAL: Automatic test code.
    %
    % PROPOSAL: Merge higherThresholdAVolt* to a struct somehow.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(GetAccess=private, SetAccess=immutable)
        % How long the sliding window should be when using CDF data.
        cwfSlidingWindowLengthSec

        % Threshold for the sample-length weighted fraction of TSF-labelled
        % samples within a sliding window (CWF) or snapshot for that entire
        % sliding window or snapshot to be labelled as saturated.
        tsfFractionThreshold

        % Higher thresholds for saturation. Sample values above these values, or
        % below the negated value, count as threshold-saturated (TSF).
        higherThresholdAVoltDcSingle
        higherThresholdAVoltDcDiff
        higherThresholdAVoltAcDiffLowGain
        higherThresholdAVoltAcDiffHighGain
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = Saturation(SETTINGS)
            obj.cwfSlidingWindowLengthSec          = SETTINGS.get_fv('PROCESSING.SATURATION.CWF_SLIDING_WINDOW_LENGTH_SEC');
            obj.tsfFractionThreshold               = SETTINGS.get_fv('PROCESSING.SATURATION.TSF_FRACTION_THRESHOLD');

            obj.higherThresholdAVoltDcSingle       = SETTINGS.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.SINGLE');
            obj.higherThresholdAVoltDcDiff         = SETTINGS.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.DIFF');
            obj.higherThresholdAVoltAcDiffLowGain  = SETTINGS.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.LOW_GAIN');
            obj.higherThresholdAVoltAcDiffHighGain = SETTINGS.get_fv('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.HIGH_GAIN');



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
            assert_positive_float(obj.higherThresholdAVoltAcDiffLowGain)
            assert_positive_float(obj.higherThresholdAVoltAcDiffHighGain)
        end



        % NOTE: 2023-12-01: FINISHED BUT CURRENTLY UNUSED, BUT IS PLANNED TO BE
        % USED.
        %
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
        %       Arbitrarily size array. May contain NaN.
        %
        %
        % RETURN VALUE
        % ============
        % tsfAr
        %       float. Same size as samplesAVolt. Whether corresponding elements
        %       in samplesAVolt are deemed to be outside the relevant
        %       thresholds. False is returned for all input elements if there
        %       are no thresholds for this kind of data (e.g. for non-ASR
        %       sources). False is returned for NaN input elements.
        %
        function tsfAr = get_TSF(obj, samplesAVolt, Ssid, biasHighGainFpa)
            % PROPOSAL: Better name.
            %   ~sample-to-TSF
            %       PRO: Can use same maing scheme for TSF-to-SWSF function (for CWF).
            %   ~threshold_saturation

            assert(isfloat(samplesAVolt))
            assert(isa(Ssid, 'bicas.proc.L1L2.SignalSourceId'))
            assert(isa(biasHighGainFpa, 'bicas.utils.FPArray') && isscalar(biasHighGainFpa))

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

                biasHighGain = biasHighGainFpa.logical2doubleNan();
                if Ssid.Asid.is_AC()
                    % CASE: AC diff
                    % -------------
                    if biasHighGain == 0
                        highThresholdAVolt = obj.higherThresholdAVoltAcDiffLowGain;
                    elseif biasHighGain == 1
                        highThresholdAVolt = obj.higherThresholdAVoltAcDiffHighGain;
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



        % NOTE: 2023-12-01: FINISHED BUT CURRENTLY UNUSED, BUT IS PLANNED TO BE
        % USED.
        %
        % Determine whether ONE snapshot should be labelled as saturated.
        %
        % ARGUMENTS
        % =========
        % samplesAVolt
        %   Snapshot samples. (iSampleInSnapshot).
        %   NOTE: Should only contain the length of the snapshot. No padding at
        %         the end of array.
        %
        % RETURN VALUE
        % ============
        % isSaturated
        %       Logical. Scalar.
        %
        function isSaturated = get_snapshot_saturation(obj, samplesAVolt, Ssid, biasHighGain)
            irf.assert.sizes(samplesAVolt, [1, NaN, 1])

            tsfAr = obj.get_TSF(samplesAVolt, Ssid, biasHighGain);

            isSaturated = (sum(tsfAr, 'all') / numel(samplesAVolt)) > obj.tsfFractionThreshold;
        end



        % NOTE: 2023-12-01: FINISHED BUT CURRENTLY UNUSED, BUT IS PLANNED TO BE
        % USED.
        %
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
                zvNValidSamplesPerRecord, zvSamplesAVolt, Ssid, biasHighGain)

            nRecs = irf.assert.sizes(...
                zvNValidSamplesPerRecord, [-1],  ...
                zvSamplesAVolt,           [-1, NaN, 1]);

            isSaturatedAr = false(nRecs, 1);
            for iRec = 1:nRecs
                isSaturatedAr(iRec) = obj.get_snapshot_saturation(...
                    zvSamplesAVolt(iRec, 1:zvNValidSamplesPerRecord(iRec)), ...
                    Ssid, biasHighGain);
            end
        end



        % NOTE: 2023-12-01: UNFINISHED AND UNUSED, BUT IS PLANNED TO BE USED.
        %
        % Given ZV-like variables, get saturation bits for quality bitmask.
        %
        % RETURN VALUE
        % ============
        % isSaturatedAr
        %       (iCdfRecords). Logical. Quality bit for saturation.
        %
        function isSaturatedAr = get_voltage_saturation_quality_bit(...
                obj, tt2000Ar, AsrSamplesAVoltSrm, zvNValidSamplesPerRecord, ...
                bdmFpa, dlrFpa, lrx, biasHighGainFpa, hasSnapshotFormat, L)

            % TODO: Test code

            % ASSERTIONS
            bicas.utils.assert_ZV_Epoch(tt2000Ar)
            assert(islogical(hasSnapshotFormat) && isscalar(hasSnapshotFormat))
            assert(strcmp(bdmFpa.mc, 'uint8'))   % Must be some integer.
            assert(strcmp(dlrFpa.mc, 'logical'))
            assert(isa(lrx, 'double'))
            nRows = irf.assert.sizes(...
                bdmFpa, [-1], ...
                dlrFpa, [-1], ...
                lrx,    [-1]);
            assert(AsrSamplesAVoltSrm.nRows == nRows)

            %===================================================================
            % (1) Find continuous subsequences of records with identical
            %     saturation thresholds.
            % (2) Detect TSFs separately for each such sequence.
            %===================================================================
            [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
                bdmFpa.int2doubleNan(), ...
                dlrFpa.logical2doubleNan(), ...
                biasHighGainFpa.logical2doubleNan(), ...
                lrx);

            L.logf('info', ...
                ['Detecting threshold saturation (voltages) -', ...
                ' One sequence of records with identical settings at a time.'])

            % Bit array, one element per CDF record
            % -------------------------------------
            % NOTE: The meaning and usage of the bits is different for different
            % values of hasSnapshotFormat, hence the neutral name.
            bitAr = false(nRows, 1);

            for iSs = 1:nSs
                iRec1 = iRec1Ar(iSs);
                iRec2 = iRec2Ar(iSs);
                ssNRows = iRec2-iRec1 + 1;

                % CV = Constant values = Values which are constant for the
                %      entire subsequence of records.
                Cv = [];
                Cv.bdmFpa          = bdmFpa(         iRec1);
                Cv.dlrFpa          = dlrFpa(         iRec1);
                Cv.biasHighGainFpa = biasHighGainFpa(iRec1);
                Cv.lrx             = lrx(            iRec1);
                % NOTE: Below variables do not vary over CDF records anyhow.
                Cv.hasSnapshotFormat = hasSnapshotFormat;

                RoutingAr = bicas.proc.L1L2.demuxer.get_routings(...
                    Cv.bdmFpa, Cv.dlrFpa);

                % Array of BLTS's for which there is data.
                iBltsAr = bicas.proc.utils.interpret_LRX(Cv.lrx);

                % Bits for the current sub-sequence.
                ssBitAr = false(ssNRows, 1);

                for iBlts = iBltsAr(:)'
                    Ssid = RoutingAr(iBlts).Ssid;

                    if Ssid.is_ASR()
                        ssBltsSamplesAVolt = AsrSamplesAVoltSrm(Ssid.Asid.s);
                        ssBltsSamplesAVolt = ssBltsSamplesAVolt(iRec1:iRec2, :);

                        % Set ssBltsBitAr=bits for the current subsequence and
                        % BLTS.
                        if hasSnapshotFormat
                            ssBltsBitAr = obj.get_snapshot_saturation_many(...
                                 zvNValidSamplesPerRecord(iRec1:iRec2), ...
                                 ssBltsSamplesAVolt, ...
                                 Ssid, Cv.biasHighGainFpa);
                        else
                            ssBltsBitAr = obj.get_TSF(...
                                ssBltsSamplesAVolt, ...
                                Ssid, Cv.biasHighGainFpa);
                        end

                        % Merge (OR) bits over BLTS's (for current subsequence).
                        ssBitAr = any([ssBitAr, ssBltsBitAr], 2);
                    end
                end    % for iBlts = ...

                bitAr(iRec1:iRec2) = ssBitAr;
            end    % for iSs = ...

            if hasSnapshotFormat
                isSaturatedAr = bitAr;
            else
                isSaturatedAr = bicas.proc.L1L2.qual.sliding_window_over_fraction(...
                    tt2000Ar, bitAr, ...
                    obj.tsfFractionThreshold, obj.cwfSlidingWindowLengthSec);
            end
        end



    end    % methods(Access=public)



end
