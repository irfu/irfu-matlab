%
% Class that collects "processing functions" as public static methods. Only
% covers processing L1/L1R-->L2.
%
% This class is not meant to be instantiated.
%
%
% CODE CONVENTIONS
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like"
%   data, use the first MATLAB array index to represent CDF records.
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% See bicas.proc.L1L2.cal and readme.txt.
%
%
% SOME INTERMEDIATE PROCESSING DATA FORMATS
% =========================================
% - PreDC = Pre-(Demuxing & Calibration) Data
%       Generic data format that can represent all forms of input datasets
%       before demuxing and calibration. Can use an arbitrary number of samples
%       per record. Some variables are therefore not used in CWF output
%       datasets.
% - PostDC = Post-(Demuxing & Calibration) Data
%       Data format that includes calibrated currents & calibrated & demuxed
%       voltages.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
classdef L1L2
%#######################################################################################################################
%
% PROPOSAL: POLICY: Include all functions which set "policy"/configure the output of datasets. -- ABANDONED?
%
% PROPOSAL: Split into smaller files.
%   NOTE: All functions 2021-05-25:
%         function HkSciTime = process_HK_CDF_to_HK_on_SCI_TIME(InSci, InHk, SETTINGS, L)
%         function [PreDc, PostDc] = process_quality_filter_L2(...
%         function CALIBRATION_TABLE_INDEX = normalize_CALIBRATION_TABLE_INDEX(...
%         function assert_PreDC(PreDc)
%         function assert_PostDC(PostDc)
%         function zvUseFillValues = get_UFV_records_from_settings(...
%         function log_UFV_records(zvEpoch, zvUfv, logHeaderStr, L)
%
% PROPOSAL: Submit zVar variable attributes.
%   PRO: Can interpret fill values.
%       Ex: Can doublecheck TDS RSWF snapshot length using fill values and compare with zVar SAMPS_PER_CH (which seems
%           to be bad).
%
% PROPOSAL: Return (to execute_sw_mode), global attributes.
%   PRO: Needed for output datasets: CALIBRATION_TABLE, CALIBRATION_VERSION
%       ~CON: CALIBRATION_VERSION refers to algorithm and should maybe be a SETTING.
%
% PROPOSAL: Classes for PreDc, PostDc, PreDc.Zv, PostDc.Zv.
%   PRO: Better documentation of formats.
%   PRO: Can abolish
%       bicas.proc.L1L2.assert_PreDC
%       bicas.proc.L1L2.assert_PostDC
% PROPOSAL: Class for HkSciTime.
% PROPOSAL: Class for DemuxerOutput.
%
%#######################################################################################################################



    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################



        % Processing function
        %
        % NOTE: Only converts relevant HK zVars to be on SCI Epoch. Later
        % (other) code decides whether to use it (mux mode).
        function HkSciTime = process_HK_CDF_to_HK_on_SCI_TIME(InSci, InHk, SETTINGS, L)
            % PROPOSAL: Separate function for the actual interpolation of data
            %           (changing time array HK-->SCI).
            
            

            % ASSERTIONS
            assert(isa(InSci, 'bicas.InputDataset'))
            assert(isa(InHk,  'bicas.InputDataset'))

            HkSciTime = [];



            %===================================================================
            % Select whether HK should use
            %   (1) (HK) Epoch, or
            %   (2) (HK) ACQUISITION_TIME (not always available).
            % ---------------------------------------------------
            % IMPLEMENTATION NOTE: Historically, there have been datasets where
            % Epoch contains errors, but ACQUISITION_TIME seems OK. This should
            % be phased out eventually.
            %===================================================================
            ACQUISITION_TIME_EPOCH_UTC = SETTINGS.get_fv('INPUT_CDF.ACQUISITION_TIME_EPOCH_UTC');
            USE_ZV_ACQUISITION_TIME_HK = SETTINGS.get_fv('PROCESSING.HK.USE_ZV_ACQUISITION_TIME');
            if USE_ZV_ACQUISITION_TIME_HK
                hkEpoch = bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
                    InHk.Zv.ACQUISITION_TIME, ...
                    ACQUISITION_TIME_EPOCH_UTC);

                L.logf('warning', 'Using HK zVar ACQUISITION_TIME instead of Epoch.')
            else
                hkEpoch = InHk.Zv.Epoch;
            end



            %==================================================================
            % Log time intervals to enable comparing available SCI and HK data
            %==================================================================
            TimeVars = [];    % Temporary struct only used for logging.
            TimeVars.HK_Epoch  = InHk.Zv.Epoch;
            TimeVars.SCI_Epoch = InSci.Zv.Epoch;
            if isfield(InHk.Zv, 'ACQUISITION_TIME')
                TimeVars.HK_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
                        InHk.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            if isfield(InSci.Zv, 'ACQUISITION_TIME') && ~isempty(InSci.Zv.ACQUISITION_TIME)
                TimeVars.SCI_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
                    InSci.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            bicas.utils.log_zVars(TimeVars, SETTINGS, L);



            %===================
            % WARNINGS / ERRORS
            %===================
            if ~issorted(hkEpoch, 'strictascend')
                % Ex: zVar ACQUISITION_TIME in test file
                % TDS___TESTDATA_RGTS_TDS_CALBA_V0.8.6/
                % solo_HK_rpw-bia_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
                % is not monotonically increasing (in fact, it is completely
                % strange).
                error(...
                    ['HK timestamps do not increase monotonically', ...
                    ' (USE_ZV_ACQUISITION_TIME_HK=%g).'], ...
                    USE_ZV_ACQUISITION_TIME_HK)
            end
            if ~irf.utils.ranges_intersect(InSci.Zv.Epoch, hkEpoch)
                %---------------------------------------
                % CASE: SCI does not overlap HK in time
                %---------------------------------------

                % NOTE: "WARNING" (rather than error) only makes sense if it is
                % possible to later meaningfully permit non-intersection.
                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY');
                bicas.default_anomaly_handling(L, ...
                    settingValue, settingKey, 'E+W+illegal', ...
                    'SCI and HK time ranges do not overlap in time.', ...
                    'BICAS:SWModeProcessing')
            elseif ~irf.utils.is_range_subset(InSci.Zv.Epoch, hkEpoch)
                %-------------------------------------------------
                % CASE: SCI does not cover a subset of HK in time
                %-------------------------------------------------
                % NOTE: This anomaly is obviously implied by the anomaly above
                % (SCI, HK do not overlap). It is therefore only meaningful to
                % detect it if the above anomaly is not detected.
                hk1RelativeSec = 1e-9 * (min(hkEpoch) - min(InSci.Zv.Epoch));
                hk2RelativeSec = 1e-9 * (max(hkEpoch) - max(InSci.Zv.Epoch));

                anomalyDescrMsg = sprintf(...
                    ['HK time range is not a superset of SCI time range.', ...
                    ' Can not reliably interpolate HK data for all of SCI.', ...
                    ' HK begins %g s AFTER SCI begins. HK ends %g s BEFORE SCI ends.'], ...
                    hk1RelativeSec, ...
                    -hk2RelativeSec);

                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.HK.TIME_NOT_SUPERSET_OF_SCI_POLICY');
                bicas.default_anomaly_handling(L, ...
                    settingValue, settingKey, 'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:DatasetFormat:SWModeProcessing')
            end



            % Derive time margin within which the nearest HK value will be used.
            % NOTE: Requires >=2 records.
            hkEpochExtrapMargin = mode(diff(hkEpoch)) / 2;
            
            

            %=============================================================
            % Derive MUX_SET
            % --------------
            % NOTE: Only obtains one MUX_SET per record
            %       ==> Can not change MUX_SET in the middle of a record.
            % NOTE: Can potentially obtain MUX_SET from LFR SCI.
            %=============================================================
            HkSciTime.MUX_SET = bicas.utils.interpolate_nearest(...
                hkEpochExtrapMargin, ...
                hkEpoch, ...
                InHk.Zv.HK_BIA_MODE_MUX_SET, ...
                InSci.Zv.Epoch);



            %==================================================================
            % Derive DIFF_GAIN
            % ----------------
            % NOTE: Not perfect handling of time when 1 snapshot/record, since
            % one should ideally use time stamps for every LFR _sample_.
            %==================================================================
            HkSciTime.DIFF_GAIN = bicas.utils.interpolate_nearest(...
                hkEpochExtrapMargin, ...
                hkEpoch, ...
                InHk.Zv.HK_BIA_DIFF_GAIN, ...
                InSci.Zv.Epoch);



            % ASSERTIONS
            irf.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
        end



        % Processing function (L1/L1R-->L2; not L2-->L3).
        %
        % Overwrite selected data in selected CDF records with fill values/NaN.
        % Modify quality zVariables.
        %
        % NOTE: Almost does not modify PreDc.
        %   Exception: Modifies PreDc.Zv.QUALITY_FLAG
        %
        % Sets
        %   PreDc.Zv.QUALITY_FLAG (modifies)
        %   PostDc.Zv.L2_QUALITY_BITMASK
        %   PostDc.Zv.DemuxerOutput
        %   PostDc.Zv.currentAAmpere
        %
        %
        % RATIONALE
        % =========
        % Does NOT want to operate on structs that mimic the input or output
        % datasets, but on struct that are as similiar as possible for all forms
        % of L1R-->L2 processing.
        %
        function [PreDc, PostDc] = process_quality_filter_L2(...
                PreDc, PostDc, NsoTable, SETTINGS, L)

            % NOTE: Adds zVar L2_QUALITY_FLAG to PostDc, technically altering the PostDc format.
            %   NOTE: Also overwrites voltage with fill values.
            %   PROPOSAL: Treat output PostDc as another format?
            %   PROPOSAL: Initialize empty L2_QUALITY_FLAG when PostDc first created.
            %   PROPOSAL: Keep as is. List as optional field in assert_PostDc
            %
            % PROPOSAL: Abolish test functionality. -- IMPLEMENTED
            %   PRO: Test functionality can lead to bugs.
            %   PRO: Can use alternative NSO table (using NSO table path override
            %        setting).
            %
            % PROPOSAL: Generalize function to be used in L3.
            %   CON: Can not be done since this function is meant to have access
            %        to arbitrary L1/L1R and L2 data to make decisions, although
            %        this is not much used yet.

            % ASSERTION
            bicas.proc.L1L2.assert_PreDC(PreDc)
            bicas.proc.L1L2.assert_PostDC(PostDc)
            nRecords = irf.assert.sizes(PreDc.Zv.Epoch, [-1]);



            % NOTE: Preallocates and adds zVar to PostDc.
            PostDc.Zv.L2_QUALITY_BITMASK = zeros(nRecords, 1, 'uint16');



            %============================================
            % Find CDF records to remove due to settings
            %============================================
            zvUfvSettings = bicas.proc.L1L2.get_UFV_records_from_settings(...
                PreDc.Zv.Epoch, PreDc.Zv.MUX_SET, PreDc.isLfr, SETTINGS, L);

            zvUfv = PreDc.Zv.useFillValues | zvUfvSettings;



            %========================================
            % Take actions based on NSO events table
            %========================================
            %testNsoIdsEnabled = SETTINGS.get_fv('PROCESSING.RCS_NSO.TEST_IDS_ENABLED');

            % Variable naming convention:
            % CDF event    = NSO event that overlaps with CDF records.
            % Global event = NSO event in global NSO event table.

            % NOTE: iCdfEventNa = CDF events as indices to global events.
            [bCdfEventRecordsCa, cdfEventNsoIdCa, iCdfEventNa] = ...
                NsoTable.get_NSO_timestamps(PreDc.Zv.Epoch);
            nCdfEvents    = numel(cdfEventNsoIdCa);
            nGlobalEvents = numel(NsoTable.evtNsoIdCa);
            L.logf('info', ...
                ['Searched non-standard operations (NSO) table.', ...
                ' Found %i relevant NSO events out of a total of %i NSO events.'], ...
                nCdfEvents, nGlobalEvents);

            % Index into LOCAL/CDF NSO events table.
            for kCdfEvent = 1:nCdfEvents

                % Index into GLOBAL NSO events table.
                iGlobalEvent = iCdfEventNa(kCdfEvent);
                eventNsoId   = cdfEventNsoIdCa{kCdfEvent};

                %===========================================================
                % Log the relevant NSO event in the GLOBAL NSO events table
                %===========================================================
                L.logf('info', '    %s -- %s %s', ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGlobalEvent)), ...
                    irf.cdf.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array( iGlobalEvent)), ...
                    eventNsoId);



                %==========================================================
                % TEST FUNCTIONALITY
                % ------------------
                % Optionally translate (selected) TEST NSO IDs into actual
                % NSO IDs.
                %==========================================================
%                 eventNsoIdTranslated = irf.utils.translate({...
%                     {bicas.constants.NSOID.TEST_THRUSTER_FIRING}, ...
%                      bicas.constants.NSOID.THRUSTER_FIRING}, ...
%                     eventNsoId, eventNsoId);
%                 if ~testNsoIdsEnabled && ~strcmp(eventNsoId, eventNsoIdTranslated)
%                     % CASE:   (1) Not test mode
%                     %       & (2) NSO ID was translated (changed).
%                     % ==> Original NSO ID was a TEST NSO ID
%                     % ==> NSO should be ignored.
%                     eventNsoIdTranslated = 'nothing';   % Local constant.
%                 end
%                 eventNsoId = eventNsoIdTranslated;
                %========================================================



                %=================================
                % Take action depending on NSO ID
                %=================================
                % Temporary shorter variable name.
                zv_QUALITY_FLAG       = PreDc.Zv.QUALITY_FLAG       (bCdfEventRecordsCa{kCdfEvent});
                zv_L2_QUALITY_BITMASK = PostDc.Zv.L2_QUALITY_BITMASK(bCdfEventRecordsCa{kCdfEvent});

                switch(eventNsoId)

                    %=====================================================
                    % TEST FUNCTIONALITY
                    % Can test the setting of QUALITY_FLAG and zvUfv.
%                     case bicas.constants.NSOID.TEST_QF0
%                         if testNsoIdsEnabled
%                             zv_QUALITY_FLAG = min(zv_QUALITY_FLAG, 0, ...
%                                 'includeNaN');
%                         end
%                     case bicas.constants.NSOID.TEST_UFV
%                         if testNsoIdsEnabled
%                             zvUfv = zvUfv | bCdfEventRecordsCa{kCdfEvent};
%                         end
                    %=====================================================

                    case bicas.constants.NSOID.PARTIAL_SATURATION
                        zv_QUALITY_FLAG       = min(zv_QUALITY_FLAG, 1, 'includeNaN');
                        zv_L2_QUALITY_BITMASK = bitor(...
                            zv_L2_QUALITY_BITMASK, ...
                            bicas.constants.L2QBM_PARTIAL_SATURATION);

                    case bicas.constants.NSOID.FULL_SATURATION
                        zv_QUALITY_FLAG       = min(zv_QUALITY_FLAG, 0, 'includeNaN');
                        zv_L2_QUALITY_BITMASK = bitor(...
                            zv_L2_QUALITY_BITMASK, ...
                            bicas.constants.L2QBM_FULL_SATURATION);
                        zv_L2_QUALITY_BITMASK = bitor(...
                            zv_L2_QUALITY_BITMASK, ...
                            bicas.constants.L2QBM_PARTIAL_SATURATION);
                        % NOTE: Also set PARTIAL saturation bit when FULL
                        % saturation. /YK 2020-10-02.

                    case bicas.constants.NSOID.THRUSTER_FIRING
                        zv_QUALITY_FLAG = min(zv_QUALITY_FLAG, 1, 'includeNaN');
                        % NOTE: There will be an L1 QUALITY_BITMASK bit for
                        % thruster firings eventually according to
                        % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
                        % Therefore(?) not setting any bit in
                        % L2_QUALITY_BITMASK. (YK 2020-11-03 did not ask for any
                        % to be set.)

%                     case 'nothing'
%                         % CASE: Do nothing.
%                         % This case is necessary so that test NSO IDs can be
%                         % translated to something harmless when tests are
%                         % disabled.

                    otherwise
                        % ASSERTION
                        % NOTE: Not perfect assertion on legal NSO IDs since
                        % code only checks those relevant for the data (time
                        % interval) currently processed. (Therefore also checks
                        % all NSO IDs when reads NSO table.)
                        error('Can not interpret RCS NSO ID "%s".', ...
                            cdfEventNsoIdCa{kCdfEvent})

                end
                PreDc.Zv.QUALITY_FLAG       (bCdfEventRecordsCa{kCdfEvent}) = zv_QUALITY_FLAG;
                PostDc.Zv.L2_QUALITY_BITMASK(bCdfEventRecordsCa{kCdfEvent}) = zv_L2_QUALITY_BITMASK;

            end    % for



            %=================================================================
            % Set zVariables for CURRENTS and VOLTAGES to NaN based on zvUfv.
            %=================================================================
            % Log
            logHeaderStr = sprintf(...
                ['All interval(s) of CDF records for which data should be set', ...
                ' to fill values (i.e. removed), regardless of reason.\n']);
            bicas.proc.L1L2.log_UFV_records(PreDc.Zv.Epoch, zvUfv, logHeaderStr, L)
            %
            PostDc.Zv.currentAAmpere(zvUfv, :) = NaN;
            %
            fnCa = fieldnames(PostDc.Zv.DemuxerOutput);
            for iFn = 1:numel(fnCa)
                PostDc.Zv.DemuxerOutput.(fnCa{iFn})(zvUfv, :, :) = NaN;
            end



            % ASSERTION
            bicas.proc.L1L2.assert_PreDC(PreDc)
            bicas.proc.L1L2.assert_PostDC(PostDc)

        end    % process_quality_filter_L2



        % Utility function to shorten code.
        %
        % NOTE: Operates on entire ZvStruct since CALIBRATION_TABLE_INDEX exists
        % for L1R, but not L1, and the corresponding field may thus be or not be
        % present.
        function CALIBRATION_TABLE_INDEX = normalize_CALIBRATION_TABLE_INDEX(...
                ZvStruct, nRecords, inputDsi)

            C = bicas.classify_BICAS_L1_L1R_to_L2_DATASET_ID(inputDsi);

            if C.isL1r
                CALIBRATION_TABLE_INDEX = ZvStruct.CALIBRATION_TABLE_INDEX;
            elseif C.isL1
                CALIBRATION_TABLE_INDEX = nan(nRecords, 2);
            else
                error(...
                    ['Can not normalize CALIBRATION_TABLE_INDEX', ...
                    ' for this DATASET_ID classification.'])
            end

            irf.assert.sizes(CALIBRATION_TABLE_INDEX, [nRecords, 2])
        end



        function assert_PreDC(PreDc)
            irf.assert.struct(PreDc, ...
                {'Zv', 'Ga', 'hasSnapshotFormat', 'isLfr', 'isTdsCwf'}, {});

            irf.assert.struct(PreDc.Zv, ...
                {'Epoch', 'samplesCaTm', 'freqHz', 'nValidSamplesPerRecord', ...
                'iLsf', 'DIFF_GAIN', ...
                'MUX_SET', 'QUALITY_BITMASK', 'QUALITY_FLAG', 'SYNCHRO_FLAG', ...
                'DELTA_PLUS_MINUS', 'CALIBRATION_TABLE_INDEX', ...
                'useFillValues', 'lfrRx'}, ...
                {'BW'});

            bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(PreDc.Zv);

            assert(isa(PreDc.Zv.freqHz, 'double'))
        end



        function assert_PostDC(PostDc)
            irf.assert.struct(PostDc, ...
                {'Zv'}, {});

            irf.assert.struct(PostDc.Zv, ...
                {'DemuxerOutput', 'currentAAmpere'}, {'L2_QUALITY_BITMASK'});

            bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(PostDc.Zv);
        end



    end    % methods(Static, Access=public)



    %##############################
    %##############################
    methods(Static, Access=private)
    %##############################
    %##############################



        % Find CDF records to remove based on settings (not data itself, almost,
        % since MUX mode is data).
        %
        % Ex: Sweeps
        %
        function zvUseFillValues = get_UFV_records_from_settings(...
                zvEpoch, zv_MUX_SET, isLfr, SETTINGS, L)
            % PROPOSAL: Only derive UFV records based on settings. Not take
            %           previously found UFV records (BW) into account. Merging UFV
            %           records from settings and BW respectively can be done
            %           outside (trivial).
            % PROPOSAL: Separate function for logging which records that should be removed.

            bicas.utils.assert_zv_Epoch(zvEpoch)
            assert(islogical(isLfr));

            %===============
            % Read settings
            %===============
            [muxModesRemove, settingMuxModesKey] = SETTINGS.get_fv(...
                'PROCESSING.L2.REMOVE_DATA.MUX_MODES');
            if     isLfr   settingMarginKey = 'PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S';    % LFR
            else           settingMarginKey = 'PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S';    % TDS
            end
            [removeMarginSec, settingMarginKey] = SETTINGS.get_fv(settingMarginKey);

            %==========================================
            % Find exact indices/CDF records to remove
            %==========================================
            zvUseFillValues = irf.utils.true_with_margin(...
                zvEpoch, ...
                ismember(zv_MUX_SET, muxModesRemove), ...
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
                settingMuxModesKey, ...
                strjoin(irf.str.sprintf_many('%g', muxModesRemove), ', '), ...
                settingMarginKey, ...
                removeMarginSec);
            bicas.proc.L1L2.log_UFV_records(zvEpoch, zvUseFillValues, logHeaderStr, L)
        end



        % Log UFV records
        %
        % NOTE: Only logs (including header) if there are records to remove.
        function log_UFV_records(zvEpoch, zvUfv, logHeaderStr, L)
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
                    utc1  = irf.cdf.TT2000_to_UTC_str(zvEpoch(iCdfRecord1));
                    utc2  = irf.cdf.TT2000_to_UTC_str(zvEpoch(iCdfRecord2));
                    L.logf(LL, '    Records %7i-%7i, %s -- %s', ...
                        iCdfRecord1, iCdfRecord2, utc1, utc2);
                end
            end

        end



    end    % methods(Static, Access=private)

end
