%
% Class that collects "processing functions" as public static methods. Only
% covers processing L2-->L3 ("23" in the name).
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
% See bicas.calib.
% ZV  : CDF zVariable, or something analogous to it. If refers to CDF:ish
%       content, then the first index corresponds to the CDF record.
% SPR : Samples Per (CDF) Record. Only refers to actual data (currents,
%       voltages), not metadata.
% UFV : Use Fill Values (refers to records which data should overwritten with
%       fill values)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
classdef proc_sub23
%##############################################################################################
%
% PROPOSAL: POLICY: Include all functions which set "policy"/configure the output of datasets.
%
% PROPOSAL: Move to here: bicas.proc_utils.downsample_Epoch()
%
%##############################################################################################



    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################

    
    
        % Processing function for processing L2-->L3.
        %
        function [OutEfield,  OutEfieldDwns, ...
                  OutScpot,   OutScpotDwns, ...
                  OutDensity, OutDensityDwns] ...
                = process_L2_to_L3(InLfrCwf, SETTINGS, L)

            % PROPOSAL: Split up in one part for non-downsampled and
            %           downsampled.
            % PROPOSAL: Split up into different parts for EFIELD, SCPOT, DENSITY
            %           (still combine non-downsampled and downsampled).
            %   CON: Slows down over all processing.
            %       PRO: Must read same L2 dataset multiple times.
            %       PRO: Must read L3 SCPOT dataset to produce L3 DENSITY dataset.
            %   CON: There is much shared functionality for 3 quality zVars.
            %       PRO: Same ~constants
            %           Ex: INPUT_DATASET_ID, BIN_LENGTH_WOLS_NS, BIN_TIMESTAMP_POS_WOLS_NS
            %       PRO: Read setting QUALITY_FLAG_MIN_FOR_USE
            %       PRO: Normalizing CWF zVar names.
            %       PRO: Preparations for downsampled.
            %           Bin locations, bundling of records,
            %           Downsampling of quality variables
            %               (QUALITY_FLAG, QUALITY_BITMASK, L2_QUALITY_BITMASK).
            %           DELTA_PLUS_MINUS_dwns
            %       --
            %       CON-PROPOSAL: Put shared functionality in function.
            %           CON: Slows down processing.
            %               CON: Probably negligible.
            %
            % PROPOSAL: Downsampled records with fewer than N samples should
            %           set voltage to fill value.
            %   NOTE: May affect QUALITY_FLAG(?)
            %   PROPOSAL: Take into account whether samples are fill values.
            %       NOTE: Leads to different behaviour for different downsampled
            %             datasets.
            %       NOTE: May be different for different "channels" (vary over
            %             non-record dimensions) within the same zVar.
            %
            % BUG: Fill values in the INPUT QUALITY_FLAG,
            % QUALITY_BITMASK, L2_QUALITY_BITMASK are not known since the
            % variables are not double.
            %   PROPOSAL: Use double also for CDF integer variables so NaN can
            %             represent fill value also for these.



            % The only acceptable input DATASET_ID.
            INPUT_DATASET_ID          = 'SOLO_L2_RPW-LFR-SURV-CWF-E';
            % Define length of bins, and relative position of corresponding
            % bin timestamps.
            BIN_LENGTH_WOLS_NS        = int64(10e9);
            BIN_TIMESTAMP_POS_WOLS_NS = int64(BIN_LENGTH_WOLS_NS / 2);

            QUALITY_FLAG_MIN_FOR_USE  = SETTINGS.get_fv(...
                'PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN');



            %======================
            % Normalize zVar names
            %======================
            [InLfrCwf.Zv, fnChangeList] = ...
                EJ_library.utils.normalize_struct_fieldnames(InLfrCwf.Zv, ...
                {{{'VDC', 'V'}, 'VDC'}}, 'Assert one matching candidate');

            bicas.proc_utils.handle_zv_name_change(...
                fnChangeList, INPUT_DATASET_ID, SETTINGS, L, 'VDC', ...
                'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')



            zv_QUALITY_FLAG = InLfrCwf.Zv.QUALITY_FLAG;



            %====================================================================
            % Calculate both
            %   (1) ELECTRIC FIELD, and
            %   (2) SPACECRAFT POTENTIALS
            % via the same BICAS-external code (inside irfu-matlab)
            % -----------------------------------------------------
            % NOTE: Needs to be careful with the units, and incompatible updates
            % to solo.vdccal without the knowledge of the BICAS author.
            % Therefore uses extra assertions to detect such changes.
            %====================================================================
            % Set some records to NaN.
            zv_VDC = InLfrCwf.Zv.VDC;
            zv_VDC(zv_QUALITY_FLAG < QUALITY_FLAG_MIN_FOR_USE, :) = NaN;
            zv_EDC = InLfrCwf.Zv.EDC;
            zv_EDC(zv_QUALITY_FLAG < QUALITY_FLAG_MIN_FOR_USE, :) = NaN;
            %
            % NOTE: Should TSeries objects really use TensorOrder=1 and
            % repres={x,y,z}?!! VDC and EDC are not time series of vectors, but
            % fo three scalars. Probably does not matter. solo.vdccal() does
            % indeed use VDC.x, EDC.x etc.
            VdcTs = TSeries(...
                EpochTT(InLfrCwf.Zv.Epoch), zv_VDC, ...
                'TensorOrder', 1, ...
                'repres', {'x', 'y', 'z'});
            EdcTs = TSeries(...
                EpochTT(InLfrCwf.Zv.Epoch), zv_EDC, ...
                'TensorOrder', 1, ...
                'repres', {'x', 'y', 'z'});
            %-----------------------------------------------------------------
            % CALL EXTERNAL CODE
            [EdcSrfTs, PspTs, ScpotTs, vdccalCodeVerStr, vdccalMatVerStr] ...
                = solo.vdccal(VdcTs, EdcTs);
            %-----------------------------------------------------------------
            EJ_library.assert.sizes(...
                InLfrCwf.Zv.Epoch, [-1, 1], ...
                EdcSrfTs.data,     [-1, 3], ...
                PspTs.data,        [-1, 1], ...
                ScpotTs.data,      [-1, 1]);
            assert(strcmp(EdcSrfTs.units,            'mV/m'))
            assert(strcmp(EdcSrfTs.coordinateSystem, 'SRF'))
            assert(strcmp(PspTs.units,               'V'))
            assert(strcmp(ScpotTs.units,             'V'))



            %===================================================================
            % Normalize E-field
            % -----------------
            % Set E_x = NaN, but ONLY if assertion deems that the corresponding
            % information is missing.
            %
            % IMPLEMENTATION NOTE: solo.vdccal set antenna 1 to be zero, if the
            % source data is non-fill value/NaN, but NaN if fill value. Must
            % therefore check for both zero and NaN.
            % Ex: Dataset 2020-08-01
            %===================================================================
            zvEdcMvpm = EdcSrfTs.data;    % MVPM = mV/m
            % IMPLEMENTATION NOTE: ismember does not work for NaN.
            assert(all(zvEdcMvpm(:, 1) == 0 | isnan(zvEdcMvpm(:, 1))), ...
                ['EDC for antenna 1 returned from', ...
                ' solo.vdccal() is not zero or NaN and can therefore not be', ...
                ' assumed to be unknown anymore.', ...
                ' BICAS needs to be updated to reflect this.'])
            zvEdcMvpm(:, 1) = NaN;
            clear TsEdc



            %====================================================================
            % Calculate DENSITY via a BICAS-external code (inside irfu-matlab)
            % ----------------------------------------------------------------
            % NOTE: Needs to be careful with the units, and incompatible updates
            % to solo.vdccal without the knowledge of the BICAS author.
            % Therefore uses extra assertions to detect such changes.
            %
            % NOTE: Empirically, some return values are NaN.
            % NOTE: "SCP" comes from the return variable name in solo.psp2ne().
            % Do not know what it means.
            %====================================================================
            %-----------------------------
            % CALL EXTERNAL CODE
            [NeScpTs, psp2neCodeVerStr] = solo.psp2ne(PspTs);
            %-----------------------------
            EJ_library.assert.sizes(...
                PspTs.data,   [-1, 1], ...
                NeScpTs.data, [-1, 1]);
            assert(all( (NeScpTs.data > 0) | isnan(NeScpTs.data)), ...
                'solo.psp2ne() returned non-positive (non-NaN) plasma density.')
            assert(strcmp(NeScpTs.units, 'cm^-3'))



            %====================================================================
            % Derive values for CDF global attribute "Misc_calibration_versions"
            %====================================================================
            % Reg.exp. equivalent to: yyyy-mm-ddThh:mm:ss
            CODE_VER_STR_REGEXP = '[0-9]{4}-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]';
            assert(isempty(vdccalMatVerStr), ...
                ['solo.vdccal() no longer returns empty vdccalMatVerStr.', ...
                ' BICAS needs to be updated.'])
            EJ_library.assert.castring_regexp(vdccalCodeVerStr, CODE_VER_STR_REGEXP)
            EJ_library.assert.castring_regexp(psp2neCodeVerStr, CODE_VER_STR_REGEXP)
            %
            % NOTE: Does not set BICAS version since this is already encoded in
            % global attribute "Software_version" (together with
            % "Software_name").
            gaEfieldScpot_Misc_calibration_versions = {};
            gaEfieldScpot_Misc_calibration_versions{end+1} = ...
                ['solo.vdccal() code version ', vdccalCodeVerStr];
            %
            gaDensity_Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            gaDensity_Misc_calibration_versions{end+1}     = ...
                ['solo.psp2ne() code version ', psp2neCodeVerStr];



            %==================
            % Shared variables
            %==================
            % Global attributes, shared between all 3x2 datasets
            Ga = struct();
%             Ga.Datetime           = InLfrCwf.Ga.Datetime;
            Ga.OBS_ID             = InLfrCwf.Ga.OBS_ID;
            Ga.SOOP_TYPE          = InLfrCwf.Ga.SOOP_TYPE;
            % zVariables, shared between all non-downsampled datasets
            InitialZv = struct();
            InitialZv.Epoch              = InLfrCwf.Zv.Epoch;
            InitialZv.QUALITY_BITMASK    = InLfrCwf.Zv.QUALITY_BITMASK;
            InitialZv.L2_QUALITY_BITMASK = InLfrCwf.Zv.L2_QUALITY_BITMASK;
            InitialZv.QUALITY_FLAG       = zv_QUALITY_FLAG;
            InitialZv.DELTA_PLUS_MINUS   = InLfrCwf.Zv.DELTA_PLUS_MINUS;



            %====================================
            % zVars for EFIELD (not downsampled)
            %====================================
            OutEfield = struct();
            OutEfield.Ga                           = Ga;
            OutEfield.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            OutEfield.Zv                           = InitialZv;
            %
            OutEfield.Zv.EDC_SRF                   = zvEdcMvpm;



            %===================================
            % zVars for SCPOT (not downsampled)
            %===================================
            OutScpot = struct();
            OutScpot.Ga                           = Ga;
            OutScpot.Ga.Misc_calibration_versions = gaEfieldScpot_Misc_calibration_versions;
            OutScpot.Zv                           = InitialZv;
            %
            OutScpot.Zv.SCPOT                     = ScpotTs.data;
            OutScpot.Zv.PSP                       = PspTs.data;



            %=====================================
            % zVars for DENSITY (not downsampled)
            %=====================================
            OutDensity = struct();
            OutDensity.Ga                           = Ga;
            OutDensity.Ga.Misc_calibration_versions = gaDensity_Misc_calibration_versions;
            OutDensity.Zv                           = InitialZv;
            %
            OutDensity.Zv.DENSITY                   = NeScpTs.data;



            %====================================================
            % Calculate values used for all downsampled datasets
            %====================================================
            % Find bin boundary reference timestamp. This is used for
            % setting the bin boundaries together with the bin length.
            v = spdfbreakdowntt2000(InLfrCwf.Zv.Epoch(1));
            % UTC subsecond (milliseconds, microseconds, nanoseconds)
            v(7:9) = 0;
            v(6)   = 5;   % UTC second
            boundaryRefTt2000 = spdfcomputett2000(v);
            % Find
            % (1) bin timestamps (downsampled timestamps), and
            % (2) which (non-downsampled) records belong to which bins
            %     (=downsampled records).
            [zvEpochDwns, iRecordsDwnsCa, binSizeArrayNs] = ...
                bicas.proc_utils.downsample_Epoch(...
                    InLfrCwf.Zv.Epoch, boundaryRefTt2000, ...
                    BIN_LENGTH_WOLS_NS,      BIN_TIMESTAMP_POS_WOLS_NS);
            nRecordsDwns = numel(zvEpochDwns);
%             for i = 1:nRecordsDwns
%                 % TODO-DEC: Bad to remove non-donwsampled bins since quality
%                 % variables depend on them ??!!!
%
%                 nSamplesPerBin = numel(iRecordsDwnsCa{i});
%                 if (1 <= nSamplesPerBin) ...
%                 &&  (nSamplesPerBin < bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN)
%                     % NOTE: Does not have to test for 1 <= nSamplesPerBin, but
%                     % it makes it possible to detect (test) whether criterion is
%                     % used.
%                     iRecordsDwnsCa{i} = [];
%                 end
%             end

            % Quality zVariables
            % ------------------
            %
            % Set zVariable-like quality variables with "thought-out" values
            % also for empty bins. Later code can then decide whether to use
            % these empty bin values or not.
            %
            % Correct zVar data types
            % -----------------------
            % "QUALITY_FLAG shall be a CDF_UINT1 flag"
            % "QUALITY_BITMASK shall be a CDF_UINT2 flag"
            % Source: SOL-SGS-TN-0009, "Metadata Definition for Solar Orbiter
            % Science Data"
            % --
            % "The optional CDF_UINT2 zVariable L2_QUALITY_BITMASK /.../"
            % Source:
            % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
            %
            QUALITY_FLAG_dwns       = zeros(nRecordsDwns, 1, 'uint8');
            QUALITY_BITMASK_dwns    = zeros(nRecordsDwns, 1, 'uint16');
            L2_QUALITY_BITMASK_dwns = zeros(nRecordsDwns, 1, 'uint16');
            for i = 1:nRecordsDwns
                k = iRecordsDwnsCa{i};

                QUALITY_FLAG_dwns(i)       = ...
                    bicas.proc_sub23.downsample_bin_QUALITY_FLAG(...
                        zv_QUALITY_FLAG( k) );

                % IMPLEMENTATION NOTE: 2020-11-23: L2 zVar "QUALITY_BITMASK" is
                % mistakenly uint8/CDF_UINT1 when it should be uint16/CDF_UINT2.
                % Must therefore TYPECAST.
                QUALITY_BITMASK_dwns(i)    = ...
                    bicas.proc_sub23.downsample_bin_L12_QUALITY_BITMASK(...
                        uint16( InLfrCwf.Zv.QUALITY_BITMASK( k )) );

                L2_QUALITY_BITMASK_dwns(i) = ...
                    bicas.proc_sub23.downsample_bin_L12_QUALITY_BITMASK(...
                        InLfrCwf.Zv.L2_QUALITY_BITMASK( k ) );
            end

            % Set DELTA_PLUS_MINUS_dwns
            % -------------------------
            % Takes leap seconds into account.
            %
            % NOTE/BUG: Not perfect since the bin timestamp is not centered for
            % leap seconds. Epoch+-DELTA_PLUS_MINUS will thus go outside/inside
            % the bin boundaries for leap seconds. The same problem exists for
            % both positive and negative leap seconds.
            DELTA_PLUS_MINUS_dwns = double(binSizeArrayNs / 2);



            %====================================================
            % Shared zVariables between all downsampled datasets
            %====================================================
            InitialDwnsZv = struct();
            InitialDwnsZv.Epoch              = zvEpochDwns;
            % Below: Pre-allocations.
            % IMPLEMENTATION NOTE: Does not set the final values here to keep it
            % open exactly how they should be set when there are no or too few
            % samples. For loops further below should decide.
            InitialDwnsZv.QUALITY_FLAG       = NaN(nRecordsDwns, 1);
            InitialDwnsZv.QUALITY_BITMASK    = NaN(nRecordsDwns, 1);
            InitialDwnsZv.L2_QUALITY_BITMASK = NaN(nRecordsDwns, 1);
            InitialDwnsZv.DELTA_PLUS_MINUS   = NaN(nRecordsDwns, 1);



            %==============================
            % zVars for EFIELD DOWNSAMPLED
            %==============================
            OutEfieldDwns = [];
            OutEfieldDwns.Ga            = OutEfield.Ga;
            OutEfieldDwns.Zv            = InitialDwnsZv;
            %
            OutEfieldDwns.Zv.EDC_SRF    = NaN(nRecordsDwns, 3);
            OutEfieldDwns.Zv.EDCSTD_SRF = NaN(nRecordsDwns, 3);

            for i = 1:nRecordsDwns
                k = iRecordsDwnsCa{i};
%                 if ~isempty(k)

                    OutEfieldDwns.Zv.QUALITY_FLAG(i)       = QUALITY_FLAG_dwns(i);
                    OutEfieldDwns.Zv.QUALITY_BITMASK(i)    = QUALITY_BITMASK_dwns(i);
                    OutEfieldDwns.Zv.L2_QUALITY_BITMASK(i) = L2_QUALITY_BITMASK_dwns(i);
                    OutEfieldDwns.Zv.DELTA_PLUS_MINUS(i)   = DELTA_PLUS_MINUS_dwns(i);

                    [edc_srf, edcStd_srf] = bicas.proc_sub23.downsample_bin_sci_values(...
                        OutEfield.Zv.EDC_SRF(k, :), ...
                        bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN);

                    OutEfieldDwns.Zv.EDC_SRF(i, :)         = edc_srf;
                    OutEfieldDwns.Zv.EDCSTD_SRF(i, :)      = edcStd_srf;
%                 end
            end



            %=============================
            % zVars for SCPOT DOWNSAMPLED
            %=============================
            OutScpotDwns = [];
            OutScpotDwns.Ga          = OutScpot.Ga;
            OutScpotDwns.Zv          = InitialDwnsZv;
            %
            OutScpotDwns.Zv.SCPOT    = NaN(nRecordsDwns, 1);
            OutScpotDwns.Zv.SCPOTSTD = NaN(nRecordsDwns, 1);
            OutScpotDwns.Zv.PSP      = NaN(nRecordsDwns, 1);
            OutScpotDwns.Zv.PSPSTD   = NaN(nRecordsDwns, 1);

            for i = 1:nRecordsDwns
                k = iRecordsDwnsCa{i};
%                 if ~isempty(k)

                    OutScpotDwns.Zv.QUALITY_FLAG(i)       = QUALITY_FLAG_dwns(i);
                    OutScpotDwns.Zv.QUALITY_BITMASK(i)    = QUALITY_BITMASK_dwns(i);
                    OutScpotDwns.Zv.L2_QUALITY_BITMASK(i) = L2_QUALITY_BITMASK_dwns(i);
                    OutScpotDwns.Zv.DELTA_PLUS_MINUS(i)   = DELTA_PLUS_MINUS_dwns(i);

                    [scpot, scpotStd] = bicas.proc_sub23.downsample_bin_sci_values(...
                        OutScpot.Zv.SCPOT(k, :), ...
                        bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN);
                    [psp, pspstd]     = bicas.proc_sub23.downsample_bin_sci_values(...
                        OutScpot.Zv.PSP(  k, :), ...
                        bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN);

                    OutScpotDwns.Zv.SCPOT(i, :)           = scpot;
                    OutScpotDwns.Zv.SCPOTSTD(i, :)        = scpotStd;
                    OutScpotDwns.Zv.PSP(i)                = psp;
                    OutScpotDwns.Zv.PSPSTD(i)             = pspstd;
%                 end
            end



            %===============================
            % zVars for DENSITY DOWNSAMPLED
            %===============================
            OutDensityDwns = [];
            OutDensityDwns.Ga            = OutDensity.Ga;
            OutDensityDwns.Zv            = InitialDwnsZv;
            %
            OutDensityDwns.Zv.DENSITY    = NaN(nRecordsDwns, 1);
            OutDensityDwns.Zv.DENSITYSTD = NaN(nRecordsDwns, 1);

            for i = 1:nRecordsDwns
                k = iRecordsDwnsCa{i};
%                 if ~isempty(k)

                    OutDensityDwns.Zv.QUALITY_FLAG(i)       = QUALITY_FLAG_dwns(i);
                    OutDensityDwns.Zv.QUALITY_BITMASK(i)    = QUALITY_BITMASK_dwns(i);
                    OutDensityDwns.Zv.L2_QUALITY_BITMASK(i) = L2_QUALITY_BITMASK_dwns(i);
                    OutDensityDwns.Zv.DELTA_PLUS_MINUS(i)   = DELTA_PLUS_MINUS_dwns(i);

                    [density, densityStd] = bicas.proc_sub23.downsample_bin_sci_values(...
                        OutDensity.Zv.DENSITY(k, :), ...
                        bicas.constants.N_MIN_SAMPLES_PER_DWNS_BIN);

                    OutDensityDwns.Zv.DENSITY(i, :)         = density;
                    OutDensityDwns.Zv.DENSITYSTD(i, :)      = densityStd;
%                 end
            end

        end    % process_L2_to_L3



        % Derive median and modified standard deviation over dimension 1. For a
        % range of CDF records in a zVariable (at most 1D/CDF record), construct
        % two zVariables for median+MSTD for the corresponding downsampled CDF
        % record.
        %
        % NOTE: Can handle zero input records.
        % NOTE: Function is only public so that automated test code can access
        % it.
        %
        %
        % ARGUMENTS
        % =========
        % zVarSegment
        %       (iCdfRecord, iChannel).
        % nMinReqSamples
        %       Minimum number of samples (fill value or not) for not returning
        %       fill value.
        %
        %
        % RETURN VALUES
        % =============
        % med  : (1, iChannel). 1xN. Median
        % msdt : (1, iChannel). 1xN. Modified STandard Deviation (MSTD).
        %
        %
        %
        function [med, mstd] = downsample_bin_sci_values(...
                zVarSegment, nMinReqSamples)

            % PROPOSAL: Move to bicas.proc_utils.
            % PROPOSAL: Argument for minimum number of samples in each bin. If
            %           number of samples per bin is below limit, then return
            %           NaN.
            %   PROPOSAL: Take NaN samples into account. Exclude them.
            %       CON: Can not do for
            %
            % PROPOSAL: Include the loop over downsampled records.
            %   PRO: Same procedure for all downsampled datasets.
            %   NOTE: Future variations in procedure could be handle using
            %   parameters.
            %   CON: ~Can/should still not eliminate setting quality zVariables
            %        in loop.
            %
            % PROPOSAL: Merge with
            %   downsample_bin_QUALITY_FLAG
            %   downsample_bin_L12_QUALITY_BITMASK
            %   PRO: Centralizes the conversion from bin to downsampled CDF
            %        record.
            %       CON: SCPOT has two different downsampled zVariables and
            %            therefore calls downsample_bin_sci_values() twice per
            %            bin.

            % ASSERTION
            % Only first two dimensions may be size non-one (with current
            % implementation).
            assert(ismatrix(zVarSegment))
            assert(isscalar(nMinReqSamples))



            nRecords = size(zVarSegment, 1);
            nSpr     = size(zVarSegment, 2);   % SPR = Samples Per (CDF) Record

            % ~NORMALIZATION
            if nRecords < nMinReqSamples
                % CASE: Too few samples. ==> Remove all samples.
                zVarSegment = zVarSegment([], :);
            end

            med  = median(zVarSegment, 1);
            mstd = NaN(1, nSpr);    % Pre-allocate.
            for i = 1:nSpr
                mstd(1, i) = bicas.utils.modif_std_deviation(...
                    zVarSegment(:, i), med(i), 1);
            end
        end
        


    end    % methods(Static, Access=public)



    %##############################
    %##############################
    methods(Static, Access=private)
    %##############################
    %##############################



        % Derive QUALITY_FLAG for one downsampled CDF record, from corresponding
        % non-downsampled records (bin).
        %
        % NOTE: Handles empty bins.
        %
        function QUALITY_FLAG = downsample_bin_QUALITY_FLAG(zv_QUALITY_FLAG_segment)
            % Return NaN or 0 for empty bin?

            % IMPLEMENTATION NOTE: Just using min([zv_QUALITY_FLAG; 0]) does not work.
            if isempty(zv_QUALITY_FLAG_segment)
                QUALITY_FLAG = 0;
            else
                QUALITY_FLAG = min(zv_QUALITY_FLAG_segment);
            end
        end



        % Derive a quality bitmask for one downsampled CDF record, from
        % corresponding non-downsampled records (bin).
        %
        % NOTE: "L12_QUALITY_BITMASK" refers to both zVariables
        %   (1) QUALITY_BITMASK (set in L1), and
        %   (2) L2_QUALITY_BITMASK.
        %
        % NOTE: Handles empty bins.
        %
        function L12_QUALITY_BITMASK = downsample_bin_L12_QUALITY_BITMASK(...
                zv_L12_QUALITY_BITMASK_segment)
            % Return NaN or 0 for empty bin?

            % IMPLEMENTATION NOTE: 2020-11-23: L2 zVar "QUALITY_BITMASK" is
            % mistakenly uint8/CDF_UINT1 when it should be uint16/CDF_UINT2.
            assert(isa(zv_L12_QUALITY_BITMASK_segment, 'uint16'))

            if isempty(zv_L12_QUALITY_BITMASK_segment)
                L12_QUALITY_BITMASK = 0;   % Appropriate?!! Fill value/NaN?
            else
                L12_QUALITY_BITMASK = bicas.utils.bitops.or(zv_L12_QUALITY_BITMASK_segment);
            end
        end

        

    end    % methods(Static, Access=private)

end
