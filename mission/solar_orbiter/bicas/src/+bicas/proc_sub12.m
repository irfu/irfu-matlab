%
% Class that collects "processing functions" as public static methods. Only
% covers processing L1/L1R-->L2 ("12" in the name).
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
classdef proc_sub12
%#######################################################################################################################
%
% PROPOSAL: POLICY: Include all functions which set "policy"/configure the output of datasets.
%
% PROPOSAL: Split into smaller files.
%   PROPOSAL: proc_LFR
%   PROPOSAL: proc_TDS
%   PROPOSAL: proc_demux_calib
%   PROPOSAL: Local utility functions are moved to bicas.proc_utils.
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
% PROPOSAL:   process_calibrate_demux
%           & calibrate_demux_voltages
%           should only accept the needed zVars and variables.
%   NOTE: Needs some way of packaging/extracting only the relevant zVars/fields
%         from struct.
%
% PROPOSAL: Redefine PostDC to something that include L2_QUALITY_BITMASK.
%
%#######################################################################################################################



    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################



        % Processing function
        function HkSciTime = process_HK_CDF_to_HK_on_SCI_TIME(InSci, InHk, SETTINGS, L)

            % ASSERTIONS
            EJ_library.assert.struct(InSci, {'Zv', 'Ga', 'filePath'}, {})
            EJ_library.assert.struct(InHk,  {'Zv', 'Ga', 'filePath'}, {})

            HkSciTime = [];



            %===================================================================
            % Select whether HK should use
            %   (1) Epoch, or
            %   (2) ACQUISITION_TIME (not always available).
            % ----------------------------------------------
            % IMPLEMENTATION NOTE: Historically, there have been datasets where
            % Epoch is contains errors, but ACQUISITION_TIME seems OK. This
            % should be phased out eventually.
            %===================================================================
            ACQUISITION_TIME_EPOCH_UTC = SETTINGS.get_fv('INPUT_CDF.ACQUISITION_TIME_EPOCH_UTC');
            USE_ZV_ACQUISITION_TIME_HK = SETTINGS.get_fv('PROCESSING.HK.USE_ZV_ACQUISITION_TIME');
            if USE_ZV_ACQUISITION_TIME_HK
                hkEpoch = bicas.proc_utils.ACQUISITION_TIME_to_TT2000(...
                    InHk.Zv.ACQUISITION_TIME, ...
                    ACQUISITION_TIME_EPOCH_UTC);

                L.logf('warning', 'Using HK zVar ACQUISITION_TIME instead of Epoch.')
            else
                hkEpoch = InHk.Zv.Epoch;
            end



            %==================================================================
            % Log time intervals to enable comparing available SCI and HK data
            %==================================================================
            TimeVars = [];
            TimeVars.HK_Epoch  = InHk.Zv.Epoch;
            TimeVars.SCI_Epoch = InSci.Zv.Epoch;
            if isfield(InHk.Zv, 'ACQUISITION_TIME')
                TimeVars.HK_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc_utils.ACQUISITION_TIME_to_TT2000(...
                        InHk.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            if isfield(InSci.Zv, 'ACQUISITION_TIME') && ~isempty(InSci.Zv.ACQUISITION_TIME)
                TimeVars.SCI_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc_utils.ACQUISITION_TIME_to_TT2000(...
                    InSci.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            bicas.proc_utils.log_zVars(TimeVars, SETTINGS, L);



            if SETTINGS.get_fv('INPUT_CDF.HK.MOVE_TIME_TO_SCI')
                L.log('warning', '===================================================================')
                L.log('warning', 'Moving/adjusting HK time to begin at the same timestamp as voltage.')
                L.log('warning', '===================================================================')
                hkEpoch = hkEpoch - hkEpoch(1) + InSci.Zv.Epoch(1);
            end



            %===================
            % WARNINGS / ERRORS
            %===================
            if ~issorted(hkEpoch, 'strictascend')
                % NOTE: ACQUISITION_TIME in test file
                % TDS___TESTDATA_RGTS_TDS_CALBA_V0.8.6/solo_HK_rpw-bia_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
                % is not monotonically increasing (in fact, it is completely
                % strange).
                error([...
                    'HK timestamps do not increase monotonically', ...
                    ' (USE_ZV_ACQUISITION_TIME_HK=%g).'], ...
                    USE_ZV_ACQUISITION_TIME_HK)
            end
            if ~EJ_library.utils.is_range_subset(InSci.Zv.Epoch, hkEpoch)
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
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
            end
            if ~EJ_library.utils.ranges_intersect(InSci.Zv.Epoch, hkEpoch)

                % NOTE: "WARNING" (rather than error) only makes sense if it is
                % possible to later meaningfully permit non-intersection.
                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY');
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    'SCI and HK time ranges do not overlap in time.', ...
                    'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
            end



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
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
        end



        function currentSAmpere = process_CUR_to_CUR_on_SCI_TIME(sciEpoch, InCur, SETTINGS, L)
            % PROPOSAL: Change function name. process_* implies converting struct-->struct.

            % ASSERTIONS
            EJ_library.assert.struct(InCur, {'Zv', 'Ga', 'filePath'}, {})



            %===================================================================
            % CDF ASSERTION: CURRENT data begins before SCI data (i.e. there is
            % enough CURRENT data).
            %===================================================================
            if ~(min(InCur.Zv.Epoch) <= min(sciEpoch))
                curRelativeSec    = 1e-9 * (min(InCur.Zv.Epoch) - min(sciEpoch));
                sciEpochUtcStr    = EJ_library.cdf.TT2000_to_UTC_str(min(sciEpoch));
                curEpochMinUtcStr = EJ_library.cdf.TT2000_to_UTC_str(min(InCur.Zv.Epoch));

                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY');

                anomalyDescrMsg = sprintf(...
                    ['Bias current data begins %g s (%s) AFTER voltage data begins (%s).', ....
                    ' Can therefore not determine currents for all voltage timestamps.'], ...
                    curRelativeSec, curEpochMinUtcStr, sciEpochUtcStr);

                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:proc_sub:SWModeProcessing')
            end



            %====================================================================
            % CDF ASSERTION: Epoch increases (not monotonically)
            % --------------------------------------------------
            % NOTE: bicas.proc_sub12.zv_TC_to_current() checks (and handles)
            % that Epoch increases monotonically, but only for each antenna
            % separately (which does not capture all cases). Therefore checks
            % that Epoch is (non-monotonically) increasing.
            % Ex: Timestamps, iAntenna = mod(iRecord,3): 1,2,3,5,4,6
            %       ==> Monotonically increasing sequences for each antenna
            %           separately, but not even increasing when combined.
            %====================================================================
            assert(issorted(InCur.Zv.Epoch), ...
                'BICAS:proc_sub:DatasetFormat', ...
                'CURRENT timestamps zVar Epoch does not increase (all antennas combined).')

            % NOTE: bicas.proc_sub12.zv_TC_to_current() checks that Epoch
            % increases monotonically.
            currentNanoSAmpere = [];
            currentNanoSAmpere(:,1) = bicas.proc_sub12.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_1, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,2) = bicas.proc_sub12.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_2, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,3) = bicas.proc_sub12.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_3, sciEpoch, L, SETTINGS);

            currentSAmpere = 1e-9 * currentNanoSAmpere;
        end



        % Processing function. Only "normalizes" data to account for technically
        % illegal input LFR datasets. This should try to:
        % ** modify L1 to look like L1R
        % ** mitigate historical bugs (in the input datasets)
        % ** mitigate for not yet implemented features (in input datasets)
        %
        function InSciNorm = process_LFR_CDF_normalize(InSci, inSciDsi, SETTINGS, L)

            % Default behaviour: Copy values, except for values which are
            % modified later
            InSciNorm = InSci;

            nRecords = EJ_library.assert.sizes(InSci.Zv.Epoch, [-1]);



            %===================================
            % Normalize CALIBRATION_TABLE_INDEX
            %===================================
            InSciNorm.Zv.CALIBRATION_TABLE_INDEX = ...
                bicas.proc_sub12.normalize_CALIBRATION_TABLE_INDEX(...
                    InSci.Zv, nRecords, inSciDsi);



            %========================
            % Normalize SYNCHRO_FLAG
            %========================
            has_SYNCHRO_FLAG      = isfield(InSci.Zv, 'SYNCHRO_FLAG');
            has_TIME_SYNCHRO_FLAG = isfield(InSci.Zv, 'TIME_SYNCHRO_FLAG');
            if      has_SYNCHRO_FLAG && ~has_TIME_SYNCHRO_FLAG

                % CASE: Everything nominal.
                InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.SYNCHRO_FLAG;

            elseif ~has_SYNCHRO_FLAG &&  has_TIME_SYNCHRO_FLAG

                % CASE: Input CDF uses wrong zVar name.
                [settingValue, settingKey] = SETTINGS.get_fv('INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY');
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    'Found zVar TIME_SYNCHRO_FLAG instead of SYNCHRO_FLAG.')
                L.log('warning', 'Using illegally named zVar TIME_SYNCHRO_FLAG as SYNCHRO_FLAG.')
                InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.TIME_SYNCHRO_FLAG;

            elseif has_SYNCHRO_FLAG &&  has_TIME_SYNCHRO_FLAG

                % CASE: Input CDF has two zVars: one with correct name, one with
                % incorrect name

                %------------------------
                % "Normal" normalization
                %------------------------
                % 2020-01-21: Based on skeletons (.skt; L1R, L2), SYNCHRO_FLAG
                % seems to be the correct zVar.
                if SETTINGS.get_fv('INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED') ...
                        && isempty(InSci.Zv.SYNCHRO_FLAG)
                    %----------------------------------------------------------
                    % Workaround: Normalize LFR data to handle variations that
                    % should not exist
                    %----------------------------------------------------------
                    % Handle that SYNCHRO_FLAG (empty) and TIME_SYNCHRO_FLAG
                    % (non-empty) may BOTH be present. "DEFINITION BUG" in
                    % definition of datasets/skeleton?
                    % Ex: LFR___TESTDATA_RGTS_LFR_CALBUT_V0.7.0/ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_4129f0b_CNE_V02.cdf /2020-03-17

                    InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.TIME_SYNCHRO_FLAG;
                else
                    error('BICAS:process_LFR_CDF_normalize:DatasetFormat', ...
                        'Input dataset has both zVar SYNCHRO_FLAG and TIME_SYNCHRO_FLAG.')
                end
            else
                error('BICAS:process_LFR_CDF_normalize:DatasetFormat', ...
                    'Input dataset does not have zVar SYNCHRO_FLAG as expected.')
            end



            %=======================================================================================================
            % Set QUALITY_BITMASK, QUALITY_FLAG:
            % Replace illegally empty data with fill values/NaN
            % ------------------------------------------------------------------
            % IMPLEMENTATION NOTE: QUALITY_BITMASK, QUALITY_FLAG have been found
            % empty in test data, but should have attribute DEPEND_0 = "Epoch"
            % ==> Should have same number of records as Epoch.
            %
            % Can not save CDF with zVar with zero records (crashes when reading
            % CDF). ==> Better create empty records.
            %
            % Examples of QUALITY_FLAG = empty:
            %  MYSTERIOUS_SIGNAL_1_2016-04-15_Run2__7729147__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-SWF_7729147_CNE_V01.cdf
            %  ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_4129f0b_CNE_V02.cdf (TESTDATA_RGTS_LFR_CALBUT_V1.1.0)
            %  ROC-SGSE_L1R_RPW-LFR-SBM2-CWF-E_6b05822_CNE_V02.cdf (TESTDATA_RGTS_LFR_CALBUT_V1.1.0)
            %=======================================================================================================
            % PROPOSAL: Move to the code that reads CDF datasets instead. Generalize to many zVariables.
            % PROPOSAL: Regard as "normalization" code. ==> Group together with other normalization code.
            %=======================================================================================================
            [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY');

            InSciNorm.Zv.QUALITY_BITMASK = bicas.proc_sub12.normalize_LFR_zVar_empty(...
                L, settingValue, settingKey, nRecords, ...
                InSci.Zv.QUALITY_BITMASK, 'QUALITY_BITMASK');

            InSciNorm.Zv.QUALITY_FLAG    = bicas.proc_sub12.normalize_LFR_zVar_empty(...
                L, settingValue, settingKey, nRecords, ...
                InSci.Zv.QUALITY_FLAG,    'QUALITY_FLAG');

            % ASSERTIONS
            EJ_library.assert.sizes(...
                InSciNorm.Zv.QUALITY_BITMASK, [nRecords, 1], ...
                InSciNorm.Zv.QUALITY_FLAG,    [nRecords, 1])

        end    % process_LFR_CDF_normalize



        % Processing function. Convert LFR CDF data to PreDC.
        %
        % IMPLEMENTATION NOTE: Does not modify InSci in an attempt to save RAM
        % (should help MATLAB's optimization). Unclear if actually works.
        %
        function PreDc = process_LFR_CDF_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
            %
            % PROBLEM: Hard-coded CDF data types (MATLAB classes).
            % MINOR PROBLEM: Still does not handle LFR zVar TYPE for determining
            % "virtual snapshot" length. Should only be relevant for
            % V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF (not V02) which should expire.

            % ASSERTIONS: VARIABLES
            EJ_library.assert.struct(InSci,     {'Zv', 'Ga', 'filePath'}, {})
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})

            % ASSERTIONS: CDF
            assert(issorted(InSci.Zv.Epoch, 'strictascend'), ...
                'BICAS:process_LFR_CDF_to_PreDC:DatasetFormat', ...
                'Voltage (science) dataset timestamps Epoch do not increase monotonously.')
            nRecords = EJ_library.assert.sizes(InSci.Zv.Epoch, [-1]);



            C = EJ_library.so.adm.classify_BICAS_L1_L1R_to_L2_DATASET_ID(inSciDsi);



            %============
            % Set iLsfZv
            %============
            if     C.isLfrSbm1   iLsfZv = ones(nRecords, 1) * 2;   % Always value "2" (F1, "FREQ = 1").
            elseif C.isLfrSbm2   iLsfZv = ones(nRecords, 1) * 3;   % Always value "3" (F2, "FREQ = 2").
            else                 iLsfZv = InSci.Zv.FREQ + 1;
                % NOTE: Translates from LFR's FREQ values (0=F0 etc) to LSF
                % index values (1=F0) used in loaded RCT data structs.
            end
            EJ_library.assert.sizes(iLsfZv, [nRecords])



            % NOTE: Needed also for 1 SPR.
            zvFreqHz = EJ_library.so.get_LFR_frequency( iLsfZv );

            % Obtain the relevant values (one per record) from zVariables R0,
            % R1, R2, and the virtual "R3".
            zv_Rx = EJ_library.so.get_LFR_Rx(...
                InSci.Zv.R0, ...
                InSci.Zv.R1, ...
                InSci.Zv.R2, ...
                iLsfZv );



            %===================================================================
            % IMPLEMENTATION NOTE: E & V must be floating-point so that values
            % can be set to NaN.
            %
            % Switch last two indices of E.
            % ==> index 2 = "snapshot" sample index, including for CWF
            %               (sample/record, "snapshots" consisting of 1 sample).
            %     index 3 = E1/E2 component
            %               NOTE: 1/2=index into array; these are diffs but not
            %               equivalent to any particular diffs).
            %===================================================================
            E = single(permute(InSci.Zv.E, [1,3,2]));

            % ASSERTIONS
            nCdfSamplesPerRecord = EJ_library.assert.sizes(...
                InSci.Zv.V, [nRecords, -1], ...
                E,          [nRecords, -1, 2]);
            if C.isLfrSurvSwf   assert(nCdfSamplesPerRecord == EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH)
            else                assert(nCdfSamplesPerRecord == 1)
            end



            PreDc = [];

            PreDc.Zv.samplesCaTm    = cell(5,1);
            PreDc.Zv.samplesCaTm{1} = single(InSci.Zv.V);
            % Copy values, except when zvRx==0 (==>NaN).
            PreDc.Zv.samplesCaTm{2} = bicas.proc_utils.filter_rows( E(:,:,1), zv_Rx==0 );
            PreDc.Zv.samplesCaTm{3} = bicas.proc_utils.filter_rows( E(:,:,2), zv_Rx==0 );
            PreDc.Zv.samplesCaTm{4} = bicas.proc_utils.filter_rows( E(:,:,1), zv_Rx==1 );
            PreDc.Zv.samplesCaTm{5} = bicas.proc_utils.filter_rows( E(:,:,2), zv_Rx==1 );

            PreDc.Zv.Epoch                   = InSci.Zv.Epoch;
            PreDc.Zv.DELTA_PLUS_MINUS        = bicas.proc_utils.derive_DELTA_PLUS_MINUS(...
                zvFreqHz, nCdfSamplesPerRecord);
            PreDc.Zv.freqHz                  = zvFreqHz;
            PreDc.Zv.nValidSamplesPerRecord  = ones(nRecords, 1) * nCdfSamplesPerRecord;
            PreDc.Zv.BW                      = InSci.Zv.BW;
            PreDc.Zv.useFillValues           = ~logical(InSci.Zv.BW);
            PreDc.Zv.DIFF_GAIN               = HkSciTime.DIFF_GAIN;
            PreDc.Zv.iLsf                    = iLsfZv;

            PreDc.Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
            PreDc.Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;

            PreDc.Zv.QUALITY_BITMASK         = InSci.Zv.QUALITY_BITMASK;
            PreDc.Zv.QUALITY_FLAG            = InSci.Zv.QUALITY_FLAG;



            %==================================================================
            % Set MUX_SET
            % -----------
            % Select which source of mux mode is used: LFR datasets or BIAS HK
            %==================================================================
            [value, key] = SETTINGS.get_fv('PROCESSING.LFR.MUX_MODE_SOURCE');
            switch(value)
                case 'BIAS_HK'
                    L.log('debug', 'Using BIAS HK mux mode.')
                    PreDc.Zv.MUX_SET = HkSciTime.MUX_SET;
                case 'LFR_SCI'
                    L.log('debug', 'Using LFR SCI mux mode.')
                    PreDc.Zv.MUX_SET = InSci.Zv.BIAS_MODE_MUX_SET;
                otherwise
                    error('BICAS:proc_sub:ConfigurationBug', ...
                        'Illegal settings value %s="%s"', key, value)
            end



%             PreDc.Ga.Datetime       = InSci.Ga.Datetime;
            PreDc.Ga.OBS_ID         = InSci.Ga.OBS_ID;
            PreDc.Ga.SOOP_TYPE      = InSci.Ga.SOOP_TYPE;

            PreDc.hasSnapshotFormat = C.isLfrSurvSwf;
            PreDc.isLfr             = true;
            PreDc.isTdsCwf          = false;



            % ASSERTIONS
            bicas.proc_sub12.assert_PreDC(PreDc)

        end    % process_LFR_CDF_to_PreDC



        % Processing function. Only "normalizes" data to account for technically
        % illegal input TDS datasets. It should try to:
        % ** modify L1 to look like L1R
        % ** mitigate historical bugs in the input datasets
        % ** mitigate for not yet implemented features in input datasets
        %
        function InSciNorm = process_TDS_CDF_normalize(InSci, inSciDsi, SETTINGS, L)

            % Default behaviour: Copy values, except for values which are
            % modified later
            InSciNorm = InSci;

            nRecords = EJ_library.assert.sizes(InSci.Zv.Epoch, [-1]);

            C = EJ_library.so.adm.classify_BICAS_L1_L1R_to_L2_DATASET_ID(inSciDsi);


            %===================================
            % Normalize CALIBRATION_TABLE_INDEX
            %===================================
            InSciNorm.Zv.CALIBRATION_TABLE_INDEX = bicas.proc_sub12.normalize_CALIBRATION_TABLE_INDEX(...
                InSci.Zv, nRecords, inSciDsi);



            %===========================================================
            % Normalize zVar name SYNCHRO_FLAG
            % --------------------------------
            % Both zVars TIME_SYNCHRO_FLAG, SYNCHRO_FLAG found in input
            % datasets. Unknown why. "DEFINITION BUG" in definition of
            % datasets/skeleton? /2020-01-05
            % Based on skeletons (.skt; L1R, L2), SYNCHRO_FLAG seems
            % to be the correct one. /2020-01-21
            %===========================================================
            [InSci.Zv, fnChangeList] = EJ_library.utils.normalize_struct_fieldnames(...
                InSci.Zv, ...
                {{{'TIME_SYNCHRO_FLAG', 'SYNCHRO_FLAG'}, 'SYNCHRO_FLAG'}}, ...
                'Assert one matching candidate');

            bicas.proc_utils.handle_zv_name_change(...
                fnChangeList, inSciDsi, SETTINGS, L, ...
                'SYNCHRO_FLAG', 'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')



            %=========================
            % Normalize SAMPLING_RATE
            %=========================
            if any(InSci.Zv.SAMPLING_RATE == 255)
                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY');
                anomalyDescrMsg = ...
                    ['Finds illegal, stated sampling frequency', ...
                    ' 255 in TDS L1/L1R LFM-RSWF dataset.'];

                if C.isTdsRswf
                    switch(settingValue)
                        case 'CORRECT'
                            %===================================================
                            % IMPLEMENTATION NOTE: Has observed test file
                            % TESTDATA_RGTS_TDS_CALBA_V0.8.5C:
                            % solo_L1R_rpw-tds-lfm-rswf-e_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
                            % to have SAMPLING_RATE == 255, which is likely a
                            % BUG in the dataset.
                            % /Erik P G Johansson 2019-12-03
                            % Is bug in TDS RCS.  /David Pisa 2019-12-03
                            % Setting it to what is probably the correct value.
                            %===================================================
                            InSciNorm.Zv.SAMPLING_RATE(InSci.Zv.SAMPLING_RATE == 255) = 32768;
                            L.logf('warning', ...
                                'Using workaround to modify instances of sampling frequency 255-->32768.')
                            bicas.default_anomaly_handling(L, ...
                                settingValue, settingKey, 'other', anomalyDescrMsg)

                        otherwise
                            bicas.default_anomaly_handling(L, ...
                                settingValue, settingKey, 'E+W+illegal', ...
                                anomalyDescrMsg, 'BICAS:process_TDS_CDF_normalize:DatasetFormat')
                    end
                else
                    error(anomalyDescrMsg)
                end
            end



            if C.isTdsRswf
                %============================================================
                % Check for and handle illegal input data, zVar SAMPS_PER_CH
                % ----------------------------------------------------------
                % NOTE: Has observed invalid SAMPS_PER_CH value 16562 in
                % ROC-SGSE_L1R_RPW-TDS-LFM-RSWF-E_73525cd_CNE_V03.CDF.
                % 2019-09-18, David Pisa: Not a flaw in TDS RCS but in the
                % source L1 dataset.
                %============================================================
                zv_SAMPS_PER_CH_corrected = round(2.^round(log2(double(InSci.Zv.SAMPS_PER_CH))));
                zv_SAMPS_PER_CH_corrected = cast(zv_SAMPS_PER_CH_corrected, class(InSci.Zv.SAMPS_PER_CH));
                zv_SAMPS_PER_CH_corrected = max( zv_SAMPS_PER_CH_corrected, EJ_library.so.constants.TDS_RSWF_SNAPSHOT_LENGTH_MIN);
                zv_SAMPS_PER_CH_corrected = min( zv_SAMPS_PER_CH_corrected, EJ_library.so.constants.TDS_RSWF_SNAPSHOT_LENGTH_MAX);

                if any(zv_SAMPS_PER_CH_corrected ~= InSci.Zv.SAMPS_PER_CH)
                    % CASE: SAMPS_PER_CH has at least one illegal value

                    SAMPS_PER_CH_badValues = unique(InSci.Zv.SAMPS_PER_CH(zv_SAMPS_PER_CH_corrected ~= InSci.Zv.SAMPS_PER_CH));

                    badValuesDisplayStr = strjoin(arrayfun(...
                        @(n) sprintf('%i', n), SAMPS_PER_CH_badValues, 'uni', false), ', ');
                    anomalyDescrMsg = sprintf(...
                        ['TDS LFM RSWF zVar SAMPS_PER_CH contains unexpected', ...
                        ' value(s) which are not on the form 2^n and in the', ...
                        ' interval %.0f to %.0f: %s'], ...
                        EJ_library.so.constants.TDS_RSWF_SNAPSHOT_LENGTH_MIN, ...
                        EJ_library.so.constants.TDS_RSWF_SNAPSHOT_LENGTH_MAX, ...
                        badValuesDisplayStr);

                    [settingValue, settingKey] = SETTINGS.get_fv(...
                        'PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY');
                    switch(settingValue)
                        case 'ROUND'
                            bicas.default_anomaly_handling(...
                                L, settingValue, settingKey, 'other', ...
                                anomalyDescrMsg, ...
                                'BICAS:proc_sub:process_TDS_CDF_normalize:Assertion:DatasetFormat')
                            % NOTE: Logging the mitigation, NOT the anomaly
                            % itself.
                            L.logf('warning', ...
                                ['Replacing TDS RSWF zVar SAMPS_PER_CH', ...
                                ' values with values, rounded to valid', ...
                                ' values due to setting %s.'], ...
                                settingKey)

                            InSciNorm.Zv.SAMPS_PER_CH = zv_SAMPS_PER_CH_corrected;

                        otherwise
                            bicas.default_anomaly_handling(L, ...
                                settingValue, settingKey, 'E+W+illegal', ...
                                anomalyDescrMsg, ...
                                'BICAS:proc_sub:process_TDS_CDF_normalize:Assertion:DatasetFormat')

                    end    % switch
                end    % if
            end    % if

        end    % process_TDS_CDF_normalize



        % Processing function. Convert TDS CDF data (PDs) to PreDC.
        function PreDc = process_TDS_CDF_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
        %
        % BUG?: Does not use CHANNEL_STATUS_INFO.
        % NOTE: BIAS output datasets do not have a variable for the length of
        % snapshots. Need to use NaN/fill value.

            % ASSERTIONS: VARIABLES
            EJ_library.assert.struct(InSci,     {'Zv', 'Ga', 'filePath'}, {})
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})

            C = EJ_library.so.adm.classify_BICAS_L1_L1R_to_L2_DATASET_ID(inSciDsi);



            % ASSERTIONS: CDF
            assert(issorted(InSci.Zv.Epoch, 'strictascend'), ...
                'BICAS:process_TDS_CDF_to_PreDC:DatasetFormat', ...
                'Voltage (science) dataset timestamps Epoch do not increase monotonously.')
            [nRecords, WAVEFORM_DATA_nChannels, nCdfSamplesPerRecord] = EJ_library.assert.sizes(...
                InSci.Zv.Epoch, [-1], ...
                InSci.Zv.WAVEFORM_DATA, [-1, -2, -3]);
            if     C.isL1r   WAVEFORM_DATA_nChannels_expected = 3;
            elseif C.isL1    WAVEFORM_DATA_nChannels_expected = 8;
            end
            assert(...
                WAVEFORM_DATA_nChannels == WAVEFORM_DATA_nChannels_expected, ...
                'BICAS:proc_sub:process_TDS_CDF_to_PreDC:Assertion:DatasetFormat', ...
                'TDS zVar WAVEFORM_DATA has an unexpected size.')
            if C.isTdsRswf   assert(nCdfSamplesPerRecord == EJ_library.so.constants.TDS_RSWF_SAMPLES_PER_RECORD)
            else             assert(nCdfSamplesPerRecord == 1)
            end



            % TODO-NI: Why convert to double? To avoid precision problems when
            % doing math with other variables?
            freqHzZv = double(InSci.Zv.SAMPLING_RATE);



            PreDc = [];

            PreDc.Zv.Epoch                   = InSci.Zv.Epoch;
            PreDc.Zv.DELTA_PLUS_MINUS        = bicas.proc_utils.derive_DELTA_PLUS_MINUS(...
                freqHzZv, nCdfSamplesPerRecord);
            PreDc.Zv.freqHz                  = freqHzZv;
            PreDc.Zv.QUALITY_BITMASK         = InSci.Zv.QUALITY_BITMASK;
            PreDc.Zv.QUALITY_FLAG            = InSci.Zv.QUALITY_FLAG;
            PreDc.Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
            PreDc.Zv.MUX_SET                 = HkSciTime.MUX_SET;
            PreDc.Zv.DIFF_GAIN               = HkSciTime.DIFF_GAIN;
            PreDc.Zv.useFillValues           = false(nRecords, 1);
            PreDc.Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;



            %=====================================
            % Set PreDc.Zv.nValidSamplesPerRecord
            %=====================================
            if C.isTdsRswf
                %================================================================
                % NOTE: This might only be appropriate for TDS's "COMMON_MODE"
                % mode. TDS also has a "FULL_BAND" mode with 2^18=262144 samples
                % per snapshot. You should never encounter FULL_BAND in any
                % dataset (even on ground), only used for calibration and
                % testing. /David Pisa & Jan Soucek in emails, 2016.
                % --
                % FULL_BAND mode has each snapshot divided into 2^15
                % samples/record * 8 records.  /Unknown source. Unclear what
                % value SAMPS_PER_CH should have for FULL_BAND mode. How does
                % Epoch work for FULL_BAND snapshots?
                %================================================================
                % Converting to double because code did so before code
                % reorganization. Reason unknown. Needed to avoid precision
                % problems when doing math with other variables?
                PreDc.Zv.nValidSamplesPerRecord = double(InSci.Zv.SAMPS_PER_CH);
            else
                PreDc.Zv.nValidSamplesPerRecord = ones(nRecords, 1) * 1;
            end
            assert(all(PreDc.Zv.nValidSamplesPerRecord <= nCdfSamplesPerRecord), ...
                'BICAS:proc_sub:process_TDS_CDF_to_PreDC:Assertion:DatasetFormat', ...
                ['Dataset indicates that the number of valid samples per CDF', ...
                ' record (max(PreDc.Zv.nValidSamplesPerRecord)=%i) is', ...
                ' NOT fewer than the number of indices per CDF record', ...
                ' (nCdfMaxSamplesPerSnapshot=%i).'], ...
                max(PreDc.Zv.nValidSamplesPerRecord), ...
                nCdfSamplesPerRecord)



            %==========================
            % Set PreDc.Zv.samplesCaTm
            %==========================
            modif_WAVEFORM_DATA = double(permute(InSci.Zv.WAVEFORM_DATA, [1,3,2]));

            PreDc.Zv.samplesCaTm    = cell(5,1);
            PreDc.Zv.samplesCaTm{1} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,1), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{2} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,2), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{3} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,3), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{4} = nan(nRecords, nCdfSamplesPerRecord);
            PreDc.Zv.samplesCaTm{5} = nan(nRecords, nCdfSamplesPerRecord);



%             PreDc.Ga.Datetime       = InSci.Ga.Datetime;
            PreDc.Ga.OBS_ID         = InSci.Ga.OBS_ID;
            PreDc.Ga.SOOP_TYPE      = InSci.Ga.SOOP_TYPE;

            PreDc.isLfr             = false;
            PreDc.isTdsCwf          = C.isTdsCwf;
            PreDc.hasSnapshotFormat = C.isTdsRswf;
            % Only set because the code shared with LFR requires it.
            PreDc.Zv.iLsf           = nan(nRecords, 1);



            % ASSERTIONS
            bicas.proc_sub12.assert_PreDC(PreDc)

        end    % process_TDS_CDF_to_PreDC



        function [OutSci] = process_PostDC_to_LFR_CDF(SciPreDc, SciPostDc, outputDsi, L)
            OutSci    = bicas.proc_sub12.process_PostDC_to_TDS_CDF(...
                SciPreDc, SciPostDc, outputDsi, L);

            OutSci.Zv.BW = SciPreDc.Zv.BW;
        end



        % Processing function. Convert PreDc+PostDC to something that
        % (1) almost an LFR dataset (the rest is done in a wrapper).
        % (2) a TDS dataset (hence the name), and
        %
        % This function only changes the data format (and selects data to send
        % to CDF).
        %
        function [OutSci] = process_PostDC_to_TDS_CDF(SciPreDc, SciPostDc, outputDsi, L)

            % ASSERTIONS
            bicas.proc_sub12.assert_PostDC(SciPostDc)



            nSamplesPerRecordChannel  = size(SciPostDc.Zv.DemuxerOutput.dcV1, 2);
            nRecords                  = size(SciPreDc.Zv.Epoch, 1);

            OutSci = [];

            OutSci.Zv.Epoch              = SciPreDc.Zv.Epoch;
            OutSci.Zv.QUALITY_BITMASK    = SciPreDc.Zv.QUALITY_BITMASK;
            OutSci.Zv.L2_QUALITY_BITMASK = SciPostDc.Zv.L2_QUALITY_BITMASK;
            OutSci.Zv.QUALITY_FLAG       = SciPreDc.Zv.QUALITY_FLAG;
            OutSci.Zv.DELTA_PLUS_MINUS   = SciPreDc.Zv.DELTA_PLUS_MINUS;
            OutSci.Zv.SYNCHRO_FLAG       = SciPreDc.Zv.SYNCHRO_FLAG;
            OutSci.Zv.SAMPLING_RATE      = SciPreDc.Zv.freqHz;

            % NOTE: Convert aampere --> nano-aampere
            OutSci.Zv.IBIAS1 = SciPostDc.Zv.currentAAmpere(:, 1) * 1e9;
            OutSci.Zv.IBIAS2 = SciPostDc.Zv.currentAAmpere(:, 2) * 1e9;
            OutSci.Zv.IBIAS3 = SciPostDc.Zv.currentAAmpere(:, 3) * 1e9;

%             OutSci.Ga.Datetime  = SciPreDc.Ga.Datetime;
            OutSci.Ga.OBS_ID    = SciPreDc.Ga.OBS_ID;
            OutSci.Ga.SOOP_TYPE = SciPreDc.Ga.SOOP_TYPE;



            C = EJ_library.so.adm.classify_BICAS_L1_L1R_to_L2_DATASET_ID(outputDsi);

            % NOTE: The two cases are different in the indexes they use for
            % OutSciZv.
            if C.isCwf

                % ASSERTIONS
                assert(nSamplesPerRecordChannel == 1, ...
                    'BICAS:proc_sub:Assertion:IllegalArgument', ...
                    ['Number of samples per CDF record is not 1, as expected.', ...
                    ' Bad input CDF?'])
                EJ_library.assert.sizes(...
                    OutSci.Zv.QUALITY_BITMASK, [nRecords, 1], ...
                    OutSci.Zv.QUALITY_FLAG,    [nRecords, 1])

                % Try to pre-allocate to save RAM/speed up.
                tempNaN = nan(nRecords, 3);
                OutSci.Zv.VDC = tempNaN;
                OutSci.Zv.EDC = tempNaN;
                OutSci.Zv.EAC = tempNaN;

                OutSci.Zv.VDC(:,1) = SciPostDc.Zv.DemuxerOutput.dcV1;
                OutSci.Zv.VDC(:,2) = SciPostDc.Zv.DemuxerOutput.dcV2;
                OutSci.Zv.VDC(:,3) = SciPostDc.Zv.DemuxerOutput.dcV3;

                OutSci.Zv.EDC(:,1) = SciPostDc.Zv.DemuxerOutput.dcV12;
                OutSci.Zv.EDC(:,2) = SciPostDc.Zv.DemuxerOutput.dcV13;
                OutSci.Zv.EDC(:,3) = SciPostDc.Zv.DemuxerOutput.dcV23;

                OutSci.Zv.EAC(:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                OutSci.Zv.EAC(:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                OutSci.Zv.EAC(:,3) = SciPostDc.Zv.DemuxerOutput.acV23;

            elseif C.isSwf

                if     C.isLfr
                    SAMPLES_PER_RECORD_CHANNEL = ...
                        EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH;
                elseif C.isTds
                    SAMPLES_PER_RECORD_CHANNEL = ...
                        EJ_library.so.constants.TDS_RSWF_SAMPLES_PER_RECORD;
                else
                    error(...
                        'BICAS:proc_sub:Assertion', ...
                        'Illegal DATASET_ID classification.')
                end

                % ASSERTION
                assert(nSamplesPerRecordChannel == SAMPLES_PER_RECORD_CHANNEL, ...
                    'BICAS:proc_sub:Assertion:IllegalArgument', ...
                    ['Number of samples per CDF record (%i) is not', ...
                    ' %i, as expected. Bad Input CDF?'], ...
                    nSamplesPerRecordChannel, ...
                    SAMPLES_PER_RECORD_CHANNEL)

                % Try to pre-allocate to save RAM/speed up.
                tempNaN = nan(nRecords, nSamplesPerRecordChannel, 3);
                OutSci.Zv.VDC = tempNaN;
                OutSci.Zv.EDC = tempNaN;
                OutSci.Zv.EAC = tempNaN;

                OutSci.Zv.VDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV1;
                OutSci.Zv.VDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV2;
                OutSci.Zv.VDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV3;

                OutSci.Zv.EDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV12;
                OutSci.Zv.EDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV13;
                OutSci.Zv.EDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV23;

                OutSci.Zv.EAC(:,:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                OutSci.Zv.EAC(:,:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                OutSci.Zv.EAC(:,:,3) = SciPostDc.Zv.DemuxerOutput.acV23;

            else
                error('BICAS:proc_sub:Assertion:IllegalArgument', ...
                    'Function can not produce outputDsi=%s.', outputDsi)
            end



            % ASSERTION
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(OutSci.Zv);
            % NOTE: Not really necessary since the list of zVars will be checked
            % against the master CDF?
            % NOTE: Includes zVar "BW" (LFR L2 only).
            EJ_library.assert.struct(OutSci.Zv, {...
                'IBIAS1', 'IBIAS2', 'IBIAS3', 'VDC', 'EDC', 'EAC', 'Epoch', ...
                'QUALITY_BITMASK', 'L2_QUALITY_BITMASK', 'QUALITY_FLAG', ...
                'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'SAMPLING_RATE'}, {})

        end    % process_PostDC_to_TDS_CDF



        % Processing function. Derive PostDC from PreDc, i.e. demux and
        % calibrate data. Function is in large part a wrapper around
        % "calibrate_demux_voltages".
        %
        % NOTE: Public function as opposed to the other demuxing/calibration
        % functions.
        %
        function PostDc = process_calibrate_demux(PreDc, InCurPd, Cal, SETTINGS, L)
            % PROPOSAL: Separate function for calibrating currents.

            tTicToc = tic();

            % ASSERTION
            bicas.proc_sub12.assert_PreDC(PreDc);



            % IMPLEMENTATION NOTE: Only copy fields PreDc-->PostDc which are
            % known to be needed in order to conserve memory (not sure if
            % meaningful).
            PostDc = [];



            %############################
            % DEMUX & CALIBRATE VOLTAGES
            %############################
            PostDc.Zv.DemuxerOutput = bicas.proc_sub12.calibrate_demux_voltages(PreDc, Cal, L);



            %#########################
            % Calibrate bias CURRENTS
            %#########################
            currentSAmpere = bicas.proc_sub12.process_CUR_to_CUR_on_SCI_TIME(...
                PreDc.Zv.Epoch, InCurPd, SETTINGS, L);
            currentTm      = bicas.calib.calibrate_current_sampere_to_TM(currentSAmpere);

            currentAAmpere = nan(size(currentSAmpere));    % Variable to fill/set.
            iCalibLZv      = Cal.get_calibration_time_L(PreDc.Zv.Epoch);
            [iFirstList, iLastList, nSubseq] = EJ_library.utils.split_by_change(iCalibLZv);
            L.logf('info', ...
                ['Calibrating currents -', ...
                ' One sequence of records with identical settings at a time.'])
            for iSubseq = 1:nSubseq
                iFirst = iFirstList(iSubseq);
                iLast  = iLastList(iSubseq);

                iRecords = iFirst:iLast;

                L.logf('info', 'Records %7i-%7i : %s -- %s', ...
                    iFirst, iLast, ...
                    bicas.proc_utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.proc_utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iLast)))

                for iAnt = 1:3
                    %--------------------
                    % CALIBRATE CURRENTS
                    %--------------------
                    currentAAmpere(iRecords, iAnt) = Cal.calibrate_current_TM_to_aampere(...
                        currentTm( iRecords, iAnt), iAnt, iCalibLZv(iRecords));
                end
            end
            PostDc.Zv.currentAAmpere = currentAAmpere;



            % ASSERTION
            bicas.proc_sub12.assert_PostDC(PostDc)

            nRecords = size(PreDc.Zv.Epoch, 1);
            bicas.log_speed_profiling(L, ...
                'bicas.proc_sub12.process_calibrate_demux', tTicToc, ...
                nRecords, 'record')
        end    % process_calibrate_demux



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

            % NOTE: Adds zVar L2_QUALITY_FLAG to PostDc, technically altering the format.
            %   NOTE: Also overwrites voltage with fill values.
            %   PROPOSAL: Treat output PostDc as another format?
            %   PROPOSAL: Initialize empty L2_QUALITY_FLAG when PostDc first created.
            %   PROPOSAL: Keep as is. List as optional field in assert_PostDc
            %
            % PROPOSAL: Abolish test functionality.
            %   PRO: Test functionality can lead to bugs.
            %   PRO: Can use alternative NSO table (using NSO table path override
            %        setting).
            %
            % PROPOSAL: Generalize function be used in L3.
            %   CON: Can not be done since this function is meant to have access
            %        to arbitrary L1/L1R and L2 data to make decisions, although
            %        this is not much used yet.

            % ASSERTION
            bicas.proc_sub12.assert_PreDC(PreDc)
            bicas.proc_sub12.assert_PostDC(PostDc)
            nRecords = EJ_library.assert.sizes(PreDc.Zv.Epoch, [-1]);



            % NOTE: Preallocates and adds zVar to PostDc.
            PostDc.Zv.L2_QUALITY_BITMASK = zeros(nRecords, 1, 'uint16');



            %============================================
            % Find CDF records to remove due to settings
            %============================================
            zvUfvSettings = bicas.proc_sub12.get_UFV_records_from_settings(...
                PreDc.Zv.Epoch, PreDc.Zv.MUX_SET, PreDc.isLfr, SETTINGS, L);

            zvUfv = PreDc.Zv.useFillValues | zvUfvSettings;



            %========================================
            % Take actions based on NSO events table
            %========================================
            testNsoIdsEnabled = SETTINGS.get_fv('PROCESSING.RCS_NSO.TEST_IDS_ENABLED');

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
                    EJ_library.cdf.TT2000_to_UTC_str(NsoTable.evtStartTt2000Array(iGlobalEvent)), ...
                    EJ_library.cdf.TT2000_to_UTC_str(NsoTable.evtStopTt2000Array(iGlobalEvent)), ...
                    eventNsoId);



                %==========================================================
                % TEST FUNCTIONALITY
                % ------------------
                % Optionally translate (selected) TEST NSO IDs into actual
                % NSO IDs.
                %==========================================================
                eventNsoIdTranslated = EJ_library.utils.translate({...
                    {bicas.constants.NSOID.TEST_THRUSTER_FIRING}, ...
                     bicas.constants.NSOID.THRUSTER_FIRING}, ...
                    eventNsoId, eventNsoId);
                if ~testNsoIdsEnabled && ~strcmp(eventNsoId, eventNsoIdTranslated)
                    % CASE:   (1) Not test mode
                    %       & (2) NSO ID was translated (changed).
                    % ==> Original NSO ID was a TEST NSO ID
                    % ==> NSO should be ignored.
                    eventNsoIdTranslated = 'nothing';   % Local constant.
                end
                eventNsoId = eventNsoIdTranslated;
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
                    case bicas.constants.NSOID.TEST_QF0
                        if testNsoIdsEnabled
                            zv_QUALITY_FLAG = min(zv_QUALITY_FLAG, 0, ...
                                'includeNaN');
                        end
                    case bicas.constants.NSOID.TEST_UFV
                        if testNsoIdsEnabled
                            zvUfv = zvUfv | bCdfEventRecordsCa{kCdfEvent};
                        end
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

                    case 'nothing'
                        % CASE: Do nothing.
                        % This case is necessary so that test NSO IDs can be
                        % translated to something harmless when tests are
                        % disabled.

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
            bicas.proc_sub12.log_UFV_records(PreDc.Zv.Epoch, zvUfv, logHeaderStr, L)
            %
            PostDc.Zv.currentAAmpere(zvUfv, :) = NaN;
            %
            fnCa = fieldnames(PostDc.Zv.DemuxerOutput);
            for iFn = 1:numel(fnCa)
                PostDc.Zv.DemuxerOutput.(fnCa{iFn})(zvUfv, :, :) = NaN;
            end



            % ASSERTION
            bicas.proc_sub12.assert_PreDC(PreDc)
            bicas.proc_sub12.assert_PostDC(PostDc)

        end    % process_quality_filter_L2



    end    % methods(Static, Access=public)



    %##############################
    %##############################
    methods(Static, Access=private)
    %##############################
    %##############################



        % Local utility function to shorten & clarify code.
        %
        % ARGUMENTS
        % =========
        % zv1 : zVar-like variabel or empty. Column vector (Nx1) or empty.
        %
        % RETURN VALUE
        % ============
        % zv2 : If zv1 is non-empty, then zv2=zv1.
        %       If zv1 is empty,     then error/mitigate.
        %
        function zv2 = normalize_LFR_zVar_empty(...
                L, settingValue, settingKey, nRecords, zv1, zvName)

            if ~isempty(zv1)
                % Do nothing (except assertion later).
                zv2 = zv1;
            else
                anomalyDescrMsg = sprintf(...
                    'zVar "%s" from the LFR SCI source dataset is empty.', zvName);
                switch(settingValue)
                    case 'USE_FILL_VALUE'
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, 'other', ...
                            anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')

                        L.logf('warning', 'Using fill values for %s.', zvName)
                        zv2 = nan(nRecords, 1);

                    otherwise
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+illegal', ...
                            anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
                end
            end

            EJ_library.assert.sizes(zv2, [NaN])
        end



        % Utility function to shorten code.
        %
        % NOTE: Operates on entire ZvStruct since CALIBRATION_TABLE_INDEX exists
        % for L1R, but not L1.
        function CALIBRATION_TABLE_INDEX = normalize_CALIBRATION_TABLE_INDEX(...
                ZvStruct, nRecords, inputDsi)

            C = EJ_library.so.adm.classify_BICAS_L1_L1R_to_L2_DATASET_ID(inputDsi);

            if C.isL1r
                CALIBRATION_TABLE_INDEX = ZvStruct.CALIBRATION_TABLE_INDEX;
            elseif C.isL1
                CALIBRATION_TABLE_INDEX = nan(nRecords, 2);
            else
                error(...
                    ['Can not normalize CALIBRATION_TABLE_INDEX', ...
                    ' for this DATASET_ID classification.'])
            end

            EJ_library.assert.sizes(CALIBRATION_TABLE_INDEX, [nRecords, 2])
        end



        % Wrapper around EJ_library.so.CURRENT_zv_to_current_interpolate for
        % anomaly handling.
        function sciZv_IBIASx = zv_TC_to_current(...
                curZv_Epoch, curZv_IBIAS_x, sciZv_Epoch, L, SETTINGS)

            % TEST:
%             curZv_Epoch(4:6) = curZv_Epoch(4:6)-(86400e9)
%             [~, iSort] = sort(curZv_Epoch)
%             curZv_Epoch   = curZv_Epoch  (iSort);
%             curZv_IBIAS_x = curZv_IBIAS_x(iSort);

            %====================
            % Calibrate currents
            %====================
            [sciZv_IBIASx, duplicateAnomaly] = EJ_library.so.CURRENT_zv_to_current_interpolate(...
                curZv_Epoch, ...
                curZv_IBIAS_x, ...
                sciZv_Epoch);



            if duplicateAnomaly
                %====================================================
                % Handle anomaly: Non-monotonically increasing Epoch
                %====================================================
                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY');
                anomalyDescriptionMsg = [...
                    'Bias current data contain duplicate settings, with identical timestamps', ...
                    ' and identical bias settings on the same antenna.'];

                switch(settingValue)
                    case 'REMOVE_DUPLICATES'
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, 'other', ...
                            anomalyDescriptionMsg)
                        L.log('warning', ...
                            ['Removed duplicated bias current settings with', ...
                            ' identical timestamps on the same antenna.'])

                    otherwise
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, 'E+illegal', ...
                            anomalyDescriptionMsg, ...
                            'BICAS:proc_sub:SWModeProcessing:DatasetFormat')
                end
            end

        end    % bicas.proc_sub12.zv_TC_to_current



        function assert_PreDC(PreDc)
            EJ_library.assert.struct(PreDc, ...
                {'Zv', 'Ga', 'hasSnapshotFormat', 'isLfr', 'isTdsCwf'}, {});

            EJ_library.assert.struct(PreDc.Zv, ...
                {'Epoch', 'samplesCaTm', 'freqHz', 'nValidSamplesPerRecord', ...
                'iLsf', 'DIFF_GAIN', ...
                'MUX_SET', 'QUALITY_BITMASK', 'QUALITY_FLAG', 'SYNCHRO_FLAG', ...
                'DELTA_PLUS_MINUS', 'CALIBRATION_TABLE_INDEX', 'useFillValues'}, ...
                {'BW'});

            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(PreDc.Zv);

            assert(isa(PreDc.Zv.freqHz, 'double'))
        end



        function assert_PostDC(PostDc)
            EJ_library.assert.struct(PostDc, ...
                {'Zv'}, {});

            EJ_library.assert.struct(PostDc.Zv, ...
                {'DemuxerOutput', 'currentAAmpere'}, {'L2_QUALITY_BITMASK'});

            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(PostDc.Zv);
        end



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

            bicas.proc_utils.assert_zv_Epoch(zvEpoch)
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
            zvUseFillValues = EJ_library.utils.true_with_margin(...
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
                strjoin(EJ_library.str.sprintf_many('%g', muxModesRemove), ', '), ...
                settingMarginKey, ...
                removeMarginSec);
            bicas.proc_sub12.log_UFV_records(zvEpoch, zvUseFillValues, logHeaderStr, L)
        end



        % Log UFV records
        %
        % NOTE: Only logs (including header) if there are records to remove.
        function log_UFV_records(zvEpoch, zvUfv, logHeaderStr, L)
            LL = 'info';    % LL = Log Level

            [i1Array, i2Array] = EJ_library.utils.split_by_false(zvUfv);
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
                    utc1  = EJ_library.cdf.TT2000_to_UTC_str(zvEpoch(iCdfRecord1));
                    utc2  = EJ_library.cdf.TT2000_to_UTC_str(zvEpoch(iCdfRecord2));
                    L.logf(LL, '    Records %7i-%7i, %s -- %s', ...
                        iCdfRecord1, iCdfRecord2, utc1, utc2);
                end
            end

        end



        % Demultiplex and calibrate voltages.
        %
        % NOTE: Can handle arrays of any size if the sizes are
        % consistent.
        %
        function AsrSamplesAVolt = calibrate_demux_voltages(PreDc, Cal, L)
        % PROPOSAL: Sequence of constant settings includes dt (for CWF)
        %   PROBLEM: Not clear how to implement it since it is a property of two records, not one.
        %       PROPOSAL: Use other utility function(s).
        %           PROPOSAL: Function that finds changes in dt.
        %           PROPOSAL: Function that further splits list of index intervals ~on the form iFirstList, iLastList.
        %           PROPOSAL: Write functions such that one can detect suspicious jumps in dt (under some threshold).
        %               PROPOSAL: Different policies/behaviours:
        %                   PROPOSAL: Assertion on expected constant dt.
        %                   PROPOSAL: Always split sequence at dt jumps.
        %                   PROPOSAL: Never  split sequence at dt jumps.
        %                   PROPOSAL: Have threshold on dt when expected constant dt.
        %                       PROPOSAL: Below dt jump threshold, never split sequence
        %                       PROPOSAL: Above dt jump threshold, split sequence
        %                       PROPOSAL: Above dt jump threshold, assert never/give error
        %
        % PROPOSAL: Sequence of constant settings includes constant NaN/non-NaN for CWF.
        %
        % PROPOSAL: Integrate into bicas.demultiplexer (as method).
        % NOTE: Calibration is really separate from the demultiplexer. Demultiplexer only needs to split into
        %       subsequences based on mux mode and latching relay, nothing else.
        %   PROPOSAL: Separate out demultiplexer. Do not call from this function.
        %
        % PROPOSAL: Function for dtSec.
        %     PROPOSAL: Some kind of assertion (assumption of) constant sampling frequency.
        %
        % PROPOSAL: Move the different conversion of CWF/SWF (one/many cell arrays) into the calibration function?!!
        %
        % PROPOSAL: Move processing of one subsequence (one for-loop iteration) into its own function.

            %tTicToc  = tic();

            % ASSERTIONS
            assert(isscalar(PreDc.hasSnapshotFormat))
            assert(iscell(  PreDc.Zv.samplesCaTm))
            EJ_library.assert.vector(PreDc.Zv.samplesCaTm)
            assert(numel(PreDc.Zv.samplesCaTm) == 5)
            bicas.proc_utils.assert_cell_array_comps_have_same_N_rows(PreDc.Zv.samplesCaTm)
            [nRecords, nSamplesPerRecordChannel] = EJ_library.assert.sizes(...
                PreDc.Zv.MUX_SET,        [-1,  1], ...
                PreDc.Zv.DIFF_GAIN,      [-1,  1], ...
                PreDc.Zv.samplesCaTm{1}, [-1, -2]);



            % Pre-allocate
            % ------------
            % IMPLEMENTATION NOTE: Very important for speeding up LFR-SWF which
            % tends to be broken into subsequences of 1 record.
            tempVoltageArray = nan(nRecords, nSamplesPerRecordChannel);
            AsrSamplesAVolt = struct(...
                'dcV1',  tempVoltageArray, ...
                'dcV2',  tempVoltageArray, ...
                'dcV3',  tempVoltageArray, ...
                'dcV12', tempVoltageArray, ...
                'dcV13', tempVoltageArray, ...
                'dcV23', tempVoltageArray, ...
                'acV12', tempVoltageArray, ...
                'acV13', tempVoltageArray, ...
                'acV23', tempVoltageArray);

            dlrUsing12zv = bicas.demultiplexer_latching_relay(PreDc.Zv.Epoch);
            iCalibLZv    = Cal.get_calibration_time_L(        PreDc.Zv.Epoch);
            iCalibHZv    = Cal.get_calibration_time_H(        PreDc.Zv.Epoch);



            %===================================================================
            % (1) Find continuous subsequences of records with identical
            %     settings.
            % (2) Process data separately for each such sequence.
            % NOTE: Just finding continuous subsequences can take a significant
            % amount of time.
            % NOTE: Empirically, this is not useful for real LFR SWF datasets
            % where the LFR sampling frequency changes in every record, meaning
            % that the subsequences are all 1 record long.
            %===================================================================
            [iFirstList, iLastList, nSubseq] = EJ_library.utils.split_by_change(...
                PreDc.Zv.MUX_SET, ...
                PreDc.Zv.DIFF_GAIN, ...
                dlrUsing12zv, ...
                PreDc.Zv.freqHz, ...
                iCalibLZv, ...
                iCalibHZv, ...
                PreDc.Zv.iLsf, ...
                PreDc.Zv.useFillValues, ...
                PreDc.Zv.CALIBRATION_TABLE_INDEX);
            L.logf('info', ...
                ['Calibrating voltages - ', ...
                ' One sequence of records with identical settings at a time.'])

            for iSubseq = 1:nSubseq

                iFirst = iFirstList(iSubseq);
                iLast  = iLastList (iSubseq);

                % Extract SCALAR settings to use for entire subsequence of
                % records.
                % SS = Subsequence (single, constant value valid for entire
                %      subsequence)
                MUX_SET_ss                 = PreDc.Zv.MUX_SET  (              iFirst);
                DIFF_GAIN_ss               = PreDc.Zv.DIFF_GAIN(              iFirst);
                dlrUsing12_ss              = dlrUsing12zv(                    iFirst);
                freqHz_ss                  = PreDc.Zv.freqHz(                 iFirst);
                iCalibL_ss                 = iCalibLZv(                       iFirst);
                iCalibH_ss                 = iCalibHZv(                       iFirst);
                iLsf_ss                    = PreDc.Zv.iLsf(                   iFirst);
                ufv_ss                     = PreDc.Zv.useFillValues(          iFirst);
                CALIBRATION_TABLE_INDEX_ss = PreDc.Zv.CALIBRATION_TABLE_INDEX(iFirst, :);

                % PROPOSAL: Make into "proper" table.
                %   NOTE: Can not use EJ_library.str.assist_print_table since
                %         it requires the entire table to pre-exist.
                %   PROPOSAL: Print after all iterations.
                % NOTE: DIFF_GAIN needs three characters two fit in "NaN".
                L.logf('info', ['Records %8i-%8i : %s -- %s', ...
                    ' MUX_SET=%i; DIFF_GAIN=%-3i; dlrUsing12=%i;', ...
                    ' freqHz=%5g; iCalibL=%i; iCalibH=%i; ufv=%i', ...
                    ' CALIBRATION_TABLE_INDEX=[%i, %i]'], ...
                    iFirst, iLast, ...
                    bicas.proc_utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.proc_utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iLast)), ...
                    MUX_SET_ss, DIFF_GAIN_ss, dlrUsing12_ss, freqHz_ss, ...
                    iCalibL_ss, iCalibH_ss, ufv_ss, ...
                    CALIBRATION_TABLE_INDEX_ss(1), ...
                    CALIBRATION_TABLE_INDEX_ss(2))

                %=======================================
                % DEMULTIPLEXER: FIND ASR-BLTS ROUTINGS
                %=======================================
                % NOTE: Call demultiplexer with no samples. Only for collecting
                % information on which BLTS channels are connected to which
                % ASRs.
                [BltsSrcAsrArray, ~] = bicas.demultiplexer.main(...
                    MUX_SET_ss, dlrUsing12_ss, {[],[],[],[],[]});



                % Extract subsequence of DATA records to "demux".
                ssSamplesTm                = bicas.proc_utils.select_row_range_from_cell_comps(...
                    PreDc.Zv.samplesCaTm, iFirst, iLast);
                % NOTE: "zVariable" (i.e. first index=record) for only the
                % current subsequence.
                ssZvNValidSamplesPerRecord = PreDc.Zv.nValidSamplesPerRecord(iFirst:iLast);
                if PreDc.hasSnapshotFormat
                    % NOTE: Vector of constant numbers (one per snapshot).
                    ssDtSec = 1 ./ PreDc.Zv.freqHz(iFirst:iLast);
                else
                    % NOTE: Scalar (one for entire sequence).
                    ssDtSec = double(...
                        PreDc.Zv.Epoch(iLast) - PreDc.Zv.Epoch(iFirst)) ...
                        / (iLast-iFirst) * 1e-9;   % TEMPORARY?
                end

                biasHighGain = DIFF_GAIN_ss;



                %=====================
                % ITERATE OVER BLTS's
                %=====================
                ssSamplesAVolt = cell(5,1);
                for iBlts = 1:5

                    if strcmp(BltsSrcAsrArray(iBlts).category, 'Unknown')
                        % ==> Calibrated data == NaN.
                        ssSamplesAVolt{iBlts} = nan(size(ssSamplesTm{iBlts}));

                    elseif ismember(BltsSrcAsrArray(iBlts).category, {'GND', '2.5V Ref'})
                        % ==> No calibration.
                        ssSamplesAVolt{iBlts} = ssSamplesTm{iBlts};

                    else
                        assert(BltsSrcAsrArray(iBlts).is_ASR())
                        % ==> Calibrate (unless explicitly stated that should
                        % not)

                        if PreDc.hasSnapshotFormat
                            ssSamplesCaTm = bicas.proc_utils.convert_matrix_to_cell_array_of_vectors(...
                                double(ssSamplesTm{iBlts}), ssZvNValidSamplesPerRecord);
                        else
                            assert(all(ssZvNValidSamplesPerRecord == 1))
                            ssSamplesCaTm = {double(ssSamplesTm{iBlts})};
                        end

                        %######################
                        %######################
                        %  CALIBRATE VOLTAGES
                        %######################
                        %######################
                        % IMPLEMENTATION NOTE: Must explicitly disable
                        % calibration for LFR zVar BW=0
                        % ==> CALIBRATION_TABLE_INDEX(1,:) illegal value.
                        % ==> Can not calibrate.
                        % Therefore uses ufv_ss to disable calibration.
                        % It is thus not enough to overwrite the values later.
                        % This incidentally also potentially speeds up the code.
                        % Ex: LFR SWF 2020-02-25, 2020-02-28.
                        CalSettings = struct();
                        CalSettings.iBlts        = iBlts;
                        CalSettings.BltsSrc      = BltsSrcAsrArray(iBlts);
                        CalSettings.biasHighGain = biasHighGain;
                        CalSettings.iCalibTimeL  = iCalibL_ss;
                        CalSettings.iCalibTimeH  = iCalibH_ss;
                        CalSettings.iLsf         = iLsf_ss;
                        %#######################################################
                        ssSamplesCaAVolt = Cal.calibrate_voltage_all(...
                            ssDtSec, ssSamplesCaTm, ...
                            PreDc.isLfr, PreDc.isTdsCwf, CalSettings, ...
                            CALIBRATION_TABLE_INDEX_ss, ufv_ss);
                        %#######################################################

                        if PreDc.hasSnapshotFormat
                            [ssSamplesAVolt{iBlts}, ~] = bicas.proc_utils.convert_cell_array_of_vectors_to_matrix(...
                                ssSamplesCaAVolt, size(ssSamplesTm{iBlts}, 2));
                        else
                            % NOTE: Must be column array.
                            ssSamplesAVolt{iBlts} = ssSamplesCaAVolt{1};
                        end
                    end
                end    % for iBlts = 1:5

                %====================================
                % DEMULTIPLEXER: DERIVE MISSING ASRs
                %====================================
                [~, SsAsrSamplesAVolt] = bicas.demultiplexer.main(...
                    MUX_SET_ss, dlrUsing12_ss, ssSamplesAVolt);

                % Add demuxed sequence to the to-be complete set of records.
                AsrSamplesAVolt = bicas.proc_utils.set_struct_field_rows(...
                    AsrSamplesAVolt, SsAsrSamplesAVolt, iFirst:iLast);

            end    % for iSubseq = 1:length(iFirstList)



            % NOTE: Assumes no "return" statement.
            %bicas.log_speed_profiling(L, 'bicas.proc_sub12.calibrate_demux_voltages', tTicToc, nRecords, 'record')
            %bicas.log_memory_profiling(L, 'bicas.proc_sub12.calibrate_demux_voltages:end')
        end    % calibrate_demux_voltages



    end    % methods(Static, Access=private)

end