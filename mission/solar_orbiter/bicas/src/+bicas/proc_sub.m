%
% Class that collects "processing functions" as public static methods.
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
% UFV : Use Fill Values
%
%
% SOME INTERMEDIATE PROCESSING DATA FORMATS
% =========================================
% - PreDC = Pre-Demuxing-Calibration Data
%       Generic data format that can represent all forms of input datasets
%       before demuxing and calibration. Can use an arbitrary number of samples
%       per record. Some variables are therefore not used in CWF output
%       datasets.
% - PostDC = Post-Demuxing-Calibration Data
%       Like PreDC but with additional fields. Tries to capture a superset of
%       the information that goes into any dataset produced by BICAS, and the
%       exact set of variables that goes into the output datasets.
% 
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
classdef proc_sub
%#######################################################################################################################
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
% TODO: add_UFV_records_from_settings should know whether output is L2 or not.
%
% PROPOSAL: Merge (most of)
%   process_PostDC_to_LFR and
%   process_PostDC_to_TDS into one function that shares most of the functionality.
%
% PROPOSAL:   process_calibrate_demux
%           & calibrate_demux_voltages
%           should only accept the needed zVars and variables.
%   NOTE: Needs some way of packaging/extracting only the relevant zVars/fields
%         from struct.
%
%#######################################################################################################################

    %#############################
    %#############################
    methods(Static, Access=public)
    %#############################
    %#############################
        
        
        
        % Processing function
        function HkSciTime = process_HK_to_HK_on_SCI_TIME(InSci, InHk, SETTINGS, L)
        
            % ASSERTIONS
            EJ_library.assert.struct(InSci, {'Zv', 'Ga'}, {})
            EJ_library.assert.struct(InHk,  {'Zv', 'Ga'}, {})

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
                hkEpoch = bicas.proc_utils.ACQUISITION_TIME_to_tt2000(...
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
                    bicas.proc_utils.ACQUISITION_TIME_to_tt2000(...
                        InHk.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            if isfield(InSci.Zv, 'ACQUISITION_TIME') && ~isempty(InSci.Zv.ACQUISITION_TIME)
                TimeVars.SCI_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc_utils.ACQUISITION_TIME_to_tt2000(...
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
                % is not monotonically increasing (in fact, it is completely strange).
                error('HK timestamps do not increase monotonically (USE_ZV_ACQUISITION_TIME_HK=%g).', USE_ZV_ACQUISITION_TIME_HK)
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
                
                [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.HK.TIME_NOT_SUPERSET_OF_SCI_POLICY');
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
            end
            if ~EJ_library.utils.ranges_intersect(InSci.Zv.Epoch, hkEpoch)
                
                % NOTE: "WARNING" (rather than error) only makes sense if it is
                % possible to later meaningfully permit non-intersection.
                [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY');
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    'SCI and HK time ranges do not overlap in time.', 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
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
            EJ_library.assert.struct(InCur, {'Zv', 'Ga'}, {})
            
            
            
            %===================================================================
            % CDF ASSERTION: CURRENT data begins before SCI data (i.e. there is
            % enough CURRENT data).
            %===================================================================
            if ~(min(InCur.Zv.Epoch) <= min(sciEpoch))
                curRelativeSec    = 1e-9 * (min(InCur.Zv.Epoch) - min(sciEpoch));
                sciEpochUtcStr    = EJ_library.cdf.tt2000_to_UTC_str(min(sciEpoch));
                curEpochMinUtcStr = EJ_library.cdf.tt2000_to_UTC_str(min(InCur.Zv.Epoch));
                
                [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY');
                
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
            % NOTE: bicas.proc_sub.interpolate_current checks (and handles) that
            % Epoch increases monotonically, but only for each antenna
            % separately (which does not capture all cases).
            % Ex: Timestamps, iAntenna = mod(iRecord,3): 1,2,3,5,4,6
            %       ==> Monotonically increasing sequences for each antenna
            %           separately, but not even increasing when combined.
            %====================================================================
            if ~issorted(InCur.Zv.Epoch)
                error('CURRENT timestamps do not increase (all antennas combined).')
            end
            
            % NOTE: bicas.proc_sub.interpolate_current checks that Epoch
            % increases monotonically.
            currentNanoSAmpere = [];
            currentNanoSAmpere(:,1) = bicas.proc_sub.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_1, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,2) = bicas.proc_sub.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_2, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,3) = bicas.proc_sub.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_3, sciEpoch, L, SETTINGS);
            
            currentSAmpere = 1e-9 * currentNanoSAmpere;
        end
        
        
        
        % Processing function. Only "normalizes" data to account for technically
        % illegal input LFR datasets. This should try to:
        % ** modify L1 to look like L1R
        % ** mitigate historical bugs (in the input datasets)
        % ** mitigate for not yet implemented features (in input datasets)
        %
        function InSciNorm = process_LFR_normalize(InSci, inSciDsi, SETTINGS, L)
            
            % Default behaviour: Copy values, except for values which are
            % modified later
            InSciNorm = InSci;
            
            nRecords = EJ_library.assert.sizes(InSci.Zv.Epoch, [-1]);
            
            
            
            %===================================
            % Normalize CALIBRATION_TABLE_INDEX
            %===================================
            InSciNorm.Zv.CALIBRATION_TABLE_INDEX = bicas.proc_sub.normalize_CALIBRATION_TABLE_INDEX(...
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
                    error('BICAS:process_LFR_normalize:DatasetFormat', 'Input dataset has both zVar SYNCHRO_FLAG and TIME_SYNCHRO_FLAG.')
                end
            else
                error('BICAS:process_LFR_normalize:DatasetFormat', 'Input dataset does not have zVar SYNCHRO_FLAG as expected.')
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
            
            InSciNorm.Zv.QUALITY_BITMASK = bicas.proc_sub.normalize_LFR_zVar_empty(...
                L, settingValue, settingKey, nRecords, ...
                InSci.Zv.QUALITY_BITMASK, 'QUALITY_BITMASK');
            
            InSciNorm.Zv.QUALITY_FLAG    = bicas.proc_sub.normalize_LFR_zVar_empty(...
                L, settingValue, settingKey, nRecords, ...
                InSci.Zv.QUALITY_FLAG,    'QUALITY_FLAG');
            
            % ASSERTIONS
            EJ_library.assert.sizes(...
                InSciNorm.Zv.QUALITY_BITMASK, [nRecords, 1], ...
                InSciNorm.Zv.QUALITY_FLAG,    [nRecords, 1])

        end    % process_LFR_normalize
        
        
        
        % Processing function. Convert LFR CDF data to PreDC.
        %
        % IMPLEMENTATION NOTE: Does not modify InSci in an attempt to save RAM
        % (should help MATLAB's optimization). Unclear if actually works.
        %
        function PreDc = process_LFR_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
            %
            % PROBLEM: Hard-coded CDF data types (MATLAB classes).
            % MINOR PROBLEM: Still does not handle LFR zVar TYPE for determining
            % "virtual snapshot" length. Should only be relevant for
            % V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF (not V02) which should expire.
            
            % ASSERTIONS
            EJ_library.assert.struct(InSci,     {'Zv', 'Ga'}, {})
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
            
            % CDF ASSERTION
            assert(issorted(InSci.Zv.Epoch, 'strictascend'), ...
                'Voltage (science) dataset timestamps do not increase.')
            
            
            
            nRecords = size(InSci.Zv.Epoch, 1);
            C = EJ_library.so.adm.classify_DATASET_ID(inSciDsi);
            


            %============
            % Set iLsfZv
            %============
            if     C.isLfrSbm1   iLsfZv = ones(nRecords, 1) * 2;   % Always value "2" (F1, "FREQ = 1").
            elseif C.isLfrSbm2   iLsfZv = ones(nRecords, 1) * 3;   % Always value "3" (F2, "FREQ = 2").
            else                 iLsfZv = InSci.Zv.FREQ + 1;
                % NOTE: Translates from LFR's FREQ values (0=F0 etc) to LSF
                % index values (1=F0) used in loaded RCT data structs.
            end
            EJ_library.assert.sizes(iLsfZv, [nRecords, 1])



            % NOTE: Needed also for 1 SPR.
            zvFreqHz = EJ_library.so.get_LFR_frequency( iLsfZv );

            % Obtain the relevant values (one per record) from zVariables R0,
            % R1, R2, and the virtual "R3".
            zv_Rx = EJ_library.so.get_LFR_Rx(...
                InSci.Zv.R0, ...
                InSci.Zv.R1, ...
                InSci.Zv.R2, ...
                iLsfZv );

            
            
            % IMPLEMENTATION NOTE: E & V must be floating-point so that values
            % can be set to NaN.
            % 
            % Switch last two indices of E.
            % ==> index 2 = "snapshot" sample index, including for CWF
            %               (sample/record, "snapshots" consisting of 1 sample).
            %     index 3 = E1/E2 component
            %               NOTE: 1/2=index into array; these are diffs but not
            %               equivalent to any particular diffs).
            E = single(permute(InSci.Zv.E, [1,3,2]));
            
            % ASSERTIONS
            nCdfSamplesPerRecord = EJ_library.assert.sizes(InSci.Zv.V, [nRecords, -1], E, [nRecords, -1, 2]);
            if C.isLfrSurvSwf   assert(nCdfSamplesPerRecord == EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH)
            else                assert(nCdfSamplesPerRecord == 1)
            end



            PreDc = [];
            
            PreDc.Zv.samplesCaTm    = cell(5,1);
            PreDc.Zv.samplesCaTm{1} = single(InSci.Zv.V);
            PreDc.Zv.samplesCaTm{2} = bicas.proc_utils.filter_rows( E(:,:,1), zv_Rx==0 );    % Copy values, except when zvRx==0 (==>NaN).
            PreDc.Zv.samplesCaTm{3} = bicas.proc_utils.filter_rows( E(:,:,2), zv_Rx==0 );
            PreDc.Zv.samplesCaTm{4} = bicas.proc_utils.filter_rows( E(:,:,1), zv_Rx==1 );
            PreDc.Zv.samplesCaTm{5} = bicas.proc_utils.filter_rows( E(:,:,2), zv_Rx==1 );
            
            PreDc.Zv.Epoch                   = InSci.Zv.Epoch;
            PreDc.Zv.DELTA_PLUS_MINUS        = bicas.proc_utils.derive_DELTA_PLUS_MINUS(zvFreqHz, nCdfSamplesPerRecord);            
            PreDc.Zv.freqHz                  = zvFreqHz;
            PreDc.Zv.nValidSamplesPerRecord  = ones(nRecords, 1) * nCdfSamplesPerRecord;
            PreDc.Zv.BW                      = InSci.Zv.BW;
            PreDc.Zv.useFillValues           = ~logical(InSci.Zv.BW);
            PreDc.Zv.DIFF_GAIN               = HkSciTime.DIFF_GAIN;
            PreDc.Zv.iLsf                    = iLsfZv;
            
            PreDc.Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
            PreDc.Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;
            
            PreDc.Zv.QUALITY_BITMASK         = InSci.Zv.QUALITY_BITMASK;
            PreDc.Zv.QUALITY_FLAG = min(...
                InSci.Zv.QUALITY_FLAG, ...
                SETTINGS.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX'), 'includeNaN');
            


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
                    error('BICAS:proc_sub:ConfigurationBug', 'Illegal settings value %s="%s"', key, value)
            end

            
            
            PreDc.hasSnapshotFormat = C.isLfrSurvSwf;
            PreDc.isLfr             = true;
            PreDc.isTdsCwf          = false;
            

            
            % ASSERTIONS
            bicas.proc_sub.assert_PreDC(PreDc)
            
        end    % process_LFR_to_PreDC
        
        
        
        % Processing function. Only "normalizes" data to account for technically
        % illegal input TDS datasets. This should try to:
        % ** modify L1 to look like L1R
        % ** mitigate historical bugs (in the input datasets)
        % ** mitigate for not yet implemented features (in input datasets)
        %
        function InSciNorm = process_TDS_normalize(InSci, inSciDsi, SETTINGS, L)
            
            % Default behaviour: Copy values, except for values which are
            % modified later
            InSciNorm = InSci;
            
            nRecords = EJ_library.assert.sizes(InSci.Zv.Epoch, [-1]);
            
            C = EJ_library.so.adm.classify_DATASET_ID(inSciDsi);
            
            
            
            %===================================
            % Normalize CALIBRATION_TABLE_INDEX
            %===================================
            InSciNorm.Zv.CALIBRATION_TABLE_INDEX = bicas.proc_sub.normalize_CALIBRATION_TABLE_INDEX(...
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
            [InSci.Zv, fnChangeList] = EJ_library.utils.normalize_struct_fieldnames(InSci.Zv, ...
                {{{'TIME_SYNCHRO_FLAG', 'SYNCHRO_FLAG'}, 'SYNCHRO_FLAG'}}, ...
                'Assert one matching candidate');
            
            bicas.proc_sub.handle_zv_name_change(...
                fnChangeList, inSciDsi, SETTINGS, L, ...
                'SYNCHRO_FLAG', 'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')

            
            
            %=========================
            % Normalize SAMPLING_RATE
            %=========================
            if any(InSci.Zv.SAMPLING_RATE == 255)
                [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY');
                anomalyDescrMsg = 'Finds illegal, stated sampling frequency 255 in TDS L1/L1R LFM-RSWF dataset.';
                
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
                                settingValue, settingKey, 'E+W+illegal', anomalyDescrMsg, 'BICAS:process_TDS_normalize:DatasetFormat')
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
                SAMPS_PER_CH_MIN_VALID    = 2^10;
                SAMPS_PER_CH_MAX_VALID    = 2^15;

                zv_SAMPS_PER_CH_corrected = round(2.^round(log2(double(InSci.Zv.SAMPS_PER_CH))));
                zv_SAMPS_PER_CH_corrected = cast(zv_SAMPS_PER_CH_corrected, class(InSci.Zv.SAMPS_PER_CH));
                zv_SAMPS_PER_CH_corrected = max( zv_SAMPS_PER_CH_corrected, SAMPS_PER_CH_MIN_VALID);
                zv_SAMPS_PER_CH_corrected = min( zv_SAMPS_PER_CH_corrected, SAMPS_PER_CH_MAX_VALID);
                
                if any(zv_SAMPS_PER_CH_corrected ~= InSci.Zv.SAMPS_PER_CH)
                    % CASE: SAMPS_PER_CH has at least one illegal value
                    
                    SAMPS_PER_CH_badValues = unique(InSci.Zv.SAMPS_PER_CH(zv_SAMPS_PER_CH_corrected ~= InSci.Zv.SAMPS_PER_CH));
                    
                    badValuesDisplayStr = strjoin(arrayfun(...
                        @(n) sprintf('%i', n), SAMPS_PER_CH_badValues, 'uni', false), ', ');
                    anomalyDescrMsg = sprintf(...
                        'TDS LFM RSWF zVar SAMPS_PER_CH contains unexpected value(s) which are not on the form 2^n and in the interval %.0f to %.0f: %s', ...
                        SAMPS_PER_CH_MIN_VALID, ...
                        SAMPS_PER_CH_MAX_VALID, ...
                        badValuesDisplayStr);
                    
                    [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY');
                    switch(settingValue)
                        case 'ROUND'
                            bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', ...
                                anomalyDescrMsg, 'BICAS:proc_sub:process_TDS_normalize:Assertion:DatasetFormat')
                            L.log('warning', ...
                                ['Replacing TDS RSWF zVar SAMPS_PER_CH values with values, rounded to valid', ...
                                ' values due to setting PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY.'])
                            
                            InSciNorm.Zv.SAMPS_PER_CH = zv_SAMPS_PER_CH_corrected;
                            
                        otherwise
                            bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                                anomalyDescrMsg, 'BICAS:proc_sub:process_TDS_normalize:Assertion:DatasetFormat')

                    end    % switch
                end    % if
            end    % if
            
        end    % process_TDS_normalize
        
        
        
        % Processing function. Convert TDS CDF data (PDs) to PreDC.
        function PreDc = process_TDS_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
        %
        % BUG?: Does not use CHANNEL_STATUS_INFO.
        % NOTE: BIAS output datasets do not have a variable for the length of
        % snapshots. Need to use NaN/fill value.

            % ASSERTIONS
            EJ_library.assert.struct(InSci,     {'Zv', 'Ga'}, {})
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})

            % CDF ASSERTION
            if ~issorted(InSci.Zv.Epoch, 'strictascend')
                error('Voltage timestamps do not increase (all antennas combined).')
            end

            C = EJ_library.so.adm.classify_DATASET_ID(inSciDsi);

            
            
            nRecords                  = size(InSci.Zv.Epoch, 1);
            % Number of samples in the zVariable, not necessarily actual data.
            nCdfMaxSamplesPerSnapshot = size(InSci.Zv.WAVEFORM_DATA, 3);

            
            
            % Why convert to double? To avoid precision problems when doing math
            % with other variables?
            freqHzZv = double(InSci.Zv.SAMPLING_RATE);
            
            
            
            PreDc = [];
            
            PreDc.Zv.Epoch                   = InSci.Zv.Epoch;
            PreDc.Zv.DELTA_PLUS_MINUS        = bicas.proc_utils.derive_DELTA_PLUS_MINUS(freqHzZv, nCdfMaxSamplesPerSnapshot);
            PreDc.Zv.freqHz                  = freqHzZv;
            PreDc.Zv.QUALITY_BITMASK         = InSci.Zv.QUALITY_BITMASK;
            PreDc.Zv.QUALITY_FLAG = min(...
                InSci.Zv.QUALITY_FLAG, ...
                SETTINGS.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX'), 'includeNaN');
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
                % samples/record * 8 records.  /Unknown source Unclear what
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
            assert(all(PreDc.Zv.nValidSamplesPerRecord <= nCdfMaxSamplesPerSnapshot), ...
                'BICAS:proc_sub:process_TDS_to_PreDC:Assertion:DatasetFormat', ...
                ['Dataset indicates that the number of valid samples per CDF', ...
                ' record (max(PreDc.Zv.nValidSamplesPerRecord)=%i) is', ...
                ' NOT fewer than the number of indices per CDF record', ...
                ' (nCdfMaxSamplesPerSnapshot=%i).'], ...
                max(PreDc.Zv.nValidSamplesPerRecord), ...
                nCdfMaxSamplesPerSnapshot)
            


            %==========================
            % Set PreDc.Zv.samplesCaTm
            %==========================
            % CDF ASSERTION
            if     C.isL1R   WAVEFORM_DATA_nChannels = 3;
            elseif C.isL1    WAVEFORM_DATA_nChannels = 8;
            end
            % NOTE: NOT using method EJ_library.assert.sizes directly in
            % order to produce customized error message instead.
            assert(...
                EJ_library.utils.sizes(...
                    InSci.Zv.WAVEFORM_DATA, ...
                    [nRecords, WAVEFORM_DATA_nChannels, nCdfMaxSamplesPerSnapshot]), ...
                'BICAS:proc_sub:process_TDS_to_PreDC:Assertion:DatasetFormat', ...
                'TDS zVar WAVEFORM_DATA has an unexpected size.')
            modif_WAVEFORM_DATA = double(permute(InSci.Zv.WAVEFORM_DATA, [1,3,2]));
            
            PreDc.Zv.samplesCaTm    = cell(5,1);
            PreDc.Zv.samplesCaTm{1} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,1), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{2} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,2), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{3} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,3), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{4} = nan(nRecords, nCdfMaxSamplesPerSnapshot);
            PreDc.Zv.samplesCaTm{5} = nan(nRecords, nCdfMaxSamplesPerSnapshot);

            
            
            PreDc.isLfr             = false;
            PreDc.isTdsCwf          = C.isTdsCwf;
            PreDc.hasSnapshotFormat = C.isTdsRswf;
            % Only set because the code shared with LFR requires it.
            PreDc.Zv.iLsf           = nan(nRecords, 1);



            % ASSERTIONS
            bicas.proc_sub.assert_PreDC(PreDc)
            
        end    % process_TDS_to_PreDC
        


        function [OutSciZv] = process_PostDC_to_LFR(SciPostDc, outputDsi, L)
            OutSciZv    = bicas.proc_sub.process_PostDC_to_LFR_TDS_main(SciPostDc, outputDsi, L);
            OutSciZv.BW = SciPostDc.Zv.BW;
        end



        % Processing function. Convert PostDC to either
        % (1) a TDS dataset, or
        % (2) almost to an LFR dataset (the rest is done in a wrapper).
        function [OutSciZv] = process_PostDC_to_LFR_TDS_main(SciPostDc, outputDsi, L)

            % ASSERTIONS
            bicas.proc_sub.assert_PostDC(SciPostDc)



            nSamplesPerRecordChannel  = size(SciPostDc.Zv.DemuxerOutput.dcV1, 2);
            nRecords                  = size(SciPostDc.Zv.Epoch, 1);

            OutSciZv = [];
            
            OutSciZv.Epoch            = SciPostDc.Zv.Epoch;
            OutSciZv.QUALITY_BITMASK  = SciPostDc.Zv.QUALITY_BITMASK;
            OutSciZv.QUALITY_FLAG     = SciPostDc.Zv.QUALITY_FLAG;
            OutSciZv.DELTA_PLUS_MINUS = SciPostDc.Zv.DELTA_PLUS_MINUS;
            OutSciZv.SYNCHRO_FLAG     = SciPostDc.Zv.SYNCHRO_FLAG;
            OutSciZv.SAMPLING_RATE    = SciPostDc.Zv.freqHz;

            % NOTE: Convert aampere --> nano-aampere
            OutSciZv.IBIAS1           = SciPostDc.Zv.currentAAmpere(:, 1) * 1e9;
            OutSciZv.IBIAS2           = SciPostDc.Zv.currentAAmpere(:, 2) * 1e9;
            OutSciZv.IBIAS3           = SciPostDc.Zv.currentAAmpere(:, 3) * 1e9;
            
            
            
            C = EJ_library.so.adm.classify_DATASET_ID(outputDsi);

            EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH;
            EJ_library.so.constants.TDS_RSWF_SAMPLES_PER_RECORD;
            
            % NOTE: The two cases are different in the indexes they use for
            % OutSciZv.
            if C.isCwf
                
                % ASSERTIONS
                assert(nSamplesPerRecordChannel == 1, ...
                    'BICAS:proc_sub:Assertion:IllegalArgument', ...
                    'Number of samples per CDF record is not 1, as expected. Bad input CDF?')
                EJ_library.assert.sizes(...
                    OutSciZv.QUALITY_BITMASK, [nRecords, 1], ...
                    OutSciZv.QUALITY_FLAG,    [nRecords, 1])
                
                % Try to pre-allocate to save RAM/speed up.
                OutSciZv.VDC = nan(nRecords, 3);
                OutSciZv.EDC = nan(nRecords, 3);
                OutSciZv.EAC = nan(nRecords, 3);
                
                OutSciZv.VDC(:,1) = SciPostDc.Zv.DemuxerOutput.dcV1;
                OutSciZv.VDC(:,2) = SciPostDc.Zv.DemuxerOutput.dcV2;
                OutSciZv.VDC(:,3) = SciPostDc.Zv.DemuxerOutput.dcV3;
                
                OutSciZv.EDC(:,1) = SciPostDc.Zv.DemuxerOutput.dcV12;
                OutSciZv.EDC(:,2) = SciPostDc.Zv.DemuxerOutput.dcV13;
                OutSciZv.EDC(:,3) = SciPostDc.Zv.DemuxerOutput.dcV23;
                
                OutSciZv.EAC(:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                OutSciZv.EAC(:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                OutSciZv.EAC(:,3) = SciPostDc.Zv.DemuxerOutput.acV23;
                
            elseif C.isSwf
                
                if     C.isLfr   SAMPLES_PER_RECORD_CHANNEL = EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH;
                elseif C.isTds   SAMPLES_PER_RECORD_CHANNEL = EJ_library.so.constants.TDS_RSWF_SAMPLES_PER_RECORD;
                else             error('BICAS:proc_sub:Assertion', 'Illegal DATASET_ID classification.')
                end
                
                % ASSERTION
                assert(nSamplesPerRecordChannel == SAMPLES_PER_RECORD_CHANNEL, ...
                    'BICAS:proc_sub:Assertion:IllegalArgument', ...
                    'Number of samples per CDF record (%i) is not %i, as expected. Bad Input CDF?', ...
                    nSamplesPerRecordChannel, ...
                    SAMPLES_PER_RECORD_CHANNEL)
                
                % Try to pre-allocate to save RAM/speed up.
                temp = nan(nRecords, nSamplesPerRecordChannel, 3);
                OutSciZv.VDC = temp;
                OutSciZv.EDC = temp;
                OutSciZv.EAC = temp;
                
                OutSciZv.VDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV1;
                OutSciZv.VDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV2;
                OutSciZv.VDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV3;
                
                OutSciZv.EDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV12;
                OutSciZv.EDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV13;
                OutSciZv.EDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV23;
                
                OutSciZv.EAC(:,:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                OutSciZv.EAC(:,:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                OutSciZv.EAC(:,:,3) = SciPostDc.Zv.DemuxerOutput.acV23;
                
            else
                error('BICAS:proc_sub:Assertion:IllegalArgument', ...
                    'Function can not produce outputDsi=%s.', outputDsi)
            end
            
            
            
            % ASSERTION
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(OutSciZv);
            % NOTE: Not really necessary since the list of zVars will be checked
            % against the master CDF?
            % NOTE: Includes zVar "BW" (LFR L2 only).
            EJ_library.assert.struct(OutSciZv, {...
                'IBIAS1', 'IBIAS2', 'IBIAS3', 'VDC', 'EDC', 'EAC', 'Epoch', ...
                'QUALITY_BITMASK', 'QUALITY_FLAG', ...
                'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'SAMPLING_RATE'}, {})
            
        end    % process_PostDC_to_LFR_TDS_main
        
        
        
        % Processing function. Converts PreDC to PostDC, i.e. demux and
        % calibrate data. Function is in large part a wrapper around
        % "calibrate_demux_voltages".
        %
        % NOTE: Public function as opposed to the other demuxing/calibration
        % functions.
        %
        function PostDc = process_calibrate_demux(PreDc, InCurPd, Cal, SETTINGS, L)
            
            tTicToc = tic();

            % ASSERTION
            bicas.proc_sub.assert_PreDC(PreDc);
            
            
            
            % IMPLEMENTATION NOTE: Only copy fields PreDc-->PostDc which are
            % known to be needed in order to conserve memory.
            PostDc = [];
            
            % Copy relevant zVars.
            PostDc.Zv.Epoch            = PreDc.Zv.Epoch;
            PostDc.Zv.QUALITY_BITMASK  = PreDc.Zv.QUALITY_BITMASK;
            PostDc.Zv.QUALITY_FLAG     = PreDc.Zv.QUALITY_FLAG;
            PostDc.Zv.DELTA_PLUS_MINUS = PreDc.Zv.DELTA_PLUS_MINUS;
            PostDc.Zv.SYNCHRO_FLAG     = PreDc.Zv.SYNCHRO_FLAG;
            PostDc.Zv.freqHz           = PreDc.Zv.freqHz;
            if isfield(PreDc.Zv, 'BW')
                PostDc.Zv.BW               = PreDc.Zv.BW;
            end
            


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DEMUX & CALIBRATE VOLTAGES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            PostDc.Zv.DemuxerOutput = bicas.proc_sub.calibrate_demux_voltages(PreDc, Cal, L);
            


            %=========================
            % Calibrate bias CURRENTS
            %=========================
            currentSAmpere = bicas.proc_sub.process_CUR_to_CUR_on_SCI_TIME(PreDc.Zv.Epoch, InCurPd, SETTINGS, L);
            currentTm      = bicas.calib.calibrate_current_sampere_to_TM(currentSAmpere);
            
            currentAAmpere = nan(size(currentSAmpere));    % Variable to fill/set.
            iCalibLZv      = Cal.get_calibration_time_L(PreDc.Zv.Epoch);
            [iFirstList, iLastList, nSubseq] = EJ_library.utils.split_by_change(iCalibLZv);
            L.logf('info', 'Calibrating currents - One sequence of records with identical settings at a time.')
            for iSubseq = 1:nSubseq
                iFirst = iFirstList(iSubseq);
                iLast  = iLastList(iSubseq);
                
                iRecords = iFirst:iLast;
                
                L.logf('info', 'Records %7i-%7i : %s -- %s', ...
                    iFirst, iLast, ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iLast)))
                
                for iAnt = 1:3
                    %%%%%%%%%%%%%%%%%%%%%
                    % CALIBRATE CURRENTS
                    %%%%%%%%%%%%%%%%%%%%%
                    currentAAmpere(iRecords, iAnt) = Cal.calibrate_current_TM_to_aampere(...
                        currentTm( iRecords, iAnt), iAnt, iCalibLZv(iRecords));
                end
            end
            PostDc.Zv.currentAAmpere = currentAAmpere;
            
            
            
            % ASSERTION
            bicas.proc_sub.assert_PostDC(PostDc)
            
            nRecords = size(PreDc.Zv.Epoch, 1);
            bicas.log_speed_profiling(L, 'bicas.proc_sub.process_calibrate_demux', tTicToc, nRecords, 'record')
        end    % process_calibrate_demux
        
        
        
        % Processing function
        %
        % Overwrite selected data in selected CDF records with fill values/NaN.
        function PostDc = process_PostDc_filter(PreDc, PostDc, SETTINGS, L)
            
            %============================================
            % Find CDF records to remove due to settings
            %============================================
            zvUfvSettings = bicas.proc_sub.get_UFV_records_from_settings(...
                PreDc.Zv.Epoch, PreDc.Zv.MUX_SET, PreDc.isLfr, SETTINGS, L);
            
            zvUfvFinal = PreDc.Zv.useFillValues | zvUfvSettings;

            % Log
            logHeaderStr = sprintf(...
                ['All interval(s) of CDF records for which data should be set', ...
                ' to fill values (i.e. removed), regardless of reason.\n']);
            bicas.proc_sub.log_UFV_records(PreDc.Zv.Epoch, zvUfvFinal, logHeaderStr, L)



            %==================================================================
            % Set CURRENTS and VOLTAGES to NaN based on PreDc.Zv.useFillValues
            %==================================================================
            PostDc.Zv.currentAAmpere(zvUfvFinal, :) = NaN;
            %
            fnCa = fieldnames(PostDc.Zv.DemuxerOutput);
            for iFn = 1:numel(fnCa)
                PostDc.Zv.DemuxerOutput.(fnCa{iFn})(zvUfvFinal, :, :) = NaN;
            end
            
        end


        
        % Wrapper around bicas.proc_sub.handle_struct_name_change to be used
        % locally.
        % NOTE: Also used in bicas.proc.process_L3. Therefore public.
        %
        % ARGUMENTS
        % =========
        % inSciDsi : Input SCI DATASET_ID which contains the zVariable.
        % varargin : Passed on to bicas.handle_struct_name_change as its
        %            varargin.
        %
        function handle_zv_name_change(fnChangeList, inSciDsi, SETTINGS, L, varargin)
            anomalyDescrMsgFunc = @(oldFieldname, newFieldname) (sprintf(...
                'Input dataset DATASET_ID=%s uses an alternative but illegal(?) zVariable name "%s" instead of "%s".', ...
                inSciDsi, oldFieldname, newFieldname));
            
            bicas.handle_struct_name_change(fnChangeList, SETTINGS, L, anomalyDescrMsgFunc, varargin{:})
        end
        
        
        
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
        function zv2 = normalize_LFR_zVar_empty(L, settingValue, settingKey, nRecords, zv1, zvName)
            
            if ~isempty(zv1)
                % Do nothing (except assertion later).
                zv2 = zv1;
            else
                anomalyDescrMsg = sprintf('zVar "%s" from the LFR SCI source dataset is empty.', zvName);
                switch(settingValue)
                    case 'USE_FILL_VALUE'
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', ...
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
        function CALIBRATION_TABLE_INDEX = normalize_CALIBRATION_TABLE_INDEX(ZvStruct, nRecords, inputDsi)
            
            C = EJ_library.so.adm.classify_DATASET_ID(inputDsi);
            
            if C.isL1R
                CALIBRATION_TABLE_INDEX = ZvStruct.CALIBRATION_TABLE_INDEX;
            elseif C.isL1
                CALIBRATION_TABLE_INDEX = nan(nRecords, 2);
            else
                error('Can not normalize CALIBRATION_TABLE_INDEX for this DATASET_ID classification.')
            end
            
            EJ_library.assert.sizes(CALIBRATION_TABLE_INDEX, [nRecords, 2])
        end

        
        
        % Wrapper around EJ_library.so.CURRENT_zv_to_current_interpolate for
        % anomaly handling.
        function sciZv_IBIASx = zv_TC_to_current(curZv_Epoch, curZv_IBIAS_x, sciZv_Epoch, L, SETTINGS)
            
            %====================
            % Calibrate currents
            %====================
            [sciZv_IBIASx, duplicateAnomaly] = EJ_library.so.CURRENT_zv_to_current_interpolate(...
                double(curZv_Epoch), ...
                curZv_IBIAS_x, ...
                sciZv_Epoch);
            
            
            
            if duplicateAnomaly
                %====================================================
                % Handle anomaly: Non-monotonically increasing Epoch
                %====================================================
                [settingValue, settingKey] = SETTINGS.get_fv('INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY');
                anomalyDescriptionMsg = [...
                    'Bias current data contain duplicate settings, with identical timestamps', ...
                    ' and identical bias settings on the same antenna.'];
                
                switch(settingValue)
                    case 'REMOVE_DUPLICATES'
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, 'other', ...
                            anomalyDescriptionMsg)
                        L.log('warning', ...
                            'Removed duplicated bias current settings with identical timestamps on the same antenna.')

                    otherwise
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, 'E+illegal', ...
                            anomalyDescriptionMsg, 'BICAS:proc_sub:SWModeProcessing:DatasetFormat')
                end
            end
            
        end    % bicas.proc_sub.zv_TC_to_current
        
        
        
        function assert_PreDC(PreDc)
            EJ_library.assert.struct(PreDc, ...
                {'Zv', 'hasSnapshotFormat', 'isLfr', 'isTdsCwf'}, {});
            
            EJ_library.assert.struct(PreDc.Zv, ...
                {'Epoch', 'samplesCaTm', 'freqHz', 'nValidSamplesPerRecord', 'iLsf', 'DIFF_GAIN', ...
                'MUX_SET', 'QUALITY_BITMASK', 'QUALITY_FLAG', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', ...
                'CALIBRATION_TABLE_INDEX', 'useFillValues'}, ...
                {'BW'});
            
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(PreDc.Zv);

            assert(isa(PreDc.Zv.freqHz, 'double'))
        end



        function assert_PostDC(PostDc)
            EJ_library.assert.struct(PostDc, ...
                {'Zv'}, {});
            
            EJ_library.assert.struct(PostDc.Zv, ...
                {'Epoch', 'freqHz', ...
                'QUALITY_BITMASK', 'QUALITY_FLAG', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', ...
                'DemuxerOutput', 'currentAAmpere'}, ...
                {'BW'});
            
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
            [muxModesRemove, settingMuxModesKey] = SETTINGS.get_fv('PROCESSING.L2.REMOVE_DATA.MUX_MODES');
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
            bicas.proc_sub.log_UFV_records(zvEpoch, zvUseFillValues, logHeaderStr, L)
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
                    utc1  = EJ_library.cdf.tt2000_to_UTC_str(zvEpoch(iCdfRecord1));
                    utc2  = EJ_library.cdf.tt2000_to_UTC_str(zvEpoch(iCdfRecord2));
                    L.logf(LL, '    Records %7i-%7i, %s -- %s', iCdfRecord1, iCdfRecord2, utc1, utc2);
                end
            end
            
        end



        % Demultiplex and calibrate voltages.
        %
        % NOTE: Can handle arrays of any size as long as the sizes are
        % consistent.
        %
        function AsrSamplesAVolt = calibrate_demux_voltages(PreDc, Cal, L)
        % PROPOSAL: Incorporate into processing function process_calibrate_demux_filter.
        % PROPOSAL: Assert same nbr of "records" for MUX_SET, DIFF_GAIN as for BIAS_x.
        %
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



            % Pre-allocate. Important for speeding up LFR-SWF which tends to be
            % broken into subsequences of 1 record.
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
                PreDc.Zv.CALIBRATION_TABLE_INDEX);
            L.logf('info', 'Calibrating voltages - One sequence of records with identical settings at a time.')
            
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
                CALIBRATION_TABLE_INDEX_ss = PreDc.Zv.CALIBRATION_TABLE_INDEX(iFirst, :);
                
                % PROPOSAL: Make into "proper" table.
                %   NOTE: Can not use EJ_library.str.assist_print_table since it requires the entire table to
                %         pre-exist.
                %   PROPOSAL: Print after all iterations.
                L.logf('info', ['Records %7i-%7i : %s -- %s', ...
                    ' MUX_SET=%i; DIFF_GAIN=%i; dlrUsing12=%i; freqHz=%5g; iCalibL=%i; iCalibH=%i;', ...
                    ' CALIBRATION_TABLE_INDEX=[%i, %i]'], ...
                    iFirst, iLast, ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iLast)), ...
                    MUX_SET_ss, DIFF_GAIN_ss, dlrUsing12_ss, freqHz_ss, iCalibL_ss, iCalibH_ss, ...
                    CALIBRATION_TABLE_INDEX_ss(1), ...
                    CALIBRATION_TABLE_INDEX_ss(2))

                %============================================
                % FIND DEMUXER ROUTING, BUT DO NOT CALIBRATE
                %============================================
                % NOTE: Call demultiplexer with no samples. Only for collecting
                % information on which BLTS channels are connected to which
                % ASRs.
                [BltsSrcAsrArray, ~] = bicas.demultiplexer.main(MUX_SET_ss, dlrUsing12_ss, {[],[],[],[],[]});



                % Extract subsequence of DATA records to "demux".
                ssSamplesTm                = bicas.proc_utils.select_row_range_from_cell_comps(PreDc.Zv.samplesCaTm, iFirst, iLast);
                % NOTE: "zVariable" (i.e. first index=record) for only the
                % current subsequence.
                ssZvNValidSamplesPerRecord = PreDc.Zv.nValidSamplesPerRecord(iFirst:iLast);
                if PreDc.hasSnapshotFormat
                    % NOTE: Vector of constant numbers (one per snapshot).
                    ssDtSec = 1 ./ PreDc.Zv.freqHz(iFirst:iLast);
                else
                    % NOTE: Scalar (one for entire sequence).
                    ssDtSec = double(PreDc.Zv.Epoch(iLast) - PreDc.Zv.Epoch(iFirst)) / (iLast-iFirst) * 1e-9;   % TEMPORARY
                end
                
                biasHighGain = DIFF_GAIN_ss;



                %===================
                % ITERATE OVER BLTS
                %===================
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
                        
                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%
                        %  CALIBRATE VOLTAGES
                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%
                        CalSettings = struct();
                        CalSettings.iBlts        = iBlts;
                        CalSettings.BltsSrc      = BltsSrcAsrArray(iBlts);
                        CalSettings.biasHighGain = biasHighGain;
                        CalSettings.iCalibTimeL  = iCalibL_ss;
                        CalSettings.iCalibTimeH  = iCalibH_ss;
                        CalSettings.iLsf         = iLsf_ss;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ssSamplesCaAVolt = Cal.calibrate_voltage_all(ssDtSec, ssSamplesCaTm, ...
                            PreDc.isLfr, PreDc.isTdsCwf, CalSettings, CALIBRATION_TABLE_INDEX_ss);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        if PreDc.hasSnapshotFormat
                            [ssSamplesAVolt{iBlts}, ~] = bicas.proc_utils.convert_cell_array_of_vectors_to_matrix(...
                                ssSamplesCaAVolt, size(ssSamplesTm{iBlts}, 2));
                        else
                            ssSamplesAVolt{iBlts} = ssSamplesCaAVolt{1};   % NOTE: Must be column array.
                        end
                    end
                end    % for iBlts = 1:5
                
                %====================
                % CALL DEMULTIPLEXER
                %====================
                [~, SsAsrSamplesAVolt] = bicas.demultiplexer.main(MUX_SET_ss, dlrUsing12_ss, ssSamplesAVolt);
                
                % Add demuxed sequence to the to-be complete set of records.
                AsrSamplesAVolt = bicas.proc_utils.set_struct_field_rows(AsrSamplesAVolt, SsAsrSamplesAVolt, iFirst:iLast);
                
            end    % for iSubseq = 1:length(iFirstList)
            
            
            
            % NOTE: Assumes no "return" statement.
            %bicas.log_speed_profiling(L, 'bicas.proc_sub.calibrate_demux_voltages', tTicToc, nRecords, 'record')
            %bicas.log_memory_profiling(L, 'bicas.proc_sub.calibrate_demux_voltages:end')
        end    % calibrate_demux_voltages



    end    % methods(Static, Access=private)
        
end
