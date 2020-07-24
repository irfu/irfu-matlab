%
% Class that collects "processing functions" as public static methods.
%
% This class is not meant to be instantiated.
%
%
% CODE CONVENTIONS
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like" data, use the first MATLAB array index to
%   represent CDF records.
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% See bicas.calib.
% ZV   : CDF zVariable, or something analogous to it. If refers to CDF:ish content, then the first index corresponds to
%        the CDF record.
% SPR  : Samples Per (CDF) Record. Only refers to actual data (currents, voltages), not metadata.
%
%
% SOME INTERMEDIATE PROCESSING DATA FORMATS
% =========================================
% - PreDC = Pre-Demuxing-Calibration Data
%       Generic data format that can represent all forms of input datasets before demuxing and calibration. Can use an
%       arbitrary number of samples per record. Some variables are therefore not used in CWF output datasets.
% - PostDC = Post-Demuxing-Calibration Data
%       Like PreDC but with additional fields. Tries to capture a superset of the information that goes into any
%       dataset produced by BICAS, and the exact set of variables that goes into the output datasets.
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
%
% PROPOSAL: Merge process_PostDC_to_LFR and process_PostDC_to_TDS.
%
% PROPOSAL: Instantiate class, use instance methods instead of static.
%   PRO: Can have SETTINGS and constants as instance variable instead of calling global variables.
%
% PROPOSAL: Submit zVar variable attributes.
%   PRO: Can interpret fill values.
%       Ex: Can doublecheck TDS RSWF snapshot length using fill values and compare with zVar SAMPS_PER_CH (which seems
%           to be bad).
%
% PROPOSAL: Return (to execute_sw_mode), global attributes.
%   PRO: Needed for output datasets: CALIBRATION_TABLE, CALIBRATION_VERSION
%       ~CON: CALIBRATION_VERSION refers to algorithm and should maybe be a SETTING.
%#######################################################################################################################

    methods(Static, Access=public)
        
        function HkSciTime = process_HK_to_HK_on_SCI_TIME(InSci, InHk, SETTINGS, L)
        % Processing function
        
            % ASSERTIONS
            EJ_library.assert.struct(InSci, {'Zv', 'Ga'}, {})
            EJ_library.assert.struct(InHk,  {'Zv', 'Ga'}, {})
            
            HkSciTime = [];
            
            
            
            %=========================================================================================================
            % Select whether HK should use
            %   (1) Epoch, or
            %   (2) ACQUISITION_TIME (not always available).
            % ----------------------------------------------
            % IMPLEMENTATION NOTE: Historically, there have been datasets where Epoch is contains errors, but
            % ACQUISITION_TIME seems OK. This should be phased out eventually.
            %=========================================================================================================
            ACQUISITION_TIME_EPOCH_UTC = SETTINGS.get_fv('INPUT_CDF.ACQUISITION_TIME_EPOCH_UTC');
            USE_ZV_ACQUISITION_TIME_HK = SETTINGS.get_fv('PROCESSING.HK.USE_ZV_ACQUISITION_TIME');
            if USE_ZV_ACQUISITION_TIME_HK
                hkEpoch = bicas.proc_utils.ACQUISITION_TIME_to_tt2000(...
                    InHk.Zv.ACQUISITION_TIME, ...
                    ACQUISITION_TIME_EPOCH_UTC);
                L.logf('warning', 'Using HK zVar ACQUISITION_TIME instead of Epoch.')
                
                % NOTE: ACQUISITION_TIME in test file
                % TDS___TESTDATA_RGTS_TDS_CALBA_V0.8.6/solo_HK_rpw-bia_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
                % is not monotonically increasing (in fact, it is completely strange).
%                 assert(issorted(hkEpoch, 'strictascend'), ...
%                     'BICAS:proc_sub:Assertion:DatasetFormat', ...
%                     'Trying to use HK zVar ACQUISITION_TIME but it is not monotonically increasing (when converting to tt2000).')
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
                    bicas.proc_utils.ACQUISITION_TIME_to_tt2000(InHk.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            end
            if isfield(InSci.Zv, 'ACQUISITION_TIME') && ~isempty(InSci.Zv.ACQUISITION_TIME)
                TimeVars.SCI_ACQUISITION_TIME_tt2000 = ...
                    bicas.proc_utils.ACQUISITION_TIME_to_tt2000(InSci.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
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
                
                % NOTE: "WARNING" (rather than error) only makes sense if it is possible to later meaningfully permit
                % non-intersection.
                [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY');
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    'SCI and HK time ranges do not overlap in time.', 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
            end
            
            hkEpochExtrapMargin = mode(diff(hkEpoch)) / 2;   % NOTE: Requires >=2 records.

            %=========================================================================================================
            % Derive MUX_SET
            % --------------
            % NOTE: Only obtains one MUX_SET per record ==> Can not change MUX_SET in the middle of a record.
            % NOTE: Can potentially obtain MUX_SET from LFR SCI.
            %=========================================================================================================            
            HkSciTime.MUX_SET = bicas.utils.interpolate_nearest(...
                hkEpochExtrapMargin, ...
                hkEpoch, ...
                InHk.Zv.HK_BIA_MODE_MUX_SET, ...
                InSci.Zv.Epoch);



            %=========================================================================================================
            % Derive DIFF_GAIN
            % ----------------
            % NOTE: Not perfect handling of time when 1 snapshot/record, since one should ideally use time stamps
            % for every LFR _sample_.
            %=========================================================================================================
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
            
            
            
            %========================================================================================
            % CDF ASSERTION: CURRENT data begins before SCI data (i.e. there is enough CURRENT data).
            %========================================================================================
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
            
            
            
            %===========================================================================================================
            % CDF ASSERTION: Epoch increases (not monotonically)
            % --------------------------------------------------
            % NOTE: bicas.proc_sub.interpolate_current checks (and handles) that Epoch increases monotonically, but only
            % for each antenna separately (which does not capture all cases).
            % Ex: Timestamps, iAntenna = mod(iRecord,3): 1,2,3,5,4,6
            %       ==> Monotonically increasing sequences for each antenna separately, but not even increasing when
            %           combined.
            %===========================================================================================================
            if ~issorted(InCur.Zv.Epoch)
                error('CURRENT timestamps do not increase (all antennas combined).')
            end
            
            % NOTE: bicas.proc_sub.interpolate_current checks that Epoch increases monotonically.
            currentNanoSAmpere = [];
            currentNanoSAmpere(:,1) = bicas.proc_sub.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_1, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,2) = bicas.proc_sub.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_2, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,3) = bicas.proc_sub.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_3, sciEpoch, L, SETTINGS);
            
            currentSAmpere = 1e-9 * currentNanoSAmpere;
        end

        
        
        % Wrapper around EJ_library.so.CURRENT_zv_to_current_interpolate for anomaly handling.
        function sciZv_IBIASx = zv_TC_to_current(curZv_Epoch, curZv_IBIAS_x, sciZv_Epoch, L, SETTINGS)
            
            [sciZv_IBIASx, duplicateAnomaly] = EJ_library.so.CURRENT_zv_to_current_interpolate(...
                double(curZv_Epoch), ...
                curZv_IBIAS_x, ...
                sciZv_Epoch);
            
            if duplicateAnomaly
                %======================================================================================================
                % CDF ASSERTION
                % Handle non-monotonically increasing Epoch
                %======================================================================================================

                [settingValue, settingKey] = SETTINGS.get_fv('INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY');
                anomalyDescriptionMsg = 'Bias current data contain dupicate settings, with identical timestamps and identical bias settings on the same antenna.';
                
                switch(settingValue)
                    case 'REMOVE_DUPLICATES'
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', ...
                            anomalyDescriptionMsg)
                        L.log('warning', ...
                            'Removed duplicated bias current settings with identical timestamps on the same antenna.')

                    otherwise
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+illegal', ...
                            anomalyDescriptionMsg, 'BICAS:proc_sub:SWModeProcessing:DatasetFormat')
                end
            end
        end    % TC_to_current
        
        
        
        function PreDc = process_LFR_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
        % Processing function. Convert LFR CDF data to PreDC.
        %
        % Keeps number of samples/record. Treats 1 samples/record "length-one snapshots".
        
        % PROBLEM: Hardcoded CDF data types (MATLAB classes).
        % MINOR PROBLEM: Still does not handle LFR zVar TYPE for determining "virtual snapshot" length.
        % Should only be relevant for V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF (not V02) which should expire.
        
            % ASSERTIONS
            EJ_library.assert.struct(InSci,     {'Zv', 'Ga'}, {})
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})



            %=================
            % Normalize zVars
            %=================
            %---------------------------------------------------------------------------
            % Workaround: Normalize LFR data to handle variations that should not exist
            %---------------------------------------------------------------------------
            % Handle that SYNCHRO_FLAG (empty) and TIME_SYNCHRO_FLAG (non-empty) may BOTH be present.
            % "DEFINITION BUG" in definition of datasets/skeleton?
            % Ex: LFR___TESTDATA_RGTS_LFR_CALBUT_V0.7.0/ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_4129f0b_CNE_V02.cdf /2020-03-17
            if SETTINGS.get_fv('INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED') ...
                    && isfield(InSci.Zv, 'SYNCHRO_FLAG') && isempty(InSci.Zv.SYNCHRO_FLAG)
                InSci.Zv = rmfield(InSci.Zv, 'SYNCHRO_FLAG');
            end
            %------------------------
            % "Normal" normalization
            %------------------------
            % 2020-01-21: Based on skeletons (.skt; L1R, L2), SYNCHRO_FLAG seems to be the correct one.
            [InSci.Zv, fnChangeList] = EJ_library.utils.normalize_struct_fieldnames(InSci.Zv, ...
                {{{'TIME_SYNCHRO_FLAG', 'SYNCHRO_FLAG'}, 'SYNCHRO_FLAG'}}, 'Assert one matching candidate');

            bicas.proc_sub.handle_zv_name_change(...
                fnChangeList, inSciDsi, SETTINGS, L, 'SYNCHRO_FLAG', 'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')



            nRecords = size(InSci.Zv.Epoch, 1);
            C = EJ_library.so.adm.classify_DATASET_ID(inSciDsi);
            
            
            
            % CDF ASSERTION
            if ~issorted(InSci.Zv.Epoch, 'strictascend')
                error('Voltage (science) dataset timestamps do not increase.')
            end
            
            V = InSci.Zv.V;
            E = permute(InSci.Zv.E, [1,3,2]);
            % Switch last two indices of E.
            % ==> index 2 = "snapshot" sample index, including for CWF (sample/record, "snapshots" consisting of 1 sample).
            %     index 3 = E1/E2 component
            %               NOTE: 1/2=index into array; these are diffs but not equivalent to any particular diffs).
            
            nCdfSamplesPerRecord = EJ_library.assert.sizes(V, [nRecords, -1], E, [nRecords, -1, 2]);
            
            % ASSERTIONS
            if C.isLfrSurvSwf   assert(nCdfSamplesPerRecord == EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH)
            else                assert(nCdfSamplesPerRecord == 1)
            end



            % Set iLsfZv.
            if     C.isLfrSbm1   iLsfZv = ones(nRecords, 1) * 2;   % Always value "2" (F1, "FREQ = 1").
            elseif C.isLfrSbm2   iLsfZv = ones(nRecords, 1) * 3;   % Always value "3" (F2, "FREQ = 2").
            else                 iLsfZv = InSci.Zv.FREQ + 1;
                % NOTE: Translates from LFR's FREQ values (0=F0 etc) to LSF index values (1=F0) used in loaded RCT data
                % structs.
            end
            EJ_library.assert.size(iLsfZv, [NaN, 1])



            % NOTE: Needed also for 1 SPR.
            zvFreqHz = EJ_library.so.get_LFR_frequency( iLsfZv );

            % Obtain the relevant values (one per record) from zVariables R0, R1, R2, and the virtual "R3".
            zvRx = EJ_library.so.get_LFR_Rx(...
                InSci.Zv.R0, ...
                InSci.Zv.R1, ...
                InSci.Zv.R2, ...
                iLsfZv );   % NOTE: Function also handles the imaginary zVar "R3".

            PreDc = [];
            PreDc.Zv.Epoch                  = InSci.Zv.Epoch;
            PreDc.Zv.DELTA_PLUS_MINUS       = bicas.proc_utils.derive_DELTA_PLUS_MINUS(zvFreqHz, nCdfSamplesPerRecord);            
            PreDc.Zv.freqHz                 = zvFreqHz;
            PreDc.Zv.nValidSamplesPerRecord = ones(nRecords, 1) * nCdfSamplesPerRecord;
            PreDc.Zv.SYNCHRO_FLAG           = InSci.Zv.SYNCHRO_FLAG;
            PreDc.Zv.BW                     = InSci.Zv.BW;
            PreDc.Zv.useFillValues          = ~logical(InSci.Zv.BW);
            if isfield(InSci.Zv, 'CALIBRATION_TABLE_INDEX')
                % NOTE: CALIBRATION_TABLE_INDEX exists for L1R, but not L1.
                PreDc.Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;
            end



            %===========================================================================================================
            % Replace illegally empty data with fill values/NaN
            % -------------------------------------------------
            % IMPLEMENTATION NOTE: QUALITY_FLAG, QUALITY_BITMASK have been found empty in test data, but should have
            % attribute DEPEND_0 = "Epoch" ==> Should have same number of records as Epoch.
            % Can not save CDF with zVar with zero records (crashes when reading CDF). ==> Better create empty records.
            %
            % Test data:
            %  MYSTERIOUS_SIGNAL_1_2016-04-15_Run2__7729147__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-SWF_7729147_CNE_V01.cdf
            %  ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_4129f0b_CNE_V02.cdf
            %  ROC-SGSE_L1R_RPW-LFR-SBM2-CWF-E_6b05822_CNE_V02.cdf
            %
            % PROPOSAL: Move to the code that reads CDF datasets instead. Generalize to many zVariables.
            %===========================================================================================================
            PreDc.Zv.QUALITY_FLAG    = InSci.Zv.QUALITY_FLAG;
            PreDc.Zv.QUALITY_BITMASK = InSci.Zv.QUALITY_BITMASK;
            
            [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY');
            
            if isempty(PreDc.Zv.QUALITY_FLAG)
                
                anomalyDescrMsg = 'QUALITY_FLAG from the LFR SCI source dataset is empty.';
                switch(settingValue)
                    case 'USE_FILL_VALUE'
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', ...
                            anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
                        L.log('warning', 'Filling with empty values.')
                        PreDc.Zv.QUALITY_FLAG = bicas.proc_utils.create_NaN_array([nRecords, 1]);
                    otherwise
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+illegal', ...
                            anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
                end
                
            end
            if isempty(PreDc.Zv.QUALITY_BITMASK)
                
                anomalyDescrMsg = 'QUALITY_BITMASK from the LFR SCI source dataset is empty.';
                switch(settingValue)
                    case 'USE_FILL_VALUE'
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', ...
                            anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
                        L.log('warning', 'Filling with empty values.')
                        PreDc.Zv.QUALITY_BITMASK = bicas.proc_utils.create_NaN_array([nRecords, 1]);               
                    otherwise
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+illegal', ...
                            anomalyDescrMsg, 'BICAS:proc_sub:DatasetFormat:SWModeProcessing')
                end
                
            end
            
            % ASSERTIONS
            % LFR QUALITY_FLAG, QUALITY_BITMASK not set yet (2019-09-17), but I presume they should have just one value
            % per record. BIAS output datasets should.
            EJ_library.assert.sizes(...
                PreDc.Zv.QUALITY_FLAG,    [nRecords, 1], ...
                PreDc.Zv.QUALITY_BITMASK, [nRecords, 1])



            % E must be floating-point so that values can be set to NaN.
            % ==> V must be floating-point so that demultiplexer can subtract/add V/E, and calibration can work.
            % bicas.proc_utils.filter_rows requires this. Variable may be integer if integer in source CDF.
            V = single(V);
            E = single(E);

            PreDc.Zv.samplesCaTm    = cell(5,1);
            PreDc.Zv.samplesCaTm{1} = V;
            PreDc.Zv.samplesCaTm{2} = bicas.proc_utils.filter_rows( E(:,:,1), zvRx==0 );    % Copy values, except when zvRx==0 (==>NaN).
            PreDc.Zv.samplesCaTm{3} = bicas.proc_utils.filter_rows( E(:,:,2), zvRx==0 );
            PreDc.Zv.samplesCaTm{4} = bicas.proc_utils.filter_rows( E(:,:,1), zvRx==1 );
            PreDc.Zv.samplesCaTm{5} = bicas.proc_utils.filter_rows( E(:,:,2), zvRx==1 );
            
            
            
            %==================================================================
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
            
            PreDc.Zv.DIFF_GAIN         = HkSciTime.DIFF_GAIN;
            PreDc.Zv.iLsf              = iLsfZv;
            
            PreDc.hasSnapshotFormat    = C.isLfrSurvSwf;
            PreDc.nRecords             = nRecords;
            PreDc.nCdfSamplesPerRecord = nCdfSamplesPerRecord;
            PreDc.isLfr                = true;
            PreDc.isTdsCwf             = false;



            % ASSERTIONS
            bicas.proc_sub.assert_PreDC(PreDc)
        end
        
        
        
        function PreDc = process_TDS_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
        % Processing function. Convert TDS CDF data (PDs) to PreDC.
        %
        % Keeps number of samples/record. Treats 1 samples/record "length-one snapshots".
        %
        % BUG?: Does not use CHANNEL_STATUS_INFO.
        % NOTE: BIAS output datasets do not have a variable for the length of snapshots. Need to use NaN/fill value.

            % ASSERTIONS
            EJ_library.assert.struct(InSci,     {'Zv', 'Ga'}, {})
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})

            C = EJ_library.so.adm.classify_DATASET_ID(inSciDsi);

            % CDF ASSERTION
            if ~issorted(InSci.Zv.Epoch, 'strictascend')
                error('Voltage timestamps do not increase (all antennas combined).')
            end
            
            %===============================================================================================
            % Normalize zVar names
            % --------------------
            % Both zVars TIME_SYNCHRO_FLAG, SYNCHRO_FLAG found in input datasets (2020-01-05). Unknown why.
            % "DEFINITION BUG" in definition of datasets/skeleton?
            % 2020-01-21: Based on skeletons (.skt; L1R, L2), SYNCHRO_FLAG seems to be the correct one.
            %===============================================================================================
            [InSci.Zv, fnChangeList] = EJ_library.utils.normalize_struct_fieldnames(InSci.Zv, ...
                {{{'TIME_SYNCHRO_FLAG', 'SYNCHRO_FLAG'}, 'SYNCHRO_FLAG'}}, 'Assert one matching candidate');
            
            bicas.proc_sub.handle_zv_name_change(...
                fnChangeList, inSciDsi, SETTINGS, L, 'SYNCHRO_FLAG', 'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')



            nRecords             = size(InSci.Zv.Epoch, 1);
            nCdfSamplesPerRecord = size(InSci.Zv.WAVEFORM_DATA, 3);    % Number of samples in the zVariable, not necessarily actual data.

            % CDF ASSERTION
            if ~issorted(InSci.Zv.Epoch, 'strictascend')
                error('Voltage timestamps do not increase (all antennas combined).')
            end            
            
            freqHzZv = double(InSci.Zv.SAMPLING_RATE);
            
            if any(freqHzZv == 255)
                [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY');
                anomalyDescrMsg = 'Finds illegal stated sampling frequency 255 in TDS L1/L1R LFM-RSWF dataset.';
                
                if C.isTdsRswf
                    switch(settingValue)
                        case 'CORRECT'
                            % IMPLEMENTATION NOTE: Has observed test file
                            % TESTDATA_RGTS_TDS_CALBA_V0.8.5C: solo_L1R_rpw-tds-lfm-rswf-e_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
                            % to have SAMPLING_RATE == 255, which is likely a BUG in the dataset. /Erik P G Johansson 2019-12-03
                            % Bug in TDS RCS.  /David Pisa 2019-12-03
                            % Setting it to what is probably the correct value.
                            freqHzZv(freqHzZv == 255) = 32768;
                            L.logf('warning', ...
                                'Using workaround to modify instances of sampling frequency 255-->32768.')
                            bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', anomalyDescrMsg)
                            
                        otherwise
                            bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', anomalyDescrMsg, 'BICAS:DatasetFormat')
                    end
                else
                    error(anomalyDescrMsg)
                end
            end
            
            PreDc = [];
            
            PreDc.Zv.Epoch            = InSci.Zv.Epoch;
            PreDc.Zv.DELTA_PLUS_MINUS = bicas.proc_utils.derive_DELTA_PLUS_MINUS(freqHzZv, nCdfSamplesPerRecord);
            PreDc.Zv.freqHz           = freqHzZv;
            PreDc.Zv.QUALITY_FLAG     = InSci.Zv.QUALITY_FLAG;
            PreDc.Zv.QUALITY_BITMASK  = InSci.Zv.QUALITY_BITMASK;
            PreDc.Zv.SYNCHRO_FLAG     = InSci.Zv.SYNCHRO_FLAG;
            PreDc.Zv.MUX_SET          = HkSciTime.MUX_SET;
            PreDc.Zv.DIFF_GAIN        = HkSciTime.DIFF_GAIN;
            PreDc.Zv.useFillValues    = false(nRecords, 1);
            if isfield(InSci.Zv, 'CALIBRATION_TABLE_INDEX')
                % NOTE: CALIBRATION_TABLE_INDEX exists for L1R, but not L1.
                PreDc.Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;
            end



            %=====================================
            % Set PreDc.Zv.nValidSamplesPerRecord
            %=====================================
            if C.isTdsRswf
                %====================================================================================================
                % Check for and handle illegal input data, zVar SAMPS_PER_CH
                % ----------------------------------------------------------
                % NOTE: Has observed invalid SAMPS_PER_CH value 16562 in
                % ROC-SGSE_L1R_RPW-TDS-LFM-RSWF-E_73525cd_CNE_V03.CDF.
                % 2019-09-18, David Pisa: Not a flaw in TDS RCS but in the source L1 dataset.
                %====================================================================================================
                SAMPS_PER_CH_MIN_VALID    = 2^10;
                SAMPS_PER_CH_MAX_VALID    = 2^15;
                zv_SAMPS_PER_CH           = double(InSci.Zv.SAMPS_PER_CH);
                zv_SAMPS_PER_CH_rounded   = round(2.^round(log2(zv_SAMPS_PER_CH)));
                zv_SAMPS_PER_CH_rounded(zv_SAMPS_PER_CH_rounded < SAMPS_PER_CH_MIN_VALID) = SAMPS_PER_CH_MIN_VALID;
                zv_SAMPS_PER_CH_rounded(zv_SAMPS_PER_CH_rounded > SAMPS_PER_CH_MAX_VALID) = SAMPS_PER_CH_MAX_VALID;
                if any(zv_SAMPS_PER_CH_rounded ~= zv_SAMPS_PER_CH)
                    SAMPS_PER_CH_badValues = unique(zv_SAMPS_PER_CH(zv_SAMPS_PER_CH_rounded ~= zv_SAMPS_PER_CH));
                    
                    badValuesDisplayStr = strjoin(arrayfun(...
                        @(n) sprintf('%i', n), SAMPS_PER_CH_badValues, 'uni', false), ', ');
                    anomalyDescrMsg = sprintf(...
                        'TDS LFM RSWF zVar SAMPS_PER_CH contains unexpected value(s), not 2^n: %s', ...
                        badValuesDisplayStr);
                    
                    [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY');
                    switch(settingValue)
                        case 'ROUND'
                            bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', ...
                                anomalyDescrMsg, 'BICAS:proc_sub:Assertion:DatasetFormat')
                            L.log('warning', ...
                                ['Replacing TDS RSWF zVar SAMPS_PER_CH values with values, rounded to valid', ...
                                ' values due to setting PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY.'])
                            
                            zv_SAMPS_PER_CH = zv_SAMPS_PER_CH_rounded;
                            
                        otherwise
                            bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                                anomalyDescrMsg, 'BICAS:proc_sub:Assertion:DatasetFormat')

                    end
                end
                
                % NOTE: This might only be appropriate for TDS's "COMMON_MODE" mode. TDS also has a "FULL_BAND" mode
                % with 2^18=262144 samples per snapshot. You should never encounter FULL_BAND in any dataset (even on
                % ground), only used for calibration and testing. /David Pisa & Jan Soucek in emails, 2016.
                % --
                % FULL_BAND mode has each snapshot divided into 2^15 samples/record * 8 records.  /Unknown source
                % Unclear what value SAMPS_PER_CH should have for FULL_BAND mode. How does Epoch work for FULL_BAND
                % snapshots?
                PreDc.Zv.nValidSamplesPerRecord = zv_SAMPS_PER_CH;
            else
                PreDc.Zv.nValidSamplesPerRecord = ones(nRecords, 1) * 1;
            end



            if C.isL1R
                assert(size(InSci.Zv.WAVEFORM_DATA, 2) == 3, ...
                    'BICAS:proc_sub:process_TDS_to_PreDC:Assertion:DatasetFormat', 'TDS zVar WAVEFORM_DATA has an unexpected size.')
            elseif C.isL1
                assert(size(InSci.Zv.WAVEFORM_DATA, 2) == 8, ...
                    'BICAS:proc_sub:process_TDS_to_PreDC:Assertion:DatasetFormat', 'TDS zVar WAVEFORM_DATA has an unexpected size.')
            end
            modif_WAVEFORM_DATA = double(permute(InSci.Zv.WAVEFORM_DATA, [1,3,2]));

            PreDc.Zv.samplesCaTm    = {};
            PreDc.Zv.samplesCaTm{1} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,1), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{2} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,2), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{3} = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,3), PreDc.Zv.nValidSamplesPerRecord );
            PreDc.Zv.samplesCaTm{4} = bicas.proc_utils.create_NaN_array([nRecords, nCdfSamplesPerRecord]);
            PreDc.Zv.samplesCaTm{5} = bicas.proc_utils.create_NaN_array([nRecords, nCdfSamplesPerRecord]);

            PreDc.isLfr                = false;
            PreDc.isTdsCwf             = C.isTdsCwf;
            PreDc.hasSnapshotFormat    = C.isTdsRswf;
            PreDc.nRecords             = nRecords;
            PreDc.nCdfSamplesPerRecord = nCdfSamplesPerRecord;
            PreDc.Zv.iLsf              = zeros(nRecords, 1) * NaN;   % Only set becuse the code shared with LFR requires it.



            % ASSERTIONS
            bicas.proc_sub.assert_PreDC(PreDc)
        end



        function assert_PreDC(PreDc)
            EJ_library.assert.struct(PreDc, ...
                {'Zv', 'hasSnapshotFormat', 'nRecords', 'nCdfSamplesPerRecord', 'isLfr', 'isTdsCwf'}, {});
            EJ_library.assert.struct(PreDc.Zv, ...
                {'Epoch', 'samplesCaTm', 'freqHz', 'nValidSamplesPerRecord', 'iLsf', 'DIFF_GAIN', ...
                'MUX_SET', 'QUALITY_FLAG', 'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'useFillValues'}, ...
                {'CALIBRATION_TABLE_INDEX', 'BW'});
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(PreDc.Zv);

            assert(isa(PreDc.Zv.freqHz, 'double'))
        end



        function assert_PostDC(PostDc)
            EJ_library.assert.struct(PostDc, ...
                {'Zv', 'hasSnapshotFormat', 'nRecords', 'nCdfSamplesPerRecord', 'isLfr', 'isTdsCwf'}, {});
            EJ_library.assert.struct(PostDc.Zv, ...
                {'Epoch', 'samplesCaTm', 'freqHz', 'nValidSamplesPerRecord', 'iLsf', 'DIFF_GAIN', ...
                'MUX_SET', 'QUALITY_FLAG', 'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'DemuxerOutput', ...
                'currentAAmpere', 'DemuxerOutput', 'useFillValues'}, ...
                {'CALIBRATION_TABLE_INDEX', 'BW'});
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(PostDc.Zv);
        end



        function [OutSciZv] = process_PostDC_to_LFR(SciPostDc, outputDsi)
        % Processing function. Convert PostDC to any one of several similar LFR dataset PDs.
        
            % ASSERTIONS
            bicas.proc_sub.assert_PostDC(SciPostDc)
            
            OutSciZv = [];
            
            nSamplesPerRecord = size(SciPostDc.Zv.DemuxerOutput.dcV1, 2);   % Samples per record.
            
            OutSciZv.Epoch            = SciPostDc.Zv.Epoch;
            OutSciZv.QUALITY_BITMASK  = SciPostDc.Zv.QUALITY_BITMASK;
            OutSciZv.QUALITY_FLAG     = SciPostDc.Zv.QUALITY_FLAG;
            OutSciZv.BW               = SciPostDc.Zv.BW;
            OutSciZv.DELTA_PLUS_MINUS = SciPostDc.Zv.DELTA_PLUS_MINUS;
            OutSciZv.SYNCHRO_FLAG     = SciPostDc.Zv.SYNCHRO_FLAG;
            OutSciZv.SAMPLING_RATE    = SciPostDc.Zv.freqHz;

            % NOTE: Convert AAmpere --> (antenna) nA
            OutSciZv.IBIAS1           = SciPostDc.Zv.currentAAmpere(:, 1) * 1e9;
            OutSciZv.IBIAS2           = SciPostDc.Zv.currentAAmpere(:, 2) * 1e9;
            OutSciZv.IBIAS3           = SciPostDc.Zv.currentAAmpere(:, 3) * 1e9;
            
            % NOTE: The two cases are different in the indexes they use for OutSciZv.
            switch(outputDsi)
                case  {'SOLO_L2_RPW-LFR-SURV-CWF-E' ...
                       'SOLO_L2_RPW-LFR-SBM1-CWF-E' ...
                       'SOLO_L2_RPW-LFR-SBM2-CWF-E'}

                    % ASSERTION
                    assert(nSamplesPerRecord == 1, ...
                        'BICAS:proc_sub:Assertion:IllegalArgument', ...
                        'Number of samples per CDF record is not 1, as expected. Bad input CDF?')
                    assert(size(OutSciZv.QUALITY_FLAG,    2) == 1)
                    assert(size(OutSciZv.QUALITY_BITMASK, 2) == 1)
                    
                    OutSciZv.VDC(:,1) = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.VDC(:,2) = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.VDC(:,3) = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.EDC(:,1) = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.EDC(:,2) = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.EDC(:,3) = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,3) = SciPostDc.Zv.DemuxerOutput.acV23;
                    
                case  {'SOLO_L2_RPW-LFR-SURV-SWF-E'}
                    
                    % ASSERTION
                    assert(nSamplesPerRecord == EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH, ...
                        'BICAS:proc_sub:Assertion:IllegalArgument', ...
                        'Number of samples per CDF record is not %i, as expected. Bad Input CDF?', ...
                        EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH)
                    
                    OutSciZv.VDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.VDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.VDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.EDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.EDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.EDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,:,3) = SciPostDc.Zv.DemuxerOutput.acV23;

                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', ...
                        'Function can not produce outputDsi=%s.', outputDsi)
            end
            
            
            
            % ASSERTION
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(OutSciZv);
            % NOTE: Not really necessary since the list of zVars will be checked against the master CDF?
            EJ_library.assert.struct(OutSciZv, {...
                'IBIAS1', 'IBIAS2', 'IBIAS3', 'VDC', 'EDC', 'EAC', 'Epoch', 'QUALITY_BITMASK', 'QUALITY_FLAG', 'BW', ...
                'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'SAMPLING_RATE'}, {})
        end   % process_PostDC_to_LFR



        function OutSciZv = process_PostDC_to_TDS(SciPostDc, outputDsi)
            
            % ASSERTIONS
            bicas.proc_sub.assert_PostDC(SciPostDc)
            
            OutSciZv = [];
            
            OutSciZv.Epoch            = SciPostDc.Zv.Epoch;
            OutSciZv.QUALITY_FLAG     = SciPostDc.Zv.QUALITY_FLAG;
            OutSciZv.QUALITY_BITMASK  = SciPostDc.Zv.QUALITY_BITMASK;
            OutSciZv.DELTA_PLUS_MINUS = SciPostDc.Zv.DELTA_PLUS_MINUS;
            OutSciZv.SYNCHRO_FLAG     = SciPostDc.Zv.SYNCHRO_FLAG;
            OutSciZv.SAMPLING_RATE    = SciPostDc.Zv.freqHz;

            % NOTE: Convert AAmpere --> (antenna) nA
            OutSciZv.IBIAS1           = SciPostDc.Zv.currentAAmpere(:, 1) * 1e9;
            OutSciZv.IBIAS2           = SciPostDc.Zv.currentAAmpere(:, 2) * 1e9;
            OutSciZv.IBIAS3           = SciPostDc.Zv.currentAAmpere(:, 3) * 1e9;
            
            % NOTE: The two cases are actually different in the indexes they use for OutSciZv.
            switch(outputDsi)
                
                case {'SOLO_L2_RPW-TDS-LFM-CWF-E'}

                    OutSciZv.VDC(:,1)   = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.VDC(:,2)   = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.VDC(:,3)   = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.EDC(:,1)   = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.EDC(:,2)   = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.EDC(:,3)   = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,1)   = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,2)   = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,3)   = SciPostDc.Zv.DemuxerOutput.acV23;
                    
                case {'SOLO_L2_RPW-TDS-LFM-RSWF-E'}
                    OutSciZv.VDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.VDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.VDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.EDC(:,:,1) = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.EDC(:,:,2) = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.EDC(:,:,3) = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,:,3) = SciPostDc.Zv.DemuxerOutput.acV23;
                    
                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', ...
                        'Function can not produce outputDsi=%s.', outputDsi)
            end



            % ASSERTION
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(OutSciZv);
            % NOTE: Not really necessary since the list of zVars will be checked against the master CDF?
            EJ_library.assert.struct(OutSciZv, {...
                'IBIAS1', 'IBIAS2', 'IBIAS3', 'VDC', 'EDC', 'EAC', 'Epoch', 'QUALITY_BITMASK', 'QUALITY_FLAG', ...
                'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'SAMPLING_RATE'}, {})
        end
        
        
        
        % Processing function. Converts PreDC to PostDC, i.e. demux and calibrate data.
        % Function is in large part a wrapper around "calibrate_demux_voltages".
        %
        % NOTE: Public function as opposed to the other demuxing/calibration functions.
        %
        function PostDc = process_calibrate_demux_filter(PreDc, InCurPd, Cal, SETTINGS, L)
        %
        % PROPOSAL: Move the setting of IBIASx (bias current) somewhere else?
        %   PRO: Unrelated to demultiplexing.
        %   CON: Related to calibration.

            % ASSERTION
            bicas.proc_sub.assert_PreDC(PreDc);
            
            
            
            %=======
            % DEMUX
            %=======
            PostDc = PreDc;    % Copy all values, to later overwrite a subset of them.
            
            nRecords = size(PostDc.Zv.Epoch, 1);
            tTicToc = tic();            
            
            PostDc.Zv.DemuxerOutput = bicas.proc_sub.calibrate_demux_voltages(PreDc, Cal, SETTINGS, L);
            
            wallTimeSec = toc(tTicToc);
            L.logf('info', ...
                'bicas.proc_sub.calibrate_demux_voltages: Time used for execution (wall time): %g [s], %g [s/record]', ...
                wallTimeSec, wallTimeSec/nRecords)



            %=========================
            % Calibrate bias currents
            %=========================
            currentSAmpere = bicas.proc_sub.process_CUR_to_CUR_on_SCI_TIME(PostDc.Zv.Epoch, InCurPd, SETTINGS, L);
            currentTm      = bicas.calib.calibrate_current_sampere_to_TM(currentSAmpere);
            
            currentAAmpere = bicas.proc_utils.create_NaN_array(size(currentSAmpere));    % Variable to fill/set.
            iCalibLZv      = Cal.get_calibration_time_L(PostDc.Zv.Epoch);
            iEdgeList      = bicas.proc_utils.find_constant_sequences(iCalibLZv, PreDc.Zv.useFillValues);
            [iFirstList, iLastList] = bicas.proc_utils.index_edges_2_first_last(iEdgeList);
            L.logf('info', 'Calibrating currents - One sequence of records with identical settings at a time.')
            for iSubseq = 1:length(iFirstList)
                iFirst = iFirstList(iSubseq);
                iLast  = iLastList(iSubseq);
                
                iRecords = iFirst:iLast;
                
                L.logf('info', 'Records %7i-%7i : %s -- %s; useFillValues=%g', ...
                    iFirst, iLast, ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iLast)), ...
                    PreDc.Zv.useFillValues(iFirst))
                
                for iAnt = 1:3
                    if PreDc.Zv.useFillValues(iFirst)
                        % Set samples to NaN.
                        currentAAmpere(iRecords, iAnt) = ones(size(currentTm(iRecords, iAnt))) * NaN;
                    else
                        currentAAmpere(iRecords, iAnt) = Cal.calibrate_current_TM_to_aampere(...
                            currentTm( iRecords, iAnt), iAnt, iCalibLZv(iRecords));
                    end
                end
            end
            
            PostDc.Zv.currentAAmpere = currentAAmpere;
            
            
            
            % Remove data because of settings.
            [...
                PostDc.Zv.currentAAmpere, ...
                PostDc.Zv.DemuxerOutput] ...
                = bicas.proc_sub.remove_PostDc_data(...
                PostDc.Zv.currentAAmpere, ...
                PostDc.Zv.DemuxerOutput, ...
                PostDc.Zv.Epoch, ...
                PostDc.Zv.MUX_SET, ...
                PostDc.isLfr, ...
                SETTINGS, L);

            
            
            % ASSERTION
            bicas.proc_sub.assert_PostDC(PostDc)
        end
        
    end    % methods(Static, Access=public)
            
    %###################################################################################################################
    
    methods(Static, Access=private)
    %methods(Static, Access=public)
    
    
    
        % Remove data because of settings
        function [zvCurrentAAmpere, DemuxerOutput] = remove_PostDc_data(...
                zvCurrentAAmpere, DemuxerOutput, zvEpoch, zv_MUX_SET, isLfr, SETTINGS, L)
            % -------------------------------
            % PROPOSAL: Separate function
            %   NOTE: PreDc as argument, or
            %   [samplesCaTm, currentAAmpere] = (isLfr, Epoch, MUX_SET, samplesCaTm, currentAAmpere, SETTINGS, L)
            %   CON: Increases memory usage?
            
            LL = 'info';    % LL = Log Level
            
            % Read settings
            % NOTE: There is no PostDc.isTds, only .isTdsCwf
            [muxModesRemove, settingMuxModesKey] = SETTINGS.get_fv('PROCESSING.L2.REMOVE_DATA.MUX_MODES');
            if     isLfr   settingMarginKey = 'PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S';
            else           settingMarginKey = 'PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S';
                %else            error('BICAS:proc_sub:Assertion', 'Neither LFR or TDS data.')
            end
            removeMarginSec = SETTINGS.get_fv(settingMarginKey);
            
            % Algorithm for finding what to remove.
            bRemoveArray = EJ_library.utils.true_with_margin(...
                zvEpoch, ...
                ismember(zv_MUX_SET, muxModesRemove), ...
                removeMarginSec * 1e9);
            
            % Remove and log
            [i1Array, i2Array] = EJ_library.utils.split_by_false(bRemoveArray);
            nRemoveIntervals = numel(i1Array);
            if nRemoveIntervals > 0
                
                %=====
                % Log
                %=====
                L.logf(LL, 'Setting science data to fill value in selected CDF records due to settings: ');
                L.logf(LL, '    Setting %s = [%s]', ...
                    settingMuxModesKey, ...
                    strjoin(EJ_library.str.sprintf_many('%g', muxModesRemove), ', '));
                L.logf(LL, '    Setting %s = %g', settingMarginKey, removeMarginSec);
                
                for iRi = 1:nRemoveIntervals
                    iCdfRecord1 = i1Array(iRi);
                    iCdfRecord2 = i2Array(iRi);
                    utc1  = EJ_library.cdf.tt2000_to_UTC_str(zvEpoch(iCdfRecord1));
                    utc2  = EJ_library.cdf.tt2000_to_UTC_str(zvEpoch(iCdfRecord2));
                    L.logf(LL, '    Records %7i-%7i, %s--%s', iCdfRecord1, iCdfRecord2, utc1, utc2);
                end
                
                %=============
                % Remove data
                %=============
                DemuxerOutput = structfun(...
                    @(x) (bicas.proc_utils.filter_rows(x, bRemoveArray)), ...
                    DemuxerOutput, 'UniformOutput', false);
                zvCurrentAAmpere(bRemoveArray, :) = NaN;     % BUG: Sets variable that does not exist!   ???
            end            
        end
    
    

        % Demultiplex and calibrate voltages.
        %
        % NOTE: Can handle arrays of any size as long as the sizes are consistent.
        %
        function AsrSamplesAVolt = calibrate_demux_voltages(PreDc, Cal, SETTINGS, L)
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
        % subsequences based on mux mode and latching relay, nothing else.
        %   PROPOSAL: Separate out demultiplexer. Do not call from this function.
        %
        % PROPOSAL: Function for dtSec.
        %     PROPOSAL: Some kind of assertion (assumption of) constant sampling frequency.
        %
        % PROPOSAL: Move the different conversion of CWF/SWF (one/many cell arrays) into the calibration function?!!
        %
        % PROPOSAL: Move processing of one subsequence (one for loop iteration) into one function.

            % ASSERTIONS
            assert(isscalar(PreDc.hasSnapshotFormat))
            assert(iscell(PreDc.Zv.samplesCaTm))
            EJ_library.assert.vector(PreDc.Zv.samplesCaTm)
            assert(numel(PreDc.Zv.samplesCaTm) == 5)
            bicas.proc_utils.assert_cell_array_comps_have_same_N_rows(PreDc.Zv.samplesCaTm)
            EJ_library.assert.all_equal([...
                size(PreDc.Zv.MUX_SET,        1), ...
                size(PreDc.Zv.DIFF_GAIN,      1), ...
                size(PreDc.Zv.samplesCaTm{1}, 1)])



            % Create empty 1x1 structure to which new array components can be added.
            % NOTE: Unit is avolt. Not including unit in the field names to keep them short.
            AsrSamplesAVolt = EJ_library.utils.empty_struct([1], ...
                'dcV1',  'dcV2',  'dcV3', ...
                'dcV12', 'dcV23', 'dcV13', ...
                'acV12', 'acV23', 'acV13');

            dlrUsing12zv = bicas.demultiplexer_latching_relay(PreDc.Zv.Epoch);
            iCalibLZv    = Cal.get_calibration_time_L(        PreDc.Zv.Epoch);
            iCalibHZv    = Cal.get_calibration_time_H(        PreDc.Zv.Epoch);

            
            
            if isfield(PreDc.Zv, 'CALIBRATION_TABLE_INDEX')
                % NOTE: CALIBRATION_TABLE_INDEX exists for L1R, but not L1.
                
                % ASSERTION
                % NOTE: Checks both LFR & TDS CDF files.
                % NOTE: Has observed breaking assertion in LFR test files "LFR___LFR_suggested_2019-01-17".
                % PROPOSAL: Abolish somehow.
                
                zv_CALIBRATION_TABLE_INDEX = PreDc.Zv.CALIBRATION_TABLE_INDEX;
%                 hasLegalCtiSize = size(zv_CALIBRATION_TABLE_INDEX, 2) == 2;
%                 if ~hasLegalCtiSize
%                     [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.L1R.ZV_CALIBRATION_TABLE_INDEX_ILLEGAL_SIZE_REPLACE');
%                     if settingValue
%                         L.log('warning', ...
%                             'Setting CALIBRATION_TABLE_INDEX to NaN due to setting %s = "%i".', settingKey, settingValue)
%                         zv_CALIBRATION_TABLE_INDEX = zeros(PreDc.nRecords, 2) * NaN;
%                     else
%                         error('BICAS:proc_sub:Assertion', 'zVar CALIBRATION_TABLE_INDEX has illegal width=%i (<>2).', ...
%                             size(zv_CALIBRATION_TABLE_INDEX, 2))
%                     end
%                 end
            else
                % NOTE: Technically, this should only happen for L1 input.
                
                % Create "empty" zv_CALIBRATION_TABLE_INDEX.
                L.log('warning', 'Creating NaN-valued CALIBRATION_TABLE_INDEX due to zVar not being present in input CDF.')
                zv_CALIBRATION_TABLE_INDEX = zeros(PreDc.nRecords, 2) * NaN;
            end
            
            
            
            %===================================================================================
            % (1) Find continuous subsequences of records with identical settings.
            % (2) Process data separately for each such sequence.
            % NOTE: Just finding continuous subsequences can take a significant amount of time.
            %===================================================================================
            [iEdgeList] = bicas.proc_utils.find_constant_sequences(...
                PreDc.Zv.MUX_SET, ...
                PreDc.Zv.DIFF_GAIN, ...
                dlrUsing12zv, ...
                PreDc.Zv.freqHz, ...
                iCalibLZv, ...
                iCalibHZv, ...
                PreDc.Zv.iLsf, ...
                zv_CALIBRATION_TABLE_INDEX, ...
                PreDc.Zv.useFillValues);
            [iFirstList, iLastList] = bicas.proc_utils.index_edges_2_first_last(iEdgeList);
            
            for iSubseq = 1:length(iFirstList)

                iFirst = iFirstList(iSubseq);
                iLast  = iLastList (iSubseq);

                % Extract SCALAR settings to use for entire subsequence of records.
                % SS = Subsequence (single, constant value valid for entire subsequence)
                MUX_SET_ss                 = PreDc.Zv.MUX_SET  (        iFirst);
                DIFF_GAIN_ss               = PreDc.Zv.DIFF_GAIN(        iFirst);
                dlrUsing12_ss              = dlrUsing12zv(              iFirst);
                freqHz_ss                  = PreDc.Zv.freqHz(           iFirst);
                iCalibL_ss                 = iCalibLZv(                 iFirst);
                iCalibH_ss                 = iCalibHZv(                 iFirst);
                iLsf_ss                    = PreDc.Zv.iLsf(             iFirst);
                useFillValues_ss           = PreDc.Zv.useFillValues(    iFirst);
                CALIBRATION_TABLE_INDEX_ss = zv_CALIBRATION_TABLE_INDEX(iFirst, :);
                
                % PROPOSAL: Make into "proper" table.
                %   NOTE: Can not use EJ_library.utils.assist_print_table since it requires the entire table to pre-exist.
                %   PROPOSAL: Print after all iterations.
                L.logf('info', ['Records %7i-%7i : %s -- %s', ...
                    ' MUX_SET=%i; DIFF_GAIN=%i; dlrUsing12=%i; freqHz=%5g; iCalibL=%i; iCalibH=%i;', ...
                    ' CALIBRATION_TABLE_INDEX=[%i, %i]; useFillValues=%g'], ...
                    iFirst, iLast, ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iLast)), ...
                    MUX_SET_ss, DIFF_GAIN_ss, dlrUsing12_ss, freqHz_ss, iCalibL_ss, iCalibH_ss, ...
                    CALIBRATION_TABLE_INDEX_ss(1), ...
                    CALIBRATION_TABLE_INDEX_ss(2), ...
                    useFillValues_ss)

                %============================================
                % FIND DEMUXER ROUTING, BUT DO NOT CALIBRATE
                %============================================
                % NOTE: Call demultiplexer with no samples. Only for collecting information on which BLTS channels are
                % connected to which ASRs.
                [BltsSrcAsrArray, ~] = bicas.demultiplexer.main(MUX_SET_ss, dlrUsing12_ss, {[],[],[],[],[]});



                % Extract subsequence of DATA records to "demux".
                ssSamplesTm                = bicas.proc_utils.select_row_range_from_cell_comps(PreDc.Zv.samplesCaTm, iFirst, iLast);
                % NOTE: "zVariable" (i.e. first index=record) for only the current subsequence.
                ssZvNValidSamplesPerRecord = PreDc.Zv.nValidSamplesPerRecord(iFirst:iLast);
                if PreDc.hasSnapshotFormat
                    % NOTE: Vector of constant numbers (one per snapshot).
                    ssDtSec = 1 ./ PreDc.Zv.freqHz(iFirst:iLast);
                else
                    % NOTE: Scalar (one for entire sequence).
                    ssDtSec = double(PreDc.Zv.Epoch(iLast) - PreDc.Zv.Epoch(iFirst)) / (iLast-iFirst) * 1e-9;   % TEMPORARY
                end
                
                biasHighGain = DIFF_GAIN_ss;    % NOTE: Not yet sure that this is correct.



                %=======================
                % ITERATE OVER CHANNELS
                %=======================
                ssSamplesAVolt = cell(5,1);
                for iBlts = 1:5

                    if strcmp(BltsSrcAsrArray(iBlts).category, 'Unknown')
                        % Calibrated data is NaN.
                        ssSamplesAVolt{iBlts} = NaN * zeros(size(ssSamplesTm{iBlts}));

                    elseif ismember(BltsSrcAsrArray(iBlts).category, {'GND', '2.5V Ref'})
                        % ==> No calibration.
                        ssSamplesAVolt{iBlts} = ssSamplesTm{iBlts};
                        
                    else
                        assert(BltsSrcAsrArray(iBlts).is_ASR())
                        % ==> Calibrate
                        
                        if PreDc.hasSnapshotFormat
                            ssSamplesCaTm = bicas.proc_utils.convert_matrix_to_cell_array_of_vectors(...
                                double(ssSamplesTm{iBlts}), ssZvNValidSamplesPerRecord);
                        else
                            assert(all(ssZvNValidSamplesPerRecord == 1))
                            ssSamplesCaTm = {double(ssSamplesTm{iBlts})};
                        end
                        
                        if useFillValues_ss
                            % Create equivalent cell array of NaN-valued sample sequences.
                            ssSamplesCaAVolt = cell(size(ssSamplesCaTm));
                            for i = 1:numel(ssSamplesCaTm)
                                ssSamplesCaAVolt{i} = ones(size(ssSamplesCaTm{i})) * NaN;
                            end
                        else
                            %%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%
                            %  CALIBRATE
                            %%%%%%%%%%%%%%
                            %%%%%%%%%%%%%%
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
                        end
                        
                        if PreDc.hasSnapshotFormat
                            [ssSamplesAVolt{iBlts}, ~] = bicas.proc_utils.convert_cell_array_of_vectors_to_matrix(...
                                ssSamplesCaAVolt, size(ssSamplesTm{iBlts}, 2));
                        else
                            ssSamplesAVolt{iBlts} = ssSamplesCaAVolt{1};   % NOTE: Must be column array.
                        end
                    end
                end
                
                %====================
                % CALL DEMULTIPLEXER
                %====================
                [~, SsAsrSamplesAVolt] = bicas.demultiplexer.main(MUX_SET_ss, dlrUsing12_ss, ssSamplesAVolt);
                
                % Add demuxed sequence to the to-be complete set of records.
                %DemuxerOutput = bicas.proc_utils.add_rows_to_struct_fields(DemuxerOutput, DemuxerOutputSubseq);
                AsrSamplesAVolt = bicas.proc_utils.add_rows_to_struct_fields(AsrSamplesAVolt, SsAsrSamplesAVolt);
                
            end
            
        end   % simple_demultiplex


        
        % Wrapper around bicas.proc_sub.handle_struct_name_change to be used locally.
        %
        % ARGUMENTS
        % =========
        % inSciDsi : Input SCI DATASET_ID which contains the zVariable.
        %
        function handle_zv_name_change(fnChangeList, inSciDsi, SETTINGS, L, varargin)
            anomalyDescrMsgFunc = @(oldFieldname, newFieldname) (sprintf(...
                'Input dataset DATASET_ID=%s uses an alternative but illegal(?) zVariable name "%s" instead of "%s".', ...
                inSciDsi, oldFieldname, newFieldname));
            
            bicas.handle_struct_name_change(fnChangeList, SETTINGS, L, anomalyDescrMsgFunc, varargin{:})
        end



    end   % methods(Static, Access=private)
        
end
