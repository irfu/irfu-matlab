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
% PROPOSAL: Use double for all numeric zVariables in the processing. Do not produce or require proper type, e.g. integers, in any
%           intermediate processing. Only convert to the proper data type/class when writing to CDF.
%   PRO: Variables can keep NaN to represent fill/pad value, also for "integers".
%   PRO: The knowledge of the dataset CDF formats is not spread out over the code.
%       Ex: Setting default values for PreDc.QUALITY_FLAG, PreDc.QUALITY_BITMASK, PreDc.DELTA_PLUS_MINUS.
%       Ex: ACQUISITION_TIME.
%   CON: Less assertions can be made in utility functions.
%       Ex: proc_utils.ACQUISITION_TIME_*, proc_utils.tt2000_* functions.
%   CON: ROUNDING ERRORS. Can not be certain that values which are copied, are actually copied.
%   --
%   NOTE: Functions may in principle require integer math to work correctly.
% --
% PROPOSAL: Derive DIFF_GAIN (from BIAS HK using time interpolation) in one code common to both LFR & TDS.
%   PROPOSAL: Function
%   PRO: Uses flag for selecting interpolation time in one place.
% PROPOSAL: Derive HK_BIA_MODE_MUX_SET (from BIAS SCI or HK using time interpolation for HK) in one code common to both LFR & TDS.
%   PROPOSAL: Function
%   PRO: Uses flag for selecting HK/SCI DIFF_GAIN in one place.
%   PRO: Uses flag for selecting interpolation time in one place.
%--
% NOTE: Both BIAS HK and LFR SURV CWF contain MUX data (only LFR has one timestamp per snapshot). True also for other input datasets?
%
% PROPOSAL: Instantiate class, use instance methods instead of static.
%   PRO: Can have SETTINGS and constants as instance variable instead of calling global variables.
%
% PROPOSAL: Submit zVar variable attributes.
%   PRO: Can interpret fill values.
%       Ex: Can doublecheck TDS RSWF snapshot length using fill values and compare with zVar SAMPS_PER_CH (which seems to be
%       bad).
% PROPOSAL: Return (to execute_sw_mode), global attributes.
%   PRO: Needed for output datasets: CALIBRATION_TABLE, CALIBRATION_VERSION
%       ~CON: CALIBRATION_VERSION refers to algorithm and should maybe be a SETTING.
%
% PROPOSAL: Separate LFR and TDS in different files.
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
            bicas.proc_utils.log_zVars(TimeVars, L);



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
                %error('BICAS:proc_sub:DatasetFormat:SWModeProcessing', 'SCI and HK time ranges do not overlap in time.')
                
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
%             HkSciTime.MUX_SET = bicas.proc_utils.interpolate_float_records(...
%                 hkEpoch, ...
%                 double(InHk.Zv.HK_BIA_MODE_MUX_SET), ...
%                 InSci.Zv.Epoch, ...
%                 'nearest');
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
%             HkSciTime.DIFF_GAIN = bicas.proc_utils.interpolate_float_records(...
%                 hkEpoch, ...
%                 double(InHk.Zv.HK_BIA_DIFF_GAIN), ...
%                 InSci.Zv.Epoch, ...
%                 'nearest');
            HkSciTime.DIFF_GAIN = bicas.utils.interpolate_nearest(...
                hkEpochExtrapMargin, ...
                hkEpoch, ...
                InHk.Zv.HK_BIA_DIFF_GAIN, ...
                InSci.Zv.Epoch);



            % ASSERTIONS
            EJ_library.assert.struct(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
        end
        
        
        
        function currentSAmpere = process_CUR_to_CUR_on_SCI_TIME(sciEpoch, InCur, SETTINGS, L)
            % ASSERTIONS
            EJ_library.assert.struct(InCur, {'Zv', 'Ga'}, {})
            
            
            
            if SETTINGS.get_fv('INPUT_CDF.CUR.PREPEND_TEST_DATA')
                L.log('warning', '==========================================================')
                L.log('warning', 'WARNING: PREPENDING MADE-UP TEST DATA TO BIAS CURRENT DATA')
                L.log('warning', '==========================================================')
                
                %==========================================
                % ADD (PREPEND) TEST DATA TO BIAS CURRENTS
                %==========================================
                j  = 0:20;
                i3 = 1:3:(numel(j)-2);   % Useful for setting components to NaN.
                % One timestamp every minute.
                % NOTE: Must have timestamps with non-NaN bias values beginning before first voltage Epoch value.
                EpochAmend = (sciEpoch(1) + int64(j *60*1e9 - 5*60*1e9))';
                IBIAS_amend_1 =  mod(j,60)';
                IBIAS_amend_2 = -mod(j,60)';
                IBIAS_amend_3 =  mod(j,60)';
                IBIAS_amend_1(i3+1) = NaN;
                IBIAS_amend_1(i3+2) = NaN;
                IBIAS_amend_2(i3+0) = NaN;
                IBIAS_amend_2(i3+2) = NaN;
                IBIAS_amend_3(i3+0) = NaN;
                IBIAS_amend_3(i3+1) = NaN;
                
                InCur.Zv.Epoch   = [EpochAmend    ; InCur.Zv.Epoch];
                InCur.Zv.IBIAS_1 = [IBIAS_amend_1 ; InCur.Zv.IBIAS_1];
                InCur.Zv.IBIAS_2 = [IBIAS_amend_2 ; InCur.Zv.IBIAS_2];
                InCur.Zv.IBIAS_3 = [IBIAS_amend_3 ; InCur.Zv.IBIAS_3];
                assert(issorted(InCur.Zv.Epoch), 'TEST CODE failed')
                EJ_library.assert.all_equal([...
                    size(InCur.Zv.Epoch,   1), ...
                    size(InCur.Zv.IBIAS_1, 1), ...
                    size(InCur.Zv.IBIAS_2, 1), ...
                    size(InCur.Zv.IBIAS_3, 1)])
            end



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

            nRecords = size(InSci.Zv.Epoch, 1);            
            C = EJ_library.so.adm.classify_DATASET_ID(inSciDsi);
            
            
            
            % CDF ASSERTION
            if ~issorted(InSci.Zv.Epoch, 'strictascend')
                error('Voltage (science) dataset timestamps do not increase.')
            end
            
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



            V = InSci.Zv.V;
            E = permute(InSci.Zv.E, [1,3,2]);
            % Switch last two indices of E.
            % ==> index 2 = "snapshot" sample index, including for CWF (sample/record, "snapshots" consisting of 1 sample).
            %     index 3 = E1/E2 component
            %               NOTE: 1/2=index into array; these are diffs but not equivalent to any particular diffs).
            
            nCdfSamplesPerRecord = size(V, 2);
            
            % ASSERTIONS
            if C.isLfrSurvSwf
                assert(nCdfSamplesPerRecord == EJ_library.so.constants.LFR_SWF_SNAPSHOT_LENGTH)
            else
                assert(nCdfSamplesPerRecord == 1)
            end
            assert(size(E, 3) == 2)



            % Set iLsfZv.
            if     C.isLfrSbm1
                iLsfZv = ones(nRecords, 1) * 2;   % Always value "2" (F1, "FREQ = 1").
            elseif C.isLfrSbm2
                iLsfZv = ones(nRecords, 1) * 3;   % Always value "3" (F2, "FREQ = 2").
            else
                % NOTE: Translates from LFR's FREQ values (0=F0 etc) to LSF index values (1=F0) used in loaded RCT data structs.
                iLsfZv = InSci.Zv.FREQ + 1;
            end
            EJ_library.assert.size(iLsfZv, [NaN, 1])



            % NOTE: Needed also for 1 SPR.
            zvFreqHz = EJ_library.so.get_LFR_frequency( iLsfZv );

            % Obtain the relevant values (one per record) from zVariables R0, R1, R2, and the virtual "R3".
            zvRx = bicas.proc_utils.get_LFR_Rx(...
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
            % Test data: MYSTERIOUS_SIGNAL_1_2016-04-15_Run2__7729147__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-SWF_7729147_CNE_V01.cdf
            %
            % PROPOSAL: Move to the code that reads CDF datasets instead. Generalize to many zVariables.
            %===========================================================================================================
            PreDc.Zv.QUALITY_FLAG    = InSci.Zv.QUALITY_FLAG;
            PreDc.Zv.QUALITY_BITMASK = InSci.Zv.QUALITY_BITMASK;
            if isempty(PreDc.Zv.QUALITY_FLAG)
                L.log('warning', 'QUALITY_FLAG from the LFR SCI source dataset is empty. Filling with empty values.')
                PreDc.Zv.QUALITY_FLAG = bicas.proc_utils.create_NaN_array([nRecords, 1]);
            end
            if isempty(PreDc.Zv.QUALITY_BITMASK)
                L.log('warning', 'QUALITY_BITMASK from the LFR SCI source dataset is empty. Filling with empty values.')
                PreDc.Zv.QUALITY_BITMASK = bicas.proc_utils.create_NaN_array([nRecords, 1]);
            end
            
            % ASSERTIONS
            % LFR QUALITY_FLAG, QUALITY_BITMASK not set yet (2019-09-17), but I presume they should have just one value
            % per record. BIAS output datasets should.
            assert(size(PreDc.Zv.QUALITY_FLAG,    2) == 1)
            assert(size(PreDc.Zv.QUALITY_BITMASK, 2) == 1)



            % E must be floating-point so that values can be set to NaN.
            % ==> V must be floating-point so that demultiplexer can subtract/add V/E, and calibration can work.
            % bicas.proc_utils.filter_rows requires this. Variable may be integer if integer in source CDF.
            V = single(V);
            E = single(E);

            PreDc.Zv.samplesCaTm    = {};
            PreDc.Zv.samplesCaTm{1} = V;
            PreDc.Zv.samplesCaTm{2} = bicas.proc_utils.filter_rows( E(:,:,1), zvRx==1 );
            PreDc.Zv.samplesCaTm{3} = bicas.proc_utils.filter_rows( E(:,:,2), zvRx==1 );
            PreDc.Zv.samplesCaTm{4} = bicas.proc_utils.filter_rows( E(:,:,1), zvRx==0 );
            PreDc.Zv.samplesCaTm{5} = bicas.proc_utils.filter_rows( E(:,:,2), zvRx==0 );
            
            
            
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
            
            freqHz = double(InSci.Zv.SAMPLING_RATE);
            
            if any(freqHz == 255)
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
                            freqHz(freqHz == 255) = 32768;
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
            PreDc.Zv.DELTA_PLUS_MINUS = bicas.proc_utils.derive_DELTA_PLUS_MINUS(freqHz, nCdfSamplesPerRecord);
            PreDc.Zv.freqHz           = freqHz;
            PreDc.Zv.QUALITY_FLAG     = InSci.Zv.QUALITY_FLAG;
            PreDc.Zv.QUALITY_BITMASK  = InSci.Zv.QUALITY_BITMASK;
            PreDc.Zv.SYNCHRO_FLAG     = InSci.Zv.SYNCHRO_FLAG;
            PreDc.Zv.MUX_SET          = HkSciTime.MUX_SET;
            PreDc.Zv.DIFF_GAIN        = HkSciTime.DIFF_GAIN;
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
                SAMPS_PER_CH_MIN_VALID = 2^10;
                SAMPS_PER_CH_MAX_VALID = 2^15;
                SAMPS_PER_CH_zv        = double(InSci.Zv.SAMPS_PER_CH);
                SAMPS_PER_CH_rounded   = round(2.^round(log2(SAMPS_PER_CH_zv)));
                SAMPS_PER_CH_rounded(SAMPS_PER_CH_rounded < SAMPS_PER_CH_MIN_VALID) = SAMPS_PER_CH_MIN_VALID;
                SAMPS_PER_CH_rounded(SAMPS_PER_CH_rounded > SAMPS_PER_CH_MAX_VALID) = SAMPS_PER_CH_MAX_VALID;
                if any(SAMPS_PER_CH_rounded ~= SAMPS_PER_CH_zv)
                    SAMPS_PER_CH_badValues = unique(SAMPS_PER_CH_zv(SAMPS_PER_CH_rounded ~= SAMPS_PER_CH_zv));
                    
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
                            
                            SAMPS_PER_CH_zv = SAMPS_PER_CH_rounded;
                            
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
                PreDc.Zv.nValidSamplesPerRecord = SAMPS_PER_CH_zv;
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
                'MUX_SET', 'QUALITY_FLAG', 'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG'}, ...
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
                'currentAAmpere', 'DemuxerOutput'}, ...
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
            OutSciZv.IBIAS1 = SciPostDc.Zv.currentAAmpere(:, 1) * 1e9;
            OutSciZv.IBIAS2 = SciPostDc.Zv.currentAAmpere(:, 2) * 1e9;
            OutSciZv.IBIAS3 = SciPostDc.Zv.currentAAmpere(:, 3) * 1e9;
            
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
                    assert(nSamplesPerRecord == 2048, ...
                        'BICAS:proc_sub:Assertion:IllegalArgument', ...
                        'Number of samples per CDF record is not 2048, as expected. Bad Input CDF?')
                    
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
            
            %nSamplesPerRecord = size(SciPostDc.Zv.DemuxerOutput.dcV1, 2);   % Samples per record.
            
            OutSciZv.Epoch            = SciPostDc.Zv.Epoch;
            OutSciZv.QUALITY_FLAG     = SciPostDc.Zv.QUALITY_FLAG;
            OutSciZv.QUALITY_BITMASK  = SciPostDc.Zv.QUALITY_BITMASK;
            OutSciZv.DELTA_PLUS_MINUS = SciPostDc.Zv.DELTA_PLUS_MINUS;
            OutSciZv.IBIAS1           = SciPostDc.Zv.currentAAmpere(:, 1);
            OutSciZv.IBIAS2           = SciPostDc.Zv.currentAAmpere(:, 2);
            OutSciZv.IBIAS3           = SciPostDc.Zv.currentAAmpere(:, 3);
            OutSciZv.SYNCHRO_FLAG     = SciPostDc.Zv.SYNCHRO_FLAG;
            OutSciZv.SAMPLING_RATE    = SciPostDc.Zv.freqHz;

            OutSciZv.IBIAS1     = SciPostDc.Zv.currentAAmpere(:, 1) * 1e9;
            OutSciZv.IBIAS2     = SciPostDc.Zv.currentAAmpere(:, 2) * 1e9;
            OutSciZv.IBIAS3     = SciPostDc.Zv.currentAAmpere(:, 3) * 1e9;
            
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
        % Function is in large part a wrapper around "simple_demultiplex".
        %
        % NOTE: Public function as opposed to the other demuxing/calibration functions.
        %
        function PostDc = process_demuxing_calibration(PreDc, InCurPd, Cal, SETTINGS, L)
        % PROPOSAL: Move the setting of IBIASx (bias current) somewhere else?
        %   PRO: Unrelated to demultiplexing.
        %   CON: Related to calibration.

            % ASSERTION
            bicas.proc_sub.assert_PreDC(PreDc);

            %=======
            % DEMUX
            %=======
            PostDc = PreDc;    % Copy all values, to later overwrite a subset of them.
            
            nRecords = size(PreDc.Zv.Epoch, 1);
            tTicToc = tic();            
            
            PostDc.Zv.DemuxerOutput = bicas.proc_sub.simple_demultiplex(PreDc, Cal, SETTINGS, L);
            
            wallTimeSec = toc(tTicToc);
            L.logf('info', 'simple_demultiplex: Time used for execution (wall time): %g [s], %g [s/record]', wallTimeSec, wallTimeSec/nRecords)
            
            
            
            %================================
            % Set (calibrated) bias currents
            %================================
            % BUG / TEMP: Set default values since the real bias current values are not available.
            currentSAmpere = bicas.proc_sub.process_CUR_to_CUR_on_SCI_TIME(PreDc.Zv.Epoch, InCurPd, SETTINGS, L);
            currentTm      = bicas.calib.calibrate_set_current_to_bias_current(currentSAmpere);
            
            currentAAmpere = bicas.proc_utils.create_NaN_array(size(currentSAmpere));    % Variable to fill/set.
            iCalibLZv      = Cal.get_calibration_time_L(PreDc.Zv.Epoch);
            iEdgeList      = bicas.proc_utils.find_constant_sequences(iCalibLZv);
            [iFirstList, iLastList] = bicas.proc_utils.index_edges_2_first_last(iEdgeList);
            for iSubseq = 1:length(iFirstList)
                iRecords = iFirstList(iSubseq) : iLastList(iSubseq);
                
                for iAnt = 1:3
                    currentAAmpere(iRecords, iAnt) = Cal.calibrate_TC_bias_TM_to_bias_current(...
                        currentTm(iRecords, iAnt), iAnt, iCalibLZv(iRecords));
                end
            end
            
            PostDc.Zv.currentAAmpere = currentAAmpere;

            % ASSERTION
            bicas.proc_sub.assert_PostDC(PostDc)
        end
        
    end    % methods(Static, Access=public)
            
    %###################################################################################################################
    
    methods(Static, Access=private)
    %methods(Static, Access=public)

        % Wrapper around "simple_demultiplex_subsequence_OLD" to be able to handle multiple CDF records with changing
        % settings.
        %
        % NOTE: Can handle arrays of any size as long as the sizes are consistent.
        %
        function AsrSamplesAVolt = simple_demultiplex(PreDc, Cal, SETTINGS, L)
        % PROPOSAL: Incorporate into processing function process_demuxing_calibration.
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
            % NOTE: Unit is AVolt. Not including in the field names to keep them short.
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
                
                CALIBRATION_TABLE_INDEX_zv = PreDc.Zv.CALIBRATION_TABLE_INDEX;
                hasLegalCtiSize = size(CALIBRATION_TABLE_INDEX_zv, 2) == 2;
                if ~hasLegalCtiSize
                    [settingValue, settingKey] = SETTINGS.get_fv('PROCESSING.L1R.ZV_CALIBRATION_TABLE_INDEX_ILLEGAL_SIZE_REPLACE');
                    if settingValue
                        L.log('warning', ...
                            'Setting CALIBRATION_TABLE_INDEX to NaN due to setting %s = "%i".', settingKey, settingValue)
                        CALIBRATION_TABLE_INDEX_zv = zeros(PreDc.nRecords, 2) * NaN;
                    else
                        error('BICAS:proc_sub:Assertion', 'zVar CALIBRATION_TABLE_INDEX has illegal width=%i (<>2).', size(CALIBRATION_TABLE_INDEX_zv, 2))
                    end
                end
            else
                % NOTE: Technically, this should only happen for L1 input.
                
                % Create "empty" CALIBRATION_TABLE_INDEX_zv.
                L.log('warning', 'Creating NaN-valued CALIBRATION_TABLE_INDEX due to zVar not being present in input CDF.')
                CALIBRATION_TABLE_INDEX_zv = zeros(PreDc.nRecords, 2) * NaN;
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
                CALIBRATION_TABLE_INDEX_zv);
            [iFirstList, iLastList] = bicas.proc_utils.index_edges_2_first_last(iEdgeList);
            
            for iSubseq = 1:length(iFirstList)

                iFirst = iFirstList(iSubseq);
                iLast  = iLastList (iSubseq);

                % Extract SCALAR settings to use for entire subsequence of records.
                % SS = Subsequence (single, constant value valid for entire subsequence)
                MUX_SET_ss                 = PreDc.Zv.MUX_SET  (iFirst);
                DIFF_GAIN_ss               = PreDc.Zv.DIFF_GAIN(iFirst);
                dlrUsing12_ss              = dlrUsing12zv(      iFirst);
                freqHz_ss                  = PreDc.Zv.freqHz(   iFirst);
                iCalibL_ss                 = iCalibLZv(         iFirst);
                iCalibH_ss                 = iCalibHZv(         iFirst);
                iLsf_ss                    = PreDc.Zv.iLsf(     iFirst);
                CALIBRATION_TABLE_INDEX_ss = CALIBRATION_TABLE_INDEX_zv(iFirst, :);
                
                % PROPOSAL: Make into "proper" table.
                %   NOTE: Can not use EJ_library.utils.assist_print_table since it requires the entire table to pre-exist.
                %   PROPOSAL: Print after all iterations.
                L.logf('info', ['Records %7i-%7i : %s -- %s ', ...
                    'MUX_SET=%i; DIFF_GAIN=%i; dlrUsing12=%i; freqHz=%5g; iCalibL=%i; iCalibH=%i; CALIBRATION_TABLE_INDEX=[%i, %i]'], ...
                    iFirst, iLast, ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.proc_utils.tt2000_to_UTC_str(PreDc.Zv.Epoch(iLast)), ...
                    MUX_SET_ss, DIFF_GAIN_ss, dlrUsing12_ss, freqHz_ss, iCalibL_ss, iCalibH_ss, CALIBRATION_TABLE_INDEX_ss(1), CALIBRATION_TABLE_INDEX_ss(2))

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

                    elseif strcmp(BltsSrcAsrArray(iBlts).category, 'GND') || strcmp(BltsSrcAsrArray(iBlts).category, '2.5V Ref')
                        % No calibration.
                        ssSamplesAVolt{iBlts} = ssSamplesTm{iBlts};
                        
                    else
                        % PROPOSAL: Check .category and add else-assertion.
                        
                        if PreDc.hasSnapshotFormat
                            ssSamplesCaTm = bicas.proc_utils.convert_matrix_to_cell_array_of_vectors(...
                                double(ssSamplesTm{iBlts}), ssZvNValidSamplesPerRecord);
                        else
                            assert(all(ssZvNValidSamplesPerRecord == 1))
                            ssSamplesCaTm = {double(ssSamplesTm{iBlts})};
                        end
                        
                        %%%%%%%%%%%%
                        %%%%%%%%%%%%
                        % CALIBRATE
                        %%%%%%%%%%%%
                        %%%%%%%%%%%%
                        CalSettings = struct();
                        CalSettings.iBlts        = iBlts;
                        CalSettings.BltsSrc      = BltsSrcAsrArray(iBlts);
                        CalSettings.biasHighGain = biasHighGain;
                        CalSettings.iCalibTimeL  = iCalibL_ss;
                        CalSettings.iCalibTimeH  = iCalibH_ss;
                        CalSettings.iLsf         = iLsf_ss;
                        ssSamplesCaAVolt = Cal.calibrate_voltage_all(ssDtSec, ssSamplesCaTm, ...
                            PreDc.isLfr, PreDc.isTdsCwf, CalSettings, CALIBRATION_TABLE_INDEX_ss);
%                         ssSamplesCaAVolt = Cal.calibrate_voltage_all(ssDtSec, ssSamplesCaTm, ...
%                             PreDc.isLfr, PreDc.isTdsCwf, iBlts, ...
%                             BltsSrcAsrArray(iBlts), biasHighGain, iCalibL_ss, iCalibH_ss, iLsf_ss, CALIBRATION_TABLE_INDEX_ss);
                        
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
