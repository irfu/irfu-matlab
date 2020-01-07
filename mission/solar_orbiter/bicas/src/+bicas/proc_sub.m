% Class that collects "processing functions" as public static methods.
%
% This class is not meant to be instantiated.
% 
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
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
        
        function HkSciTime = process_HK_to_HK_on_SCI_TIME(InSci, InHk, SETTINGS)
        % Processing function
        
            % ASSERTIONS
            EJ_library.utils.assert.struct2(InSci, {'Zv', 'Ga'}, {})
            EJ_library.utils.assert.struct2(InHk,  {'Zv', 'Ga'}, {})
            
            HkSciTime = [];
            
            
            
            % Define local convenience variables. AT = ACQUISITION_TIME
            ACQUISITION_TIME_EPOCH_UTC = SETTINGS.get_fv('PROCESSING.ACQUISITION_TIME_EPOCH_UTC');
            
            hkAtTt2000  = bicas.proc_utils.ACQUISITION_TIME_to_tt2000(  InHk.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            sciAtTt2000 = bicas.proc_utils.ACQUISITION_TIME_to_tt2000( InSci.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            hkEpoch     = InHk.Zv.Epoch;
            sciEpoch    = InSci.Zv.Epoch;
            
            %==================================================================
            % Log time intervals to enable comparing available SCI and HK data
            %==================================================================
            bicas.proc_utils.log_tt2000_array('HK  ACQUISITION_TIME', hkAtTt2000)
            bicas.proc_utils.log_tt2000_array('SCI ACQUISITION_TIME', sciAtTt2000)
            bicas.proc_utils.log_tt2000_array('HK  Epoch           ', hkEpoch)
            bicas.proc_utils.log_tt2000_array('SCI Epoch           ', sciEpoch)
            
            
            
            % WARNINGS
            if ~(bicas.proc_utils.ranges_overlap(hkAtTt2000, sciAtTt2000))
                bicas.log('warning', 'zVar ACQUSITION_TIME in HK and SCI input datasets do not overlap in time.')
            end
            if ~(bicas.proc_utils.ranges_overlap(hkEpoch, sciEpoch))
                bicas.log('warning', 'zVar Epoch in HK and SCI input datasets do not overlap in time.')
            end



            %=========================================================================================================
            % 1) Convert time to something linear in time that can be used for processing (not storing time to file).
            % 2) Effectively also chooses which time to use for the purpose of processing:
            %       (a) ACQUISITION_TIME, or
            %       (b) Epoch.
            %=========================================================================================================
            if SETTINGS.get_fv('PROCESSING.USE_ZV_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION')
                bicas.log('info', 'Using HK & SCI zVariable ACQUISITION_TIME (not Epoch) for interpolating HK dataset data to SCI dataset time, due to setting PROCESSING.USE_ZV_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION.')
                hkInterpolationTimeTt2000  = hkAtTt2000;
                sciInterpolationTimeTt2000 = sciAtTt2000;
            else
                bicas.log('info', 'Using HK & SCI zVariable Epoch (not ACQUISITION_TIME) for interpolating HK dataset data to SCI dataset time.')
                hkInterpolationTimeTt2000  = hkEpoch;
                sciInterpolationTimeTt2000 = sciEpoch;
            end
            clear hkAtTt2000 sciAtTt2000
            clear hkEpoch    sciEpoch

            % ASSERTION
            if ~(bicas.proc_utils.ranges_overlap(hkInterpolationTimeTt2000, sciInterpolationTimeTt2000))
                error('BICAS:proc_sub:Assertion', ...
                'Time zVariables (Epoch or ACQUISITION_TIME) in HK and SCI input datasets, used for converting HK data to SCI timestamps, do not overlap in time.')
            end
            

            
            %=========================================================================================================
            % Derive MUX_SET
            % --------------
            % NOTE: Only obtains one MUX_SET per record ==> Can not change MUX_SET in the middle of a record.
            % NOTE: Can potentially obtain MUX_SET from LFR SCI.
            %=========================================================================================================            
            HkSciTime.MUX_SET = bicas.proc_utils.nearest_interpolate_float_records(...
                double(InHk.Zv.HK_BIA_MODE_MUX_SET), ...
                hkInterpolationTimeTt2000, ...
                sciInterpolationTimeTt2000);   % Use BIAS HK.
            %PreDc.MUX_SET = LFR_cdf.BIAS_MODE_MUX_SET;    % Use LFR SCI. NOTE: Only possible for ___LFR___.



            %=========================================================================================================
            % Derive DIFF_GAIN
            % ----------------
            % NOTE: Not perfect handling of time when 1 snapshot/record, since one should ideally use time stamps
            % for every LFR _sample_.
            %=========================================================================================================
            HkSciTime.DIFF_GAIN = bicas.proc_utils.nearest_interpolate_float_records(...
                double(InHk.Zv.HK_BIA_DIFF_GAIN), hkInterpolationTimeTt2000, sciInterpolationTimeTt2000);



            % ASSERTIONS
            EJ_library.utils.assert.struct2(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
        end



        function PreDc = process_LFR_to_PreDC(InSci, inSciDsi, HkSciTime)
        % Processing function. Convert LFR CDF data to PreDC.
        %
        % Keeps number of samples/record. Treats 1 samples/record "length-one snapshots".
        
        % PROBLEM: Hardcoded CDF data types (MATLAB classes).
        % MINOR PROBLEM: Still does not handle LFR zVar TYPE for determining "virtual snapshot" length.
        % Should only be relevant for V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF (not V02) which should expire.
        
            LFR_SWF_SNAPSHOT_LENGTH = 2048;
        
            % ASSERTIONS
            EJ_library.utils.assert.struct2(InSci,  {'Zv', 'Ga'}, {})
            EJ_library.utils.assert.struct2(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
            
            nRecords = size(InSci.Zv.Epoch, 1);            
            C = bicas.proc_utils.classify_DATASET_ID(inSciDsi);
           
            % Both zVars TIME_SYNCHRO_FLAG, SYNCHRO_FLAG found in input datasets (2020-01-05). Unknown why.
            % "DEFINITION BUG" in definition of datasets/skeleton?
            InSci.Zv = bicas.utils.normalize_struct_fieldnames(InSci.Zv, ...
                {{{'TIME_SYNCHRO_FLAG', 'SYNCHRO_FLAG'}, 'TIME_SYNCHRO_FLAG'}});
            
            V = InSci.Zv.V;
            E = permute(InSci.Zv.E, [1,3,2]);
            % Switch last two indices of E.
            % ==> index 2 = "snapshot" sample index, including for CWF (sample/record, "snapshots" consisting of 1 sample).
            %     index 3 = E1/E2 component
            %               NOTE: 1/2=index into array; these are diffs but not equivalent to any particular diffs).
            
            nCdfSamplesPerRecord = size(V, 2);
            
            % ASSERTIONS
            if C.isLfrSwf
                assert(nCdfSamplesPerRecord == LFR_SWF_SNAPSHOT_LENGTH)
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
            EJ_library.utils.assert.size(iLsfZv, [NaN, 1])
            %assert(size(iLsfZv, 2) == 1)
            
            
            
            zvFreqHz = bicas.proc_utils.get_LFR_frequency( iLsfZv );   % NOTE: Needed also for 1 SPR.

            % Obtain the relevant values (one per record) from zVariables R0, R1, R2, and the virtual "R3".
            zvRx = bicas.proc_utils.get_LFR_Rx( ...
                InSci.Zv.R0, ...
                InSci.Zv.R1, ...
                InSci.Zv.R2, ...
                iLsfZv );   % NOTE: Function also handles the imaginary zVar "R3".

            PreDc = [];
            PreDc.Zv.Epoch                  = InSci.Zv.Epoch;
            PreDc.Zv.ACQUISITION_TIME       = InSci.Zv.ACQUISITION_TIME;
            PreDc.Zv.DELTA_PLUS_MINUS       = bicas.proc_utils.derive_DELTA_PLUS_MINUS(zvFreqHz, nCdfSamplesPerRecord);            
            PreDc.Zv.freqHz                 = zvFreqHz;
            PreDc.Zv.nValidSamplesPerRecord = ones(nRecords, 1) * nCdfSamplesPerRecord;
            PreDc.Zv.SYNCHRO_FLAG           = InSci.Zv.TIME_SYNCHRO_FLAG;   % NOTE: Different zVar name in input and output datasets.
            if isfield(InSci.Zv, 'CALIBRATION_TABLE_INDEX')
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
                bicas.log('warning', 'QUALITY_FLAG from the LFR SCI source dataset is empty. Filling with empty values.')
                PreDc.Zv.QUALITY_FLAG = bicas.proc_utils.create_NaN_array([nRecords, 1]);
            end
            if isempty(PreDc.Zv.QUALITY_BITMASK)
                bicas.log('warning', 'QUALITY_BITMASK from the LFR SCI source dataset is empty. Filling with empty values.')
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
            
            PreDc.Zv.MUX_SET        = HkSciTime.MUX_SET;
            PreDc.Zv.DIFF_GAIN      = HkSciTime.DIFF_GAIN;            
            
            PreDc.hasSnapshotFormat    = C.isLfrSwf;
            PreDc.nRecords             = nRecords;
            PreDc.nCdfSamplesPerRecord = nCdfSamplesPerRecord;
            PreDc.isLfr                = true;
            PreDc.isTdsCwf             = false;
            PreDc.Zv.iLsf              = iLsfZv;



            % ASSERTIONS
            bicas.proc_sub.assert_PreDC(PreDc)
        end
        
        
        
        function PreDc = process_TDS_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS)
        % Processing function. Convert TDS CDF data (PDs) to PreDC.
        %
        % Keeps number of samples/record. Treats 1 samples/record "length-one snapshots".
        %
        % BUG?: Does not use CHANNEL_STATUS_INFO.
        % NOTE: BIAS output datasets do not have a variable for the length of snapshots. Need to use NaN/fill value.

            % ASSERTIONS
            EJ_library.utils.assert.struct2(InSci,     {'Zv', 'Ga'}, {})
            EJ_library.utils.assert.struct2(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})

            C = bicas.proc_utils.classify_DATASET_ID(inSciDsi);

            % Both zVars TIME_SYNCHRO_FLAG, SYNCHRO_FLAG found in input datasets (2019-12-xx). Unknown why.
            % "DEFINITION BUG" in definition of datasets/skeleton?
            InSci.Zv = bicas.utils.normalize_struct_fieldnames(InSci.Zv, ...
                {{{'TIME_SYNCHRO_FLAG', 'SYNCHRO_FLAG'}, 'TIME_SYNCHRO_FLAG'}});

            nRecords             = size(InSci.Zv.Epoch, 1);
            nCdfSamplesPerRecord = size(InSci.Zv.WAVEFORM_DATA, 3);    % Number of samples in the zVariable, not necessarily actual data.

            freqHz = double(InSci.Zv.SAMPLING_RATE);
            if C.isL1R && C.isTdsRswf && SETTINGS.get_fv('PROCESSING.L1R.TDS.RSWF_L1R_ZV_SAMPLING_RATE_DATASET_BUGFIX_ENABLED')
                % TEMPORARY
                % IMPLEMENTATION NOTE: Has observed test file
                % TESTDATA_RGTS_TDS_CALBA_V0.8.5C: solo_L1R_rpw-tds-lfm-rswf-e_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
                % to have SAMPLING_RATE == 255, which is likely a BUG in the dataset. /Erik P G Johansson 2019-12-03
                % Bug in TDS RCS.  /David Pisa 2019-12-03
                % Setting it to what is probably the correct value.
                freqHz(freqHz == 255) = 32768;
                bicas.logf('warning', 'Correcting presumed bug in TDS L1R LFM-RSWF dataset due to setting PROCESSING.L1R.TDS.RSWF_L1R_ZV_SAMPLING_RATE_DATASET_BUGFIX_ENABLED. Modifying the frequency 255-->32768.')
            end
            
            PreDc = [];
            
            PreDc.Zv.Epoch            = InSci.Zv.Epoch;
            PreDc.Zv.ACQUISITION_TIME = InSci.Zv.ACQUISITION_TIME;
            PreDc.Zv.DELTA_PLUS_MINUS = bicas.proc_utils.derive_DELTA_PLUS_MINUS(freqHz, nCdfSamplesPerRecord);
            PreDc.Zv.freqHz           = freqHz;
            PreDc.Zv.QUALITY_FLAG     = InSci.Zv.QUALITY_FLAG;
            PreDc.Zv.QUALITY_BITMASK  = InSci.Zv.QUALITY_BITMASK;
            PreDc.Zv.SYNCHRO_FLAG     = InSci.Zv.TIME_SYNCHRO_FLAG;   % NOTE: Different zVar name in input and output datasets.
            PreDc.Zv.MUX_SET          = HkSciTime.MUX_SET;
            PreDc.Zv.DIFF_GAIN        = HkSciTime.DIFF_GAIN;
            if isfield(InSci.Zv, 'CALIBRATION_TABLE_INDEX')
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
                    badValuesDisplayStr = strjoin(arrayfun(@(n) sprintf('%i', n), SAMPS_PER_CH_badValues, 'uni', false), ', ');                    
                    logErrorMsg = sprintf('TDS LFM RSWF zVar SAMPS_PER_CH contains unexpected value(s), not 2^n: %s', badValuesDisplayStr);
                    
                    actionSettingValue = SETTINGS.get_fv('PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_ACTION');
                    switch(actionSettingValue)
                        case 'ERROR'
                            error('BICAS:proc_sub:Assertion:DatasetFormat', logErrorMsg)
                        case 'PERMIT'
                            bicas.logf('warning', [logErrorMsg, 'Permitting due to setting PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_ACTION.'])
                            % Do nothing
                        case 'ROUND'
                            bicas.logf('warning', logErrorMsg)
                            bicas.log('warning', 'Replacing TDS RSWF zVar SAMPS_PER_CH values with values, rounded to valid values due to setting PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_ACTION.')
                            SAMPS_PER_CH_zv = SAMPS_PER_CH_rounded;
                        otherwise
                            error('BICAS:proc_sub:Assertion:ConfigurationBug', ...
                                'Illegal value PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_ACTION=%s', actionSettingValue)
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
            EJ_library.utils.assert.struct2(PreDc, {'Zv', 'hasSnapshotFormat', 'nRecords', 'nCdfSamplesPerRecord', 'isLfr', 'isTdsCwf'}, {});
            EJ_library.utils.assert.struct2(PreDc.Zv, {...
                'Epoch', 'ACQUISITION_TIME', 'samplesCaTm', 'freqHz', 'nValidSamplesPerRecord', 'iLsf', 'DIFF_GAIN', 'MUX_SET', 'QUALITY_FLAG', ...
                'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG'}, {'CALIBRATION_TABLE_INDEX'});
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(PreDc.Zv);

            assert(isa(PreDc.Zv.freqHz, 'double'))
        end



        function assert_PostDC(PostDc)
            EJ_library.utils.assert.struct2(PostDc, {'Zv', 'hasSnapshotFormat', 'nRecords', 'nCdfSamplesPerRecord', 'isLfr', 'isTdsCwf'}, {});
            EJ_library.utils.assert.struct2(PostDc.Zv, {...
                'Epoch', 'ACQUISITION_TIME', 'samplesCaTm', 'freqHz', 'nValidSamplesPerRecord', 'iLsf', 'DIFF_GAIN', 'MUX_SET', 'QUALITY_FLAG', ...
                'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'DemuxerOutput', 'IBIAS1', 'IBIAS2', 'IBIAS3', 'DemuxerOutput'}, {'CALIBRATION_TABLE_INDEX'});
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(PostDc.Zv);
        end



        function [OutSciZv] = process_PostDC_to_LFR(SciPostDc, outputDsi, outputVersion)
        % Processing function. Convert PostDC to any one of several similar LFR dataset PDs.
        
            % ASSERTIONS
            bicas.proc_sub.assert_PostDC(SciPostDc)
            
            OutSciZv = [];
            
            nSamplesPerRecord = size(SciPostDc.Zv.DemuxerOutput.dcV1, 2);   % Samples per record.
            
            outputDvid = bicas.construct_DVID(outputDsi, outputVersion);
            
            OutSciZv.Epoch            = SciPostDc.Zv.Epoch;
            OutSciZv.ACQUISITION_TIME = SciPostDc.Zv.ACQUISITION_TIME;
            OutSciZv.QUALITY_BITMASK  = SciPostDc.Zv.QUALITY_BITMASK;
            OutSciZv.QUALITY_FLAG     = SciPostDc.Zv.QUALITY_FLAG;
            OutSciZv.DELTA_PLUS_MINUS = SciPostDc.Zv.DELTA_PLUS_MINUS;
            OutSciZv.IBIAS1           = SciPostDc.Zv.IBIAS1;
            OutSciZv.IBIAS2           = SciPostDc.Zv.IBIAS2;
            OutSciZv.IBIAS3           = SciPostDc.Zv.IBIAS3;
            OutSciZv.SYNCHRO_FLAG     = SciPostDc.Zv.SYNCHRO_FLAG;
                    
            % NOTE: The two cases are different in the indexes they use for OutSciZv.
            switch(outputDvid)
                case  {'V05_SOLO_L2_RPW-LFR-SURV-CWF-E' ...
                       'V05_SOLO_L2_RPW-LFR-SBM1-CWF-E' ...
                       'V05_SOLO_L2_RPW-LFR-SBM2-CWF-E'}

                    % ASSERTION
                    assert(nSamplesPerRecord == 1, 'BICAS:proc_sub:Assertion:IllegalArgument', 'Number of samples per CDF record is not 1, as expected. Bad input CDF?')
                    assert(size(OutSciZv.QUALITY_FLAG,    2) == 1)
                    assert(size(OutSciZv.QUALITY_BITMASK, 2) == 1)

                    assert(size(OutSciZv.IBIAS1, 2) == 1)
                    assert(size(OutSciZv.IBIAS2, 2) == 1)
                    assert(size(OutSciZv.IBIAS3, 2) == 1)
                    
                    OutSciZv.V(:,1)   = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.V(:,2)   = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.V(:,3)   = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.E(:,1)   = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.E(:,2)   = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.E(:,3)   = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,3) = SciPostDc.Zv.DemuxerOutput.acV23;
                    
                case  {'V05_SOLO_L2_RPW-LFR-SURV-SWF-E'}
                    
                    % ASSERTION
                    assert(nSamplesPerRecord == 2048, 'BICAS:proc_sub:Assertion:IllegalArgument', 'Number of samples per CDF record is not 2048, as expected. Bad Input CDF?')
                    
                    OutSciZv.V(:,:,1)   = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.V(:,:,2)   = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.V(:,:,3)   = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.E(:,:,1)   = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.E(:,:,2)   = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.E(:,:,3)   = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,:,3) = SciPostDc.Zv.DemuxerOutput.acV23;

                    % Only in LFR SWF (not CWF): F_SAMPLE, SAMP_DTIME
                    OutSciZv.F_SAMPLE   = SciPostDc.Zv.freqHz;
                    
                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', 'Function can not produce outputDvid=%s.', outputDvid)
            end


            
            % ASSERTION
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(OutSciZv);
            % NOTE: Not really necessary since the list of zVars will be checked against the master CDF?
            EJ_library.utils.assert.struct2(OutSciZv, {'IBIAS1', 'IBIAS2', 'IBIAS3', 'V', 'E', 'EAC', 'Epoch', ...
                'QUALITY_BITMASK', 'QUALITY_FLAG', 'DELTA_PLUS_MINUS', 'ACQUISITION_TIME', 'SYNCHRO_FLAG'}, {'F_SAMPLE'})
        end   % process_PostDC_to_LFR



        function OutSciZv = process_PostDC_to_TDS(SciPostDc, outputDsi, outputVersion)
            
            % ASSERTIONS
            bicas.proc_sub.assert_PostDC(SciPostDc)
            
            OutSciZv = [];
            
            outputDvid = bicas.construct_DVID(outputDsi, outputVersion);

            OutSciZv.Epoch            = SciPostDc.Zv.Epoch;
            OutSciZv.ACQUISITION_TIME = SciPostDc.Zv.ACQUISITION_TIME;
            OutSciZv.QUALITY_FLAG     = SciPostDc.Zv.QUALITY_FLAG;
            OutSciZv.QUALITY_BITMASK  = SciPostDc.Zv.QUALITY_BITMASK;
            OutSciZv.DELTA_PLUS_MINUS = SciPostDc.Zv.DELTA_PLUS_MINUS;
            OutSciZv.IBIAS1           = SciPostDc.Zv.IBIAS1;
            OutSciZv.IBIAS2           = SciPostDc.Zv.IBIAS2;
            OutSciZv.IBIAS3           = SciPostDc.Zv.IBIAS3;
            OutSciZv.SYNCHRO_FLAG     = SciPostDc.Zv.SYNCHRO_FLAG;

            % NOTE: The two cases are actually different in the indexes they use for OutSciZv.
            switch(outputDvid)
                
                case {'V05_SOLO_L2_RPW-TDS-LFM-CWF-E'}
                    OutSciZv.V(:,1)     = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.V(:,2)     = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.V(:,3)     = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.E(:,1)     = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.E(:,2)     = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.E(:,3)     = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,1)   = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,2)   = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,3)   = SciPostDc.Zv.DemuxerOutput.acV23;
                    
                case {'V05_SOLO_L2_RPW-TDS-LFM-RSWF-E'}
                    OutSciZv.V(:,:,1)   = SciPostDc.Zv.DemuxerOutput.dcV1;
                    OutSciZv.V(:,:,2)   = SciPostDc.Zv.DemuxerOutput.dcV2;
                    OutSciZv.V(:,:,3)   = SciPostDc.Zv.DemuxerOutput.dcV3;
                    OutSciZv.E(:,:,1)   = SciPostDc.Zv.DemuxerOutput.dcV12;
                    OutSciZv.E(:,:,2)   = SciPostDc.Zv.DemuxerOutput.dcV13;
                    OutSciZv.E(:,:,3)   = SciPostDc.Zv.DemuxerOutput.dcV23;
                    OutSciZv.EAC(:,:,1) = SciPostDc.Zv.DemuxerOutput.acV12;
                    OutSciZv.EAC(:,:,2) = SciPostDc.Zv.DemuxerOutput.acV13;
                    OutSciZv.EAC(:,:,3) = SciPostDc.Zv.DemuxerOutput.acV23;
                    
                    OutSciZv.F_SAMPLE = SciPostDc.Zv.freqHz;
                    
                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', 'Function can not produce outputDvid=%s.', outputDvid)
            end



            % ASSERTION
            bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(OutSciZv);
            % NOTE: Not really necessary since the list of zVars will be checked against the master CDF?
            EJ_library.utils.assert.struct2(OutSciZv, {'IBIAS1', 'IBIAS2', 'IBIAS3', 'V', 'E', 'EAC', 'Epoch', ...
                'QUALITY_BITMASK', 'QUALITY_FLAG', 'DELTA_PLUS_MINUS', 'ACQUISITION_TIME', 'SYNCHRO_FLAG'}, {'F_SAMPLE'})
        end
        
        

        % Processing function. Converts PreDC to PostDC, i.e. demux and calibrate data.
        % Function is in large part a wrapper around "simple_demultiplex".
        %
        % NOTE: Public function as opposed to the other demuxing/calibration functions.
        %
        function PostDc = process_demuxing_calibration(PreDc, Cal, SETTINGS)
        % PROPOSAL: Move the setting of IBIASx (bias current) somewhere else?
        %   PRO: Unrelated to demultiplexing.
        %   CON: Related to calibration.
        % PROPOSAL: Change name. Will not calibrate measured samples here, only currents, maybe.

            % ASSERTION
            bicas.proc_sub.assert_PreDC(PreDc);

            %=======
            % DEMUX
            %=======
            PostDc = PreDc;    % Copy all values, to later overwrite a subset of them.
            PostDc.Zv.DemuxerOutput = bicas.proc_sub.simple_demultiplex(PreDc, Cal, SETTINGS);

            %================================
            % Set (calibrated) bias currents
            %================================
            % BUG / TEMP: Set default values since the real bias current values are not available.
            PostDc.Zv.IBIAS1 = bicas.proc_utils.create_NaN_array([PostDc.nRecords, PostDc.nCdfSamplesPerRecord]);
            PostDc.Zv.IBIAS2 = bicas.proc_utils.create_NaN_array([PostDc.nRecords, PostDc.nCdfSamplesPerRecord]);
            PostDc.Zv.IBIAS3 = bicas.proc_utils.create_NaN_array([PostDc.nRecords, PostDc.nCdfSamplesPerRecord]);
            
            % ASSERTION
            bicas.proc_sub.assert_PostDC(PostDc)
        end
        
    end   % methods(Static, Access=public)
            
    %###################################################################################################################
    
    methods(Static, Access=private)
    %methods(Static, Access=public)
        
        % Wrapper around "simple_demultiplex_subsequence_OLD" to be able to handle multiple CDF records with changing
        % settings.
        %
        % NOTE: Can handle arrays of any size as long as the sizes are consistent.
        %
        function AsrSamplesAVolt = simple_demultiplex(PreDc, Cal, SETTINGS)
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
            EJ_library.utils.assert.vector(PreDc.Zv.samplesCaTm)
            assert(numel(PreDc.Zv.samplesCaTm) == 5)
            bicas.proc_utils.assert_cell_array_comps_have_same_N_rows(PreDc.Zv.samplesCaTm)
            EJ_library.utils.assert.all_equal([...
                size(PreDc.Zv.MUX_SET,        1), ...
                size(PreDc.Zv.DIFF_GAIN,      1), ...
                size(PreDc.Zv.samplesCaTm{1}, 1)])



            % Create empty structure to which new array components can be added.
            % NOTE: Unit is AVolt. Not including in the field names to keep them short.
            AsrSamplesAVolt = struct(...
                'dcV1',  [], 'dcV2',  [], 'dcV3',  [], ...
                'dcV12', [], 'dcV23', [], 'dcV13', [], ...
                'acV12', [], 'acV23', [], 'acV13', []);

            dlrUsing12zv = bicas.demultiplexer_latching_relay(PreDc.Zv.Epoch);
            iCalibLZv    = Cal.get_calibration_time_L(        PreDc.Zv.Epoch);
            iCalibHZv    = Cal.get_calibration_time_H(        PreDc.Zv.Epoch);

            
            
            if isfield(PreDc.Zv, 'CALIBRATION_TABLE_INDEX')
                % NOTE: Technically, this should only happen for L1R input.
                
                % ASSERTION
                % NOTE: Checks both LFR & TDS CDF files.
                % NOTE: Has observed breaking assertion in LFR test files "LFR___LFR_suggested_2019-01-17".
                
                CALIBRATION_TABLE_INDEX_zv = PreDc.Zv.CALIBRATION_TABLE_INDEX;
                legalCtiSize = size(CALIBRATION_TABLE_INDEX_zv, 2) == 2;
                if ~legalCtiSize
                    if SETTINGS.get_fv('PROCESSING.L1R.ZV_CALIBRATION_TABLE_INDEX_ILLEGAL_SIZE_REPLACE')
                        bicas.log('warning', 'Setting CALIBRATION_TABLE_INDEX to NaN due to setting PROCESSING.L1R.ZV_CALIBRATION_TABLE_INDEX_ILLEGAL_SIZE_REPLACE.')
                        CALIBRATION_TABLE_INDEX_zv = zeros(PreDc.nRecords, 2) * NaN;
                    else
                        error('BICAS:proc_sub:Assertion', 'zVar CALIBRATION_TABLE_INDEX has illegal width (<>2).')
                    end
                end
            else
                % NOTE: Technically, this should only happen for L1 input.
                
                % Create "empty" CALIBRATION_TABLE_INDEX_zv.
                bicas.log('warning', 'Creating NaN-valued CALIBRATION_TABLE_INDEX due to zVar not being present in input CDF.')
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
                
                bicas.logf('info', ['Records %5i-%5i : ', ...
                    'MUX_SET=%i; DIFF_GAIN=%i; dlrUsing12=%i; freqHz=%5g; iCalibL=%i; iCalibH=%i; CALIBRATION_TABLE_INDEX=[%i, %i]'], ...
                    iFirst, iLast, ...
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
                        ssSamplesCaAVolt = Cal.calibrate_voltage_all(ssDtSec, ssSamplesCaTm, ...
                            PreDc.isLfr, PreDc.isTdsCwf, iBlts, ...
                            BltsSrcAsrArray(iBlts), biasHighGain, iCalibL_ss, iCalibH_ss, iLsf_ss, CALIBRATION_TABLE_INDEX_ss);
                        
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



    end   % methods(Static, Access=private)
        
end
