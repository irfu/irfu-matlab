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
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like" data, using the first MATLAB array index to
%   represent CDF records.
%
%
% SOME INTERMEDIATE PROCESSING DATA FORMATS
% =========================================
% - PreDC = Pre-Demuxing-Calibration Data
%       Generic data format that can represent all forms of input datasets before demuxing and calibration. Can use an
%       arbitrary number of samples per record. Some variables are therefore not used in CWF output datasets.
%       Consists of struct with fields:
%           .Epoch
%           .ACQUISITION_TIME
%           .DemuxerInput : struct with fields.
%               BIAS_1 to .BIAS_5  : NxM arrays, where M may be 1 (1 sample/record) or >1.
%           .freqHz                : Snapshot frequency in Hz. Unimportant for one sample/record data.
%           .DIFF_GAIN
%           .MUX_SET
%           QUALITY_FLAG
%           QUALITY_BITMASK
%           DELTA_PLUS_MINUS
%           % SAMP_DTIME          % Only important for SWF. - Abolished?
%       Fields are "CDF-like": rows=records, all have same number of rows.
% - PostDC = Post-Demuxing-Calibration Data
%       Like PreDC but with additional fields. Tries to capture a superset of the information that goes into any
%       dataset produced by BICAS.
%       Has extra fields:
%           .DemuxerOutput   : struct with fields.
%               V1, V2, V3,   V12, V13, V23,   V12_AC, V13_AC, V23_AC.
%           .IBIAS1
%           .IBIAS2
%           .IBIAS3
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
% PROPOSAL: Every processing function should use a special function for asserting and retrieving the right set of
%           InputsMap keys and values.
%   NOTE: Current convention/scheme only checks the existence of required keys, not absence of non-required keys.
%   PRO: More assertions.
%   PRO: Clearer dependencies.
%
% PROPOSAL: Assertions after every switch statement that differentiates different processing data/dataset versions.
%           Describe what they should all "converge" on, and make sure they actually do.
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
%#######################################################################################################################
    
    methods(Static, Access=public)
        
        function HkSciTime = process_HK_to_HK_on_SCI_TIME(Sci, Hk)
        % Processing function
        
            global SETTINGS
            
            % ASSERTIONS
            EJ_library.utils.assert.struct2(Sci, {'ZVars', 'Ga'}, {})
            EJ_library.utils.assert.struct2(Hk,  {'ZVars', 'Ga'}, {})
            
            HkSciTime = [];
            
            
            
            % Define local convenience variables. AT = ACQUISITION_TIME
            ACQUISITION_TIME_EPOCH_UTC = SETTINGS.get_fv('PROCESSING.ACQUISITION_TIME_EPOCH_UTC');
            
            hkAtTt2000  = bicas.proc_utils.ACQUISITION_TIME_to_tt2000(  Hk.ZVars.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            sciAtTt2000 = bicas.proc_utils.ACQUISITION_TIME_to_tt2000( Sci.ZVars.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
            hkEpoch     = Hk.ZVars.Epoch;
            sciEpoch    = Sci.ZVars.Epoch;
            
            %==================================================================
            % Log time intervals to enable comparing available SCI and HK data
            %==================================================================
            bicas.proc_utils.log_tt2000_array('HK  ACQUISITION_TIME', hkAtTt2000)
            bicas.proc_utils.log_tt2000_array('SCI ACQUISITION_TIME', sciAtTt2000)
            bicas.proc_utils.log_tt2000_array('HK  Epoch           ', hkEpoch)
            bicas.proc_utils.log_tt2000_array('SCI Epoch           ', sciEpoch)
            
            %=========================================================================================================
            % 1) Convert time to something linear in time that can be used for processing (not storing time to file).
            % 2) Effectively also chooses which time to use for the purpose of processing:
            %       (a) ACQUISITION_TIME, or
            %       (b) Epoch.
            %=========================================================================================================
            if SETTINGS.get_fv('PROCESSING.USE_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION')
                bicas.log('info', 'Using HK & SCI zVariable ACQUISITION_TIME (not Epoch) for interpolating HK dataset data to SCI dataset time.')
                hkInterpolationTimeTt2000  = hkAtTt2000;
                sciInterpolationTimeTt2000 = sciAtTt2000;
            else
                bicas.log('info', 'Using HK & SCI zVariable Epoch (not ACQUISITION_TIME) for interpolating HK dataset data to SCI dataset time.')
                hkInterpolationTimeTt2000  = hkEpoch;
                sciInterpolationTimeTt2000 = sciEpoch;
            end
            clear hkAtTt2000 sciAtTt2000
            clear hkEpoch    sciEpoch



            %=========================================================================================================
            % Derive MUX_SET
            % --------------
            % NOTE: Only obtains one MUX_SET per record ==> Can not change MUX_SET in the middle of a record.
            % NOTE: Can potentially obtain MUX_SET from LFR SCI.
            %=========================================================================================================            
            HkSciTime.MUX_SET = bicas.proc_utils.nearest_interpolate_float_records(...
                double(Hk.ZVars.HK_BIA_MODE_MUX_SET), ...
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
                double(Hk.ZVars.HK_BIA_DIFF_GAIN), hkInterpolationTimeTt2000, sciInterpolationTimeTt2000);



            % ASSERTIONS
            EJ_library.utils.assert.struct2(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
        end        
        
        

        function PreDc = process_LFR_to_PreDC(Sci, inputSciDsi, HkSciTime)
        % Processing function. Convert LFR CDF data to PreDC.
        %
        % Keeps number of samples/record. Treats 1 samples/record "length-one snapshots".
        
        % PROBLEM: Hardcoded CDF data types (MATLAB classes).
        % MINOR PROBLEM: Still does not handle LFR zVar TYPE for determining "virtual snapshot" length.
        % Should only be relevant for V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF (not V02) which should expire.
        
            LFR_SWF_SNAPSHOT_LENGTH = 2048;
        
            % ASSERTIONS
            EJ_library.utils.assert.struct2(Sci,       {'ZVars', 'Ga'}, {})
            EJ_library.utils.assert.struct2(HkSciTime, {'MUX_SET', 'DIFF_GAIN'}, {})
            
            nRecords = size(Sci.ZVars.Epoch, 1);            
            C = bicas.proc_utils.classify_DATASET_ID(inputSciDsi);
           
            V = Sci.ZVars.V;
            E = permute(Sci.ZVars.E, [1,3,2]);
            % Switch last two indices of E.
            % ==> index 2 = "snapshot" sample index, including for CWF (sample/record, "snapshots" consisting of 1 sample).
            %     index 3 = E1/E2 component.
            
            nSamplesPerRecord = size(V, 2);
            
            % ASSERTIONS
            if C.isLfrSwf
                assert(nSamplesPerRecord == LFR_SWF_SNAPSHOT_LENGTH)
            else
                assert(nSamplesPerRecord == 1)
            end
            assert(size(E, 3) == 2)
            

            
            if     C.isLfrSbm1
                FREQ = ones(nRecords, 1) * 1;   % Always value "1" (F1).
            elseif C.isLfrSbm2
                FREQ = ones(nRecords, 1) * 2;   % Always value "2" (F2).
            else
                FREQ = Sci.ZVars.FREQ;
            end
            assert(size(FREQ, 2) == 1)
            
            
            
            freqHz = bicas.proc_utils.get_LFR_frequency( FREQ );   % NOTE: Needed also for 1 SPR.
            
            % Obtain the relevant values (one per record) from zVariables R0, R1, R2, "R3".
            Rx = bicas.proc_utils.get_LFR_Rx( ...
                Sci.ZVars.R0, ...
                Sci.ZVars.R1, ...
                Sci.ZVars.R2, ...
                FREQ );   % NOTE: Function also handles the imaginary zVar "R3".

            PreDc = [];
            PreDc.Epoch                  = Sci.ZVars.Epoch;
            PreDc.ACQUISITION_TIME       = Sci.ZVars.ACQUISITION_TIME;
            PreDc.DELTA_PLUS_MINUS       = bicas.proc_utils.derive_DELTA_PLUS_MINUS(freqHz, nSamplesPerRecord);            
            PreDc.freqHz                 = freqHz;
            PreDc.nValidSamplesPerRecord = ones(nRecords, nSamplesPerRecord);
            PreDc.SYNCHRO_FLAG           = Sci.ZVars.TIME_SYNCHRO_FLAG;   % NOTE: Different zVar name in input and output datasets.

            
            
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
            PreDc.QUALITY_FLAG    = Sci.ZVars.QUALITY_FLAG;
            PreDc.QUALITY_BITMASK = Sci.ZVars.QUALITY_BITMASK;
            if isempty(PreDc.QUALITY_FLAG)
                bicas.log('warning', 'QUALITY_FLAG from the LFR SCI source dataset is empty. Filling with empty values.')
                PreDc.QUALITY_FLAG = bicas.proc_utils.create_NaN_array([nRecords, 1]);
            end
            if isempty(PreDc.QUALITY_BITMASK)
                bicas.log('warning', 'QUALITY_BITMASK from the LFR SCI source dataset is empty. Filling with empty values.')
                PreDc.QUALITY_BITMASK = bicas.proc_utils.create_NaN_array([nRecords, 1]);
            end
            % LFR QUALITY_FLAG, QUALITY_BITMASK not set yet (2019-09-17), but I presume they should have just one value
            % per record. BIAS output datasets should.
            assert(size(PreDc.QUALITY_FLAG,    2) == 1)
            assert(size(PreDc.QUALITY_BITMASK, 2) == 1)



            % E must be floating-point so that values can be set to NaN.
            % bicas.proc_utils.filter_rows requires this. Variable may be integer if integer in source CDF.
            E = single(E);

            PreDc.DemuxerInput        = [];
            PreDc.DemuxerInput.BIAS_1 = V;
            PreDc.DemuxerInput.BIAS_2 = bicas.proc_utils.filter_rows( E(:,:,1), Rx==1 );
            PreDc.DemuxerInput.BIAS_3 = bicas.proc_utils.filter_rows( E(:,:,2), Rx==1 );
            PreDc.DemuxerInput.BIAS_4 = bicas.proc_utils.filter_rows( E(:,:,1), Rx==0 );
            PreDc.DemuxerInput.BIAS_5 = bicas.proc_utils.filter_rows( E(:,:,2), Rx==0 );

            PreDc.MUX_SET   = HkSciTime.MUX_SET;
            PreDc.DIFF_GAIN = HkSciTime.DIFF_GAIN;



            % ASSERTIONS
            bicas.proc_sub.assert_PreDC(PreDc)
        end
        
        
        
        function PreDc = process_TDS_to_PreDC(Sci, inputSciDsi, HkSciTime)
        % Processing function. Convert TDS CDF data (PDs) to PreDC.
        %
        % Keeps number of samples/record. Treats 1 samples/record "length-one snapshots".
        %
        % BUG?: Does not use CHANNEL_STATUS_INFO.
        % NOTE: BIAS output datasets do not have a variable for the length of snapshots. Need to use NaN/fill value.

            % ASSERTIONS
            EJ_library.utils.assert.struct2(Sci,        {'ZVars', 'Ga'}, {})
            EJ_library.utils.assert.struct2(HkSciTime,  {'MUX_SET', 'DIFF_GAIN'}, {})
            
            C = bicas.proc_utils.classify_DATASET_ID(inputSciDsi);
            
            nRecords                  = size(Sci.ZVars.Epoch, 1);
            nVariableSamplesPerRecord = size(Sci.ZVars.WAVEFORM_DATA, 3);   % Number of samples in the variable, not necessarily actual data.
            
            freqHz = double(Sci.ZVars.SAMPLING_RATE);
            
            PreDc = [];
            
            PreDc.Epoch            = Sci.ZVars.Epoch;
            PreDc.ACQUISITION_TIME = Sci.ZVars.ACQUISITION_TIME;
            PreDc.DELTA_PLUS_MINUS = bicas.proc_utils.derive_DELTA_PLUS_MINUS(freqHz, nVariableSamplesPerRecord);
            PreDc.freqHz           = freqHz;
            if C.isTdsRswf
                
                %====================================================================================================
                % ASSERTION WARNING: Check zVar SAMPS_PER_CH for invalid values
                %
                % NOTE: Has observed invalid SAMPS_PER_CH value 16562 in
                % ROC-SGSE_L1R_RPW-TDS-LFM-RSWF-E_73525cd_CNE_V03.CDF.
                % 2019-09-18, David Pisa: Not a flaw in TDS RCS but in the source L1 dataset.
                %====================================================================================================
                SAMPS_PER_CH_MIN_VALID = 2^10;
                SAMPS_PER_CH_MAX_VALID = 2^15;
                SAMPS_PER_CH         = double(Sci.ZVars.SAMPS_PER_CH);
                SAMPS_PER_CH_rounded = round(2.^round(log2(SAMPS_PER_CH)));
                SAMPS_PER_CH_rounded(SAMPS_PER_CH_rounded < SAMPS_PER_CH_MIN_VALID) = SAMPS_PER_CH_MIN_VALID;
                SAMPS_PER_CH_rounded(SAMPS_PER_CH_rounded > SAMPS_PER_CH_MAX_VALID) = SAMPS_PER_CH_MAX_VALID;
                if any(SAMPS_PER_CH_rounded ~= SAMPS_PER_CH)
                    SAMPS_PER_CH_badValues = unique(SAMPS_PER_CH(SAMPS_PER_CH_rounded ~= SAMPS_PER_CH));
                    badValuesDisplayStr = strjoin(arrayfun(@(n) sprintf('%i', n), SAMPS_PER_CH_badValues, 'uni', false), ', ');                    
                    bicas.logf('warning', 'TDS LFM RSWF zVar SAMPS_PER_CH contains unexpected value(s), not 2^n: %s', badValuesDisplayStr)
                    
                    % NOTE: Unclear if this is the appropriate action.
                    %bicas.log('warning', 'Replacing TDS RSWF zVar SAMPS_PER_CH values with values, rounded to valid values.')
                    %SAMPS_PER_CH = SAMPS_PER_CH_rounded;
                end
                
                % NOTE: This might only be appropriate for TDS's "COMMON_MODE" mode. TDS also has a "FULL_BAND" mode
                % with 2^18=262144 samples per snapshot. You should never encounter FULL_BAND in any dataset (even on
                % ground), only used for calibration and testing. /David Pisa & Jan Soucek in emails, 2016.
                % --
                % FULL_BAND mode has each snapshot divided into 2^15 samples/record * 8 records.  /Unknown source
                % Unclear what value SAMPS_PER_CH should have for FULL_BAND mode. How does Epoch work for FULL_BAND
                % snapshots?
                PreDc.nValidSamplesPerRecord = SAMPS_PER_CH;                
                
            else
                PreDc.nValidSamplesPerRecord = ones(nRecords, 1) * 1;
            end

            PreDc.QUALITY_FLAG    = Sci.ZVars.QUALITY_FLAG;
            PreDc.QUALITY_BITMASK = Sci.ZVars.QUALITY_BITMASK;
            PreDc.SYNCHRO_FLAG    = Sci.ZVars.TIME_SYNCHRO_FLAG;   % NOTE: Different zVar name in input and output datasets.
            
            modif_WAVEFORM_DATA = double(permute(Sci.ZVars.WAVEFORM_DATA, [1,3,2]));
            
            PreDc.DemuxerInput        = [];
            PreDc.DemuxerInput.BIAS_1 = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,1), PreDc.nValidSamplesPerRecord );
            PreDc.DemuxerInput.BIAS_2 = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,2), PreDc.nValidSamplesPerRecord );
            PreDc.DemuxerInput.BIAS_3 = bicas.proc_utils.set_NaN_after_snapshots_end( modif_WAVEFORM_DATA(:,:,3), PreDc.nValidSamplesPerRecord );
            PreDc.DemuxerInput.BIAS_4 = bicas.proc_utils.create_NaN_array([nRecords, nVariableSamplesPerRecord]);
            PreDc.DemuxerInput.BIAS_5 = bicas.proc_utils.create_NaN_array([nRecords, nVariableSamplesPerRecord]);
            
            PreDc.MUX_SET   = HkSciTime.MUX_SET;
            PreDc.DIFF_GAIN = HkSciTime.DIFF_GAIN;
            
            
            
            % ASSERTIONS
            bicas.proc_sub.assert_PreDC(PreDc)
        end



        function assert_PreDC(PreDc)
            EJ_library.utils.assert.struct2(PreDc, {...
                'Epoch', 'ACQUISITION_TIME', 'DemuxerInput', 'freqHz', 'nValidSamplesPerRecord', 'DIFF_GAIN', 'MUX_SET', 'QUALITY_FLAG', ...
                'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG'}, {});
            bicas.proc_utils.assert_unvaried_N_rows(PreDc);
            bicas.proc_utils.assert_unvaried_N_rows(PreDc.DemuxerInput);
            
            assert(isa(PreDc.freqHz, 'double'))
        end
        
        
        
        function assert_PostDC(PostDc)
            EJ_library.utils.assert.struct2(PostDc, {...
                'Epoch', 'ACQUISITION_TIME', 'DemuxerInput', 'freqHz', 'nValidSamplesPerRecord', 'DIFF_GAIN', 'MUX_SET', 'QUALITY_FLAG', ...
                'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'DemuxerOutput', 'IBIAS1', 'IBIAS2', 'IBIAS3'}, {});
            bicas.proc_utils.assert_unvaried_N_rows(PostDc);
            bicas.proc_utils.assert_unvaried_N_rows(PostDc.DemuxerOutput);
        end
        

        
        function OutSciZVars = process_PostDC_to_LFR(SciPostDc, outputDsi, outputVersion)
        % Processing function. Convert PostDC to any one of several similar LFR dataset PDs.
        
            % ASSERTIONS
            bicas.proc_sub.assert_PostDC(SciPostDc)
            
            OutSciZVars = [];
            
            nSamplesPerRecord = size(SciPostDc.DemuxerOutput.V1, 2);   % Samples per record.
            
            outputDvid = bicas.construct_DVID(outputDsi, outputVersion);
            ZVAR_FN_LIST = {'IBIAS1', 'IBIAS2', 'IBIAS3', 'V', 'E', 'EAC', 'Epoch', ...
                'QUALITY_BITMASK', 'QUALITY_FLAG', 'DELTA_PLUS_MINUS', 'ACQUISITION_TIME'};
            
            OutSciZVars.Epoch = SciPostDc.Epoch;
            OutSciZVars.ACQUISITION_TIME = SciPostDc.ACQUISITION_TIME;
            OutSciZVars.QUALITY_BITMASK  = SciPostDc.QUALITY_BITMASK;
            OutSciZVars.QUALITY_FLAG     = SciPostDc.QUALITY_FLAG;
            OutSciZVars.DELTA_PLUS_MINUS = SciPostDc.DELTA_PLUS_MINUS;
            
            switch(outputDvid)
                case  {'V04_ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E' ...
                       'V04_ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E' ...
                       'V04_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E' ...
                            'V04_SOLO_L2_RPW-LFR-SURV-CWF-E' ...
                            'V04_SOLO_L2_RPW-LFR-SBM1-CWF-E' ...
                            'V04_SOLO_L2_RPW-LFR-SBM2-CWF-E' ...
                       }
                   
                    % ASSERTION
                    assert(nSamplesPerRecord == 1, 'BICAS:proc_sub:Assertion:IllegalArgument', 'Number of samples per CDF record is not 1, as expected. Bad input CDF?')
                    assert(size(OutSciZVars.QUALITY_FLAG,    2) == 1)
                    assert(size(OutSciZVars.QUALITY_BITMASK, 2) == 1)

                    OutSciZVars.IBIAS1 = SciPostDc.IBIAS1;
                    OutSciZVars.IBIAS2 = SciPostDc.IBIAS2;
                    OutSciZVars.IBIAS3 = SciPostDc.IBIAS3;
                    assert(size(OutSciZVars.IBIAS1, 2) == 1)
                    assert(size(OutSciZVars.IBIAS2, 2) == 1)
                    assert(size(OutSciZVars.IBIAS3, 2) == 1)
                    
                    OutSciZVars.V(:,1)           = SciPostDc.DemuxerOutput.V1;
                    OutSciZVars.V(:,2)           = SciPostDc.DemuxerOutput.V2;
                    OutSciZVars.V(:,3)           = SciPostDc.DemuxerOutput.V3;
                    OutSciZVars.E(:,1)           = SciPostDc.DemuxerOutput.V12;
                    OutSciZVars.E(:,2)           = SciPostDc.DemuxerOutput.V13;
                    OutSciZVars.E(:,3)           = SciPostDc.DemuxerOutput.V23;
                    OutSciZVars.EAC(:,1)         = SciPostDc.DemuxerOutput.V12_AC;
                    OutSciZVars.EAC(:,2)         = SciPostDc.DemuxerOutput.V13_AC;
                    OutSciZVars.EAC(:,3)         = SciPostDc.DemuxerOutput.V23_AC;
                    
                case  {'V04_ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E',  ...
                            'V04_SOLO_L2_RPW-LFR-SURV-SWF-E'}
                    
                    % ASSERTION
                    assert(nSamplesPerRecord == 2048, 'BICAS:proc_sub:Assertion:IllegalArgument', 'Number of samples per CDF record is not 2048, as expected. Bad Input CDF?')
                    
                    OutSciZVars.IBIAS1           = SciPostDc.IBIAS1;
                    OutSciZVars.IBIAS2           = SciPostDc.IBIAS2;
                    OutSciZVars.IBIAS3           = SciPostDc.IBIAS3;
                    OutSciZVars.V(:,:,1)         = SciPostDc.DemuxerOutput.V1;
                    OutSciZVars.V(:,:,2)         = SciPostDc.DemuxerOutput.V2;
                    OutSciZVars.V(:,:,3)         = SciPostDc.DemuxerOutput.V3;
                    OutSciZVars.E(:,:,1)         = SciPostDc.DemuxerOutput.V12;
                    OutSciZVars.E(:,:,2)         = SciPostDc.DemuxerOutput.V13;
                    OutSciZVars.E(:,:,3)         = SciPostDc.DemuxerOutput.V23;
                    OutSciZVars.EAC(:,:,1)       = SciPostDc.DemuxerOutput.V12_AC;
                    OutSciZVars.EAC(:,:,2)       = SciPostDc.DemuxerOutput.V13_AC;
                    OutSciZVars.EAC(:,:,3)       = SciPostDc.DemuxerOutput.V23_AC;

                    % Only in LFR SWF (not CWF): F_SAMPLE, SAMP_DTIME
                    OutSciZVars.F_SAMPLE         = SciPostDc.freqHz;
                    ZVAR_FN_LIST{end+1} = 'F_SAMPLE';
                    
                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', 'Function can not produce outputDvid=%s.', outputDvid)
            end



            % PROPOSAL: Use classify_DATASET_ID flag.
            switch(outputDvid)
                case  {'V04_ROC-SGSE_L2S_RPW-LFR-SBM1-CWF-E' ...
                       'V04_ROC-SGSE_L2S_RPW-LFR-SBM2-CWF-E' ...
                       'V04_ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E' ...
                       'V04_ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E'}
                    OutSciZVars.SYNCHRO_FLAG = SciPostDc.SYNCHRO_FLAG;
                    ZVAR_FN_LIST{end+1} = 'SYNCHRO_FLAG';

                case  {'V04_SOLO_L2_RPW-LFR-SURV-CWF-E' ...
                       'V04_SOLO_L2_RPW-LFR-SURV-SWF-E' ...
                       'V04_SOLO_L2_RPW-LFR-SBM1-CWF-E' ...
                       'V04_SOLO_L2_RPW-LFR-SBM2-CWF-E'}
                    OutSciZVars.SYNCHRO_FLAG = SciPostDc.SYNCHRO_FLAG;
                    ZVAR_FN_LIST{end+1} = 'SYNCHRO_FLAG';
                    
                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', 'Function can not produce outputDvid=%s.', outputDvid)
            end
            
            
            
            % ASSERTION
            bicas.proc_utils.assert_unvaried_N_rows(OutSciZVars);
            EJ_library.utils.assert.struct2(OutSciZVars, ZVAR_FN_LIST, {})
        end   % process_PostDC_to_LFR



        function OutSciZVars = process_PostDC_to_TDS(SciPostDc, outputDsi, outputVersion)
            
            % ASSERTIONS
            bicas.proc_sub.assert_PostDC(SciPostDc)
            
            OutSciZVars = [];
            
            outputDvid = bicas.construct_DVID(outputDsi, outputVersion);
            ZVAR_FN_LIST = {'IBIAS1', 'IBIAS2', 'IBIAS3', 'V', 'E', 'EAC', 'Epoch', ...
                'QUALITY_BITMASK', 'QUALITY_FLAG', 'DELTA_PLUS_MINUS', 'ACQUISITION_TIME'};
            
            switch(outputDvid)
                case {'V04_ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E' ...
                           'V04_SOLO_L2_RPW-TDS-LFM-CWF-E'}
                    OutSciZVars.V(:,1)     = SciPostDc.DemuxerOutput.V1;
                    OutSciZVars.V(:,2)     = SciPostDc.DemuxerOutput.V2;
                    OutSciZVars.V(:,3)     = SciPostDc.DemuxerOutput.V3;
                    OutSciZVars.E(:,1)     = SciPostDc.DemuxerOutput.V12;
                    OutSciZVars.E(:,2)     = SciPostDc.DemuxerOutput.V13;
                    OutSciZVars.E(:,3)     = SciPostDc.DemuxerOutput.V23;
                    OutSciZVars.EAC(:,1)   = SciPostDc.DemuxerOutput.V12_AC;
                    OutSciZVars.EAC(:,2)   = SciPostDc.DemuxerOutput.V13_AC;
                    OutSciZVars.EAC(:,3)   = SciPostDc.DemuxerOutput.V23_AC;
                    
                case {'V04_ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E' ...
                           'V04_SOLO_L2_RPW-TDS-LFM-RSWF-E'}
                       
                    OutSciZVars.V(:,:,1)   = SciPostDc.DemuxerOutput.V1;
                    OutSciZVars.V(:,:,2)   = SciPostDc.DemuxerOutput.V2;
                    OutSciZVars.V(:,:,3)   = SciPostDc.DemuxerOutput.V3;
                    OutSciZVars.E(:,:,1)   = SciPostDc.DemuxerOutput.V12;
                    OutSciZVars.E(:,:,2)   = SciPostDc.DemuxerOutput.V13;
                    OutSciZVars.E(:,:,3)   = SciPostDc.DemuxerOutput.V23;
                    OutSciZVars.EAC(:,:,1) = SciPostDc.DemuxerOutput.V12_AC;
                    OutSciZVars.EAC(:,:,2) = SciPostDc.DemuxerOutput.V13_AC;
                    OutSciZVars.EAC(:,:,3) = SciPostDc.DemuxerOutput.V23_AC;
                    
                    OutSciZVars.F_SAMPLE       = SciPostDc.freqHz;
                    ZVAR_FN_LIST{end+1} = 'F_SAMPLE';

                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', 'Function can not produce outputDvid=%s.', outputDvid)
            end

            OutSciZVars.Epoch = SciPostDc.Epoch;
            OutSciZVars.ACQUISITION_TIME = SciPostDc.ACQUISITION_TIME;
            OutSciZVars.QUALITY_FLAG     = SciPostDc.QUALITY_FLAG;
            OutSciZVars.QUALITY_BITMASK  = SciPostDc.QUALITY_BITMASK;
            OutSciZVars.DELTA_PLUS_MINUS = SciPostDc.DELTA_PLUS_MINUS;
            OutSciZVars.IBIAS1           = SciPostDc.IBIAS1;
            OutSciZVars.IBIAS2           = SciPostDc.IBIAS2;
            OutSciZVars.IBIAS3           = SciPostDc.IBIAS3;
            
            switch(outputDvid)
                case  {'V04_ROC-SGSE_L2S_RPW-TDS-LFM-CWF-E' ...
                       'V04_ROC-SGSE_L2S_RPW-TDS-LFM-RSWF-E'}
                    OutSciZVars.SYNCHRO_FLAG = SciPostDc.SYNCHRO_FLAG;
                    ZVAR_FN_LIST{end+1} = 'SYNCHRO_FLAG';

                case  {'V04_SOLO_L2_RPW-TDS-LFM-CWF-E' ...
                       'V04_SOLO_L2_RPW-TDS-LFM-RSWF-E'}
                    OutSciZVars.SYNCHRO_FLAG = SciPostDc.SYNCHRO_FLAG;
                    ZVAR_FN_LIST{end+1} = 'SYNCHRO_FLAG';
                    
                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument', 'Function can not produce outputDvid=%s.', outputDvid)
            end

            % ASSERTION
            bicas.proc_utils.assert_unvaried_N_rows(OutSciZVars);
            EJ_library.utils.assert.struct2(OutSciZVars, ZVAR_FN_LIST, {})
        end



        % Processing function. Converts PreDC to PostDC, i.e. demux and calibrate data.
        %
        % Is in large part a wrapper around "simple_demultiplex".
        % NOTE: Public function as opposed to the other demuxing/calibration functions.
        function PostDc = process_demuxing_calibration(PreDc)
        % PROPOSAL: Move setting of IBIASx (bias current) somewhere else?
        %   PRO: Unrelated to demultiplexing.
        %   CON: Related to calibration.

            % ASSERTION
            bicas.proc_sub.assert_PreDC(PreDc);

            %=======
            % DEMUX
            %=======
            PostDc = PreDc;    % Copy all values, to later overwrite a subset of them.
            PostDc.DemuxerOutput = bicas.proc_sub.simple_demultiplex(...
                PreDc.DemuxerInput, ...
                PreDc.MUX_SET, ...
                PreDc.DIFF_GAIN);

            %================================
            % Set (calibrated) bias currents
            %================================
            % BUG / TEMP: Set default values since the real bias current values are not available.
            PostDc.IBIAS1 = bicas.proc_utils.create_NaN_array(size(PostDc.DemuxerOutput.V1));
            PostDc.IBIAS2 = bicas.proc_utils.create_NaN_array(size(PostDc.DemuxerOutput.V2));
            PostDc.IBIAS3 = bicas.proc_utils.create_NaN_array(size(PostDc.DemuxerOutput.V3));
            
            % ASSERTION
            bicas.proc_sub.assert_PostDC(PostDc)
        end
        
    end   % methods(Static, Access=public)
            
    %###################################################################################################################
    
    methods(Static, Access=private)
    %methods(Static, Access=public)
        
        % Wrapper around "simple_demultiplex_subsequence_OLD" to be able to handle multiple CDF records with changing
        % settings (mux_set, diff_gain).
        %
        % NOTE: NOT a processing function (does not derive a PDV).
        %
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % DemuxerInput = Struct with fields BIAS_1 to BIAS_5.
        % MUX_SET      = Column vector. Numbers identifying the MUX/DEMUX mode. 
        % DIFF_GAIN    = Column vector. Gains for differential measurements. 0 = Low gain, 1 = High gain.
        %
        %
        % NOTE: Can handle arrays of any size as long as the sizes are consistent.
        function DemuxerOutput = simple_demultiplex(DemuxerInput, MUX_SET, DIFF_GAIN)
        % PROPOSAL: Incorporate into processing function process_demuxing_calibration.
        % PROPOSAL: Assert same nbr of "records" for MUX_SET, DIFF_GAIN as for BIAS_x.
        
            % ASSERTIONS
            EJ_library.utils.assert.struct2(DemuxerInput, {'BIAS_1', 'BIAS_2', 'BIAS_3', 'BIAS_4', 'BIAS_5'}, {})
            bicas.proc_utils.assert_unvaried_N_rows(DemuxerInput)
            EJ_library.utils.assert.all_equal([...
                size(MUX_SET,             1), ...
                size(DIFF_GAIN,           1), ...
                size(DemuxerInput.BIAS_1, 1)])



            % Create empty structure to which new components can be added.
            DemuxerOutput = struct(...
                'V1',     [], 'V2',     [], 'V3',     [], ...
                'V12',    [], 'V23',    [], 'V13',    [], ...
                'V12_AC', [], 'V23_AC', [], 'V13_AC', []);



            %====================================================================
            % Find continuous sequences of records with identical settings, then
            % process data separately (one iteration) for those sequences.
            %====================================================================
            [iFirstList, iLastList] = bicas.proc_utils.find_sequences(MUX_SET, DIFF_GAIN);            
            for iSequence = 1:length(iFirstList)
                
                iFirst = iFirstList(iSequence);
                iLast  = iLastList (iSequence);
                
                % Extract SCALAR settings to use for entire subsequence of records.
                MUX_SET_value   = MUX_SET  (iFirst);
                DIFF_GAIN_value = DIFF_GAIN(iFirst);
                bicas.logf('info', 'Records %2i-%2i : Demultiplexing; MUX_SET=%-3i; DIFF_GAIN=%-3i', ...
                    iFirst, iLast, MUX_SET_value, DIFF_GAIN_value)    % "%-3" since value might be NaN.

                % Extract subsequence of DATA records to "demux".
                DemuxerInputSubseq = bicas.proc_utils.select_row_range_from_struct_fields(DemuxerInput, iFirst, iLast);
                
                %=================================================
                % CALL DEMUXER - See method/function for comments
                %=================================================
                DemuxerOutputSubseq = bicas.proc_sub.simple_demultiplex_subsequence_OLD(...
                    DemuxerInputSubseq, MUX_SET_value, DIFF_GAIN_value);
                
                % Add demuxed sequence to the to-be complete set of records.
                DemuxerOutput = bicas.proc_utils.add_rows_to_struct_fields(DemuxerOutput, DemuxerOutputSubseq);
                
            end
            
        end   % simple_demultiplex



        % NOTE: FUNCTION IS PLANNED TO BE PHASED OUT/OBSOLETED.
        %
        %
        % Demultiplex, with only constant factors for calibration (no transfer functions, no offsets) and exactly one
        % setting for MUX_SET and DIFF_GAIN respectively.
        %
        % This function implements Table 3 and Table 4 in "RPW-SYS-MEB-BIA-SPC-00001-IRF", iss1rev16.
        % Variable names are chosen according to these tables.
        %
        % NOTE/BUG: Does not handle latching relay.
        %
        % NOTE: Conceptually, this function does both (a) demuxing and (b) calibration which could be separated.
        % - Demuxing is done on individual samples at a specific point in time.
        % - Calibration (with transfer functions) is made on a time series (presumably of one variable, but could be several).
        %
        % NOTE: NOT a processing function (does not derive a PDV).
        %
        % NOTE: Function is intended for development/testing until there is proper code for using transfer functions.
        % NOTE: "input"/"output" refers to input/output for the function, which is (approximately) the opposite of
        % the physical signals in the BIAS hardware.
        %
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % Input     : Struct with fields BIAS_1 to BIAS_5.
        % MUX_SET   : Scalar number identifying the MUX/DEMUX mode.
        % DIFF_GAIN : Scalar gain for differential measurements. 0 = Low gain, 1 = High gain.
        % Output    : Struct with fields V1, V2, V3,   V12, V13, V23,   V12_AC, V13_AC, V23_AC.
        % --
        % NOTE: Will tolerate values of NaN for MUX_SET, DIFF_GAIN. The effect is NaN in the corresponding output values.
        % NOTE: Can handle any arrays of any size as long as the sizes are consistent.
        %
        function Output = simple_demultiplex_subsequence_OLD(Input, MUX_SET, DIFF_GAIN)
        %==========================================================================================================
        % QUESTION: How to structure the demuxing?
        % --
        % QUESTION: How split by record? How put together again? How do in a way which
        %           works for real transfer functions? How handle the many non-indexed outputs?
        % QUESTION: How handle changing values of diff_gain, mux_set, bias-dependent calibration offsets?
        % NOTE: LFR data can be either 1 sample/record or 1 snapshot/record.
        % PROPOSAL: Work with some subset of in- and out-values of each type?
        %   PROPOSAL: Work with exactly one value of each type?
        %       CON: Slow.
        %           CON: Only temporary implementation.
        %       PRO: Quick to implement.
        %   PROPOSAL: Work with only some arbitrary subset specified by array of indices.
        %   PROPOSAL: Work with only one row?
        %   PROPOSAL: Work with a continuous sequence of rows/records?
        %   PROPOSAL: Submit all values, and return structure. Only read and set subset specified by indices.
        %
        %
        % PROPOSAL: Could, maybe, be used for demuxing if the caller has already applied the
        %           transfer function calibration on the BIAS signals.
        % PROPOSAL: Validate with some "multiplexer" function?!
        % QUESTION: Does it make sense to have BIAS values as cell array? Struct fields?!
        %   PRO: Needed for caller's for loop to split up by record.
        %
        % QUESTION: Is there some better implementation than giant switch statement?! Something more similar to BIAS
        % specification Table 3-4?
        %
        % QUESTION: MUX modes 1-3 are overdetermined if we always have BIAS1-3?
        %           If so, how select what to calculate?! What if results disagree/are inconsistent? Check for it?
        %
        % PROPOSAL: Separate the multiplication with factor in other function.
        %   PRO: Can use function together with TFs.
        %
        % TODO: Implement demuxing latching relay.
        %==========================================================================================================
            
            global SETTINGS
            
            % ASSERTIONS
            EJ_library.utils.assert.struct2(Input, {'BIAS_1', 'BIAS_2', 'BIAS_3', 'BIAS_4', 'BIAS_5'}, {})
            bicas.proc_utils.assert_unvaried_N_rows(Input)
            assert(isscalar(MUX_SET))
            assert(isscalar(DIFF_GAIN))



            ALPHA = SETTINGS.get_fv('PROCESSING.CALIBRATION.SCALAR.ALPHA');
            BETA  = SETTINGS.get_fv('PROCESSING.CALIBRATION.SCALAR.BETA');
            GAMMA = bicas.proc_utils.get_simple_demuxer_gamma(DIFF_GAIN);   % NOTE: GAMMA can be NaN iff DIFF_GAIN is.
            
            % Set default values which will be returned for
            % variables which are not set by the demuxer.
            NAN_VALUES = ones(size(Input.BIAS_1)) * NaN;
            V1_LF     = NAN_VALUES;
            V2_LF     = NAN_VALUES;
            V3_LF     = NAN_VALUES;
            V12_LF    = NAN_VALUES;
            V13_LF    = NAN_VALUES;
            V23_LF    = NAN_VALUES;
            V12_LF_AC = NAN_VALUES;
            V13_LF_AC = NAN_VALUES;
            V23_LF_AC = NAN_VALUES;
            
            % IMPLEMENTATION NOTE: Avoid getting integer - single ==> error.
            Input.BIAS_1 = double(Input.BIAS_1);
            Input.BIAS_2 = double(Input.BIAS_2);
            Input.BIAS_3 = double(Input.BIAS_3);
            Input.BIAS_4 = double(Input.BIAS_4);
            Input.BIAS_5 = double(Input.BIAS_5);

            switch(MUX_SET)
                case 0   % "Standard operation" : We have all information.

                    % Summarize the INPUT DATA we have.
                    V1_DC  = Input.BIAS_1;
                    V12_DC = Input.BIAS_2;
                    V23_DC = Input.BIAS_3;
                    V12_AC = Input.BIAS_4;
                    V23_AC = Input.BIAS_5;
                    % Derive the OUTPUT DATA which are trivial.
                    V1_LF     = V1_DC  / ALPHA;
                    V12_LF    = V12_DC / BETA;
                    V23_LF    = V23_DC / BETA;
                    V12_LF_AC = V12_AC / GAMMA;
                    V23_LF_AC = V23_AC / GAMMA;
                    % Derive the OUTPUT DATA which are less trivial.
                    V13_LF    = V12_LF    + V23_LF;
                    V2_LF     = V1_LF     - V12_LF;
                    V3_LF     = V2_LF     - V23_LF;
                    V13_LF_AC = V12_LF_AC + V23_LF_AC;
                    
                case 1   % Probe 1 fails
                    
                    V2_LF     = Input.BIAS_1 / ALPHA;
                    V3_LF     = Input.BIAS_2 / ALPHA;
                    V23_LF    = Input.BIAS_3 / BETA;
                    % Input.BIAS_4 unavailable.
                    V23_LF_AC = Input.BIAS_5 / GAMMA;
                    
                case 2   % Probe 2 fails
                    
                    V1_LF     = Input.BIAS_1 / ALPHA;
                    V3_LF     = Input.BIAS_2 / ALPHA;
                    V13_LF    = Input.BIAS_3 / BETA;
                    V13_LF_AC = Input.BIAS_4 / GAMMA;
                    % Input.BIAS_5 unavailable.
                    
                case 3   % Probe 3 fails
                    
                    V1_LF     = Input.BIAS_1 / ALPHA;
                    V2_LF     = Input.BIAS_2 / ALPHA;
                    V12_LF    = Input.BIAS_3 / BETA;
                    V12_LF_AC = Input.BIAS_4 / GAMMA;
                    % Input.BIAS_5 unavailable.
                    
                case 4   % Calibration mode 0
                    
                    % Summarize the INPUT DATA we have.
                    V1_DC  = Input.BIAS_1;
                    V2_DC  = Input.BIAS_2;
                    V3_DC  = Input.BIAS_3;
                    V12_AC = Input.BIAS_4;
                    V23_AC = Input.BIAS_5;
                    % Derive the OUTPUT DATA which are trivial.
                    V1_LF     = V1_DC / ALPHA;
                    V2_LF     = V2_DC / ALPHA;
                    V3_LF     = V3_DC / ALPHA;
                    V12_LF_AC = V12_AC / GAMMA;
                    V23_LF_AC = V23_AC / GAMMA;
                    % Derive the OUTPUT DATA which are less trivial.
                    V12_LF    = V1_LF     - V2_LF;
                    V13_LF    = V1_LF     - V3_LF;
                    V23_LF    = V2_LF     - V3_LF;
                    V13_LF_AC = V12_LF_AC + V23_LF_AC;

                case {5,6,7}   % Calibration mode 1/2/3
                    
                    % Summarize the INPUT DATA we have.
                    V12_AC = Input.BIAS_4;
                    V23_AC = Input.BIAS_5;
                    % Derive the OUTPUT DATA which are trivial.
                    V12_LF_AC = V12_AC / GAMMA;
                    V23_LF_AC = V23_AC / GAMMA;
                    % Derive the OUTPUT DATA which are less trivial.
                    V13_LF_AC = V12_LF_AC + V23_LF_AC;
                    
                otherwise
                    if isnan(MUX_SET)
                        ;   % Do nothing. Allow the default values (NaN) to be returned.
                    else
                        error('BICAS:proc_sub:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value for mux_set.')
                    end
            end   % switch
            
            % Create structure to return. (Removes the "_LF" suffix.)
            Output = [];
            Output.V1     = V1_LF;
            Output.V2     = V2_LF;
            Output.V3     = V3_LF;
            Output.V12    = V12_LF;
            Output.V13    = V13_LF;
            Output.V23    = V23_LF;
            Output.V12_AC = V12_LF_AC;
            Output.V13_AC = V13_LF_AC;
            Output.V23_AC = V23_LF_AC;
            
        end  % simple_demultiplex_subsequence_OLD

        
        
        % NEW FUNCTION. NOT USED YET BUT MEANT TO REPLACE OLD FUNCTION "simple_demultiplex_subsequence_OLD".
        %
        % (1) Return the information needed for how to calibrate a BIAS-LFR/TDS signal (BIAS_i) that is a function of the demultiplexer mode,
        % (2) Derives as much as possible of all the antenna singles and diffs from the available BIAS-LFR/TDS signals
        % (BIAS_i), except the calibration (i.e. only addition and subtraction).
        %
        % Meant to be called in two different ways, typically twice for any time period with samples.
        % (1) To obtain signal type info needed for how to calibrate every BIAS-LFR/TDS signal (BIAS_i) signal given any demux mode. 
        % (2) To derive the complete set of ASR samples from the given BLTS samples.
        %
        % RATIONALE: Meant to collect all hardcoded information about the demultiplexer routing of signals.
        % NOTE: Does not perform any calibration. The closest is to calculate diffs and singles from diffs and singles.
        % 
        % 
        % ARGUMENTS
        % =========
        % MUX_SET            : Scalar value. Demultiplexer mode.
        % dlrUsing12         : 0/1, true/false. DLR = Demultiplexer Latching Relay.
        %                       False=0 = Using diffs V13_DC, V13_AC
        %                       True =1 = Using diffs V12_DC, V12_AC
        % BltsSamplesCalibVolt : Cell array of matrices, length 5. {iBlts} = Vector with sample values for that channel.
        %                        BIAS calibrated volts.
        % --
        % NOTE: No argument for diff gain since this function does not calibrate.
        %
        %
        % RETURN VALUES
        % =============
        % BltsAsrType : Struct array. (iBlts) = Number representing ASR type of the BLTS data, which depends on the mux mode.
        %               Has fields
        %                   .antennas = Numeric vector of length 0, 1 or 2.
        %                           Either [] (no signal, e.g. BIAS_4/5 for TDS), [iAnt] (single), or [iAnt1, iAnt2] (diff).
        %                           NOTE: iAnt1 < iAnt2. iAnt/iAnt1/iAnt2 = {1,2,3}.
        %                           Represents the current routing of signals.
        %                   .category = String constant representing the category/type of signal on the channel.
        %                           DC single, DC diff, AC low-gain, AC high-gain, no signal
        % AsrSamplesVolt
        %             : All representations of antenna signals which can possibly be derived from the BLTS (BIAS_i).
        %               Struct with fields named as in the BIAS specification: .Vi_LF, .Vij_LF, .Vij_LF_AC
        %               NOTE: Calibration signals GND and 2.5V Ref which are generated internally by BIAS are also
        %               stored in these variables although they are technically not ASRs. See implementation.
        %
        %
        % DEFINITIONS, NAMING CONVENTIONS
        % ===============================
        % See bicas.calib.
        %
        function [BltsAsrType, AsrSamplesVolt] ...
                = demultiplexer(MUX_SET, dlrUsing12, BltsSamplesCalibVolt)
            % PROPOSAL: Function name that implies constant settings (MUX_SET at least; DIFF_GAIN?!).
            % PROPOSAL: Convention for separating actual signal data/samples from signal "type".
            %   PROPOSAL: "samples" vs "type"
            % PROPOSAL/NOTE: BIAS calibrated volts = ASR volts (automatically for those ASR for which there is BLTS data)
            % TODO-DECISION: How handle calibration modes with fixed, constant BIAS-LFR/TDS signals?
            %
            % PROBLEM: BltsAsrType.category for AC can not include low-gain/high-gain which leads to different set of
            % alternatives than used for selecting transfer functions.
            % PROPOSAL: "Assertion" for using good combination of mux mode and latching relay. Log warning if assertion
            %           fails.
            % PROPOSAL: Use string constants (calib.m?).
            % PROPOSAL: Assertions for returned string constants (only legal).
            
            % ASSERTIONS
            assert(isscalar(MUX_SET))
            assert(isscalar(dlrUsing12))
            assert(iscell(BltsSamplesCalibVolt))
            EJ_library.utils.assert.vector(BltsSamplesCalibVolt)
            assert(numel(BltsSamplesCalibVolt)==5)
            
            % Cv = (BIAS) Calibrated (BLT) volt
            BIAS_1_Cv = BltsSamplesCalibVolt{1};
            BIAS_2_Cv = BltsSamplesCalibVolt{2};
            BIAS_3_Cv = BltsSamplesCalibVolt{3};
            BIAS_4_Cv = BltsSamplesCalibVolt{4};
            BIAS_5_Cv = BltsSamplesCalibVolt{5};
            
            NAN_VALUES = ones(size(BIAS_1_Cv)) * NaN;
            As.V1_LF     = NAN_VALUES;
            As.V2_LF     = NAN_VALUES;
            As.V3_LF     = NAN_VALUES;
            As.V12_LF    = NAN_VALUES;
            As.V13_LF    = NAN_VALUES;
            As.V23_LF    = NAN_VALUES;
            As.V12_LF_AC = NAN_VALUES;
            As.V13_LF_AC = NAN_VALUES;
            As.V23_LF_AC = NAN_VALUES;

            

            if dlrUsing12;   iAntB = 2;
            else             iAntB = 3;
            end
           
            import bicas.proc_sub.routing
            
            % NOTE: BLTS 5 = V23_LF_AC for all modes, but has written it out anyway for completeness.
            switch(MUX_SET)
                case 0   % "Standard operation" : We have all information.

                    % Summarize the routing.
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [1,iAntB], 'DC diff',   BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [2,3],     'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % Derive the ASR:s not in the BLTS.
                    if dlrUsing12
                        As.V13_LF    = As.V12_LF    + As.V23_LF;
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF    = As.V13_LF    - As.V23_LF;
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end
                    As.V2_LF     = As.V1_LF     - As.V12_LF;
                    As.V3_LF     = As.V2_LF     - As.V23_LF;
                    
                case 1   % Probe 1 fails

                    [BltsAsrType(1), As] = routing(As, [2],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [3],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [2,3],     'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % NOTE: Can not derive anything for DC. BLTS 1-3 contain redundant data.
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end
                    
                case 2   % Probe 2 fails
                    
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [3],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [1,iAntB], 'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % NOTE: Can not derive anything for DC. BLTS 1-3 contain redundant data.
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end
                    
                case 3   % Probe 3 fails
                    
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [2],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [1,iAntB], 'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % NOTE: Can not derive anything for DC. BLTS 1-3 contain redundant data.
                    if dlrUsing12
                        As.V13_LF_AC = V12_LF_AC + V23_LF_AC;
                    else
                        As.V12_LF_AC = V13_LF_AC - V23_LF_AC;
                    end
                    
                case 4   % Calibration mode 0
                    
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [2],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [3],       'DC single', BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    As.V12_LF    = As.V1_LF    - As.V2_LF;
                    As.V13_LF    = As.V1_LF    - As.V3_LF;
                    As.V23_LF    = As.V2_LF    - As.V3_LF;
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end

                case {5,6,7}   % Calibration mode 1/2/3
                    
                    switch(MUX_SET)
                        case 5
                            signalTypeCategory = '2.5V Ref';
                        case {6,7}
                            signalTypeCategory = 'GND';
                    end
                    
                    % NOTE: It is in principle arbitrary (probably) how the GND and 2.5V Ref signals, which are
                    % generated by the instrument, should be represented in the datasets, since the datasets assume that
                    % only assumes signals from the antennas. The implementation classifies them as antennas, including
                    % for diffs, but the signalTypeCategory specifies that they should be calibrated differently.
                    [BltsAsrType(1), As] = routing(As, [1],       signalTypeCategory, BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [2],       signalTypeCategory, BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [3],       signalTypeCategory, BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',               BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',               BIAS_5_Cv);

                    As.V12_LF    = As.V1_LF    - As.V2_LF;
                    As.V13_LF    = As.V1_LF    - As.V3_LF;
                    As.V23_LF    = As.V2_LF    - As.V3_LF;
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end

                otherwise
%                     if isnan(MUX_SET)
%                         % Do nothing. Allow the default values (NaN) to be returned.
%                     else
                        error('BICAS:proc_sub:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value for mux_set.')
%                     end
            end   % switch
            
            AsrSamplesVolt = As;
            
            assert(numel(BltsAsrType) == 5)
        end
        
        
        
        % Utility function for "demultiplexer".
        function [BltsAsrType, AsrSamples] = routing(AsrSamples, antennas, category, BltsSamples)
            
            % Normalize vector to row vector since "isequal" is sensitive to row/column vectors.
            antennas = antennas(:)';
            
            % Assign BltsType.
            BltsAsrType.antennas = antennas;
            BltsAsrType.category = category;
            
            % Modify AsrSamples (and assertion on arguments).
            if     isequal(antennas, [1])   && strcmp(category, 'DC single')   AsrSamples.V1_LF     = BltsSamples;
            elseif isequal(antennas, [2])   && strcmp(category, 'DC single')   AsrSamples.V2_LF     = BltsSamples;
            elseif isequal(antennas, [3])   && strcmp(category, 'DC single')   AsrSamples.V3_LF     = BltsSamples;
            elseif isequal(antennas, [1,2]) && strcmp(category, 'DC diff')     AsrSamples.V12_LF    = BltsSamples;
            elseif isequal(antennas, [1,3]) && strcmp(category, 'DC diff')     AsrSamples.V13_LF    = BltsSamples;
            elseif isequal(antennas, [2,3]) && strcmp(category, 'DC diff')     AsrSamples.V23_LF    = BltsSamples;
            elseif isequal(antennas, [1,2]) && strcmp(category, 'AC')          AsrSamples.V12_LF_AC = BltsSamples;
            elseif isequal(antennas, [1,3]) && strcmp(category, 'AC')          AsrSamples.V13_LF_AC = BltsSamples;
            elseif isequal(antennas, [2,3]) && strcmp(category, 'AC')          AsrSamples.V23_LF_AC = BltsSamples;
            else
                error('BICAS:proc_SUB:Assertion:IllegalArgument', 'Illegal combination of arguments antennas and category.')
            end
        end
        
        
        
        % Automatic test code.
        %
        % Very basic tests at this stage. Could be improved but unsure how much is meaningful.
        function demultiplexer___ATEST            
            
            new_test = @(inputs, outputs) (EJ_library.atest.CompareFuncResult(@bicas.proc_sub.demultiplexer, inputs, outputs));
            tl = {};
            
            V1   = 10;
            V2   = 11;
            V3   = 12;
            V12  = V1-V2;
            V13  = V1-V3;
            V23  = V2-V3;
            V12a = 45-56;
            V13a = 45-67;
            V23a = 56-67;

            function AsrSamplesVolt = ASR_samples(varargin)
                assert(nargin == 9)
                AsrSamplesVolt = struct(...
                    'V1_LF',     as(varargin{1}, V1), ...
                    'V2_LF',     as(varargin{2}, V2), ...
                    'V3_LF',     as(varargin{3}, V3), ...
                    'V12_LF',    as(varargin{4}, V12), ...
                    'V13_LF',    as(varargin{5}, V13), ...
                    'V23_LF',    as(varargin{6}, V23), ...
                    'V12_LF_AC', as(varargin{7}, V12a), ...
                    'V13_LF_AC', as(varargin{8}, V13a), ...
                    'V23_LF_AC', as(varargin{9}, V23a));
                
                function V = as(v,V)    % as = assign. Effectively implements ~ternary operator + constant (NaN).
                    if v; V = V;
                    else  V = NaN;
                    end
                end
            end
            
            if 1
                tl{end+1} = new_test({0, true, {V1, V12, V23, V12a, V23a}}, ...
                    {struct(...
                    'antennas', {[1], [1 2], [2 3], [1 2], [2 3]}, ...
                    'category', {'DC single', 'DC diff', 'DC diff', 'AC', 'AC'}), ...
                    ASR_samples(1,1,1, 1,1,1, 1,1,1)});
            end
            
            if 1
                tl{end+1} = new_test({1, false, {V2, V3, V23, V13a, V23a}}, ...
                    {struct(...
                    'antennas', {[2], [3], [2 3], [1 3], [2 3]}, ...
                    'category', {'DC single', 'DC single', 'DC diff', 'AC', 'AC'}), ...
                    ASR_samples(0,1,1, 0,0,1, 1,1,1)});
            end
            
            EJ_library.atest.run_tests(tl)
        end



        % Add probe signals that can be derived from already known probe signals.
        %
        % ARGUMENTS
        % =========
        % probeSignals             : Struct with an arbitrary subset of the fields ... . Fields must have the same array
        % sizes.
        % complementedProbeSignals
%         function ComplementedProbeSignals = complement_probe_signals(ProbeSignals)
%             % TODO: Lägg till helt tomma signaler, NaN.
% %                     % Derive the OUTPUT DATA which are less trivial.
% %                     V13_LF    = V12_LF    + V23_LF;
% %                     V2_LF     = V1_LF     - V12_LF;
% %                     V3_LF     = V2_LF     - V23_LF;
% %                     V13_LF_AC = V12_LF_AC + V23_LF_AC;
% %                     % Derive the OUTPUT DATA which are less trivial.
% %                     V12_LF    = V1_LF     - V2_LF;
% %                     V13_LF    = V1_LF     - V3_LF;
% %                     V23_LF    = V2_LF     - V3_LF;
% %                     V13_LF_AC = V12_LF_AC + V23_LF_AC;
% %                     % Derive the OUTPUT DATA which are less trivial.
% %                     V13_LF_AC = V12_LF_AC + V23_LF_AC;
%             
% %             if can_derive_signal(ProbeSignals, 'V13', 'V12', 'V23')
% %                 ProbeSignals.V13 = 
% %             elseif
% 
%             ProbeSignals = derive_signal_if_possible(ProbeSignals, 'V13', 'V12', 'V23', @(x1,x2) (x1+x2));
%             ProbeSignals = derive_signal_if_possible(ProbeSignals, 'V2',  'V1',  'V12', @(x1,x2) (x1-x2));
%             ComplementedProbeSignals = ProbeSignals;
%             
%             function ProbeSignals = derive_signal_if_possible(ProbeSignals, outputFieldName, inputFieldName1, inputFieldName2, funcPtr)
%                 if isfield(ProbeSignals, inputFieldName1) ...
%                     &&  isfield(ProbeSignals, inputFieldName2) ...
%                     && ~isfield(ProbeSignals, outputFieldName);
%                     ProbeSignals.(outputFieldName) = funcPtr(ProbeSignals.(inputFieldName1), ProbeSignals.(inputFieldName2));
%                 else
%                     ;   % Do nothing. Return the same ProbeSignals.
%                 end
%             end
%         end

    end   % methods(Static, Access=private)
        
end
