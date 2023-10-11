%
% Collection of LFR-related processing functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-25, from reorganized older code.
%
classdef lfr    
    % PROPOSAL: Automatic test code.

    
    
    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Processing function. Only "normalizes" data to account for technically
        % illegal input LFR datasets. It should try to:
        % ** modify L1 to look like L1R
        % ** mitigate historical bugs in input datasets
        % ** mitigate for not yet implemented features in input datasets
        %
        function InSciNorm = process_normalize_CDF(InSci, inSciDsi, SETTINGS, L)

            % Default behaviour: Copy values, except for values which are
            % modified later
            InSciNorm = InSci;

            nRecords = irf.assert.sizes(InSci.Zv.Epoch, [-1]);



            %===================================
            % Normalize CALIBRATION_TABLE_INDEX
            %===================================
            InSciNorm.Zv.CALIBRATION_TABLE_INDEX = ...
                bicas.proc.L1L2.normalize_CALIBRATION_TABLE_INDEX(...
                    InSci.Zv, nRecords, inSciDsi);



            %========================
            % Normalize SYNCHRO_FLAG
            %========================
            has_SYNCHRO_FLAG      = isfield(InSci.Zv, 'SYNCHRO_FLAG');
            has_TIME_SYNCHRO_FLAG = isfield(InSci.Zv, 'TIME_SYNCHRO_FLAG');
            if      has_SYNCHRO_FLAG && ~has_TIME_SYNCHRO_FLAG

                % CASE: Everything nominal.
                InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.SYNCHRO_FLAG;

            elseif ~has_SYNCHRO_FLAG && has_TIME_SYNCHRO_FLAG

                % CASE: Input CDF uses wrong zVar name.
                [settingValue, settingKey] = ...
                    SETTINGS.get_fv('INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY');
                bicas.default_anomaly_handling(L, ...
                    settingValue, settingKey, 'E+W+illegal', ...
                    'Found zVar TIME_SYNCHRO_FLAG instead of SYNCHRO_FLAG.')
                L.log('warning', ...
                    'Using illegally named zVar TIME_SYNCHRO_FLAG as SYNCHRO_FLAG.')
                InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.TIME_SYNCHRO_FLAG;

            elseif has_SYNCHRO_FLAG && has_TIME_SYNCHRO_FLAG

                % CASE: Input CDF has two ZVs: one with correct name, one with
                % incorrect name

                %------------------------
                % "Normal" normalization
                %------------------------
                % 2020-01-21: Based on skeletons (.skt; L1R, L2), SYNCHRO_FLAG
                % seems to be the correct zVar.
                if SETTINGS.get_fv(...
                        'INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED') ...
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
                    error('BICAS:DatasetFormat', ...
                        ['Input dataset has both zVar SYNCHRO_FLAG and', ...
                        ' TIME_SYNCHRO_FLAG.'])
                end
            else
                error('BICAS:DatasetFormat', ...
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
            [settingValue, settingKey] = SETTINGS.get_fv(...
                'PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY');

            InSciNorm.ZvFpa.QUALITY_BITMASK = bicas.proc.L1L2.lfr.normalize_ZV_empty(...
                L, settingValue, settingKey, nRecords, ...
                InSci.ZvFpa.QUALITY_BITMASK, 'QUALITY_BITMASK');

            InSciNorm.ZvFpa.QUALITY_FLAG    = bicas.proc.L1L2.lfr.normalize_ZV_empty(...
                L, settingValue, settingKey, nRecords, ...
                InSci.ZvFpa.QUALITY_FLAG,    'QUALITY_FLAG');
            
            % NOTE: Very old L1R LFR-SBM1/2 (test) data datasets have been
            % observed to have QUALITY_BITMASK with illegal data type
            % CDF_UINT1/uint8 while newer ones do not. This is not mitigated in
            % the code since it does not seem to be needed.

            % ASSERTIONS
            irf.assert.sizes(...
                InSciNorm.ZvFpa.QUALITY_BITMASK, [nRecords, 1], ...
                InSciNorm.ZvFpa.QUALITY_FLAG,    [nRecords, 1])
            
        end    % process_normalize_CDF



        % Processing function. Convert LFR CDF data to PreDC.
        %
        % IMPLEMENTATION NOTE: Does not modify InSci in an attempt to save RAM
        % (should help MATLAB's optimization). Unclear if actually works.
        %
        function PreDc = process_CDF_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
            %
            % PROBLEM: Hard-coded CDF data types (MATLAB classes).
            % MINOR PROBLEM: Still does not handle LFR zVar TYPE for determining
            % "virtual snapshot" length. Should only be relevant for
            % V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF (not V02) which should expire.

            % ASSERTIONS: VARIABLES
            assert(isa(InSci, 'bicas.InputDataset'))
            irf.assert.struct(HkSciTime, {'bdmFpa', 'biasHighGainFpa', 'dlrFpa'}, {})
            
            % ASSERTIONS: CDF
            bicas.proc.utils.assert_increasing(...
                InSci.Zv.Epoch, true, 'BICAS:DatasetFormat', ...
                ['Voltage (science) dataset timestamps Epoch do not', ...
                ' increase monotonously.']...
            )
            nRecords = irf.assert.sizes(InSci.Zv.Epoch, [-1]);



            C = bicas.classify_BICAS_L1_L1R_to_L2_DSI(inSciDsi);



            %============
            % Set iLsfZv
            %============
            if     C.isLfrSbm1   iLsfZv = ones(nRecords, 1) * 2;   % Always value "2" (F1, "FREQ = 1").
            elseif C.isLfrSbm2   iLsfZv = ones(nRecords, 1) * 3;   % Always value "3" (F2, "FREQ = 2").
            else                 iLsfZv = InSci.Zv.FREQ + 1;
                % NOTE: Translates from LFR's FREQ values (0=F0 etc) to LSF
                % index values (1=F0) used in loaded RCT data structs.
            end
            irf.assert.sizes(iLsfZv, [nRecords])



            % NOTE: Needed also for 1 SPR.
            zvFreqHz = solo.hwzv.get_LFR_frequency( iLsfZv );

            % Obtain the relevant values (one per record) from zVariables R0,
            % R1, R2, and the virtual "R3".
            zvRx = solo.hwzv.get_LFR_Rx(...
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
            nCdfSamplesPerRecord = irf.assert.sizes(...
                InSci.Zv.V, [nRecords, -1], ...
                E,          [nRecords, -1, 2]);
            if C.isLfrSurvSwf   assert(nCdfSamplesPerRecord == solo.hwzv.const.LFR_SWF_SNAPSHOT_LENGTH)
            else                assert(nCdfSamplesPerRecord == 1)
            end



            Zv    = [];

%             Zv.bltsSamplesTmCa    = cell(5,1);
%             Zv.bltsSamplesTmCa{1} = single(InSci.Zv.V);
%             % Copy values, except when zvRx==0 (==>NaN).
%             Zv.bltsSamplesTmCa{2} = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvRx==0 );
%             Zv.bltsSamplesTmCa{3} = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvRx==0 );
%             Zv.bltsSamplesTmCa{4} = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvRx==1 );
%             Zv.bltsSamplesTmCa{5} = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvRx==1 );
            Zv.bltsSamplesTm(:, :, 1) = single(InSci.Zv.V);
            % Copy values, except when zvRx==0 (==>NaN).
            Zv.bltsSamplesTm(:, :, 2) = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvRx==0 );
            Zv.bltsSamplesTm(:, :, 3) = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvRx==0 );
            Zv.bltsSamplesTm(:, :, 4) = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvRx==1 );
            Zv.bltsSamplesTm(:, :, 5) = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvRx==1 );
            
            Zv.Epoch                   = InSci.Zv.Epoch;
            Zv.DELTA_PLUS_MINUS        = bicas.proc.utils.derive_DELTA_PLUS_MINUS(...
                zvFreqHz, nCdfSamplesPerRecord);
            Zv.freqHz                  = zvFreqHz;
            Zv.nValidSamplesPerRecord  = ones(nRecords, 1) * nCdfSamplesPerRecord;
            Zv.BW                      = InSci.Zv.BW;
            Zv.ufv                     = ~logical(InSci.Zv.BW);
            Zv.biasHighGainFpa         = HkSciTime.biasHighGainFpa;
            Zv.dlrFpa                  = HkSciTime.dlrFpa;
            %Zv.dlrFpa                  = bicas.utils.FPArray(false(size(InSci.Zv.Epoch)), 'NO_FILL_POSITIONS');   % TEST: Always DLR = 0.
            Zv.iLsf                    = iLsfZv;

            Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
            Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;

            Zv.QUALITY_BITMASK         = InSci.ZvFpa.QUALITY_BITMASK;
            Zv.QUALITY_FLAG            = InSci.ZvFpa.QUALITY_FLAG;

            Zv.lfrRx                   = zvRx;



            %==========================================
            % Set BDM
            % -------
            % Select which source of mux mode is used.
            %==========================================
            [bdmSrcSettingValue, bdmSrcSettingKey] = SETTINGS.get_fv('PROCESSING.LFR.MUX_MODE_SOURCE');
            switch(bdmSrcSettingValue)
                case 'BIAS_HK'
                    L.log('debug', 'Using BIAS HK mux mode.')
                    bdmFpa = HkSciTime.bdmFpa;

                case 'LFR_SCI'
                    L.log('debug', 'Using LFR SCI mux mode.')
                    bdmFpa = InSci.ZvFpa.BIAS_MODE_MUX_SET;

                case 'BIAS_HK_LFR_SCI'
                    L.log('debug', ...
                        ['Using mux mode from BIAS HK when available, and', ...
                        ' from LFR SCI when the former is not available.'])

                    bdmFpa = HkSciTime.bdmFpa.complement(InSci.ZvFpa.BIAS_MODE_MUX_SET);

                otherwise
                    error('BICAS:ConfigurationBug', ...
                        'Illegal settings value %s="%s"', bdmSrcSettingKey, bdmSrcSettingValue)
            end
            Zv.bdmFpa = bdmFpa;



            Ga = [];
            Ga.OBS_ID    = InSci.Ga.OBS_ID;
            Ga.SOOP_TYPE = InSci.Ga.SOOP_TYPE;
            
            PreDc = bicas.proc.L1L2.PreDc(Zv, Ga, C.isLfrSurvSwf, true, false);

        end    % process_CDF_to_PreDC



        function [OutSci] = process_PostDC_to_CDF(SciPreDc, SciPostDc, outputDsi, L)
            % NOTE: Using __TDS__ function.
            OutSci = bicas.proc.L1L2.tds.process_PostDC_to_CDF(...
                SciPreDc, SciPostDc, outputDsi, L);

            OutSci.Zv.BW = SciPreDc.Zv.BW;
        end
        
        
        
    end    % methods(Static)
    
    
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Local utility function to shorten & clarify code.
        %
        % ARGUMENTS
        % =========
        % ZvFpa1
        %   zVar-like FPA. Column vector (Nx1) or empty.
        function ZvFpa2 = normalize_ZV_empty(...
                L, settingValue, settingKey, nRecords, ZvFpa1, zvName)

            if ~isempty(ZvFpa1)
                % Do nothing (except assertion later).
                ZvFpa2 = ZvFpa1;
            else
                anomalyDescrMsg = sprintf(...
                    'zVar "%s" from the LFR SCI source dataset is empty.', ...
                    zvName);

                switch(settingValue)
                    case 'USE_FILL_VALUE'
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, 'other', ...
                            anomalyDescrMsg, ...
                            'BICAS:DatasetFormat:SWMProcessing')

                        L.logf('warning', 'Using fill values for %s.', zvName)
                        ZvFpa2 = bicas.utils.FPArray(...
                            zeros(nRecords, 1, ZvFpa1.mc), 'ONLY_FILL_POSITIONS');

                    otherwise
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, 'E+illegal', ...
                            anomalyDescrMsg, ...
                            'BICAS:DatasetFormat:SWMProcessing')
                end
            end

            irf.assert.sizes(ZvFpa2, [NaN])
        end
        
        
        
    end    % methods(Static, Access=private)

end
