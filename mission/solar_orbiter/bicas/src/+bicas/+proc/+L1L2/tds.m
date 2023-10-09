%
% Collection of TDS-related processing functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-25
%
classdef tds    
    % PROPOSAL: Automatic test code.



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Processing function. Only "normalizes" data to account for technically
        % illegal input TDS datasets. It should try to:
        % ** modify L1 to look like L1R
        % ** mitigate historical bugs in the input datasets
        % ** mitigate for not yet implemented features in input datasets
        %
        function InSciNorm = process_normalize_CDF(InSci, inSciDsi, SETTINGS, L)

            % Default behaviour: Copy values, except for values which are
            % modified later
            InSciNorm = InSci;

            nRecords = irf.assert.sizes(InSci.Zv.Epoch, [-1]);

            C = bicas.classify_BICAS_L1_L1R_to_L2_DSI(inSciDsi);



            %===================================
            % Normalize CALIBRATION_TABLE_INDEX
            %===================================
            InSciNorm.Zv.CALIBRATION_TABLE_INDEX = bicas.proc.L1L2.normalize_CALIBRATION_TABLE_INDEX(...
                InSci.Zv, nRecords, inSciDsi);



            %===========================================================
            % Normalize zVar name SYNCHRO_FLAG
            % --------------------------------
            % Both ZVs TIME_SYNCHRO_FLAG, SYNCHRO_FLAG found in input
            % datasets. Unknown why. "DEFINITION BUG" in definition of
            % datasets/skeleton? /2020-01-05
            % Based on skeletons (.skt; L1R, L2), SYNCHRO_FLAG seems
            % to be the correct one. /2020-01-21
            %===========================================================
            [InSci.Zv, fnChangeList] = irf.ds.normalize_struct_fieldnames(...
                InSci.Zv, ...
                {{{'TIME_SYNCHRO_FLAG', 'SYNCHRO_FLAG'}, 'SYNCHRO_FLAG'}}, ...
                'Assert one matching candidate');

            bicas.proc.utils.handle_ZV_name_change(...
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
                                anomalyDescrMsg, ...
                                'BICAS:DatasetFormat')
                    end
                else
                    error('BICAS:DatasetFormat', anomalyDescrMsg)
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
                zv_SAMPS_PER_CH_corrected = max( zv_SAMPS_PER_CH_corrected, solo.hwzv.const.TDS_RSWF_SNAPSHOT_LENGTH_MIN);
                zv_SAMPS_PER_CH_corrected = min( zv_SAMPS_PER_CH_corrected, solo.hwzv.const.TDS_RSWF_SNAPSHOT_LENGTH_MAX);

                if any(zv_SAMPS_PER_CH_corrected ~= InSci.Zv.SAMPS_PER_CH)
                    % CASE: SAMPS_PER_CH has at least one illegal value

                    SAMPS_PER_CH_badValues = unique(InSci.Zv.SAMPS_PER_CH(...
                        zv_SAMPS_PER_CH_corrected ~= InSci.Zv.SAMPS_PER_CH));

                    badValuesDisplayStr = strjoin(arrayfun(...
                        @(n) sprintf('%i', n), SAMPS_PER_CH_badValues, 'uni', false), ', ');
                    anomalyDescrMsg = sprintf(...
                        ['TDS LFM RSWF zVar SAMPS_PER_CH contains unexpected', ...
                        ' value(s) which are not on the form 2^n and in the', ...
                        ' interval %.0f to %.0f: %s'], ...
                        solo.hwzv.const.TDS_RSWF_SNAPSHOT_LENGTH_MIN, ...
                        solo.hwzv.const.TDS_RSWF_SNAPSHOT_LENGTH_MAX, ...
                        badValuesDisplayStr);

                    [settingValue, settingKey] = SETTINGS.get_fv(...
                        'PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY');
                    switch(settingValue)
                        
                        case 'ROUND'
                            bicas.default_anomaly_handling(...
                                L, settingValue, settingKey, 'other', ...
                                anomalyDescrMsg, ...
                                'BICAS:Assertion:DatasetFormat')
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
                                'BICAS:Assertion:DatasetFormat')

                    end    % switch
                end    % if
            end    % if

        end    % process_normalize_CDF



        % Processing function. Convert TDS CDF data (PDs) to PreDC.
        function PreDc = process_CDF_to_PreDC(InSci, inSciDsi, HkSciTime, SETTINGS, L)
        %
        % BUG?: Does not use CHANNEL_STATUS_INFO.
        % NOTE: BIAS output datasets do not have a variable for the length of
        % snapshots. Need to use NaN/fill value.

            % ASSERTIONS: VARIABLES
            assert(isa(InSci, 'bicas.InputDataset'))
            irf.assert.struct(HkSciTime, {'bdmFpa', 'biasHighGainFpa', 'dlrFpa'}, {})

            C = bicas.classify_BICAS_L1_L1R_to_L2_DSI(inSciDsi);



            % ASSERTIONS: CDF
            bicas.proc.utils.assert_increasing(...
                InSci.Zv.Epoch, true, 'BICAS:DatasetFormat', ...
                ['Voltage (science) dataset timestamps Epoch do not', ...
                ' increase monotonously.']...
            )
            [nRecords, WAVEFORM_DATA_nChannels, nCdfSamplesPerRecord] = irf.assert.sizes(...
                InSci.Zv.Epoch,         [-1], ...
                InSci.Zv.WAVEFORM_DATA, [-1, -2, -3]);
            if     C.isL1r   WAVEFORM_DATA_nChannels_expected = 3;
            elseif C.isL1    WAVEFORM_DATA_nChannels_expected = 8;
            end
            assert(...
                WAVEFORM_DATA_nChannels == WAVEFORM_DATA_nChannels_expected, ...
                'BICAS:Assertion:DatasetFormat', ...
                'TDS zVar WAVEFORM_DATA has an unexpected size.')

            % IMPORTANT NOTE
            % ==============
            % Empirically, this value varies between L1R and L1 datasets. Since
            % BICAS does not officially support reading L1 datasets, the code
            % should use the value for L1R. The code is also (probably) not
            % adapted to having another number of samples/record. This means
            % that BICAS CURRENTLY (2023-10-09) CAN NOT READ L1
            % SOLO_L1_RPW-TDS-LFM-RSWF datasets.
            % Ex:
            %     solo_L1R_rpw-tds-lfm-rswf-e-cdag_20200409_V12.cdf: 32768 samples/record
            %     solo_L1_rpw-tds-lfm-rswf-cdag_20200409_V09.cdf   : 16384 samples/record
            if C.isTdsRswf
                assert(...
                    nCdfSamplesPerRecord == solo.hwzv.const.TDS_RSWF_L1R_SAMPLES_PER_RECORD, ...
                    'Unexpected number of samples per CDF record (%i). Expected %i.', ...
                    nCdfSamplesPerRecord, solo.hwzv.const.TDS_RSWF_L1R_SAMPLES_PER_RECORD)
            else
                assert(nCdfSamplesPerRecord == 1)
            end



            % TODO-NI: Why convert to double? To avoid precision problems when
            % doing math with other variables?
            freqHzZv = double(InSci.Zv.SAMPLING_RATE);



            Zv    = [];

            Zv.Epoch                   = InSci.Zv.Epoch;
            Zv.DELTA_PLUS_MINUS        = bicas.proc.utils.derive_DELTA_PLUS_MINUS(...
                freqHzZv, nCdfSamplesPerRecord);
            Zv.freqHz                  = freqHzZv;
            Zv.QUALITY_BITMASK         = InSci.ZvFpa.QUALITY_BITMASK;
            Zv.QUALITY_FLAG            = InSci.ZvFpa.QUALITY_FLAG;
            Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
            Zv.bdmFpa                  = HkSciTime.bdmFpa;
            Zv.biasHighGainFpa         = HkSciTime.biasHighGainFpa;
            Zv.dlrFpa                  = HkSciTime.dlrFpa;
            Zv.ufv                     = false(nRecords, 1);
            Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;



            %=====================================
            % Set Zv.nValidSamplesPerRecord
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
                Zv.nValidSamplesPerRecord = double(InSci.Zv.SAMPS_PER_CH);
            else
                Zv.nValidSamplesPerRecord = ones(nRecords, 1) * 1;
            end
            assert(all(Zv.nValidSamplesPerRecord <= nCdfSamplesPerRecord), ...
                'BICAS:Assertion:DatasetFormat', ...
                ['Dataset indicates that the number of valid samples per CDF', ...
                ' record (max(Zv.nValidSamplesPerRecord)=%i) is', ...
                ' NOT fewer than the number of indices per CDF record', ...
                ' (nCdfMaxSamplesPerSnapshot=%i).'], ...
                max(Zv.nValidSamplesPerRecord), ...
                nCdfSamplesPerRecord)



            %======================
            % Set Zv.bltsSamplesTm
            %======================
            zv_WAVEFORM_DATA_modif = double(permute(InSci.Zv.WAVEFORM_DATA, [1,3,2]));

%             Zv.bltsSamplesTmCa    = cell(5,1);
%             Zv.bltsSamplesTmCa{1} = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,1), Zv.nValidSamplesPerRecord );
%             Zv.bltsSamplesTmCa{2} = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,2), Zv.nValidSamplesPerRecord );
%             Zv.bltsSamplesTmCa{3} = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,3), Zv.nValidSamplesPerRecord );
%             Zv.bltsSamplesTmCa{4} = nan(nRecords, nCdfSamplesPerRecord);
%             Zv.bltsSamplesTmCa{5} = nan(nRecords, nCdfSamplesPerRecord);
            Zv.bltsSamplesTm(:, :, 1) = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,1), Zv.nValidSamplesPerRecord );
            Zv.bltsSamplesTm(:, :, 2) = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,2), Zv.nValidSamplesPerRecord );
            Zv.bltsSamplesTm(:, :, 3) = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,3), Zv.nValidSamplesPerRecord );
            Zv.bltsSamplesTm(:, :, 4) = nan(nRecords, nCdfSamplesPerRecord);
            Zv.bltsSamplesTm(:, :, 5) = nan(nRecords, nCdfSamplesPerRecord);



            Ga = [];
            Ga.OBS_ID         = InSci.Ga.OBS_ID;
            Ga.SOOP_TYPE      = InSci.Ga.SOOP_TYPE;

            % Only set because the code shared with LFR requires it.
            Zv.iLsf           = nan(nRecords, 1);
            Zv.lfrRx          = ones(nRecords, 1);
            
            PreDc = bicas.proc.L1L2.PreDc(Zv, Ga, C.isTdsRswf, false, C.isTdsCwf);
        end    % process_CDF_to_PreDC



        % Processing function. Convert PreDc+PostDC to something that 
        % (1) represents a TDS dataset (hence the name), and
        % (2) ALMOST REPRESENTS an LFR dataset (the rest is done in a wrapper).
        %
        % This function only changes the data format (and selects data to send
        % to CDF).
        %
        function [OutSci] = process_PostDC_to_CDF(SciPreDc, SciPostDc, outputDsi, L)
            % PROPOSAL: Rename to something shared between LFR and TDS, then use
            %           two wrappers.
            %   PROPOSAL: process_PostDC_to_LFR_TDS_CDF_core
            %   TODO-DEC: Put in which future file?

            % ASSERTIONS
            assert(isa(SciPreDc,  'bicas.proc.L1L2.PreDc'))
            assert(isa(SciPostDc, 'bicas.proc.L1L2.PostDc'))



            nRecords                 = size(SciPreDc.Zv.Epoch, 1);
            nSamplesPerRecordChannel = size(SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V1'), 2);

            OutSci = [];

            OutSci.Zv.Epoch              = SciPreDc.Zv.Epoch;
            OutSci.Zv.QUALITY_BITMASK    = SciPreDc.Zv.QUALITY_BITMASK;
            OutSci.Zv.L2_QUALITY_BITMASK = SciPostDc.Zv.L2_QUALITY_BITMASK;
            OutSci.Zv.QUALITY_FLAG       = SciPostDc.Zv.QUALITY_FLAG;
            OutSci.Zv.DELTA_PLUS_MINUS   = SciPreDc.Zv.DELTA_PLUS_MINUS;
            OutSci.Zv.SYNCHRO_FLAG       = SciPreDc.Zv.SYNCHRO_FLAG;
            OutSci.Zv.SAMPLING_RATE      = SciPreDc.Zv.freqHz;

            % NOTE: Convert aampere --> nano-aampere
            OutSci.Zv.IBIAS1 = SciPostDc.Zv.currentAAmpere(:, 1) * 1e9;
            OutSci.Zv.IBIAS2 = SciPostDc.Zv.currentAAmpere(:, 2) * 1e9;
            OutSci.Zv.IBIAS3 = SciPostDc.Zv.currentAAmpere(:, 3) * 1e9;

            OutSci.Ga.OBS_ID    = SciPreDc.Ga.OBS_ID;
            OutSci.Ga.SOOP_TYPE = SciPreDc.Ga.SOOP_TYPE;



            C = bicas.classify_BICAS_L1_L1R_to_L2_DSI(outputDsi);

            % NOTE: The two cases are different in the indexes they use for
            % OutSciZv.
            if C.isCwf

                % ASSERTIONS
                assert(nSamplesPerRecordChannel == 1, ...
                    'BICAS:Assertion:IllegalArgument', ...
                    ['Number of samples per CDF record is not 1, as expected.', ...
                    ' Bad input CDF?'])
                irf.assert.sizes(...
                    OutSci.Zv.QUALITY_BITMASK, [nRecords, 1], ...
                    OutSci.Zv.QUALITY_FLAG,    [nRecords, 1])

                % Try to pre-allocate to save RAM/speed up.
                tempNaN = nan(nRecords, 3);
                OutSci.Zv.VDC = tempNaN;
                OutSci.Zv.EDC = tempNaN;
                OutSci.Zv.EAC = tempNaN;

                OutSci.Zv.VDC(:,1) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V1');
                OutSci.Zv.VDC(:,2) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V2');
                OutSci.Zv.VDC(:,3) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V3');

                OutSci.Zv.EDC(:,1) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V12');
                OutSci.Zv.EDC(:,2) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V13');
                OutSci.Zv.EDC(:,3) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V23');

                OutSci.Zv.EAC(:,1) = SciPostDc.Zv.AsrSamplesAVoltSrm('AC_V12');
                OutSci.Zv.EAC(:,2) = SciPostDc.Zv.AsrSamplesAVoltSrm('AC_V13');
                OutSci.Zv.EAC(:,3) = SciPostDc.Zv.AsrSamplesAVoltSrm('AC_V23');

            elseif C.isSwf

                if     C.isLfr
                    SAMPLES_PER_RECORD_CHANNEL = ...
                        solo.hwzv.const.LFR_SWF_SNAPSHOT_LENGTH;
                elseif C.isTds
                    SAMPLES_PER_RECORD_CHANNEL = ...
                        solo.hwzv.const.TDS_RSWF_L1R_SAMPLES_PER_RECORD;
                else
                    error(...
                        'BICAS:Assertion', ...
                        'Illegal DSI classification.')
                end

                % ASSERTION
                assert(nSamplesPerRecordChannel == SAMPLES_PER_RECORD_CHANNEL, ...
                    'BICAS:Assertion:IllegalArgument', ...
                    ['Number of samples per CDF record (%i) is not', ...
                    ' %i, as expected. Bad Input CDF?'], ...
                    nSamplesPerRecordChannel, ...
                    SAMPLES_PER_RECORD_CHANNEL)

                % Try to pre-allocate to save RAM/speed up.
                tempNaN = nan(nRecords, nSamplesPerRecordChannel, 3);
                OutSci.Zv.VDC = tempNaN;
                OutSci.Zv.EDC = tempNaN;
                OutSci.Zv.EAC = tempNaN;

                OutSci.Zv.VDC(:,:,1) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V1');
                OutSci.Zv.VDC(:,:,2) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V2');
                OutSci.Zv.VDC(:,:,3) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V3');

                OutSci.Zv.EDC(:,:,1) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V12');
                OutSci.Zv.EDC(:,:,2) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V13');
                OutSci.Zv.EDC(:,:,3) = SciPostDc.Zv.AsrSamplesAVoltSrm('DC_V23');

                OutSci.Zv.EAC(:,:,1) = SciPostDc.Zv.AsrSamplesAVoltSrm('AC_V12');
                OutSci.Zv.EAC(:,:,2) = SciPostDc.Zv.AsrSamplesAVoltSrm('AC_V13');
                OutSci.Zv.EAC(:,:,3) = SciPostDc.Zv.AsrSamplesAVoltSrm('AC_V23');

            else
                error('BICAS:Assertion:IllegalArgument', ...
                    'Function can not produce outputDsi=%s.', outputDsi)
            end



            % ASSERTION
            bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(OutSci.Zv);
            % NOTE: Not really necessary since the list of ZVs will be checked
            % against the master CDF?
            irf.assert.struct(OutSci.Zv, {...
                'IBIAS1', 'IBIAS2', 'IBIAS3', 'VDC', 'EDC', 'EAC', 'Epoch', ...
                'QUALITY_BITMASK', 'L2_QUALITY_BITMASK', 'QUALITY_FLAG', ...
                'DELTA_PLUS_MINUS', 'SYNCHRO_FLAG', 'SAMPLING_RATE'}, {})

        end    % process_PostDC_to_CDF
        
        
        
    end    % methods(Static)



end
