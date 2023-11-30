%
% Collection of processing function for demultiplexing and calibrating (DC), and
% related code (except bicas.proc.L1L2.demuxer).
%
% DC = Demux (demultiplex) & Calibrate
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-05-25
%
classdef dc
    % PROPOSAL: Better name.
    %   PRO: Processing includes quality "processing" which is not in "DC".
    %   PROPOSAL: dcq = Demux, Calibrate, Quality
    
    % PROPOSAL: Automatic test code.
    %
    % PROPOSAL:   process_calibrate_demux()
    %           & calibrate_demux_voltages()
    %           should only accept the needed ZVs and variables.
    %   NOTE: Needs some way of packaging/extracting only the relevant ZVs/fields
    %         from struct.



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)
        
        
        
        % Processing function. Derive PostDC from PreDc, i.e. demux, calibrate
        % dat, and set quality variables.
        %
        % NOTE: Public function as opposed to the other demuxing/calibration
        % functions.
        %
        function PostDc = process_calibrate_demux(PreDc, InCurPd, Cal, NsoTable, SETTINGS, L)

            tTicToc = tic();

            % ASSERTION
            assert(isa(PreDc, 'bicas.proc.L1L2.PreDc'));



            %############################
            % DEMUX & CALIBRATE VOLTAGES
            %############################
            AsrSamplesAVoltSrm = ...
                bicas.proc.L1L2.dc.calibrate_demux_voltages(PreDc, Cal, L);



            %#########################
            % Calibrate bias CURRENTS
            %#########################
            currentSAmpere = bicas.proc.L1L2.dc.convert_CUR_to_CUR_on_SCI_TIME(...
                PreDc.Zv.Epoch, InCurPd, SETTINGS, L);
            currentTm      = bicas.proc.L1L2.cal.Cal.calibrate_current_sampere_to_TM(currentSAmpere);

            currentAAmpere = nan(size(currentSAmpere));    % Variable to fill/set.
            iCalibLZv      = Cal.get_BIAS_calibration_time_L(PreDc.Zv.Epoch);
            [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(iCalibLZv);
            L.logf('info', ...
                ['Calibrating currents -', ...
                ' One sequence of records with identical settings at a time.'])
            for iSs = 1:nSs
                iRec1 = iRec1Ar(iSs);
                iRec2 = iRec2Ar(iSs);

                iRecords = iRec1:iRec2;

                L.logf('info', 'Records %7i-%7i : %s -- %s', ...
                    iRec1, iRec2, ...
                    bicas.utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iRec1)), ...
                    bicas.utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iRec2)))

                for iAnt = 1:3
                    %--------------------
                    % CALIBRATE CURRENTS
                    %--------------------
                    currentAAmpere(iRecords, iAnt) = Cal.calibrate_current_TM_to_aampere(...
                        currentTm( iRecords, iAnt), iAnt, iCalibLZv(iRecords));
                end
            end



            % ##############################################
            % Set quality variables, and apply UFV (to data)
            % ##############################################
            ZvIn = struct(...
                'Epoch',            PreDc.Zv.Epoch, ...
                'ufv',              PreDc.Zv.ufv, ...
                'bdmFpa',           PreDc.Zv.bdmFpa, ...
                'QUALITY_FLAG_Fpa', PreDc.Zv.QUALITY_FLAG);
            [zvUfv, Zv.QUALITY_FLAG, Zv.L2_QUALITY_BITMASK] = ...
                bicas.proc.L1L2.qual.modify_quality_filter(ZvIn, PreDc.isLfr, NsoTable, SETTINGS, L);
            clear ZvIn
            
            % NOTE: Function modifies AsrSamplesAVoltSrm handle object!
            Zv.currentAAmpere = bicas.proc.L1L2.qual.set_voltage_current_FV(...
                PreDc.Zv.Epoch, AsrSamplesAVoltSrm, currentAAmpere, zvUfv, L);
            Zv.AsrSamplesAVoltSrm = AsrSamplesAVoltSrm;
            clear AsrSamplesAVoltSrm
            

            
            % ############
            % END FUNCTION
            % ############
            PostDc = bicas.proc.L1L2.PostDc(Zv);

            nRecords = size(PreDc.Zv.Epoch, 1);
            bicas.log_speed_profiling(L, ...
                'bicas.proc.L1L2.dc.process_calibrate_demux', tTicToc, ...
                nRecords, 'record')
        end    % process_calibrate_demux



    end    % methods(Static)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % Demultiplex and calibrate VOLTAGES (not e.g. currents).
        %
        % NOTE: Can handle arrays of any size if the sizes are
        % consistent.
        %
        function AsrSamplesAVoltSrm = calibrate_demux_voltages(PreDc, Cal, L)
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
        % PROPOSAL: Integrate into bicas.proc.L1L2.demuxer (as method).
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
            assert(isnumeric(  PreDc.Zv.bltsSamplesTm))
%             bicas.proc.utils.assert_cell_array_comps_have_same_N_rows(...
%                 PreDc.Zv)
            [nRecords, nSamplesPerRecordChannel] = irf.assert.sizes(...
                PreDc.Zv.bdmFpa,          [-1,  1], ...
                PreDc.Zv.biasHighGainFpa, [-1,  1], ...
                PreDc.Zv.bltsSamplesTm,   [-1, -2, bicas.const.N_BLTS]);



            % Pre-allocate
            % ------------
            % IMPLEMENTATION NOTE: Preallocation is very important for speeding
            % up LFR-SWF which tends to be broken into subsequences of 1 record.
            AsrSamplesAVoltSrm = bicas.utils.SameRowsMap(...
                'char', nRecords, 'CONSTANT', ...
                nan(nRecords, nSamplesPerRecordChannel), ...
                bicas.proc.L1L2.AntennaSignalId.C.ALL_ASID_NAMES_CA);

            % Set dlrFpa
            % ----------
            if 0
                % Use hardcoded values
                % --------------------
                % IMPLEMENTATION NOTE: This exists for historical reasons and
                % should eventually be abolished, together with the function
                % called below.
                dlrFpa = bicas.proc.L1L2.demuxer_latching_relay(PreDc.Zv.Epoch);
            else
                % Use HK
                % ------
                dlrFpa = PreDc.Zv.dlrFpa;
            end
            
            iCalibLZv = Cal.get_BIAS_calibration_time_L(PreDc.Zv.Epoch);
            iCalibHZv = Cal.get_BIAS_calibration_time_H(PreDc.Zv.Epoch);



            %===================================================================
            % (1) Find continuous subsequences of records with identical
            %     settings.
            % (2) Process data separately for each such sequence.
            % ----------------------------------------------------------
            % NOTE: Just finding continuous subsequences can take a significant
            %       amount of time.
            % NOTE: Empirically, this is not useful for real LFR SWF datasets
            %       where the LFR sampling frequency changes in every record,
            %       meaning that the subsequences are all 1 record long.
            %
            % SS = Subsequence
            %===================================================================
            [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
                PreDc.Zv.bdmFpa.int2doubleNan(), ...
                PreDc.Zv.biasHighGainFpa.logical2doubleNan(), ...
                PreDc.Zv.freqHz, ...
                PreDc.Zv.iLsf, ...
                PreDc.Zv.CALIBRATION_TABLE_INDEX, ...
                PreDc.Zv.ufv, ...
                PreDc.Zv.lfrRx, ...
                dlrFpa.logical2doubleNan(), ...
                iCalibLZv, ...
                iCalibHZv);

            L.logf('info', ...
                ['Calibrating voltages -', ...
                ' One sequence of records with identical settings at a time.'])
            for iSs = 1:nSs
                iRec1 = iRec1Ar(iSs);
                iRec2 = iRec2Ar(iSs);

                SsAsrSamplesAVoltSrm = bicas.proc.L1L2.dc.calibrate_demux_subsequence(...
                    PreDc, dlrFpa, iCalibLZv, iCalibHZv, Cal, iRec1, iRec2, L);

                % Add demuxed sequence to the to-be complete set of records.
                AsrSamplesAVoltSrm.setRows(SsAsrSamplesAVoltSrm, [iRec1:iRec2]');
            end



            % NOTE: Assumes no "return" statement.
            %bicas.log_speed_profiling(L, 'bicas.proc.L1L2.dc.calibrate_demux_voltages', tTicToc, nRecords, 'record')
            %bicas.log_memory_profiling(L, 'bicas.proc.L1L2.dc.calibrate_demux_voltages:end')
        end    % calibrate_demux_voltages



        % Calibrate and demux all BLTS channels for one subsequence.
        function SsAsrSamplesAVoltSrm = calibrate_demux_subsequence(...
                PreDc, dlrFpa, iCalibLZv, iCalibHZv, Cal, iFirst, iLast, L)
            % PROPOSAL: Rename "subsequence".
            %   ~time interval
            %   ~constant settings time interval
            %
            % PROPOSAL: Move indexing outside function.
            %   CON: Must then submit all the PreDc.Zv.* variables separately
            %        (six variables). ==> More arguments.
            %       CON: Eliminates indexing arguments iFirst & iLast.
            %       CON-PROPOSAL: Use bicas.utils.SameRowsMap and reduce
            %                     the indexing to one operation (method call).
            %           PRO: Could include the additional non-PreDc.Zv
            %                variables in the same object.
            
            % IMPLEMENTATION NOTE: Function created to make loop in
            % calibrate_demux_voltages() smaller and more easy-to-understand.

            % Extract SCALAR settings to use for entire subsequence of
            % records.
            bdmFpa_ss                  = PreDc.Zv.bdmFpa(                 iFirst);
            bdm_ss                     = bdmFpa_ss.int2doubleNan();
            biasHighGainFpa_ss         = PreDc.Zv.biasHighGainFpa(        iFirst);    % Temporary variable. Should be eliminated eventually.
            biasHighGain_ss            = biasHighGainFpa_ss.logical2doubleNan();      % Convert FPA --> double-NaN
            clear biasHighGainFpa_ss
            freqHz_ss                  = PreDc.Zv.freqHz(                 iFirst);
            iLsf_ss                    = PreDc.Zv.iLsf(                   iFirst);
            CALIBRATION_TABLE_INDEX_ss = PreDc.Zv.CALIBRATION_TABLE_INDEX(iFirst, :);
            ufv_ss                     = PreDc.Zv.ufv(                    iFirst);
            dlrFpa_ss                  = dlrFpa(                          iFirst);
            iCalibL_ss                 = iCalibLZv(                       iFirst);
            iCalibH_ss                 = iCalibHZv(                       iFirst);
            % Extract subsequence of DATA records to "demux".
            ssBltsSamplesTm            = PreDc.Zv.bltsSamplesTm(          iFirst:iLast, :, :);
            % NOTE: "zVariable" (i.e. first index=record) for only the
            % current subsequence.
            ssZvNValidSamplesPerRecord = PreDc.Zv.nValidSamplesPerRecord( iFirst:iLast);

            if ~(PreDc.hasSnapshotFormat && PreDc.isLfr)
                % IMPLEMENTATION NOTE: Do not log for LFR SWF since it
                % produces unnecessarily many log messages since sampling
                % frequencies change for every CDF record.
                %
                % PROPOSAL: Make into "proper" table with top rows with column names.
                %   NOTE: Can not use irf.str.assist_print_table() since
                %         it requires the entire table to pre-exist before execution.
                %   PROPOSAL: Print after all iterations.
                %
                % NOTE: DIFF_GAIN needs three characters to fit in "NaN".
                L.logf('info', ['Records %8i-%8i : %s -- %s', ...
                    ' bdm/HK_BIA_MODE_MUX_SET=%i;', ...
                    ' biasHighGain/DIFF_GAIN=%-3i;', ...
                    ' dlr/HK_BIA_MODE_DIFF_PROBE=%i;', ...
                    ' freqHz=%5g; iCalibL=%i; iCalibH=%i; ufv=%i', ...
                    ' CALIBRATION_TABLE_INDEX=[%i, %i]'], ...
                    iFirst, iLast, ...
                    bicas.utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iFirst)), ...
                    bicas.utils.TT2000_to_UTC_str(PreDc.Zv.Epoch(iLast)), ...
                    bdm_ss, biasHighGain_ss, ...
                    dlrFpa_ss.logical2doubleNan(), ...
                    freqHz_ss, ...
                    iCalibL_ss, iCalibH_ss, ufv_ss, ...
                    CALIBRATION_TABLE_INDEX_ss(1), ...
                    CALIBRATION_TABLE_INDEX_ss(2))
            end

            %=======================================
            % DEMULTIPLEXER: FIND ASR-BLTS ROUTINGS
            %=======================================
            DemuxerRoutingArray = bicas.proc.L1L2.demuxer.get_routings(...
                bdmFpa_ss, dlrFpa_ss);



            if PreDc.hasSnapshotFormat
                % NOTE: Vector of constant numbers (one per snapshot).
                ssDtSec = 1 ./ PreDc.Zv.freqHz(iFirst:iLast);
            else
                % NOTE: Scalar (one for entire sequence).
                ssDtSec = double(...
                    PreDc.Zv.Epoch(iLast) - PreDc.Zv.Epoch(iFirst)) ...
                    / (iLast-iFirst) * 1e-9;   % TEMPORARY?
            end



            %=====================
            % ITERATE OVER BLTS's
            %=====================
            ssBltsSamplesAVolt = [];
            for iBlts = 1:bicas.const.N_BLTS
                ssBltsSamplesAVolt(:, :, iBlts) = bicas.proc.L1L2.dc.calibrate_BLTS(...
                    DemuxerRoutingArray(iBlts).ssid, ...
                    ssBltsSamplesTm(:, :, iBlts), ...
                    iBlts, ...
                    PreDc.hasSnapshotFormat, ...
                    ssZvNValidSamplesPerRecord, ...
                    biasHighGain_ss, ...         % Should use FPA (but does not).
                    iCalibL_ss, ...
                    iCalibH_ss, ...
                    iLsf_ss, ...
                    ssDtSec, ...
                    PreDc.isLfr, PreDc.isTdsCwf, ...
                    CALIBRATION_TABLE_INDEX_ss, ufv_ss, ...
                    Cal);
            end

            %====================================
            % DEMULTIPLEXER: DERIVE MISSING ASRs
            %====================================
            SsAsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_ASRs(...
                [DemuxerRoutingArray.sdid], ssBltsSamplesAVolt);
        end    % calibrate_demux_subsequence



        % Calibrate one BLTS channel.
        function samplesAVolt = calibrate_BLTS(...
                Ssid, samplesTm, iBlts, ...
                hasSnapshotFormat, ...
                zvNValidSamplesPerRecord, biasHighGain, ...
                iCalibL, iCalibH, iLsf, dtSec, ...
                isLfr, isTdsCwf, ...
                CALIBRATION_TABLE_INDEX, ufv, ...
                Cal)
            % IMPLEMENTATION NOTE: It is ugly to have this many parameters (15!),
            % but the original code made calibrate_demux_voltages() to large and
            % unwieldy. It also highlights the dependencies.
            %
            % PROPOSAL: CalSettings as parameter.
            %   PRO: Reduces number of parameters.
            %   PROPOSAL: Add values to CalSettings: isLfr, isTdsCwf, CALIBRATION_TABLE_INDEX
            %       CON: cal does not seem to use more values.
            % PROPOSAL: Reorder arguments to group them.
            %   PROPOSAL: Group arguments from PreDc.
            
            assert(isnumeric(samplesTm))
            irf.assert.sizes(samplesTm, [-1, -2])   % One BLTS channel.

            if isequaln(Ssid, bicas.proc.L1L2.SignalSourceId.C.UNKNOWN)
                % ==> Calibrated data set to NaN.
                samplesAVolt = nan(size(samplesTm));

            elseif isequaln(Ssid, bicas.proc.L1L2.SignalSourceId.C.GND) || ...
                    isequaln(Ssid, bicas.proc.L1L2.SignalSourceId.C.REF25V)
                % ==> No calibration.
                samplesAVolt = ssSamplesTm;

            else
                assert(Ssid.is_ASR())
                % ==> Calibrate (unless explicitly stated that should not)

                if hasSnapshotFormat
                    bltsSamplesTmCa = ...
                        bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(...
                            double(samplesTm), zvNValidSamplesPerRecord);
                else
                    assert(all(zvNValidSamplesPerRecord == 1))
                    bltsSamplesTmCa = {double(samplesTm)};
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
                CalSettings.Ssid         = Ssid;
                CalSettings.biasHighGain = biasHighGain;
                CalSettings.iCalibTimeL  = iCalibL;
                CalSettings.iCalibTimeH  = iCalibH;
                CalSettings.iLsf         = iLsf;
                %#######################################################
                ssBltsSamplesAVoltCa = Cal.calibrate_voltage_all(...
                    dtSec, bltsSamplesTmCa, ...
                    isLfr, isTdsCwf, CalSettings, ...
                    CALIBRATION_TABLE_INDEX, ufv);
                %#######################################################

                if hasSnapshotFormat
                    [samplesAVolt, ~] = ...
                        bicas.proc.utils.convert_cell_array_of_vectors_to_matrix(...
                            ssBltsSamplesAVoltCa, ...
                            size(samplesTm, 2));
                else
                    % Scalar cell array since not snapshot.
                    assert(isscalar(ssBltsSamplesAVoltCa))
                    % NOTE: Cell content must be column array.
                    samplesAVolt = ssBltsSamplesAVoltCa{1};
                end
            end
        end    % calibrate_BLTS



        function currentSAmpere = convert_CUR_to_CUR_on_SCI_TIME(...
                sciEpoch, InCur, SETTINGS, L)
            
            % PROPOSAL: Change function name. process_* implies converting struct-->struct.

            % ASSERTIONS
            assert(isa(InCur, 'bicas.InputDataset'))



            %===================================================================
            % CDF ASSERTION: CURRENT data begins before SCI data (i.e. there is
            % enough CURRENT data).
            %===================================================================
            if ~(min(InCur.Zv.Epoch) <= min(sciEpoch))
                curRelativeSec    = 1e-9 * (min(InCur.Zv.Epoch) - min(sciEpoch));
                sciEpochUtcStr    = irf.cdf.TT2000_to_UTC_str(min(sciEpoch));
                curEpochMinUtcStr = irf.cdf.TT2000_to_UTC_str(min(InCur.Zv.Epoch));

                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY');

                anomalyDescrMsg = sprintf(...
                    ['Bias current data begins %g s (%s) AFTER voltage data begins (%s).', ....
                    ' Can therefore not determine currents for all voltage timestamps.'], ...
                    curRelativeSec, curEpochMinUtcStr, sciEpochUtcStr);

                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:SWMProcessing')
            end



            %====================================================================
            % CDF ASSERTION: Epoch increases (not monotonically)
            % --------------------------------------------------
            % NOTE: bicas.proc.L1L2.dc.zv_TC_to_current() checks (and handles)
            % that Epoch increases monotonically, but only for each antenna
            % separately (which does not capture all cases). Therefore checks
            % that Epoch is (non-monotonically) increasing.
            % Ex: Timestamps, iAntenna = mod(iRecord,3): 1,2,3,5,4,6
            %       ==> Monotonically increasing sequences for each antenna
            %           separately, but not even increasing when combined.
            %====================================================================
            assert(issorted(InCur.Zv.Epoch), ...
                'BICAS:DatasetFormat', ...
                'CURRENT timestamps zVar Epoch does not increase (all antennas combined).')

            % NOTE: bicas.proc.L1L2.dc.zv_TC_to_current() checks that Epoch
            % increases monotonically.
            currentNanoSAmpere = [];
            currentNanoSAmpere(:,1) = bicas.proc.L1L2.dc.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_1, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,2) = bicas.proc.L1L2.dc.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_2, sciEpoch, L, SETTINGS);
            currentNanoSAmpere(:,3) = bicas.proc.L1L2.dc.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_3, sciEpoch, L, SETTINGS);

            currentSAmpere = 1e-9 * currentNanoSAmpere;
        end
        
        
        
        % Wrapper around solo.hwzv.CURRENT_ZV_to_current_interpolate for
        % anomaly handling.
        function sciZv_IBIASx = zv_TC_to_current(...
                curZv_Epoch, curZv_IBIAS_x, sciZv_Epoch, L, SETTINGS)



            %====================
            % Calibrate currents
            %====================
            [sciZv_IBIASx, duplicateAnomaly] = ...
                solo.hwzv.CURRENT_ZV_to_current_interpolate(...
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
                    'Bias current data contain duplicate settings, with', ...
                    ' identical timestamps', ...
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
                            'BICAS:SWMProcessing:DatasetFormat')
                end
            end

        end    % bicas.proc.L1L2.dc.zv_TC_to_current



    end    % methods(Static, Access=private)

end
