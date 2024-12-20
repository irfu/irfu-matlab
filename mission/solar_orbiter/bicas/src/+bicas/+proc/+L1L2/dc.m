%
% Collection of processing functions for demultiplexing and calibrating, and
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
  %
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



    % Derive DCOP from DCIP, i.e.
    % (1) demux (demultiplex),
    % (2) calibrate data, and
    % (3) set quality variables.
    %
    % NOTE: Public function as opposed to the other demuxing/calibration
    % functions.
    %
    function Dcop = process_calibrate_demux(Dcip, InCurPd, Cal, NsoTable, Bso, L)

      Tmk = bicas.utils.Timekeeper('bicas.proc.L1L2.dc.process_calibrate_demux', L);

      % ASSERTION
      assert(isa(Dcip, 'bicas.proc.L1L2.DemultiplexingCalibrationInput'));



      %############################
      % DEMUX & CALIBRATE VOLTAGES
      %############################
      AsrSamplesAVoltSrm = ...
        bicas.proc.L1L2.dc.calibrate_demux_voltages(Dcip, Cal, L);



      %#########################
      % Calibrate bias CURRENTS
      %#########################
      currentSAmpere = bicas.proc.L1L2.dc.convert_CUR_to_CUR_on_SCI_TIME(...
        Dcip.Zv.Epoch, InCurPd, Bso, L);
      currentTm      = bicas.proc.L1L2.cal.Cal.calibrate_current_sampere_to_TM(currentSAmpere);

      currentAAmpere          = nan(size(currentSAmpere));    % Preallocate.
      iCalibLZv               = Cal.get_BIAS_calibration_time_L(Dcip.Zv.Epoch);
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(iCalibLZv);
      L.logf('info', ...
        ['Calibrating currents -', ...
        ' One sequence of records with identical settings at a time.'])
      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        iRecords = iRec1:iRec2;

        L.logf('info', 'Records %8i-%8i : %s -- %s', ...
          iRec1, iRec2, ...
          bicas.utils.TT2000_to_UTC_str(Dcip.Zv.Epoch(iRec1), 9), ...
          bicas.utils.TT2000_to_UTC_str(Dcip.Zv.Epoch(iRec2), 9))

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
      % AUTODETECT SATURATION.
      Sat = bicas.proc.L1L2.Saturation(Bso);
      isAutodetectedSaturation = Sat.get_voltage_saturation_quality_bit(...
        Dcip.Zv.Epoch, ...
        AsrSamplesAVoltSrm, ...
        Dcip.Zv.nValidSamplesPerRecord, ...
        Dcip.Zv.bdmFpa, Dcip.Zv.dlrFpa, ...
        Dcip.Zv.lrx,    Dcip.Zv.isAchgFpa, ...
        Dcip.hasSwfFormat, L);

      ZvIn = struct(...
        'Epoch',            Dcip.Zv.Epoch, ...
        'bdmFpa',           Dcip.Zv.bdmFpa, ...
        'isFullSaturation', isAutodetectedSaturation);
      [zvUfv, QUALITY_FLAG, Zv.L2_QUALITY_BITMASK] = ...
        bicas.proc.L1L2.qual.get_UFV_quality_ZVs(...
        ZvIn, Dcip.isLfr, NsoTable, Bso, L);
      clear ZvIn

      Zv.QUALITY_FLAG = Dcip.Zv.QUALITY_FLAG.min(bicas.utils.FPArray(QUALITY_FLAG));
      zvUfv           = Dcip.Zv.ufv | zvUfv;

      % NOTE: Function modifies AsrSamplesAVoltSrm handle object!
      Zv.currentAAmpere = bicas.proc.L1L2.qual.set_voltage_current_FV(...
        Dcip.Zv.Epoch, AsrSamplesAVoltSrm, currentAAmpere, zvUfv, L);
      Zv.AsrSamplesAVoltSrm = AsrSamplesAVoltSrm;
      clear AsrSamplesAVoltSrm



      % ############
      % END FUNCTION
      % ############
      Dcop = bicas.proc.L1L2.DemultiplexingCalibrationOutput(Zv);

      nRecords = size(Dcip.Zv.Epoch, 1);
      Tmk.stop_log(nRecords, 'record')
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
    % NOTE: Can handle arrays of any size if the sizes are consistent.
    %
    function AsrSamplesAVoltSrm = calibrate_demux_voltages(Dcip, Cal, L)
      % PROPOSAL: Sequence of constant settings includes constant NaN/non-NaN
      %           for CWF.
      %
      % PROPOSAL: Integrate into bicas.proc.L1L2.demuxer (as method).
      % NOTE: Calibration is really separate from the demultiplexer.
      %       Demultiplexer only needs to split into subsequences based on BDM
      %       and DLR, nothing else.
      %   PROPOSAL: Separate out demultiplexer. Do not call from this function.
      %
      % PROPOSAL: Move the different conversion of CWF/SWF (one/many cell
      %           arrays) into the calibration function?!!

      % ASSERTIONS
      assert(isscalar( Dcip.hasSwfFormat))
      assert(isnumeric(Dcip.Zv.bltsSamplesTm))
      [nRecords, nSamplesPerRecordChannel] = irf.assert.sizes(...
        Dcip.Zv.bdmFpa,        [-1,  1], ...
        Dcip.Zv.isAchgFpa,     [-1,  1], ...
        Dcip.Zv.bltsSamplesTm, [-1, -2, bicas.const.N_BLTS]);



      % Pre-allocate AsrSamplesAVoltSrm: All (ASID) channels, all records
      % -----------------------------------------------------------------
      % IMPLEMENTATION NOTE: Preallocation is very important for speeding
      % up LFR-SWF which tends to be broken into subsequences of 1 record.
      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap(...
        "bicas.proc.L1L2.AntennaSignalId", nRecords, 'CONSTANT', ...
        nan(nRecords, nSamplesPerRecordChannel), ...
        bicas.proc.L1L2.AntennaSignalId.ALL_ARRAY);

      iCalibLZv = Cal.get_BIAS_calibration_time_L(Dcip.Zv.Epoch);
      iCalibHZv = Cal.get_BIAS_calibration_time_H(Dcip.Zv.Epoch);



      %======================================================================
      % (1) Find continuous subsequences of records with identical settings.
      % (2) Process data separately for each such sequence.
      % ----------------------------------------------------------
      % NOTE: Just finding continuous subsequences can take a significant
      %       amount of time.
      % NOTE: Empirically, this is not useful for real LFR SWF datasets where
      %       the LFR sampling frequency changes in every record, meaning that
      %       the subsequences are all 1 record long.
      %
      % SS = Subsequence
      %======================================================================
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        Dcip.Zv.bdmFpa.int2doubleNan(), ...
        Dcip.Zv.isAchgFpa.logical2doubleNan(), ...
        Dcip.Zv.freqHz, ...
        Dcip.Zv.iLsf, ...
        Dcip.Zv.CALIBRATION_TABLE_INDEX, ...
        Dcip.Zv.ufv, ...
        Dcip.Zv.dlrFpa.logical2doubleNan(), ...
        Dcip.Zv.lrx, ...
        iCalibLZv, ...
        iCalibHZv);

      L.logf('info', ...
        ['Calibrating voltages -', ...
        ' One sequence of records with identical settings at a time.'])
      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        % ==============================================================
        % IMPLEMENTATION NOTE: Below extraction of data from Dcip etc.
        % may seem awkward but actually clarifies the code associated
        % with bicas.proc.L1L2.dc.calibrate_demux_voltages_subsequence()
        % as compared to earlier version before refactoring.
        %
        % PRO: Clearly divides the variables/arguments into (a) constant
        %      and (2) varying variables.
        % PRO: Clarifies the information which the function needs.
        %      = Minimizes the amount of information that goes into the
        %      function.
        % PRO: Prevents the function from having to do the same.
        % PRO: Prevents the function from simultaneously having
        %      variables version for (a) entire interval of time and (b)
        %      the selected interval of time.
        % ==============================================================

        % CV = Constant Values = Values which are constant for the
        %      entire subsequence of records.
        Cv = [];
        Cv.bdmFpa    = Dcip.Zv.bdmFpa(                 iRec1);
        Cv.isAchgFpa = Dcip.Zv.isAchgFpa(              iRec1);
        Cv.freqHz    = Dcip.Zv.freqHz(                 iRec1);
        Cv.iLsf      = Dcip.Zv.iLsf(                   iRec1);
        Cv.zvcti     = Dcip.Zv.CALIBRATION_TABLE_INDEX(iRec1, :);
        Cv.ufv       = Dcip.Zv.ufv(                    iRec1);
        Cv.dlrFpa    = Dcip.Zv.dlrFpa(                 iRec1);
        % NOTE: Excluding Dcip.Zv.lrx since it is only need for
        %       splitting time/CDF record intervals, not for calibration
        %       since calibration can handle sequences of only NaN.
        Cv.iCalibL                 = iCalibLZv(                       iRec1);
        Cv.iCalibH                 = iCalibHZv(                       iRec1);
        % NOTE: Below variables do not vary over CDF records anyhow.
        Cv.hasSwfFormat            = Dcip.hasSwfFormat;
        Cv.isLfr                   = Dcip.isLfr;
        Cv.isTdsCwf                = Dcip.isTdsCwf;

        % VV = (Record-)Varying Values
        Vv = [];
        Vv.Epoch                    = Dcip.Zv.Epoch(                  iRec1:iRec2);
        Vv.bltsSamplesTm            = Dcip.Zv.bltsSamplesTm(          iRec1:iRec2, :, :);
        Vv.zvNValidSamplesPerRecord = Dcip.Zv.nValidSamplesPerRecord( iRec1:iRec2);

        if ~(Cv.hasSwfFormat && Cv.isLfr)
          % IMPLEMENTATION NOTE: Do not log for LFR SWF since it produces
          % unnecessarily many log messages since sampling frequencies change
          % for every CDF record.
          %
          % PROPOSAL: Make into "proper" table with top rows with column names.
          %   NOTE: Can not use irf.str.assist_print_table() since
          %         it requires the entire table to pre-exist before execution.
          %   PROPOSAL: Print after all iterations.
          %
          % NOTE: DIFF_GAIN needs three characters to print the string "NaN".
          L.logf('info', ['Records %8i-%8i : %s -- %s', ...
            ' bdm/HK_BIA_MODE_MUX_SET=%i;', ...
            ' isAchg/DIFF_GAIN=%-3i;', ...
            ' dlr/HK_BIA_MODE_DIFF_PROBE=%i;', ...
            ' freqHz=%5g; iCalibL=%i; iCalibH=%i; ufv=%i', ...
            ' CALIBRATION_TABLE_INDEX=[%i, %i]'], ...
            iRec1, iRec2, ...
            bicas.utils.TT2000_to_UTC_str(Vv.Epoch(1),   9), ...
            bicas.utils.TT2000_to_UTC_str(Vv.Epoch(end), 9), ...
            Cv.bdmFpa.int2doubleNan(), ...
            Cv.isAchgFpa.logical2doubleNan(), ...
            Cv.dlrFpa.logical2doubleNan(), ...
            Cv.freqHz, ...
            Cv.iCalibL, Cv.iCalibH, Cv.ufv, ...
            Cv.zvcti(1), ...
            Cv.zvcti(2))
        end

        SsAsrSamplesAVoltSrm = bicas.proc.L1L2.dc.calibrate_demux_voltages_subsequence(...
          Cv, Vv, Cal);

        % Add demuxed sequence to the to-be complete set of records.
        AsrSamplesAVoltSrm.set_rows(SsAsrSamplesAVoltSrm, [iRec1:iRec2]');
      end    % for

    end    % calibrate_demux_voltages



    % Calibrate and demux all BLTS channels for one subsequence with various
    % constant settings/values.
    function AsrSamplesAVoltSrm = calibrate_demux_voltages_subsequence(Cv, Vv, Cal)
      % PROPOSAL: Rename "subsequence".
      %   ~time interval
      %   ~constant settings time interval

      nRows = numel(Vv.Epoch);

      %=======================================
      % DEMULTIPLEXER: FIND ASR-BLTS ROUTINGS
      %=======================================
      DemuxerRoutingArray = bicas.proc.L1L2.demuxer.get_routings(...
        Cv.bdmFpa, Cv.dlrFpa);

      if Cv.hasSwfFormat
        % NOTE: Vector of constant numbers (one per snapshot).
        dtSec = ones(nRows, 1) / Cv.freqHz;
      else
        % NOTE: Scalar (one for entire sequence).
        dtSec = double( Vv.Epoch(end) - Vv.Epoch(1) ) / (nRows-1) * 1e-9;
      end

      %====================
      % CALIBRATE VOLTAGES
      %====================
      ssBltsSamplesAVolt = [];
      for iBlts = 1:bicas.const.N_BLTS
        ssBltsSamplesAVolt(:, :, iBlts) = bicas.proc.L1L2.dc.calibrate_BLTS(...
          DemuxerRoutingArray(iBlts).Ssid, ...
          Vv.bltsSamplesTm(:, :, iBlts), ...
          iBlts, ...
          Cv.hasSwfFormat, ...
          Vv.zvNValidSamplesPerRecord, ...
          Cv.isAchgFpa.logical2doubleNan(), ...
          Cv.iCalibL, ...
          Cv.iCalibH, ...
          Cv.iLsf, ...
          dtSec, ...
          Cv.isLfr, Cv.isTdsCwf, ...
          Cv.zvcti, Cv.ufv, ...
          Cal);
      end

      %========================================================
      % DEMULTIPLEXER: DERIVE AND COMPLEMENT WITH MISSING ASRs
      %========================================================
      AsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_all_ASRs(...
        [DemuxerRoutingArray.Sdid], ssBltsSamplesAVolt);
    end    % calibrate_demux_voltages_subsequence



    % Calibrate one BLTS channel.
    function samplesAVolt = calibrate_BLTS(...
        Ssid, samplesTm, iBlts, ...
        hasSwfFormat, ...
        zvNValidSamplesPerRecord, isAchg, ...
        iCalibL, iCalibH, iLsf, dtSec, ...
        isLfr, isTdsCwf, ...
        zvcti, ufv, ...
        Cal)
      % IMPLEMENTATION NOTE: It is ugly to have this many parameters (15!),
      % but the original code made calibrate_demux_voltages() to large and
      % unwieldy. It also highlights the dependencies.
      %
      % PROPOSAL: CalSettings as parameter.
      %   PRO: Reduces number of parameters.
      %   PROPOSAL: Add values to CalSettings: isLfr, isTdsCwf, zvcti
      %       CON: cal does not seem to use more values.
      % PROPOSAL: Reorder arguments to group them.
      %   PROPOSAL: Group arguments from DCIP.

      % IMPLEMENTATION NOTE: It seems that data processing submits
      % different types of floats for LFR and TDS. This difference in
      % processing is unintended and should probably ideally be
      % eliminated.
      % NOTE: Storing TM units with floats!
      if isLfr
        assert(isa(samplesTm, 'single'))
      else
        assert(isa(samplesTm, 'double'))
      end
      irf.assert.sizes(samplesTm, [-1, -2])   % One BLTS channel.

      if isequaln(Ssid, bicas.proc.L1L2.SignalSourceId.C.UNKNOWN)
        % ==> Calibrated data set to NaN.
        samplesAVolt = nan(size(samplesTm));

      elseif isequaln(Ssid, bicas.proc.L1L2.SignalSourceId.C.GND) || ...
          isequaln(Ssid, bicas.proc.L1L2.SignalSourceId.C.REF25V)
        % ==> No calibration.
        % NOTE: samplesTm stores TM units using float!
        samplesAVolt = samplesTm;

      else
        assert(Ssid.is_ASR())
        % ==> Calibrate (unless explicitly stated that should not)

        if hasSwfFormat
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
        % IMPLEMENTATION NOTE: LFR zVar BW=0
        % ==> CALIBRATION_TABLE_INDEX(1,:) illegal value.
        % ==> Can not calibrate.
        % ==> Must explicitly disable calibration.
        % Therefore uses ufv_ss to disable calibration.
        % It is thus not enough to overwrite the values later.
        % This incidentally also potentially speeds up the code.
        % Ex: LFR SWF 2020-02-25, 2020-02-28.
        CalSettings = bicas.proc.L1L2.CalibrationSettings(...
          iBlts, Ssid, isAchg, iCalibL, iCalibH, iLsf);
        %#######################################################
        ssBltsSamplesAVoltCa = Cal.calibrate_voltage_all(...
          dtSec, bltsSamplesTmCa, ...
          isLfr, isTdsCwf, CalSettings, ...
          zvcti, ufv);
        %#######################################################

        if hasSwfFormat
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
        sciEpoch, InCur, Bso, L)

      % PROPOSAL: Change function name. process_* implies converting struct-->struct.

      % ASSERTIONS
      assert(isa(InCur, 'bicas.InputDataset'))



      %===================================================================
      % CDF ASSERTION: CURRENT data begins before SCI data (i.e. there is
      % enough CURRENT data).
      %===================================================================
      if ~(min(InCur.Zv.Epoch) <= min(sciEpoch))
        curRelativeSec    = 1e-9 * (min(InCur.Zv.Epoch) - min(sciEpoch));
        sciEpochUtcStr    = bicas.utils.TT2000_to_UTC_str(min(sciEpoch),       9);
        curEpochMinUtcStr = bicas.utils.TT2000_to_UTC_str(min(InCur.Zv.Epoch), 9);

        [settingValue, settingKey] = Bso.get_fv(...
          'PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY');

        anomalyDescrMsg = sprintf(...
          ['Bias current data begins %g s (%s) AFTER voltage data begins (%s).', ....
          ' Can therefore not determine currents for all voltage timestamps.'], ...
          curRelativeSec, curEpochMinUtcStr, sciEpochUtcStr);

        bicas.default_anomaly_handling(L, settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
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
      currentNanoSAmpere(:,1) = bicas.proc.L1L2.dc.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_1, sciEpoch, L, Bso);
      currentNanoSAmpere(:,2) = bicas.proc.L1L2.dc.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_2, sciEpoch, L, Bso);
      currentNanoSAmpere(:,3) = bicas.proc.L1L2.dc.zv_TC_to_current(InCur.Zv.Epoch, InCur.Zv.IBIAS_3, sciEpoch, L, Bso);

      currentSAmpere = 1e-9 * currentNanoSAmpere;
    end



    % Wrapper around solo.hwzv.CURRENT_ZV_to_current_interpolate for
    % anomaly handling.
    function sciZv_IBIASx = zv_TC_to_current(...
        curZv_Epoch, curZv_IBIAS_x, sciZv_Epoch, L, Bso)



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
        [settingValue, settingKey] = Bso.get_fv(...
          'INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY');
        anomalyDescriptionMsg = [...
          'Bias current data contain duplicate settings, with', ...
          ' identical timestamps', ...
          ' and identical bias settings on the same antenna.'];

        switch(settingValue)
          case 'REMOVE_DUPLICATES'
            bicas.default_anomaly_handling(L, ...
              settingValue, settingKey, 'OTHER', ...
              anomalyDescriptionMsg)
            L.log('warning', ...
              ['Removed duplicated bias current settings with', ...
              ' identical timestamps on the same antenna.'])

          otherwise
            bicas.default_anomaly_handling(L, ...
              settingValue, settingKey, 'ERROR_ILLEGAL_SETTING', ...
              anomalyDescriptionMsg, ...
              'BICAS:SWMProcessing:DatasetFormat')
        end
      end

    end    % bicas.proc.L1L2.dc.zv_TC_to_current



  end    % methods(Static, Access=private)

end
