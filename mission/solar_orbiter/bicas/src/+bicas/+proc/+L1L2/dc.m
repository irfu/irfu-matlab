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
  %   process()
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

      bicas.proc.L1L2.dc.log_input_calibration_settings(Dcip, Cal, L)



      %#################################################################
      % Obtain "demultiplexer" "routings": SSID and SDID for every BLTS
      %#################################################################
      [bltsKSsidArray, bltsKSdidArray] = bicas.proc.L1L2.dc.get_KSSID_KSDID_arrays(...
        Dcip.Zv.bdmFpa, Dcip.Zv.dlrFpa);



      %############################
      % DEMUX & CALIBRATE VOLTAGES
      %############################
      AsrSamplesAVoltSrm = bicas.proc.L1L2.dc.calibrate_demux_voltages(...
        Dcip, bltsKSsidArray, bltsKSdidArray, Cal, L);



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
        bltsKSsidArray, ...
        Dcip.Zv.isAchgFpa, ...
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



    % Obtain kSSID and kSDID arrays for arrays of BDM and DLR.
    function [bltsKSsidArray, bltsKSdidArray] = get_KSSID_KSDID_arrays(bdmFpa, dlrFpa)
      nRecTot = irf.assert.sizes(...
        bdmFpa, [-1], ...
        dlrFpa, [-1]);

      [iRec1Array, iRec2Array, nSs] = irf.utils.split_by_change( ...
        bdmFpa.int2doubleNan(), dlrFpa.logical2doubleNan());

      % Preallocate
      % -----------
      % NOTE: No need for bicas.utils.FPArray since SSIDs and DSIDs handle all
      % special cases including unknown source and destination.
      bltsKSsidArray = zeros(nRecTot, bicas.const.N_BLTS, 'uint8');
      bltsKSdidArray = bltsKSsidArray;

      for iSs = 1:nSs
        iRecSs1 = iRec1Array(iSs);
        iRecSs  = iRec1Array(iSs):iRec2Array(iSs);
        nRecSs  = numel(iRecSs);

        DemuxerRoutingArray = bicas.proc.L1L2.demuxer.get_routings(...
          bdmFpa(iRecSs1), dlrFpa(iRecSs1));

        SsidArray  = [DemuxerRoutingArray.Ssid];
        kSsidArray = bicas.sconst.C.SSID_K_DICT(SsidArray);
        SdidArray  = [DemuxerRoutingArray.Sdid];
        kSdidArray = bicas.sconst.C.SDID_K_DICT(SdidArray);

        bltsKSsidArray(iRecSs, :) = repmat(kSsidArray, nRecSs, 1);
        bltsKSdidArray(iRecSs, :) = repmat(kSdidArray, nRecSs, 1);
      end
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function log_input_calibration_settings(Dcip, Cal, L)
      % IMPLEMENTATION NOTE: Implemented separately from processing functions
      % since:
      % (1) removes dependence on logger object,
      % (2) can split data into subsequences based on arbitrary choice of
      %     variables,
      % (3) reduces size of processing function where logging would otherwise
      %     be, and
      % (4) can potentially turn table into proper table with column headers.

      iCalibLZv = Cal.get_BIAS_calibration_time_L(Dcip.Zv.Epoch);
      iCalibHZv = Cal.get_BIAS_calibration_time_H(Dcip.Zv.Epoch);

      % IMPLEMENTATION NOTE: Do not log for LFR SWF since it produces
      % unnecessarily many log messages since sampling frequencies change
      % for every CDF record.
      if Dcip.hasSwfFormat && Dcip.isLfr
        return
      end

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

      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        Cv.bdmFpa    = Dcip.Zv.bdmFpa(   iRec1);
        Cv.isAchgFpa = Dcip.Zv.isAchgFpa(iRec1);
        Cv.dlrFpa    = Dcip.Zv.dlrFpa(   iRec1);

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
          bicas.utils.TT2000_to_UTC_str(Dcip.Zv.Epoch(iRec1), 9), ...
          bicas.utils.TT2000_to_UTC_str(Dcip.Zv.Epoch(iRec2), 9), ...
          Cv.bdmFpa.int2doubleNan(), ...
          Cv.isAchgFpa.logical2doubleNan(), ...
          Cv.dlrFpa.logical2doubleNan(), ...
          Dcip.Zv.freqHz(                  iRec1), ...
          iCalibLZv(                       iRec1), ...
          iCalibHZv(                       iRec1), ...
          Dcip.Zv.ufv(                     iRec1), ...
          Dcip.Zv.CALIBRATION_TABLE_INDEX( iRec1, 1), ...
          Dcip.Zv.CALIBRATION_TABLE_INDEX( iRec1, 2))
      end    % for

    end



    % Demultiplex and calibrate VOLTAGES (not e.g. currents).
    %
    % NOTE: Can handle arrays of any size if the sizes are consistent.
    %
    function AsrSamplesAVoltSrm = calibrate_demux_voltages(...
        Dcip, bltsKSsidArray, bltsKSdidArray, Cal, L)
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
      %
      % PROPOSAL: Replace Rcip with name-value arguments.
      % arguments
      %   A.hasSwfFormat
      %   Zv.Epoch
      %   Zv.bltsSamplesTm
      %   A.bltsKSsidArray
      %   Zv.isAchgFpa
      %   CALIBRATION_TABLE_INDEX
      %   ufv
      %   nValidSamplesPerRecord
      % end

      % ASSERTIONS
      assert(isscalar( Dcip.hasSwfFormat))
      assert(isnumeric(Dcip.Zv.bltsSamplesTm))
      assert(isa(bltsKSsidArray, 'uint8'))
      assert(isa(bltsKSdidArray, 'uint8'))
      [nRecords, nSamplesPerRecordChannel] = irf.assert.sizes(...
        Dcip.Zv.isAchgFpa,     [-1,  1], ...
        bltsKSsidArray,        [-1,  bicas.const.N_BLTS], ...
        bltsKSsidArray,        [-1,  bicas.const.N_BLTS], ...
        Dcip.Zv.bltsSamplesTm, [-1, -2, bicas.const.N_BLTS]);



      % Pre-allocate AsrSamplesAVoltSrm: All (ASID) channels, all records
      % -----------------------------------------------------------------
      % IMPLEMENTATION NOTE: Preallocation is very important for speeding
      % up LFR-SWF which tends to be broken into subsequences of 1 record.
      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap(...
        "bicas.proc.L1L2.AntennaSignalId", nRecords, 'CONSTANT', ...
        nan(nRecords, nSamplesPerRecordChannel), ...
        bicas.sconst.C.S_ASID_DICT.values);

      iCalibLZv = Cal.get_BIAS_calibration_time_L(Dcip.Zv.Epoch);
      iCalibHZv = Cal.get_BIAS_calibration_time_H(Dcip.Zv.Epoch);



      %========================================================================
      % (1) Find continuous subsequences of records with identical settings.
      % (2) Process data separately for each such sequence.
      % --------------------------------------------------------------------
      % NOTE: Just finding continuous subsequences can take a significant
      %       amount of time.
      % NOTE: Empirically, this is not useful for real LFR SWF datasets where
      %       the LFR sampling frequency changes in every record, meaning that
      %       the subsequences are all 1 record long.
      %
      % SS = Subsequence
      %========================================================================
      % PROPOSAL: Do not use irf.utils.split_by_change() for SWF data. It is
      %           enough to group by identical values (not use continuous
      %           blocks of CDF records).
      %   PRO: Faster
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        Dcip.Zv.isAchgFpa.logical2doubleNan(), ...
        Dcip.Zv.freqHz, ...
        Dcip.Zv.iLsf, ...
        Dcip.Zv.CALIBRATION_TABLE_INDEX, ...
        Dcip.Zv.ufv, ...
        Dcip.Zv.lrx, ...
        bltsKSsidArray, ...
        bltsKSsidArray, ...
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
        SsAsrSamplesAVoltSrm = bicas.proc.L1L2.dc.calibrate_demux_voltages_subsequence( ...
          Cal, ...
          ... % ===============================================================
          ... % NOTE: Variables which do VARY over CDF records.
          isAchgFpa      = Dcip.Zv.isAchgFpa(              iRec1), ...
          freqHz         = Dcip.Zv.freqHz(                 iRec1), ...
          iLsf           = Dcip.Zv.iLsf(                   iRec1), ...
          zvcti          = Dcip.Zv.CALIBRATION_TABLE_INDEX(iRec1, :), ...
          ufv            = Dcip.Zv.ufv(                    iRec1), ...
          bltsKSsidArray = bltsKSsidArray(                 iRec1, :), ...
          bltsKSdidArray = bltsKSdidArray(                 iRec1, :), ...
          ... % NOTE: Excluding Dcip.Zv.lrx since it is only need for
          ... %       splitting time/CDF record intervals, not for calibration
          ... %       since calibration can handle sequences of only NaN.
          iCalibL        = iCalibLZv(iRec1), ...
          iCalibH        = iCalibHZv(iRec1), ...
          ... % ===============================================================
          ... % NOTE: Variables which do not vary over CDF records.
          hasSwfFormat   = Dcip.hasSwfFormat, ...
          isLfr          = Dcip.isLfr, ...
          isTdsCwf       = Dcip.isTdsCwf, ...
          ...   % Variables which vary by CDF records.
          Epoch                    = Dcip.Zv.Epoch(                 iRec1:iRec2), ...
          bltsSamplesTm            = Dcip.Zv.bltsSamplesTm(         iRec1:iRec2, :, :), ...
          zvNValidSamplesPerRecord = Dcip.Zv.nValidSamplesPerRecord(iRec1:iRec2));

        % Add demuxed sequence to the to-be complete set of records.
        AsrSamplesAVoltSrm.set_rows(SsAsrSamplesAVoltSrm, [iRec1:iRec2]');
      end    % for

    end    % calibrate_demux_voltages



    % Calibrate and demux all BLTS channels for one subsequence with various
    % constant settings/values.
    %
    % ARGUMENTS
    % =========
    % Cv
    %       Constant values. Scalar values which do NOT VARY by CDF record.
    % Vv
    %       Varying values. Struct with values which DO VARY by CDF record.
    function AsrSamplesAVoltSrm = calibrate_demux_voltages_subsequence(Cal, Cv, Vv)
      arguments
        Cal
        %
        % NOTE: Excluding LRX since it is only need for splitting time/CDF
        %       record intervals, not for calibration since calibration can
        %       handle sequences of only NaN.
        Cv.isAchgFpa
        Cv.freqHz
        Cv.iLsf
        Cv.zvcti
        Cv.ufv
        Cv.bltsKSsidArray
        Cv.bltsKSdidArray
        Cv.iCalibL
        Cv.iCalibH
        % NOTE: Below variables do not vary over CDF records anyhow.
        Cv.hasSwfFormat
        Cv.isLfr
        Cv.isTdsCwf

        Vv.Epoch
        Vv.bltsSamplesTm
        Vv.zvNValidSamplesPerRecord
      end

      nRows = numel(Vv.Epoch);

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
          Ssid                    =bicas.sconst.C.K_SSID_DICT(Cv.bltsKSsidArray(iBlts)), ...
          samplesTm               =Vv.bltsSamplesTm(:, :, iBlts), ...
          iBlts                   =iBlts, ...
          hasSwfFormat            =Cv.hasSwfFormat, ...
          zvNValidSamplesPerRecord=Vv.zvNValidSamplesPerRecord, ...
          isAchg                  =Cv.isAchgFpa.logical2doubleNan(), ...
          iCalibL                 =Cv.iCalibL, ...
          iCalibH                 =Cv.iCalibH, ...
          iLsf                    =Cv.iLsf, ...
          dtSec                   =dtSec, ...
          isLfr                   =Cv.isLfr, ...
          isTdsCwf                =Cv.isTdsCwf, ...
          zvcti                   =Cv.zvcti, ...
          ufv                     =Cv.ufv, ...
          Cal                     =Cal);
      end

      %========================================================
      % DEMULTIPLEXER: DERIVE AND COMPLEMENT WITH MISSING ASRs
      %========================================================
      AsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_all_ASRs(...
        bicas.sconst.C.K_SDID_DICT(Cv.bltsKSdidArray), ...
        ssBltsSamplesAVolt);
    end    % calibrate_demux_voltages_subsequence



    % Calibrate one BLTS channel.
    function samplesAVolt = calibrate_BLTS(A)
      arguments
        A.Ssid
        A.samplesTm
        A.iBlts
        A.hasSwfFormat
        A.zvNValidSamplesPerRecord
        A.isAchg
        A.iCalibL
        A.iCalibH
        A.iLsf
        A.dtSec
        A.isLfr
        A.isTdsCwf
        A.zvcti
        A.ufv
        A.Cal
      end
      % IMPLEMENTATION NOTE: It is ugly to have this many parameters (15!),
      % but the original code made calibrate_demux_voltages() to large and
      % unwieldy. Having many arguments also highlights the exact dependencies.
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
      % eliminated. Can use integers or bicas.utils.FPArray?
      % NOTE: Storing TM units with floats!
      if A.isLfr
        assert(isa(A.samplesTm, 'single'))
      else
        assert(isa(A.samplesTm, 'double'))
      end
      irf.assert.sizes(A.samplesTm, [-1, -2])   % One BLTS channel.

      if isequaln(A.Ssid, bicas.sconst.C.S_SSID_DICT("UNKNOWN"))
        % ==> Calibrated data set to NaN.
        samplesAVolt = nan(size(A.samplesTm));

      elseif isequaln(A.Ssid, bicas.sconst.C.S_SSID_DICT("GND")) || ...
          isequaln(A.Ssid, bicas.sconst.C.S_SSID_DICT("REF25V"))
        % ==> No calibration.
        % NOTE: samplesTm stores TM units using float!
        samplesAVolt = A.samplesTm;

      else
        assert(A.Ssid.is_ASR())
        % ==> Calibrate (unless explicitly stated that should not)

        if A.hasSwfFormat
          bltsSamplesTmCa = ...
            bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(...
            double(A.samplesTm), A.zvNValidSamplesPerRecord);
        else
          assert(all(A.zvNValidSamplesPerRecord == 1))
          bltsSamplesTmCa = {double(A.samplesTm)};
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
          A.iBlts, A.Ssid, A.isAchg, A.iCalibL, A.iCalibH, A.iLsf);
        %#######################################################
        ssBltsSamplesAVoltCa = A.Cal.calibrate_voltage_all(...
          A.dtSec, bltsSamplesTmCa, ...
          A.isLfr, A.isTdsCwf, CalSettings, ...
          A.zvcti, A.ufv);
        %#######################################################

        if A.hasSwfFormat
          [samplesAVolt, ~] = ...
            bicas.proc.utils.convert_cell_array_of_vectors_to_matrix(...
            ssBltsSamplesAVoltCa, ...
            size(A.samplesTm, 2));
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
