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
  %   calibration, demuxing, reconstruction, quality
  %   cdr = calibration+demuxing+reconstruction (same order as execution)
  %   cdrq = calibration+demuxing+reconstruction+quality
  %     CON: Setting quality variables is not necessarily last in the execution.
  %      Ex: Blanking data due to failed antenna should be done in TM channels
  %          (planned implementation).
  %       CON: Minor. Quality variables are still set somewhere in this file.
  %
  % PROPOSAL: Automatic test code.
  %
  % PROPOSAL:   process_calibrate_demux()
  %           & calibrate_voltages_5xBLTS()
  %           should only accept the needed ZVs and variables.
  %   NOTE: Needs some way of packaging/extracting only the relevant ZVs/fields
  %         from struct.
  %
  % PROPOSAL: Reorg. code to
  %   * Consist of more isolated/modular/generic separate steps.
  %   * Be more easily testable.
  %   * Be easier to understand.
  %   * Not require splitting up in subsequences with constant "settings" (values
  %     for specific zVariables not varying as a function of CDF record), in
  %     particular not require many constant "settings".
  %   * Use more vector operations.
  %   * More natural to implement.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Derive DCOP from DCIP:
    % * Calibrate bias currents.
    % * Voltages:
    %   (1) calibrate samples
    %   (2) demux (demultiplex): Relabel samples from BLTS to SDID.
    %   (3) reconstruct (derive) samples for missing channels from calibrated
    %       samples (e.g. DC_V12 := DC_V1 - DC_V2)
    % * Set quality variables.
    %
    function Dcop = process_calibrate_demux(Dcip, InCurPd, Cal, NsoTable, Bso, L)

      Tmk = bicas.utils.Timekeeper('bicas.proc.L1L2.dc.process_calibrate_demux', L);

      % ASSERTION
      assert(isa(Dcip, 'bicas.proc.L1L2.DemultiplexingCalibrationInput'));

      bicas.proc.L1L2.dc.log_input_calibration_settings(Dcip, Cal, L)



      %#########################
      % Calibrate bias CURRENTS
      %#########################
      currentAAmpere = bicas.proc.L1L2.cur.calibrate_bias_currents(...
        Dcip.Zv.Epoch, InCurPd, Cal, Bso, L);



      %######################################################################
      % Obtain "demultiplexer" "routings" in the form of SSID-SDID pairs for
      % every BLTS (and CDF record)
      %######################################################################
      [bltsSsidArray, bltsSdidArray] = bicas.proc.L1L2.dc.get_SSID_SDID_arrays(...
        Dcip.Zv.bdmFpa, ...
        Dcip.Zv.dlrFpa);



      %##################################################################
      % CALIBRATE VOLTAGES: 5x CHANNELS LABELLED BY SSID/BLTS (not SDID)
      %##################################################################
      % NOTE: Takes most of the time LFR-SWF.
      bltsSamplesAVolt = bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS(...
        Epoch                   = Dcip.Zv.Epoch, ...
        bltsSamplesTm           = Dcip.Zv.bltsSamplesTm, ...
        isAchgFpa               = Dcip.Zv.isAchgFpa, ...
        CALIBRATION_TABLE_INDEX = Dcip.Zv.CALIBRATION_TABLE_INDEX, ...
        freqHz                  = Dcip.Zv.freqHz, ...
        iLsf                    = Dcip.Zv.iLsf, ...
        ufv                     = Dcip.Zv.ufv, ...
        bltsSsidArray           = bltsSsidArray, ...
        isTdsCwf                = Dcip.isTdsCwf, ...
        isLfr                   = Dcip.isLfr, ...
        hasSwfFormat            = Dcip.hasSwfFormat, ...
        uspr                    = Dcip.Zv.uspr, ...
        Cal                     = Cal, ...
        L                       = L);



      %#####################################
      % Get VSIB for BLTS-labelled channels
      %#####################################
      bltsVsibAr = bicas.proc.L1L2.dc.get_VSIB_5xBLTS_NEW(...
        bltsSamplesAVolt, Dcip.Zv.uspr, ...
        bltsSsidArray, Dcip.Zv.isAchgFpa, Bso, L);



      %######################################################################
      % ~"DEMUX" VOLTAGES:
      % INPUT:  5x SIGNALS LABELLED BY SSID/BLTS
      % OUTPUT: 9x SIGNALS LABELLED BY SDID + RECONSTRUCTING MISSING SIGNALS
      %######################################################################
      % NOTE: Needs VSIB for propagating VSIB to reconstructed channels.
      SchdZvm = bicas.proc.L1L2.dc.convert_samples_5xBLTS_to_9xASR_NEW(...
        bltsSamplesAVolt, bltsVsibAr, bltsSdidArray, L);
      bicas.proc.L1L2.demuxer.reconstruct_ASR_samples_NEW(SchdZvm);



      % Convert SchdZvm --> SamplesZvm + VsibZvm
      SamplesZvm = bicas.utils.ZvMap(SchdZvm.nRecords);
      VsibZvm    = bicas.utils.ZvMap(SchdZvm.nRecords);
      for keyCa = SchdZvm.keyCa'
        SamplesZvm.add(keyCa{1}, SchdZvm.get(keyCa{1}).samplesAr);
        VsibZvm.add(   keyCa{1}, SchdZvm.get(keyCa{1}).vsibAr);
      end
      clear SchdZvm



      %####################
      % Obtain quality ZVs
      %####################
      SatSettings = bicas.proc.L1L2.sat.from_BSO_extract_saturation_settings(Bso);
      [QUALITY_FLAG, L2_QUALITY_BITMASK] = ...
        bicas.proc.L1L2.qual.get_quality_ZVs(...
        Dcip.Zv.Epoch, NsoTable, ...
        string(Bso.get_fv('PROCESSING.SATURATION.QUALITY_SCHEME')), ...
        VsibZvm, Dcip.hasSwfFormat, ...
        SatSettings.vstbFractionThreshold, ...
        SatSettings.cwfSlidingWindowLengthSec, L);



      %#########################
      % Set UFV and "final" ZVs
      %#########################
      Zv = struct();
      zvUfv                 = Dcip.Zv.ufv;
      Zv.QUALITY_FLAG       = Dcip.Zv.QUALITY_FLAG.min(...
        bicas.utils.FPArray(QUALITY_FLAG));
      Zv.L2_QUALITY_BITMASK = L2_QUALITY_BITMASK;

      % NOTE: Function modifies ZVM handle object in-place (handle object)!
      Zv.currentAAmpere     = bicas.proc.L1L2.qual.set_voltage_current_FV(...
        Dcip.Zv.Epoch, SamplesZvm, currentAAmpere, zvUfv, L);
      Zv.SamplesZvm         = SamplesZvm;



      %##############
      % END FUNCTION
      %##############
      Dcop = bicas.proc.L1L2.DemultiplexingCalibrationOutput(Zv);

      nRecords = size(Dcip.Zv.Epoch, 1);
      Tmk.stop_log(nRecords, 'record')
    end    % process_calibrate_demux



    % Obtain SSID and SDID arrays for arrays of BDM and DLR.
    %
    function [bltsSsidArray, bltsSdidArray] = get_SSID_SDID_arrays(bdmFpa, dlrFpa)
      nRecTot = irf.assert.sizes(...
        bdmFpa, [-1], ...
        dlrFpa, [-1]);

      [iRec1Array, iRec2Array, nSs] = irf.utils.split_by_change( ...
        bdmFpa.int2doubleNan(), ...
        dlrFpa.logical2doubleNan());

      % Preallocate
      % -----------
      % NOTE: No need for bicas.utils.FPArray since SSIDs and SDIDs handle all
      % special cases including unknown source and destination.
      bltsSsidArray = zeros(nRecTot, bicas.const.N_BLTS, 'uint8') + 255;
      bltsSdidArray = zeros(nRecTot, bicas.const.N_BLTS, 'uint8') + 255;

      for iSs = 1:nSs
        iRecSs1 = iRec1Array(iSs);
        iRecSs  = iRec1Array(iSs):iRec2Array(iSs);
        nRecSs  = numel(iRecSs);

        DemuxerRoutingArray = bicas.proc.L1L2.demuxer.get_routings(...
          bdmFpa(iRecSs1), ...
          dlrFpa(iRecSs1));

        ssidArray = [DemuxerRoutingArray.ssid];
        sdidArray = [DemuxerRoutingArray.sdid];

        bltsSsidArray(iRecSs, :) = repmat(ssidArray, nRecSs, 1);
        bltsSdidArray(iRecSs, :) = repmat(sdidArray, nRecSs, 1);
      end
    end



    function log_input_calibration_settings(Dcip, Cal, L)
      % IMPLEMENTATION NOTE: Implemented separately from processing functions
      % since:
      % (1) removes dependence on logger object,
      % (2) can split data into subsequences based on arbitrary choice of
      %     variables,
      % (3) reduces size of processing function where logging would otherwise
      %     be, and
      % (4) can potentially turn table into proper table with column headers.

      iCalibLZv = Cal.get_BIAS_calibration_time_index_L(Dcip.Zv.Epoch);
      iCalibHZv = Cal.get_BIAS_calibration_time_index_H(Dcip.Zv.Epoch);

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
        % NOTE: HK_BIA_DIFF_GAIN needs three characters to print the string "NaN".
        L.logf('info', ['Records %8i-%8i : %s -- %s', ...
          ' bdm/HK_BIA_MODE_MUX_SET=%i;', ...
          ' isAchg/HK_BIA_DIFF_GAIN=%-3i;', ...
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



    % Demultiplex and calibrate VOLTAGES (not e.g. currents). Processes all 5x
    % BLTS --> 5x BLTS in the same call.
    %
    % NOTE: Can handle arrays of any size if the sizes are consistent.
    %
    function bltsSamplesAVolt = calibrate_voltages_5xBLTS(Zv, A)
      % PROPOSAL: Sequence of constant settings includes constant NaN/non-NaN
      %           for CWF.
      %
      % PROPOSAL: Reorg.
      %   FROM
      %     (1) Iterate over subsequences for all channels.
      %         (bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS(); this function)
      %     (2) Iterate over BLTSs (for the same subsequence).
      %         (bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS_subsequence(); sub-function)
      %     (3) Calibrate one BLTS for one subsequence.
      %   TO
      %     (1) Iterate over BLTSs (for all subsequences).
      %     (2) Iterate over subsequences for one BLTS.
      %     (3) Calibrate one BLTS for one subsequence (same as before).
      %   --
      %   PRO: Less code which sees subsequences, more code which only sees
      %        arrays covering entire datasets.
      %   PRO: Can potentially have different sets of groups/subsequences for
      %        different BLTSs.
      %     Ex: ACHG, SSID
      %     PRO: Split depending on fill value/NaN data gaps (instead of UFV).
      %       Might be required if using different blanking for different
      %       BLTS channels.
      %       PRO: Can abolish UFV in bicas.proc.L1L2.Cal.
      %   PRO: Probably simpler automated tests(?)
      %   PRO: Easier to implement/support separate calibration functions for
      %        different channels?!
      %   CON: Splits into subsequences multiple times.
      %     PRO: Potentially slower (in particular SWF).
      %       PRO: Potentially slower wrt. finding subsequences/groups.
      %       PRO: Potentially slower wrt. indexing (which may be what grouping
      %            mitigates for LFR-SWF?).
      %         CON: bltsSamplesTm is the largest array (for SWF) but it is
      %              subdivided by BLTSs first, and again for
      %              subsequences/groups in this scheme.

      arguments
        Zv.Epoch
        Zv.bltsSamplesTm
        Zv.isAchgFpa
        Zv.CALIBRATION_TABLE_INDEX
        Zv.freqHz
        Zv.iLsf
        Zv.ufv
        Zv.bltsSsidArray
        Zv.uspr
        A.isTdsCwf
        A.isLfr
        A.hasSwfFormat
        A.Cal
        A.L
      end

      % ASSERTIONS
      assert(isscalar(A.hasSwfFormat))
      assert(isnumeric(Zv.bltsSamplesTm))
      assert(isa(Zv.bltsSsidArray, 'uint8'))
      [nRecords, aspr] = irf.assert.sizes(...
        Zv.isAchgFpa,     [-1], ...
        Zv.bltsSsidArray, [-1,     bicas.const.N_BLTS], ...
        Zv.bltsSamplesTm, [-1, -2, bicas.const.N_BLTS]);



      % Pre-allocate calibrated array. All (BLTS) channels, all records
      % ---------------------------------------------------------------
      % IMPLEMENTATION NOTE: Preallocation is very important for speeding
      % up LFR-SWF which tends to be broken into subsequences of 1 record.
      bltsSamplesAVolt = nan(nRecords, aspr, bicas.const.N_BLTS);

      iCalibLZv = A.Cal.get_BIAS_calibration_time_index_L(Zv.Epoch);
      iCalibHZv = A.Cal.get_BIAS_calibration_time_index_H(Zv.Epoch);



      %==================================================================
      % (1) Find groups/subsequences of records with identical settings.
      % (2) Process data separately for each such sequence.
      % ----------------------------------------------------------------
      % SS = Subsequence
      %==================================================================
      Tmk = bicas.utils.Timekeeper(...
        'bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS:irf.utils.split_by_change', A.L);
      settingsCa = {...
        Zv.isAchgFpa.logical2doubleNan(), ...
        Zv.freqHz, ...
        Zv.iLsf, ...
        Zv.CALIBRATION_TABLE_INDEX, ...
        Zv.ufv, ...
        Zv.bltsSsidArray, ...
        iCalibLZv, ...
        iCalibHZv};
      if A.hasSwfFormat
        % IMPLEMENTATION NOTE: SWF data is calibrated in a way where consecutive
        % records (snapshots) do NOT affect each other. One can therefore group
        % CDF records in groups of non-consecutive CDF records. This is
        % important for LFR-SWF data which tends to change calibration-relevant
        % settings with every new CDF record and which then becomes a
        % performance problem. Grouping non-consecutive CDF records for LFR-SWF
        % data makes a significant speedup!!!
        % Ex:
        %     solo_L1R_rpw-lfr-surv-swf-e-cdag_20200213_V10.cdf for the relevant
        %     code section: ~29 s --> ~4 s (843 CDF records)
        iGroupArCa = bicas.utils.group_unique_rows(settingsCa{:});
      else
        % IMPLEMENTATION NOTE: CWF data is calibrated in a way where consecutive
        % records affect each other. Must therefore divide CDF records in groups
        % (subsequences) of continuous CDF records.
        iGroupArCa = bicas.utils.group_by_change(settingsCa{:});
      end
      nGroups = numel(iGroupArCa);
      Tmk.stop_log(nRecords, 'record', nGroups, 'subsequence-group')

      %===========
      % Calibrate
      %===========
      Tmk = bicas.utils.Timekeeper(...
        'bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS:Calibrating voltages', A.L);
      for iGroup = 1:nGroups
        iGroupAr = iGroupArCa{iGroup};
        iRec1    = iGroupAr(1);

        ssBltsSamplesAVolt = bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS_subsequence( ...
          Cal           = A.Cal, ...
          ... % ===============================================================
          ... % Variables which do VARY over CDF records, but not over
          ... % the subsequence/group.
          isAchgFpa     = Zv.isAchgFpa(              iRec1), ...
          freqHz        = Zv.freqHz(                 iRec1), ...
          iLsf          = Zv.iLsf(                   iRec1), ...
          zvcti         = Zv.CALIBRATION_TABLE_INDEX(iRec1, :), ...
          ufv           = Zv.ufv(                    iRec1), ...
          bltsSsidArray = Zv.bltsSsidArray(          iRec1, :), ...
          iCalibL       = iCalibLZv(                 iRec1), ...
          iCalibH       = iCalibHZv(                 iRec1), ...
          ... % ===============================================================
          ... % Variables which do not vary over CDF records ("constants").
          hasSwfFormat  = A.hasSwfFormat, ...
          isLfr         = A.isLfr, ...
          isTdsCwf      = A.isTdsCwf, ...
          aspr          = aspr, ...
          ... % ===============================================================
          ... % Variables which vary by CDF records, and which vary over
          ....% the subsequence/group.
          Epoch         = Zv.Epoch(                 iGroupAr), ...
          bltsSamplesTm = Zv.bltsSamplesTm(         iGroupAr, :, :), ...
          uspr          = Zv.uspr(iGroupAr));

        % Add subsequence/group signals to the global array (all records).
        bltsSamplesAVolt(iGroupAr, :, :) = ssBltsSamplesAVolt;
      end    % for
      Tmk.stop_log(nRecords, 'record', nGroups, 'subsequence-group')

    end    % calibrate_voltages_5xBLTS



    % Calibrate and demux all 5x BLTS channels for one subsequence with various
    % constant settings/values.
    %
    % ARGUMENTS
    % =========
    % Cv
    %       Constant values. Scalar values which do NOT VARY by CDF record.
    % Vv
    %       Varying values. Struct with values which DO VARY by CDF record.
    %
    function ssBltsSamplesAVolt = calibrate_voltages_5xBLTS_subsequence(A, Cv, Vv)
      arguments
        A.Cal
        %
        % NOTE: Excluding LRX since it is only needed for splitting time/CDF
        %       record intervals, not for calibration since calibration can
        %       handle sequences of only NaN.
        Cv.isAchgFpa
        Cv.freqHz
        Cv.iLsf
        Cv.zvcti
        Cv.ufv
        Cv.bltsSsidArray
        Cv.iCalibL
        Cv.iCalibH
        % Below variables do not vary over CDF records at all.
        Cv.aspr
        Cv.hasSwfFormat
        Cv.isLfr
        Cv.isTdsCwf
        % Below variables vary over CDF records in subsequence.
        Vv.Epoch
        Vv.bltsSamplesTm
        Vv.uspr
      end
      irf.assert.sizes( ...
        Cv.bltsSsidArray, [ 1,          bicas.const.N_BLTS], ...
        Vv.Epoch,         [-1], ...
        Vv.bltsSamplesTm, [-1, Cv.aspr, bicas.const.N_BLTS], ...
        Vv.uspr,          [-1])

      nRows = numel(Vv.Epoch);

      if Cv.hasSwfFormat
        % CASE: SWF
        % NOTE: Column vector of (identical) numbers (one per snapshot, which
        %       will all be processed separately).
        dtSec = ones(nRows, 1) / Cv.freqHz;
      else
        % CASE: CWF
        % NOTE: Scalar (since all data will be processed in one session).
        % BUG RISK: Has the possibility of data gaps been excluded at this point
        % in the code, since the calculation seems to rely on it?
        dtSec = double( Vv.Epoch(end) - Vv.Epoch(1) ) / (nRows-1) * 1e-9;
      end

      %====================
      % CALIBRATE VOLTAGES
      %====================
      ssBltsSamplesAVolt = nan(nRows, Cv.aspr, bicas.const.N_BLTS);
      for iBlts = 1:bicas.const.N_BLTS
        ssBltsSamplesAVolt(:, :, iBlts) = bicas.proc.L1L2.dc.calibrate_1xBLTS_subsequence(...
          ... % Below variables vary within the subsequence.
          samplesTm    = Vv.bltsSamplesTm(:, :, iBlts), ...
          uspr         = Vv.uspr, ...
          ... % Below variables DO NOT vary within the subsequence (almost).
          ssid         = Cv.bltsSsidArray(iBlts), ...
          iBlts        = iBlts, ...
          hasSwfFormat = Cv.hasSwfFormat, ...
          isAchg       = Cv.isAchgFpa.logical2doubleNan(), ...
          iCalibL      = Cv.iCalibL, ...
          iCalibH      = Cv.iCalibH, ...
          iLsf         = Cv.iLsf, ...
          dtSec        = dtSec, ...
          isLfr        = Cv.isLfr, ...
          isTdsCwf     = Cv.isTdsCwf, ...
          zvcti        = Cv.zvcti, ...
          ufv          = Cv.ufv, ...
          Cal          = A.Cal);
      end
    end    % calibrate_voltages_5xBLTS_subsequence



    % Calibrate one BLTS channel.
    function samplesAVolt = calibrate_1xBLTS_subsequence(A)
      arguments
        A.ssid
        A.samplesTm
        A.iBlts
        A.hasSwfFormat
        A.uspr
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
      % but the original code made calibrate_voltages_5xBLTS() to large and
      % unwieldy. Having many arguments also highlights the exact dependencies.

      % PROPOSAL: CalSettings as parameter.
      %   PRO: Reduces number of parameters.
      %   PROPOSAL: Add values to CalSettings: isLfr, isTdsCwf, zvcti
      %       CON: cal does not seem to use more values.
      % PROPOSAL: Reorder arguments to group them.
      %   PROPOSAL: Group arguments from DCIP.

      % ASSERTIONS
      % ----------
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
      irf.assert.sizes(A.samplesTm, [-1, -2])   % One BLTS channel. CWF/SWF.
      assert(isscalar(A.ufv))



      if isequaln(A.ssid, bicas.proc.L1L2.const.C.SSID_DICT("UNKNOWN"))
        % ==> Calibrated data set to NaN.
        samplesAVolt = nan(size(A.samplesTm));

      elseif ismember(...
          A.ssid, bicas.proc.L1L2.const.C.SSID_DICT(["GND", "REF25V"]))
        % ==> No calibration.
        % NOTE: samplesTm stores TM units using float!
        samplesAVolt = A.samplesTm;

      else
        assert(bicas.proc.L1L2.const.SSID_is_ASR(A.ssid))
        % ==> Calibrate (unless explicitly stated that should not)

        if A.hasSwfFormat
          bltsSamplesTmCa = ...
            bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(...
            double(A.samplesTm), A.uspr);
        else
          assert(all(A.uspr == 1))
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
          A.iBlts, A.ssid, A.isAchg, A.iCalibL, A.iCalibH, A.iLsf);
        %#######################################################
        ssBltsSamplesAVoltCa = A.Cal.calibrate_voltage_all(...
          A.dtSec, bltsSamplesTmCa, ...
          A.isLfr, A.isTdsCwf, CalSettings, ...
          A.zvcti, A.ufv);
        %#######################################################

        if A.hasSwfFormat
          % CASE: SWF
          [samplesAVolt, ~] = ...
            bicas.proc.utils.convert_cell_array_of_vectors_to_matrix(...
            ssBltsSamplesAVoltCa, ...
            size(A.samplesTm, 2));
        else
          % CASE: CWF
          % NOTE: Scalar, since not snapshot.
          assert(isscalar(ssBltsSamplesAVoltCa))
          samplesAVolt = ssBltsSamplesAVoltCa{1};
        end
      end
    end    % calibrate_1xBLTS_subsequence



    % function AsrSamplesAVoltSrm = relabel_reconstruct_samples_5xBLTS_to_9xASR_OLD(...
    %     bltsSamplesAvolt, bltsSdidArray, L)
    %   % PROPOSAL: Automated tests.
    %
    %   Tmk = bicas.utils.Timekeeper('bicas.proc.L1L2.dc.relabel_reconstruct_samples_5xBLTS_to_9xASR_OLD', L);
    %
    %   [nRecTot, aspr] = irf.assert.sizes(...
    %     bltsSamplesAvolt, [-1, -2, bicas.const.N_BLTS], ...
    %     bltsSdidArray,    [-1,     bicas.const.N_BLTS]);
    %
    %   % -----------------------------------------------------------------
    %   % Pre-allocate AsrSamplesAVoltSrm: All (ASID) channels, all records
    %   % -----------------------------------------------------------------
    %   % IMPLEMENTATION NOTE: Preallocation is very important for speeding up
    %   % LFR-SWF which tends to be broken into subsequences of 1 record due to
    %   % changing sampling rate.
    %   AsrSamplesAVoltSrm = bicas.utils.SameRowsMap(...
    %     "uint8", nRecTot, 'CONSTANT', ...
    %     NaN(nRecTot, aspr), ...
    %     bicas.proc.L1L2.const.C.ASID_DICT.values);
    %
    %
    %
    %   [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(bltsSdidArray);
    %   for iSs = 1:nSs
    %     iRec1 = iRec1Ar(iSs);
    %     iRec2 = iRec2Ar(iSs);
    %
    %     SsAsrSamplesAVoltSrm = ...
    %       bicas.proc.L1L2.demuxer.relabel_reconstruct_samples_5xBLTS_to_9xASR_subsequence_OLD(...
    %       bltsSdidArray(   iRec1,          :), ...
    %       bltsSamplesAvolt(iRec1:iRec2, :, :));
    %
    %     % Set demuxed subsequence signals (some records/indices) in the
    %     % pre-existing global arrays (all records).
    %     AsrSamplesAVoltSrm.set_rows(SsAsrSamplesAVoltSrm, [iRec1:iRec2]');
    %   end
    %
    %   Tmk.stop_log(nRecTot, 'record', nSs, 'subsequence')
    % end



    % Derive VSIB from BLTS samples. Vectorized.
    %
    % RETURN VALUE
    % ============
    % bltsVsibAr
    %       N x 5. SWF: Set if at least one bit is set for any sample within
    %       a snapshot.
    function bltsVsibAr = get_VSIB_5xBLTS_NEW(...
        bltsSamplesAVoltAr, uspr, bltsSsidAr, isAchgFpa, Bso, L)
      % PROPOSAL: SatSettings as argument.
      %   PRO: Simpler test code.
      % PROPOSAL: Replace bltsSamplesAVoltAr --> 5x SCHD
      %   PRO: Simpler handling of dimensions.
      %   PRO: Less risk of memory problems if iterates manually over BLTSs.
      %   PRO: SCHD could be modified to store USPR and
      %        handle padded snapshots explicitly.
      %   CON: SCHD contains 1 VSIB/record. This function is intended for
      %        setting that VSIB array.
      %   CON-PROPOSAL: Custom class representing jagged array (jagged in one
      %                 dimension, dim.=2 for BICAS) for any generic MATLAB
      %                 class.
      %     PRO: SCHD could use it in its implementation.
      %     CON: Can not handle fill values like FPAs.
      %       PRO: Can not replace with FPAs in the long term unless
      %            (1) the class itself uses FPAs, or
      %            (2) implements fill positions itself.
      %     Ex: Store for VSTB for SWF data.
      %     Ex: Store samples for SWF (in particular TDS-LFM SWF).

      Tmk = bicas.utils.Timekeeper('get_VSIB_5xBLTS_NEW', L);

      % size(bltsSamplesAVoltAr), size(bltsSsidAr), size(isAchgFpa)
      [nRec, aspr] = irf.assert.sizes(...
        bltsSamplesAVoltAr, [-1, -2, bicas.const.N_BLTS], ...
        uspr,               [-1], ...
        bltsSsidAr,         [-1,     bicas.const.N_BLTS], ...
        isAchgFpa,          [-1]);
      assert(bicas.proc.L1L2.const.is_SSID(bltsSsidAr))

      % Expand variables to be of the same size as bltsSamplesAVoltAr
      % -------------------------------------------------------------
      % Needed for submitting arguments to bicas.proc.L1L2.sat.get_VSTB_NEW()
      % NOTE: This could possibly lead to memory problems, which could be
      % mitigated by e.g. calling bicas.proc.L1L2.sat.get_VSTB_NEW()
      % once per BLTS.
      isAchgFpa   = repmat(        isAchgFpa,              [1, aspr, bicas.const.N_BLTS]);
      bltsSsidAr  = repmat(permute(bltsSsidAr, [1, 3, 2]), [1, aspr, 1                 ]);

      % Expand variable to size nRec x 1 x N_BLTS, needed for later comparison
      uspr        = repmat(uspr, [1, 1, 5]);

      SatSettings = bicas.proc.L1L2.sat.from_BSO_extract_saturation_settings(Bso);



      bltsVstbAr  = bicas.proc.L1L2.sat.get_VSTB_NEW(...
        SatSettings, bltsSamplesAVoltAr, bltsSsidAr, isAchgFpa);

      % Normalize CWF/SWF data to one array format (one VSIB per record & BLTS).
      % N x M x 5 --> N x 1 x 5 --> N x 5
      % VSTB --> VSIB
      if aspr >= 2
        % CASE: SWF

        % POTENTIAL BUG: Relies on that VSTB=false for padded samples/elements
        % at the end of TDS-RSWF snapshots (the arrays are padded not account
        % for varying-length snapshots).
        bltsVsibAr = (sum(bltsVstbAr, 2) ./ uspr) >= SatSettings.vstbFractionThreshold;
      else
        % CASE: CWF
        bltsVsibAr = bltsVstbAr;
      end
      bltsVsibAr = permute(bltsVsibAr, [1, 3, 2]);

      irf.assert.sizes(bltsVsibAr, [nRec, bicas.const.N_BLTS])

      Tmk.stop_log(nRec, 'record')
    end



    % Convert samples stored as 5x BLTSs to 9x ASRs (without reconstructing
    % missing data).
    %
    function SchdZvm = convert_samples_5xBLTS_to_9xASR_NEW( ...
        bltsSamplesAVoltAr, bltsVsibAr, bltsSdidAr, L)

      % Tmk = bicas.utils.Timekeeper(...
      %   'bicas.proc.L1L2.dc.convert_samples_5xBLTS_to_9xASR_NEW', L);

      [nRec, aspr] = irf.assert.sizes(...
        bltsSamplesAVoltAr, [-1, -2, bicas.const.N_BLTS], ...
        bltsSdidAr,         [-1,     bicas.const.N_BLTS], ...
        bltsVsibAr,         [-1,     bicas.const.N_BLTS]);

      %====================================
      % Convert BLTS arrays --> bltsSchdCa
      %====================================
      bltsSchdCa = cell(5, 1);
      for iBlts = 1:bicas.const.N_BLTS
        bltsSchdCa{iBlts, 1} = bicas.proc.L1L2.SingleChannelData( ...
          bltsSamplesAVoltAr(:, :, iBlts), ...
          bltsVsibAr(        :,    iBlts));
      end

      %========================================================================
      % Convert bltsSchdCa --> SchdZvm
      % ------------------------------
      % Copy values from 5x BLTSs into 9x SCHD (without reconstructing missing
      % values)
      %========================================================================
      SchdZvm = bicas.utils.ZvMap(nRec);
      for asrSdid = bicas.proc.L1L2.const.C.SDID_ASR_AR'

        % ~Preallocate empty SCHD for the current ASR/SDID.
        AsrSchd = bicas.proc.L1L2.SingleChannelData(...
          nan(  nRec, aspr), ...
          false(nRec, 1));

        for iBlts = 1:bicas.const.N_BLTS
          bRecCopy          = ( bltsSdidAr(:, iBlts) == asrSdid );

          % Copy BLTS samples into selected records/rows in previously created
          % SCHD for ASR.
          BltsSchd          = bltsSchdCa{iBlts};
          AsrSchd(bRecCopy) = BltsSchd(bRecCopy);
        end
        SchdZvm.add(asrSdid, AsrSchd)
      end

      % Tmk.stop_log(nRec, 'record')
    end



  end    % methods(Static)



end
