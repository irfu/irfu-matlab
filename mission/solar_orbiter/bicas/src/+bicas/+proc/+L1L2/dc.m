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
  %   cdr  = calibration+demuxing+reconstruction (same order as execution)
  %   cdrq = calibration+demuxing+reconstruction+quality
  %     CON: Setting quality variables is not necessarily last in the execution.
  %      Ex: Blanking data due to failed antenna should be done in TM channels
  %          (planned implementation).
  %       CON: Minor. Quality variables are still set somewhere in this file.
  %
  % PROPOSAL: More automatic test code.
  %
  % PROPOSAL: process_calibrate_demux() should only accept the necessary ZVs and
  %           variables and nothing else. Currently receives Dcip which covers
  %           all ZVs.
  %   TODO-NI: Does function actually receive input it does not use? It might
  %            actually use all input ZVs.
  %     PROPOSAL: Assertion on Dcip.Zv fields.
  %
  % PROPOSAL: Reorg. code to
  %   * Consist of more isolated/modular/generic separate steps.
  %   * Be more easily testable.
  %   * Be easier to understand.
  %   * Not require (as much) splitting up CDF records into subsequences with
  %     constant "settings" (values for specific zVariables not varying as
  %     a function of CDF record), in particular not require many constant
  %     "settings".
  %   * Be more natural to implement.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Derive DCOP from DCIP:
    % * Calibrate bias currents.
    % * Voltages:
    %   (1) calibrate current and voltage samples
    %   (2) demux (demultiplex): Relabel voltage samples from BLTS to SDID.
    %   (3) reconstruct (derive) voltage samples for missing channels from
    %       calibrated voltage samples (e.g. DC_V12 := DC_V1 - DC_V2)
    % * Set quality variables.
    %
    function Dcop = process_calibrate_demux(...
        Dcip, InCurPd, outputDsi, Vcal, Ccal, NsoTable, Bso, L)

      Tmk = bicas.utils.Timekeeper('bicas.proc.L1L2.dc.process_calibrate_demux', L);

      % ASSERTION
      assert(ischar(outputDsi))
      assert(isa(Vcal, "bicas.proc.L1L2.cal.VoltageCalibration"))
      assert(isa(Ccal, "bicas.proc.L1L2.cal.CurrentCalibrationAbstract"))
      assert(isa(Dcip, "bicas.proc.L1L2.DemultiplexingCalibrationInput"));

      bicas.proc.L1L2.dc.log_input_calibration_settings(Dcip, Vcal, L)



      %##################################
      % NSO table, L1/L1R ~ZVs --> QRCBs
      %##################################
      % Read NSO table into QRCBs ONCE, so that it does not need to be done
      % later.
      L2Qrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        bicas.const.qrc.Q.L2_QRCSM.qrcidAr, NsoTable, Dcip.Zv.Epoch, L);
      clear NsoTable
      % Convert information about BIAS ON/OFF and sweeps into QRCBs.
      L2Qrcbm.set("BIAS_HW_OFF", Dcip.Zv.biasOffQrcb);
      L2Qrcbm.set("SWEEP",       Dcip.Zv.sweepQrcb);
      % PROPOSAL: Clear Dcip.Zv.biasOffQrcb & Dcip.Zv.sweepQrcb since they
      %           should not be used after this point.



      %#########################
      % Calibrate bias CURRENTS
      %#########################
      currentAampere = bicas.proc.L1L2.cur.calibrate_bias_currents(...
        InCurPd.Zv.Epoch, ...
        [InCurPd.Zv.IBIAS_1, ...
        InCurPd.Zv.IBIAS_2, ...
        InCurPd.Zv.IBIAS_3], ...
        Dcip.Zv.Epoch, Ccal, Bso, L);



      %######################################################################
      % Obtain "demultiplexer" "routings" in the form of SSID-SDID pairs for
      % every BLTS (and CDF record)
      %######################################################################
      [bltsSsidArray, bltsSdidArray] = bicas.proc.L1L2.dc.get_SSID_SDID_arrays(...
        Dcip.Zv.bdmFpa, ...
        Dcip.Zv.dlrFpa);



      %#################################################
      % Set 5xBLTS channels to NaN/FV according to QRCs
      %#################################################
      % BUG: Setting voltage NaN/FVs from QRCs *BEFORE* calibration which means
      % that calibration of non-blanked samples is *NOT* influenced by samples
      % which are later blanked.
      %   TODO-NI: Why is this labelled as a bug?!! This seems OK!
      %            /EJ 2025-08-01
      aspr          = size(Dcip.Zv.bltsVoltageTm, 2);
      btlsSsidAr2   = repmat(permute(bltsSsidArray, [1 3 2]), [1, aspr, 1]);
      bltsVoltageTm = bicas.proc.L1L2.qrc.set_5xBLTS_voltage_samples_FV(...
        Dcip.Zv.bltsVoltageTm, btlsSsidAr2, L2Qrcbm, bicas.const.qrc.Q.L2_QRCSM);



      %##################################################################
      % CALIBRATE VOLTAGES: 5x CHANNELS LABELLED BY SSID/BLTS (not SDID)
      %##################################################################
      % NOTE: Takes most of the time LFR-SWF.
      bltsVoltageAvolt = bicas.proc.L1L2.dc.calibrate_voltage_5xBLTS(...
        tt2000       = Dcip.Zv.Epoch, ...
        voltageTm    = bltsVoltageTm, ...       % Partially blanked by QRCs.
        isAchgFpa    = Dcip.Zv.isAchgFpa, ...
        NbriFpa      = Dcip.Zv.NbriFpa, ...
        NbciFpa      = Dcip.Zv.NbciFpa, ...
        freqHz       = Dcip.Zv.freqHz, ...
        iLsf         = Dcip.Zv.iLsf, ...
        ssid         = bltsSsidArray, ...
        isTdsCwf     = Dcip.isTdsCwf, ...
        isLfr        = Dcip.isLfr, ...
        hasSwfFormat = Dcip.hasSwfFormat, ...
        uspr         = Dcip.Zv.uspr, ...
        Vcal         = Vcal, ...
        L            = L);



      %#####################################
      % Get VSIB for BLTS-labelled channels
      %#####################################
      SatSettings = bicas.proc.L1L2.sat.from_BSO_extract_saturation_settings(Bso);
      bltsVsibAr  = bicas.proc.L1L2.dc.get_VSIB_5xBLTS(...
        bltsVoltageAvolt, Dcip.hasSwfFormat, Dcip.Zv.uspr, ...
        bltsSsidArray, Dcip.Zv.isAchgFpa, SatSettings, L);



      %######################################################################
      % ~"DEMUX" VOLTAGES:
      % INPUT:  5x SIGNALS LABELLED BY SSID/BLTS
      % OUTPUT: 9x SIGNALS LABELLED BY SDID + RECONSTRUCTING MISSING SIGNALS
      %######################################################################
      % NOTE: Needs VSIB for propagating VSIB to reconstructed channels.
      SchdZvm = bicas.proc.L1L2.dc.convert_voltage_5xBLTS_to_9xASR(...
        bltsVoltageAvolt, bltsVsibAr, bltsSdidArray, L);
      bicas.proc.L1L2.demuxer.reconstruct_ASR_voltage_channels(SchdZvm);



      % Replace/split variable: SchdZvm --> VoltageZvm + VsibZvm
      VoltageZvm = bicas.utils.ZvMap(SchdZvm.nRecords);
      VsibZvm    = bicas.utils.ZvMap(SchdZvm.nRecords);
      for keyCa = SchdZvm.keyCa'
        VoltageZvm.add(keyCa{1}, SchdZvm.get(keyCa{1}).samplesAr);
        VsibZvm.add(   keyCa{1}, SchdZvm.get(keyCa{1}).vsibAr);
      end
      clear SchdZvm



      %####################
      % Obtain quality ZVs
      %####################
      SatSettings = bicas.proc.L1L2.sat.from_BSO_extract_saturation_settings(Bso);
      % NOTE: Whether voltage samples have already been blanked (set to NaN/FV)
      % or not based on QRCs, affects the saturation detection. Can not
      % autodetect saturation in blanked data.
      SaturationQrcbm = bicas.proc.L1L2.qrc.get_saturation_QRCBs(...
        Dcip.Zv.Epoch,  ...
        string(Bso.get_fv('PROCESSING.SATURATION.QUALITY_SCHEME')), ...
        VsibZvm, Dcip.hasSwfFormat, ...
        SatSettings.vstbFractionThreshold, ...
        SatSettings.cwfSlidingWindowLengthSec);
      L2Qrcbm.union(SaturationQrcbm)
      % --
      [QUALITY_FLAG, L2_QUALITY_BITMASK] = ...
        bicas.proc.qrc.QRCB_arrays_to_quality_ZVs(...
        L2Qrcbm, bicas.const.qrc.Q.L2_QRCSM, "L2_QUALITY_BITMASK");



      %#################
      % Set "final" ZVs
      %#################
      Zv = struct();
      Zv.QUALITY_FLAG       = Dcip.Zv.QUALITY_FLAG.min(...
        bicas.utils.FPArray(QUALITY_FLAG));
      Zv.L2_QUALITY_BITMASK = L2_QUALITY_BITMASK;

      % NOTE: Function modifies VoltageZvm handle object in-place!
      Zv.currentAampere     = bicas.proc.L1L2.qrc.set_current_samples_FV(...
        currentAampere, L2Qrcbm, bicas.const.qrc.Q.L2_QRCSM);
      Zv.VoltageZvm         = VoltageZvm;



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



    function log_input_calibration_settings(Dcip, Vcal, L)
      % IMPLEMENTATION NOTE: Implemented separately from processing functions
      % since:
      % (1) removes dependence on logger object,
      % (2) can split data into subsequences based on arbitrary choice of
      %     variables,
      % (3) reduces size of processing function where logging would otherwise
      %     be, and
      % (4) can potentially turn table into proper table with column headers.
      %
      % PROPOSAL: Move to where processing splits data.
      %   NOTE: Splits once per BLTS ==> More output.

      iCalibL = Vcal.get_BIAS_calibration_time_index_L(Dcip.Zv.Epoch);
      iCalibH = Vcal.get_BIAS_calibration_time_index_H(Dcip.Zv.Epoch);

      % IMPLEMENTATION NOTE: Do not log for LFR SWF since it produces
      % unnecessarily many log messages since sampling frequencies change
      % for every CDF record.
      if Dcip.hasSwfFormat && Dcip.isLfr
        return
      end

      % IMPLEMENTATION NOTE: Indexing FPA separately, since
      % (1) MATLAB does not permit first indexing and then calling a method on
      %     the result in one command.
      % (2) irf.utils.split_by_change() is significantly faster if using
      %     builtin arrays rather than FPAs.
      bdm    = Dcip.Zv.bdmFpa   .int2doubleNan();
      isAchg = Dcip.Zv.isAchgFpa.logical2doubleNan();
      dlr    = Dcip.Zv.dlrFpa   .logical2doubleNan();
      nbri   = Dcip.Zv.NbriFpa  .int2doubleNan();
      nbci   = Dcip.Zv.NbciFpa  .int2doubleNan();

      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        bdm, ...
        isAchg, ...
        Dcip.Zv.freqHz, ...
        Dcip.Zv.iLsf, ...
        nbri, ...
        nbci, ...
        dlr, ...
        Dcip.Zv.lrx, ...
        iCalibL, ...
        iCalibH);

      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        % PROPOSAL: Make into "proper" table with top rows with column names.
        %   NOTE: Can not use irf.str.assist_print_table() since
        %         it requires the entire table to pre-exist before execution.
        %   PROPOSAL: Print after all iterations.
        %
        % NOTE: HK_BIA_DIFF_GAIN needs three characters to print the string
        %       "NaN".

        L.logf('info', ['Records %8i-%8i : %s -- %s', ...
          ' bdm/HK_BIA_MODE_MUX_SET=%i;', ...
          ' isAchg/HK_BIA_DIFF_GAIN=%-3i;', ...
          ' dlr/HK_BIA_MODE_DIFF_PROBE=%i;', ...
          ' freqHz=%5g; iCalibL=%i; iCalibH=%i', ...
          ' ~CALIBRATION_TABLE_INDEX=[%i, %i]'], ...
          iRec1, iRec2, ...
          bicas.utils.TT2000_to_UTC_str(Dcip.Zv.Epoch(iRec1), 9), ...
          bicas.utils.TT2000_to_UTC_str(Dcip.Zv.Epoch(iRec2), 9), ...
          bdm(           iRec1), ...
          isAchg(        iRec1), ...
          dlr(           iRec1), ...
          Dcip.Zv.freqHz(iRec1), ...
          iCalibL(       iRec1), ...
          iCalibH(       iRec1), ...
          nbri(          iRec1), ...
          nbci(          iRec1))
      end    % for

    end



    % Demultiplex and calibrate VOLTAGES (not e.g. currents). Processes all 5x
    % BLTS --> 5x BLTS in the same call.
    %
    % NOTE: Can handle arrays of any size if the sizes are consistent.
    %
    function voltageAvolt = calibrate_voltage_5xBLTS(Cv, Zv)
      arguments
        % Variables which DO NOT VARY over CDF records at all.
        Cv.Vcal
        Cv.L
        Cv.hasSwfFormat
        Cv.isTdsCwf
        Cv.isLfr
        % Variables which VARY over CDF records.
        Zv.tt2000
        Zv.voltageTm
        Zv.ssid
        Zv.uspr
        Zv.freqHz
        Zv.isAchgFpa
        Zv.iLsf
        Zv.NbriFpa
        Zv.NbciFpa
      end

      % ASSERTIONS
      assert(isa(Cv.Vcal, "bicas.proc.L1L2.cal.VoltageCalibration"))
      assert(isscalar(Cv.hasSwfFormat))
      assert(isnumeric(Zv.voltageTm))
      assert(isa(Zv.ssid, 'uint8'))
      [nRecords, aspr] = irf.assert.sizes(...
        Zv.isAchgFpa, [-1], ...
        Zv.ssid,      [-1,     bicas.const.N_BLTS], ...
        Zv.voltageTm, [-1, -2, bicas.const.N_BLTS]);



      % Pre-allocate array of calibrated voltage samples: All (BLTS) channels
      % ---------------------------------------------------------------------
      % IMPLEMENTATION NOTE: Preallocation is very important for speeding up
      % LFR-SWF which tends to be broken into subsequences of 1 record.
      voltageAvolt = nan(nRecords, aspr, bicas.const.N_BLTS);



      %===========
      % Calibrate
      %===========
      Tmk = bicas.utils.Timekeeper(...
        'bicas.proc.L1L2.dc.calibrate_voltage_5xBLTS:Calibrating voltages', Cv.L);
      % IMPLEMENTATION NOTE: Iterating over BLTSs before iterating over
      % subsequences/groups of CDF records in subfunction because
      % PRO: Less code sees samples combined over BLTSs.
      %      ==> Simpler automated tests (Easier to hardcode data).
      % PRO: Can split CDF records into subsequences/groups differently
      %      depending on BLTS/channel. This is needed for splitting based on
      %      NaN/FV/data gap (planned; 2025-07-14).
      %   Ex: ACHG, SSID
      % PRO: Can potentially parallelize processing over different BLTSs.
      for iBlts = 1:bicas.const.N_BLTS

        bltsVoltageAvolt2 = bicas.proc.L1L2.dc.calibrate_voltage_1xBLTS( ...
          ... % ===============================================================
          ... % Variables which DO NOT VARY over CDF records at all.
          Vcal         = Cv.Vcal, ...
          L            = Cv.L, ...
          iBlts        = iBlts, ...
          hasSwfFormat = Cv.hasSwfFormat, ...
          isLfr        = Cv.isLfr, ...
          isTdsCwf     = Cv.isTdsCwf, ...
          ... % ===============================================================
          ... % Variables which VARY over CDF records.
          tt2000       = Zv.tt2000, ...
          ssid         = Zv.ssid(     :,    iBlts), ...
          voltageTm    = Zv.voltageTm(:, :, iBlts), ...
          uspr         = Zv.uspr, ...
          freqHz       = Zv.freqHz, ...
          iLsf         = Zv.iLsf, ...
          isAchgFpa    = Zv.isAchgFpa, ...
          NbriFpa      = Zv.NbriFpa, ...
          NbciFpa      = Zv.NbciFpa);

        % Add subsequence/group signals to the global array (all records).
        voltageAvolt(:, :, iBlts) = bltsVoltageAvolt2;
      end    % for
      Tmk.stop_log(nRecords, 'record')

    end    % calibrate_voltage_5xBLTS



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
    function voltageAvolt = calibrate_voltage_1xBLTS(Cv, Zv)
      % PROPOSAL: Also split sequence based constant isnan() (isfinite()?)
      %           for CWF (not SWF).
      %   NOTE: bicas.tf.apply_TF() can split based on isfinite() (not isnan()).
      % PROPOSAL: Use parfor for iterating over calls to
      %           bicas.proc.L1L2.dc.calibrate_voltage_1xBLTS_subsequence().
      %   PRO: No logging inside loop (by default).
      %   PROBLEM: MATLAB error when using parfor:
      %     Exception.identifier = "MATLAB:mir_error_parfor_bad_temporary_variable"
      %     Exception.message    = "Error: File: dc.m Line: 449 Column: 14
      %     Unable to classify the variable 'voltageAvolt' in the body of the parfor-loop. For more information, see Parallel for Loops in MATLAB, "Solve Variable Classification Issues in parfor-Loops"."
      %
      % BUG/INEFFICIENCY: Groups by ACHG also for DC data. ==> Can divide into
      % unnecesssarily small groups.
      %   PROPOSAL: Only split on ACHG for AC data. (Check SSID.)
      %     CON: SSID changes per record.
      %       CON: In practice, SSID should be either all DC or all AC per BLTS.
      %     CON: ACHG probably does not change for DC data. ==> No effect.

      arguments
        % NOTE: Excluding LRX since it is only needed for splitting time/CDF
        %       record intervals (but ssid should handle that), not for
        %       calibration since calibration can handle sequences of only NaN.
        %
        % Variables which DO NOT VARY over CDF records at all.
        Cv.Vcal
        Cv.L
        Cv.iBlts
        Cv.hasSwfFormat
        Cv.isLfr
        Cv.isTdsCwf
        % Variables which VARY over CDF records.
        Zv.tt2000
        Zv.ssid
        Zv.voltageTm
        Zv.uspr
        Zv.freqHz
        Zv.iLsf
        Zv.isAchgFpa
        Zv.NbriFpa
        Zv.NbciFpa
      end
      [nRecords, aspr] = irf.assert.sizes( ...
        Zv.ssid,      [-1,     1], ...
        Zv.tt2000,    [-1], ...
        Zv.voltageTm, [-1, -2, 1], ...
        Zv.uspr,      [-1]);

      iCalibL = Cv.Vcal.get_BIAS_calibration_time_index_L(Zv.tt2000);
      iCalibH = Cv.Vcal.get_BIAS_calibration_time_index_H(Zv.tt2000);



      %=======================================================================
      % Find groups/subsequences of records with identical settings and which
      % can be processed separately.
      %=======================================================================
      TmkGrouping = bicas.utils.Timekeeper(...
        'bicas.proc.L1L2.dc.calibrate_voltage_1xBLTS:grouping algo.', Cv.L);
      % IMPLEMENTATION NOTE: Important to convert FPAs to builtin arrays to
      % significantly speed up both bicas.utils.group_unique_rows(), and
      % bicas.utils.group_by_change().

      settingsCa = {...
        Zv.freqHz, ...
        Zv.iLsf, ...
        Zv.NbriFpa.int2doubleNan(), ...
        Zv.NbciFpa.int2doubleNan(), ...
        Zv.ssid, ...
        iCalibL, ...
        iCalibH};

      if any(bicas.proc.L1L2.const.SSID_is_AC(Zv.ssid))
        % IMPLEMENTATION NOTE: Non-ASR SSIDs count as false. Using any() is
        % conservative in the sense that it will use the variable for splitting
        % in case of uncertainty.
        settingsCa{end+1} = Zv.isAchgFpa.logical2doubleNan();
      end

      if Cv.hasSwfFormat
        %================
        % CASE: SWF/RSWF
        %================
        % IMPLEMENTATION NOTE: SWF data is by its nature calibrated in a way
        % where consecutive CDF records (snapshots) do NOT affect each other.
        % ==> (1) group SWF CDF records in groups of non-consecutive CDF
        %         records i.e. use bicas.utils.group_unique_rows()
        %     (2) does not need to look for data gaps (between
        %         snapshots/records).
        %
        % This is important for LFR-SWF data which tends to change
        % calibration-relevant settings with every new CDF record and which
        % then becomes a performance problem. Grouping non-consecutive CDF
        % records for LFR-SWF data leads to a significant speedup!!!
        % Ex:
        %     solo_L1R_rpw-lfr-surv-swf-e-cdag_20200213_V10.cdf for the
        %     relevant code section: ~29 s --> ~4 s (843 CDF records)
        %
        % NOTE: It is possible that bicas.utils.group_unique_rows() becomes
        % slow (t~O(n^2)) for large LFR-SWF CDFs with mode changes at the end
        % due to the algorithm. May need to optimize the algorithm if this is
        % observed.
        iGroupArCa = bicas.utils.group_unique_rows(settingsCa{:});
      else
        %===========
        % CASE: CWF
        %===========
        % IMPLEMENTATION NOTE: CWF data is calibrated in a way where
        % consecutive records affect each other. Must therefore divide CDF
        % records in groups (subsequences) of *CONTINUOUS* CDF records.
        % ==> (1) Use bicas.utils.group_by_change(),
        %     (2) look for data gaps.

        %--------------------------
        % Split based on data gaps
        %--------------------------
        % NOTE: Empirically, this can lead to a lot of splitting for
        % LFR-SBM2-CWF.
        % Ex: solo_L1R_rpw-lfr-sbm2-cwf-e-cdag_20220922T134335-20220922T154536_V01.cdf"
        iSegmentAr = bicas.proc.utils.find_data_gaps(...
          Zv.tt2000, Zv.freqHz, bicas.const.MAX_SAMPLE_GAP_RATIO);
        settingsCa{end+1} = iSegmentAr;

        iGroupArCa = bicas.utils.group_by_change(settingsCa{:});
      end
      nGroups = numel(iGroupArCa);
      TmkGrouping.stop_log(nRecords, 'record', nGroups, 'subsequence-group')



      %====================
      % CALIBRATE VOLTAGES
      %====================
      voltageAvolt = nan(nRecords, aspr);

      TmkCal = bicas.utils.Timekeeper(...
        'bicas.proc.L1L2.dc.calibrate_voltage_1xBLTS:calibration', Cv.L);
      for iGroup = 1:nGroups
        iGroupAr = iGroupArCa{iGroup};
        iRec1    = iGroupAr(1);

        if 0
          % DEBUG (only? make permanent?)
          % NOTE: There may be multiple isAchg values within a
          % group/subsequence when it is not used for grouping.
          iRec2     = iGroupAr(end);
          NbriFpa   = Zv.NbriFpa(  iRec1);
          NbciFpa   = Zv.NbciFpa(  iRec1);
          isAchgFpa = Zv.isAchgFpa(iRec1);
          Cv.L.logf('debug', ...
            ['Calibrating group %s -- %s (%5i-%5i: %5i): ', ...
            'iBlts=%i, ssid=%i, freqHz=%5d, isAchg=%3d, nbri=%d, nbci=%d'], ...
            bicas.utils.TT2000_to_UTC_str(Zv.tt2000(iRec1), 0), ...
            bicas.utils.TT2000_to_UTC_str(Zv.tt2000(iRec2), 0), ...
            iRec1, iRec2, iRec2-iRec1, ...
            Cv.iBlts, Zv.ssid(iRec1), Zv.freqHz(iRec1), isAchgFpa.logical2doubleNan(), ...
            NbriFpa.int2doubleNan(), ...
            NbciFpa.int2doubleNan())
        end

        voltageAvolt(iGroupAr, :) = bicas.proc.L1L2.dc.calibrate_voltage_1xBLTS_subsequence(...
          ... % ===============================================================
          ... % Variables which DO NOT VARY over CDF records at all.
          Vcal         = Cv.Vcal, ...
          iBlts        = Cv.iBlts, ...
          hasSwfFormat = Cv.hasSwfFormat, ...
          isLfr        = Cv.isLfr, ...
          isTdsCwf     = Cv.isTdsCwf, ...
          ... % ===============================================================
          ... % Variables which DO NOT VARY over the subsequence/group.
          ssid         = Zv.ssid(     iRec1), ...
          freqHz       = Zv.freqHz(   iRec1), ...
          isAchgFpa    = Zv.isAchgFpa(iRec1), ...
          iLsf         = Zv.iLsf(     iRec1), ...
          NbriFpa      = Zv.NbriFpa(  iRec1), ...
          NbciFpa      = Zv.NbciFpa(  iRec1), ...
          iCalibL      = iCalibL(     iRec1), ...
          iCalibH      = iCalibH(     iRec1), ...
          ... % ===============================================================
          ... % Variables which VARY within the subsequence/group.
          tt2000       = Zv.tt2000(   iGroupAr), ...
          voltageTm    = Zv.voltageTm(iGroupAr, :), ...
          uspr         = Zv.uspr(     iGroupAr) ...
          );
      end
      TmkCal.stop_log(nRecords, 'record', nGroups, 'subsequence-group')

    end    % calibrate_voltage_1xBLTS



    % Calibrate one BLTS channel.
    function voltageAvolt = calibrate_voltage_1xBLTS_subsequence(Cv, Zv)
      arguments
        Cv.Vcal
        Cv.iBlts
        Cv.ssid
        Cv.hasSwfFormat
        Cv.isAchgFpa
        Cv.iCalibL
        Cv.iCalibH
        Cv.freqHz
        Cv.iLsf
        Cv.isLfr
        Cv.isTdsCwf
        Cv.NbriFpa
        Cv.NbciFpa
        Zv.tt2000
        Zv.voltageTm
        Zv.uspr
      end
      % IMPLEMENTATION NOTE: It is ugly to have this many parameters (15!),
      % but the original code made calibrate_voltage_5xBLTS() to large and
      % unwieldy. Having many arguments also highlights the exact dependencies.

      % PROPOSAL: CalSettings as parameter.
      %   PRO: Reduces number of parameters.
      %   PROPOSAL: Add values to CalSettings: isLfr, isTdsCwf, nbri, nbci
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
      assert(isa(Cv.Vcal, "bicas.proc.L1L2.cal.VoltageCalibration"))
      if Cv.isLfr
        assert(isa(Zv.voltageTm, 'single'))
      else
        assert(isa(Zv.voltageTm, 'double'))
      end
      [nRecords, aspr] = irf.assert.sizes(Zv.voltageTm, [-1, -2, 1]);
      assert(isscalar(Cv.freqHz))
      assert(isscalar(Cv.iLsf))



      %--------------
      % Derive dtSec
      %--------------
      % NOTE: Different interpretation of sampling rate for CWF (samples
      % between records) and SWF (samples within snapshot). Always one dtSec
      % value per cell element in voltageTmCa.
      if Cv.hasSwfFormat
        % CASE: SWF
        % NOTE: Time difference between samples within CDF record (snapshot).
        dtSec = 1 / Cv.freqHz;
      else
        % CASE: CWF
        % NOTE: Time difference between CDF records.
        dtSec = 1 / Cv.freqHz;
      end



      if isequaln(Cv.ssid, bicas.proc.L1L2.const.C.SSID_DICT("UNKNOWN"))
        %######################
        % CASE: SSID = UNKNOWN
        %######################
        % ==> Calibrated data set to NaN.

        voltageAvolt = nan(size(Zv.voltageTm));

      elseif ismember(...
          Cv.ssid, bicas.proc.L1L2.const.C.SSID_DICT(["GND", "REF25V"]))
        %##############################
        % CASE: SSID = GND or REF 2.5V
        %##############################
        % ==> No calibration.

        % NOTE: voltageTm stores TM units using float!
        voltageAvolt = Zv.voltageTm;

      else
        %###########################
        % CASE: Actual science data
        %###########################
        assert(bicas.proc.L1L2.const.SSID_is_ASR(Cv.ssid))
        % ==> Calibrate (unless explicitly stated that should not)

        if Cv.hasSwfFormat
          voltageTmCa = ...
            bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(...
            double(Zv.voltageTm), Zv.uspr);
        else
          assert(all(Zv.uspr == 1))
          voltageTmCa = {double(Zv.voltageTm)};
        end

        %######################
        %######################
        %  CALIBRATE VOLTAGES
        %######################
        %######################
        CalSettings = bicas.proc.L1L2.cal.CalibrationSettings(...
          Cv.iBlts, Cv.ssid, Cv.isAchgFpa.logical2doubleNan(), ...
          Cv.iCalibL, Cv.iCalibH, Cv.iLsf);
        %#######################################################
        voltageAvoltCa = Cv.Vcal.calibrate_voltage_TM_to_avolt(...
          dtSec, voltageTmCa, ...
          Cv.isLfr, Cv.isTdsCwf, CalSettings, ...
          Cv.NbriFpa, Cv.NbciFpa);
        %#######################################################

        if Cv.hasSwfFormat
          % CASE: SWF
          [voltageAvolt, ~] = ...
            bicas.proc.utils.convert_cell_array_of_vectors_to_matrix(...
            voltageAvoltCa, aspr...
            );
        else
          % CASE: CWF
          % NOTE: Scalar CA, since not snapshot.
          assert(isscalar(voltageAvoltCa))

          voltageAvolt = voltageAvoltCa{1};
        end
      end
    end    % calibrate_voltage_1xBLTS_subsequence



    % Derive VSIB from BLTS samples. Vectorized.
    %
    % RETURN VALUE
    % ============
    % vsibAr
    %       N x 5. SWF: Set if at least one bit is set for any sample within
    %       a snapshot.
    %
    function vsibAr = get_VSIB_5xBLTS(...
        samplesAvoltAr, hasSwfFormat, uspr, ssidAr, isAchgFpa, SatSettings, L)

      Tmk = bicas.utils.Timekeeper('get_VSIB_5xBLTS', L);

      [nRec, aspr] = irf.assert.sizes(...
        samplesAvoltAr, [-1, -2, bicas.const.N_BLTS], ...
        uspr,           [-1], ...
        ssidAr,         [-1,     bicas.const.N_BLTS], ...
        isAchgFpa,      [-1]);
      assert(isscalar(hasSwfFormat) && islogical(hasSwfFormat))
      assert(bicas.proc.L1L2.const.is_SSID(ssidAr))

      % Expand variables to be of the same size as bltsSamplesAvoltAr
      % -------------------------------------------------------------
      % Needed for submitting arguments to bicas.proc.L1L2.sat.get_VSTB()
      % NOTE: This could possibly lead to memory problems, which could be
      % mitigated by e.g. calling bicas.proc.L1L2.sat.get_VSTB()
      % once per BLTS.
      isAchgFpa = repmat(        isAchgFpa,          [1, aspr, bicas.const.N_BLTS]);
      ssidAr    = repmat(permute(ssidAr, [1, 3, 2]), [1, aspr, 1                 ]);

      % Expand variable to size nRec x 1 x N_BLTS, needed for later comparison
      uspr        = repmat(uspr, [1, 1, 5]);



      vstbAr = bicas.proc.L1L2.sat.get_VSTB(...
        SatSettings, samplesAvoltAr, ssidAr, isAchgFpa);

      % Normalize CWF/SWF data to one array format (one VSIB per record & BLTS).
      % N x M x 5 --> N x 1 x 5 --> N x 5
      % VSTB --> VSIB
      if hasSwfFormat
        % CASE: SWF

        % POTENTIAL BUG: Relies on that VSTB=false for padded samples/elements
        % at the end of TDS-RSWF snapshots (the arrays are padded not account
        % for varying-length snapshots).
        vsibAr = (sum(vstbAr, 2) ./ uspr) >= SatSettings.vstbFractionThreshold;
      else
        % CASE: CWF
        vsibAr = vstbAr;
      end
      vsibAr = permute(vsibAr, [1, 3, 2]);

      irf.assert.sizes(vsibAr, [nRec, bicas.const.N_BLTS])

      Tmk.stop_log(nRec, 'record')
    end



    % Convert samples stored as 5x BLTSs to 9x ASRs (without reconstructing
    % missing data).
    %
    function SchdZvm = convert_voltage_5xBLTS_to_9xASR( ...
        bltsVoltageAvoltAr, bltsVsibAr, bltsSdidAr, L)

      % Tmk = bicas.utils.Timekeeper(...
      %   'bicas.proc.L1L2.dc.convert_voltage_5xBLTS_to_9xASR', L);

      [nRec, aspr] = irf.assert.sizes(...
        bltsVoltageAvoltAr, [-1, -2, bicas.const.N_BLTS], ...
        bltsSdidAr,         [-1,     bicas.const.N_BLTS], ...
        bltsVsibAr,         [-1,     bicas.const.N_BLTS]);

      %====================================
      % Convert BLTS arrays --> bltsSchdCa
      %====================================
      bltsSchdCa = cell(5, 1);
      for iBlts = 1:bicas.const.N_BLTS
        bltsSchdCa{iBlts, 1} = bicas.proc.L1L2.SingleChannelData( ...
          bltsVoltageAvoltAr(:, :, iBlts), ...
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
