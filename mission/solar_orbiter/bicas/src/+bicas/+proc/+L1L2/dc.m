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
  % PROPOSAL: More automatic test code.
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
  %   * Not require (as much) splitting up CDF records into subsequences with
  %     constant "settings" (values for specific zVariables not varying as
  %     a function of CDF record), in particular not require many constant
  %     "settings".
  %   * Use more vector operations.
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
    %   (1) calibrate samples
    %   (2) demux (demultiplex): Relabel samples from BLTS to SDID.
    %   (3) reconstruct (derive) samples for missing channels from calibrated
    %       samples (e.g. DC_V12 := DC_V1 - DC_V2)
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
      currentAAmpere = bicas.proc.L1L2.cur.calibrate_bias_currents(...
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
      aspr          = size(Dcip.Zv.bltsSamplesTm, 2);
      btlsSsidAr2   = repmat(permute(bltsSsidArray, [1 3 2]), [1, aspr, 1]);
      bltsSamplesTm = bicas.proc.L1L2.qrc.set_5xBLTS_voltage_samples_FV(...
        Dcip.Zv.bltsSamplesTm, btlsSsidAr2, L2Qrcbm, bicas.const.qrc.Q.L2_QRCSM);



      %##################################################################
      % CALIBRATE VOLTAGES: 5x CHANNELS LABELLED BY SSID/BLTS (not SDID)
      %##################################################################
      % NOTE: Takes most of the time LFR-SWF.
      bltsSamplesAVolt = bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS(...
        tt2000       = Dcip.Zv.Epoch, ...
        samplesTm    = bltsSamplesTm, ...              % Blanked by QRCs.
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
        bltsSamplesAVolt, Dcip.hasSwfFormat, Dcip.Zv.uspr, ...
        bltsSsidArray, Dcip.Zv.isAchgFpa, SatSettings, L);



      %######################################################################
      % ~"DEMUX" VOLTAGES:
      % INPUT:  5x SIGNALS LABELLED BY SSID/BLTS
      % OUTPUT: 9x SIGNALS LABELLED BY SDID + RECONSTRUCTING MISSING SIGNALS
      %######################################################################
      % NOTE: Needs VSIB for propagating VSIB to reconstructed channels.
      SchdZvm = bicas.proc.L1L2.dc.convert_samples_5xBLTS_to_9xASR(...
        bltsSamplesAVolt, bltsVsibAr, bltsSdidArray, L);
      bicas.proc.L1L2.demuxer.reconstruct_ASR_samples(SchdZvm);



      % Replace/split variable: SchdZvm --> SamplesZvm + VsibZvm
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

      % NOTE: Function modifies SamplesZvm handle object in-place!
      Zv.currentAAmpere     = bicas.proc.L1L2.qrc.set_current_samples_FV(...
        currentAAmpere, L2Qrcbm, bicas.const.qrc.Q.L2_QRCSM);
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



    function log_input_calibration_settings(Dcip, Vcal, L)
      % IMPLEMENTATION NOTE: Implemented separately from processing functions
      % since:
      % (1) removes dependence on logger object,
      % (2) can split data into subsequences based on arbitrary choice of
      %     variables,
      % (3) reduces size of processing function where logging would otherwise
      %     be, and
      % (4) can potentially turn table into proper table with column headers.

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
    function samplesAVolt = calibrate_voltages_5xBLTS(Cv, Zv)
      arguments
        % Variables which DO NOT VARY over CDF records at all.
        Cv.Vcal
        Cv.L
        Cv.hasSwfFormat
        Cv.isTdsCwf
        Cv.isLfr
        % Variables which VARY over CDF records.
        Zv.tt2000
        Zv.samplesTm
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
      assert(isnumeric(Zv.samplesTm))
      assert(isa(Zv.ssid, 'uint8'))
      [nRecords, aspr] = irf.assert.sizes(...
        Zv.isAchgFpa, [-1], ...
        Zv.ssid,      [-1,     bicas.const.N_BLTS], ...
        Zv.samplesTm, [-1, -2, bicas.const.N_BLTS]);



      % Pre-allocate array of calibrated voltage samples: All (BLTS) channels
      % ---------------------------------------------------------------------
      % IMPLEMENTATION NOTE: Preallocation is very important for speeding up
      % LFR-SWF which tends to be broken into subsequences of 1 record.
      samplesAVolt = nan(nRecords, aspr, bicas.const.N_BLTS);



      %===========
      % Calibrate
      %===========
      Tmk = bicas.utils.Timekeeper(...
        'bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS:Calibrating voltages', Cv.L);
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

        bltsSamplesAVolt2 = bicas.proc.L1L2.dc.calibrate_voltages_1xBLTS( ...
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
          samplesTm    = Zv.samplesTm(:, :, iBlts), ...
          uspr         = Zv.uspr, ...
          freqHz       = Zv.freqHz, ...
          iLsf         = Zv.iLsf, ...
          isAchgFpa    = Zv.isAchgFpa, ...
          NbriFpa      = Zv.NbriFpa, ...
          NbciFpa      = Zv.NbciFpa);

        % Add subsequence/group signals to the global array (all records).
        samplesAVolt(:, :, iBlts) = bltsSamplesAVolt2;
      end    % for
      Tmk.stop_log(nRecords, 'record')

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
    function samplesAVolt = calibrate_voltages_1xBLTS(Cv, Zv)
      % PROPOSAL: Also split sequence based constant isnan() (isfinite()?)
      %           for CWF (not SWF).
      %   NOTE: bicas.tf.apply_TF() can split based on isfinite() (not isnan()).

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
        Zv.samplesTm
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
        Zv.samplesTm, [-1, -2, 1], ...
        Zv.uspr,      [-1]);



      iCalibL = Cv.Vcal.get_BIAS_calibration_time_index_L(Zv.tt2000);
      iCalibH = Cv.Vcal.get_BIAS_calibration_time_index_H(Zv.tt2000);

      %==================================================================
      % (1) Find groups/subsequences of records with identical settings.
      % (2) Process data separately for each such sequence.
      %==================================================================
      Tmk = bicas.utils.Timekeeper(...
        'bicas.proc.L1L2.dc.calibrate_voltages_1xBLTS:grouping algo.', Cv.L);
      % IMPLEMENTATION NOTE: Important to convert FPAs to builtin arrays to
      % significantly speed up both bicas.utils.group_unique_rows(), and
      % bicas.utils.group_by_change().
      settingsCa = {...
        Zv.isAchgFpa.logical2doubleNan(), ...
        Zv.freqHz, ...
        Zv.iLsf, ...
        Zv.NbriFpa.int2doubleNan(), ...
        Zv.NbciFpa.int2doubleNan(), ...
        Zv.ssid, ...
        iCalibL, ...
        iCalibH};
      if Cv.hasSwfFormat
        % IMPLEMENTATION NOTE: SWF data is calibrated in a way where
        % consecutive records (snapshots) do NOT affect each other. One can
        % therefore group SWF CDF records in groups of non-consecutive CDF
        % records. This is important for LFR-SWF data which tends to change
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
        % IMPLEMENTATION NOTE: CWF data is calibrated in a way where
        % consecutive records affect each other. Must therefore divide CDF
        % records in groups (subsequences) of continuous CDF records.
        iGroupArCa = bicas.utils.group_by_change(settingsCa{:});
      end
      nGroups = numel(iGroupArCa);
      Tmk.stop_log(nRecords, 'record', nGroups, 'subsequence-group')

      %====================
      % CALIBRATE VOLTAGES
      %====================
      samplesAVolt = nan(nRecords, aspr);

      for iGroup = 1:nGroups
        iGroupAr = iGroupArCa{iGroup};
        iRec1    = iGroupAr(1);

        samplesAVolt(iGroupAr, :) = bicas.proc.L1L2.dc.calibrate_1xBLTS_subsequence(...
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
          samplesTm    = Zv.samplesTm(iGroupAr, :), ...
          uspr         = Zv.uspr(     iGroupAr) ...
          );
      end
    end    % calibrate_voltages_1xBLTS



    % Calibrate one BLTS channel.
    function samplesAVolt = calibrate_1xBLTS_subsequence(Cv, Zv)
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
        Zv.samplesTm
        Zv.uspr
      end
      % IMPLEMENTATION NOTE: It is ugly to have this many parameters (15!),
      % but the original code made calibrate_voltages_5xBLTS() to large and
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
        assert(isa(Zv.samplesTm, 'single'))
      else
        assert(isa(Zv.samplesTm, 'double'))
      end
      [nRecords, aspr] = irf.assert.sizes(Zv.samplesTm, [-1, -2, 1]);
      assert(isscalar(Cv.freqHz))
      assert(isscalar(Cv.iLsf))



      % Derive dtSec
      % ------------
      % NOTE: Different number of rows depending on CWF/SWF! One element per
      %       cell element in samplesTmCa.
      if Cv.hasSwfFormat
        % CASE: SWF
        % NOTE: Column vector of (identical) numbers (one per snapshot, which
        %       will all be processed separately).
        dtSec = ones(nRecords, 1) ./ Cv.freqHz;
      else
        % CASE: CWF
        % NOTE: Scalar (since all data will be processed in one session).
        % BUG RISK: Has the possibility of data gaps been excluded at this point
        % in the code, since the calculation seems to rely on it?
        dtSec = double( Zv.tt2000(end) - Zv.tt2000(1) ) / (nRecords-1) * 1e-9;
      end



      if isequaln(Cv.ssid, bicas.proc.L1L2.const.C.SSID_DICT("UNKNOWN"))
        % ==> Calibrated data set to NaN.
        samplesAVolt = nan(size(Zv.samplesTm));

      elseif ismember(...
          Cv.ssid, bicas.proc.L1L2.const.C.SSID_DICT(["GND", "REF25V"]))
        % ==> No calibration.
        % NOTE: samplesTm stores TM units using float!
        samplesAVolt = Zv.samplesTm;

      else
        assert(bicas.proc.L1L2.const.SSID_is_ASR(Cv.ssid))
        % ==> Calibrate (unless explicitly stated that should not)

        if Cv.hasSwfFormat
          samplesTmCa = ...
            bicas.proc.utils.convert_matrix_to_cell_array_of_vectors(...
            double(Zv.samplesTm), Zv.uspr);
        else
          assert(all(Zv.uspr == 1))
          samplesTmCa = {double(Zv.samplesTm)};
        end

        %######################
        %######################
        %  CALIBRATE VOLTAGES
        %######################
        %######################
        CalSettings = bicas.proc.L1L2.CalibrationSettings(...
          Cv.iBlts, Cv.ssid, Cv.isAchgFpa.logical2doubleNan(), ...
          Cv.iCalibL, Cv.iCalibH, Cv.iLsf);
        %#######################################################
        samplesAVoltCa = Cv.Vcal.calibrate_voltage_TM_to_avolt(...
          dtSec, samplesTmCa, ...
          Cv.isLfr, Cv.isTdsCwf, CalSettings, ...
          Cv.NbriFpa, Cv.NbciFpa);
        %#######################################################

        if Cv.hasSwfFormat
          % CASE: SWF
          [samplesAVolt, ~] = ...
            bicas.proc.utils.convert_cell_array_of_vectors_to_matrix(...
            samplesAVoltCa, aspr...
            );
        else
          % CASE: CWF
          % NOTE: Scalar, since not snapshot.
          assert(isscalar(samplesAVoltCa))

          samplesAVolt = samplesAVoltCa{1};
        end
      end
    end    % calibrate_1xBLTS_subsequence



    % Derive VSIB from BLTS samples. Vectorized.
    %
    % RETURN VALUE
    % ============
    % vsibAr
    %       N x 5. SWF: Set if at least one bit is set for any sample within
    %       a snapshot.
    %
    function vsibAr = get_VSIB_5xBLTS(...
        samplesAVoltAr, hasSwfFormat, uspr, ssidAr, isAchgFpa, SatSettings, L)

      Tmk = bicas.utils.Timekeeper('get_VSIB_5xBLTS', L);

      [nRec, aspr] = irf.assert.sizes(...
        samplesAVoltAr, [-1, -2, bicas.const.N_BLTS], ...
        uspr,           [-1], ...
        ssidAr,         [-1,     bicas.const.N_BLTS], ...
        isAchgFpa,      [-1]);
      assert(isscalar(hasSwfFormat) && islogical(hasSwfFormat))
      assert(bicas.proc.L1L2.const.is_SSID(ssidAr))

      % Expand variables to be of the same size as bltsSamplesAVoltAr
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
        SatSettings, samplesAVoltAr, ssidAr, isAchgFpa);

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
    function SchdZvm = convert_samples_5xBLTS_to_9xASR( ...
        bltsSamplesAVoltAr, bltsVsibAr, bltsSdidAr, L)

      % Tmk = bicas.utils.Timekeeper(...
      %   'bicas.proc.L1L2.dc.convert_samples_5xBLTS_to_9xASR', L);

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
