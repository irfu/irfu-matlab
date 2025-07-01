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
%       CON: Minor. Quality variables are still set somewhere in this file.
%   --
%   process()
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
        nValidSamplesPerRecord  = Dcip.Zv.nValidSamplesPerRecord, ...
        Cal                     = Cal, ...
        L                       = L);



      %######################################################################
      % ~"DEMUX" VOLTAGES:
      % INPUT:  5x SIGNALS LABELLED BY SDID/BLTS (not SSID)
      % OUTPUT: 9x SIGNALS LABELLED BY SDID + RECONSTRUCTING MISSING SIGNALS
      %######################################################################
      AsrSamplesAVoltSrm = bicas.proc.L1L2.dc.relabel_reconstruct_samples_5xBLTS_to_9xASR(...
        bltsSamplesAVolt, bltsSdidArray, L);

      if 1
        % ############
        % EXPERIMENTAL
        % ############
        TmkNewImpl = bicas.utils.Timekeeper(...
          'bicas.proc.L1L2.dc.process_calibrate_demux: NEW IMPLEMENTATION', L);

        % 5x SIGNALS LABELLED BY SSID/BLTS.
        % NOTE: Incomplete detection of VSQB.
        % NOTE: No SCHD as input (though as output).
        % NOTE: Is QUITE SLOW, at least for LFR SWF (?).
        bltsVsibAr = bicas.proc.L1L2.dc.get_VSIB_5xBLTS_NEW(...
          Bso, bltsSamplesAVolt, bltsSsidArray, Dcip.Zv.isAchgFpa, L);

        Cdac = bicas.proc.L1L2.dc.convert_samples_5xBLTS_to_9xASR_NEW(...
          bltsSamplesAVolt, bltsVsibAr, bltsSdidArray, L);

        Cdac        = bicas.proc.L1L2.demuxer.reconstruct_ASR_samples_NEW(Cdac);

        SatSettings = bicas.proc.L1L2.sat.from_BSO_extract_saturation_settings(Bso);
        [QUALITY_FLAG, L2_QUALITY_BITMASK] = bicas.proc.L1L2.qual.get_quality_ZVs_channel_saturation(...
          Cdac, Dcip.Zv.Epoch, Dcip.hasSwfFormat, ...
          SatSettings.vsibFractionThreshold, SatSettings.cwfSlidingWindowLengthSec);

        TmkNewImpl.stop_log()

        % TODO: Extract Cdac VSIBs and set L2_QUALITY_BITMASK.
        %   Separate function for extracting these bits specifically.
        %   TODO-DEC: Which file should function be in?! qual, sat?
        % TODO: Convert Cdac samples to AsrSamplesAVoltSrm (or at least use
        %       it).

        % PROPOSAL: Compare SDID-separated VSIBs combined into one with old
        %           VSQBs.

        % ======================================================================
        % DEBUG: Check that samples derived using OLD and NEW code are identical
        % ======================================================================
        % NOTE: Check may fail if new code splits samples into subsequences in
        %       new way, e.g. BLTS per BLTS.
        % PROPOSAL: Permit approximative equality.
        % PROPOSAL: Store max difference.

        % maxDiff = 0;
        for asrSdid = bicas.proc.L1L2.const.C.SDID_ASR_AR'
          asid = bicas.proc.L1L2.const.C.SDID_ASID_DICT(asrSdid);
          oldImplSamplesAVolt = AsrSamplesAVoltSrm(asid);               % Samples from OLD impl.
          newImplSamplesAvolt = Cdac.get_channel(asrSdid).samplesAr;    % Samples from NEW impl.

          % NOTE: Treat NaN as equal itself.
          if ~isequaln(oldImplSamplesAVolt, newImplSamplesAvolt)
            L.logf('debug', 'Samples are different for asid = %g', asid)
            if ~isequal(isnan(oldImplSamplesAVolt), isnan(newImplSamplesAvolt))
              L.logf('debug', 'isnan() is different.')
            end
            error('Samples differ.')
          end

          % NOTE: max() ignores NaN.
          % d = max(abs(A - B), [], 'ALL');
          % maxDiff = max(maxDiff, d);
        end
        % if logical(maxDiff)
        %   L.logf('debug', 'maxDiff = %g', maxDiff)
        % end
      end    % if EXPERIMENTAL



      %######################################################
      % ~Derive quality variables, and UFV from quality QRCs
      %######################################################
      % AUTODETECT SATURATION.
      % NOTE: Derives *ONE* saturation bit for *ALL* signals combined (per CDF
      %       record).
      autodetectedVsqb = bicas.proc.L1L2.sat.get_VSQB(...
        Bso, ...
        Dcip.Zv.Epoch, ...
        AsrSamplesAVoltSrm, ...
        Dcip.Zv.nValidSamplesPerRecord, ...
        bltsSsidArray, ...
        Dcip.Zv.isAchgFpa, ...
        Dcip.hasSwfFormat, L);

      ZvIn = struct(...
        'Epoch',            Dcip.Zv.Epoch, ...
        'bdmFpa',           Dcip.Zv.bdmFpa, ...
        'autodetectedVsqb', autodetectedVsqb);
      [zvUfv, QUALITY_FLAG_max, L2_QUALITY_BITMASK] = ...
        bicas.proc.L1L2.qual.get_UFV_quality_ZVs(...
        ZvIn, Dcip.isLfr, NsoTable, Bso, L);
      clear ZvIn



      %#############################
      % Set UFV, "final" zVariables
      %#############################
      Zv = struct();
      zvUfv                 = Dcip.Zv.ufv | zvUfv;
      Zv.QUALITY_FLAG       = Dcip.Zv.QUALITY_FLAG.min(...
        bicas.utils.FPArray(QUALITY_FLAG_max));
      Zv.L2_QUALITY_BITMASK = L2_QUALITY_BITMASK;

      % NOTE: Function modifies AsrSamplesAVoltSrm handle object!
      Zv.currentAAmpere = bicas.proc.L1L2.qual.set_voltage_current_FV(...
        Dcip.Zv.Epoch, AsrSamplesAVoltSrm, currentAAmpere, zvUfv, L);
      Zv.AsrSamplesAVoltSrm = AsrSamplesAVoltSrm;
      clear AsrSamplesAVoltSrm



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
      % PROPOSAL: Integrate into bicas.proc.L1L2.demuxer (as method).
      % NOTE: Calibration is really separate from the demultiplexer.
      %       Demultiplexer only needs to split into subsequences based on BDM
      %       and DLR, nothing else.
      %   PROPOSAL: Separate out demultiplexer. Do not call from this function.
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
      %   PRO: Can potentially have different sets subsequences for different
      %        BLTSs.
      %     Ex: ACHG, SSID
      %   PRO: Probably simpler tests(?)
      %   PRO: Easier to implement/support separate calibration functions for
      %        different channels?!
      %   CON: Splits into subsequences multiple times.
      %     PRO: Potentially slower (in particular SWF).

      arguments
        Zv.Epoch
        Zv.bltsSamplesTm
        Zv.isAchgFpa
        Zv.CALIBRATION_TABLE_INDEX
        Zv.freqHz
        Zv.iLsf
        Zv.ufv
        Zv.bltsSsidArray
        Zv.nValidSamplesPerRecord
        A.isTdsCwf
        A.isLfr
        A.hasSwfFormat
        A.Cal
        A.L
      end

      % ASSERTIONS
      assert(isscalar( A.hasSwfFormat))
      assert(isnumeric(Zv.bltsSamplesTm))
      assert(isa(Zv.bltsSsidArray, 'uint8'))
      [nRecords, nSamplesPerRecordChannel] = irf.assert.sizes(...
        Zv.isAchgFpa,     [-1], ...
        Zv.bltsSsidArray, [-1,     bicas.const.N_BLTS], ...
        Zv.bltsSamplesTm, [-1, -2, bicas.const.N_BLTS]);



      % Pre-allocate calibrated array. All (BLTS) channels, all records
      % ---------------------------------------------------------------
      % IMPLEMENTATION NOTE: Preallocation is very important for speeding
      % up LFR-SWF which tends to be broken into subsequences of 1 record.
      bltsSamplesAVolt = nan(nRecords, nSamplesPerRecordChannel, bicas.const.N_BLTS);

      iCalibLZv = A.Cal.get_BIAS_calibration_time_index_L(Zv.Epoch);
      iCalibHZv = A.Cal.get_BIAS_calibration_time_index_H(Zv.Epoch);



      %========================================================================
      % (1) Find continuous subsequences of records with identical settings.
      % (2) Process data separately for each such sequence.
      % --------------------------------------------------------------------
      % SS = Subsequence
      %
      % NOTE: Empirically, splitting by changing settings is not really useful
      %       for real LFR SWF datasets where the LSF changes in every record
      %       anyway, meaning that the subsequences are all 1 record long.
      %       Suspect that this implies that the indexing of arrays becomes a
      %       performance problem since removing bltsSsidArray argument (which
      %       was mistakenly added despite being unused) speed up BICAS
      %       significantly for SWF data: about half of the time added after
      %       adding bltsSsidArray+bltsSdidArray.
      %       /Erik P G Johansson, 2024-10-29
      %========================================================================
      % PROPOSAL: Do not use irf.utils.split_by_change() for SWF data. It is
      %           enough to group by identical values (not use continuous
      %           blocks of CDF records).
      %   PRO: Faster
      %   CON: Uses knowledge of how the calibration works: that CDF records
      %        are calibrated separately.
      %   PROPOSAL: Create similar function which finds groups of unique
      %             combinations.
      %           [iMembersCa, iHeadsCa] = group_unique_combinations(...)
      %     PROPOSAL: Implement using recursion: Find unique values (rows) for
      %               one argument at a time.
      %========================================================================
      %Tmk = bicas.utils.Timekeeper('bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS:irf.utils.split_by_change', A.L);
      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(...
        Zv.isAchgFpa.logical2doubleNan(), ...
        Zv.freqHz, ...
        Zv.iLsf, ...
        Zv.CALIBRATION_TABLE_INDEX, ...
        Zv.ufv, ...
        Zv.bltsSsidArray, ...
        iCalibLZv, ...
        iCalibHZv);
      %Tmk.stop_log(nRecords, 'record', nSs, 'subsequence')

      Tmk = bicas.utils.Timekeeper('bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS:Calibrating voltages', A.L);
      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        ssBltsSamplesAVolt = bicas.proc.L1L2.dc.calibrate_voltages_5xBLTS_subsequence( ...
          Cal                      = A.Cal, ...
          ... % ===============================================================
          ... % NOTE: Variables which do VARY over CDF records.
          isAchgFpa                = Zv.isAchgFpa(              iRec1), ...
          freqHz                   = Zv.freqHz(                 iRec1), ...
          iLsf                     = Zv.iLsf(                   iRec1), ...
          zvcti                    = Zv.CALIBRATION_TABLE_INDEX(iRec1, :), ...
          ufv                      = Zv.ufv(                    iRec1), ...
          bltsSsidArray            = Zv.bltsSsidArray(          iRec1, :), ...
          iCalibL                  = iCalibLZv(                 iRec1), ...
          iCalibH                  = iCalibHZv(                 iRec1), ...
          ... % ===============================================================
          ... % NOTE: Variables which do not vary over CDF records.
          hasSwfFormat             = A.hasSwfFormat, ...
          isLfr                    = A.isLfr, ...
          isTdsCwf                 = A.isTdsCwf, ...
          nSamplesPerRecordChannel = nSamplesPerRecordChannel, ...
          ...   % Variables which vary by CDF records.
          Epoch                    = Zv.Epoch(                 iRec1:iRec2), ...
          bltsSamplesTm            = Zv.bltsSamplesTm(         iRec1:iRec2, :, :), ...
          zvNValidSamplesPerRecord = Zv.nValidSamplesPerRecord(iRec1:iRec2));

        % Add subsequence signals to the global array (all records).
        bltsSamplesAVolt(iRec1:iRec2, :, :) = ssBltsSamplesAVolt;
      end    % for
      Tmk.stop_log(nRecords, 'record', nSs, 'subsequence')

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
        % NOTE: Below variables do not vary over CDF records at all.
        Cv.nSamplesPerRecordChannel
        Cv.hasSwfFormat
        Cv.isLfr
        Cv.isTdsCwf
        %
        Vv.Epoch
        Vv.bltsSamplesTm
        Vv.zvNValidSamplesPerRecord
      end
      irf.assert.sizes( ...
        Vv.Epoch,                    [-1], ...
        Vv.bltsSamplesTm,            [-1, Cv.nSamplesPerRecordChannel, bicas.const.N_BLTS], ...
        Vv.zvNValidSamplesPerRecord, [-1])

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
      ssBltsSamplesAVolt = nan(nRows, Cv.nSamplesPerRecordChannel, bicas.const.N_BLTS);
      for iBlts = 1:bicas.const.N_BLTS
        ssBltsSamplesAVolt(:, :, iBlts) = bicas.proc.L1L2.dc.calibrate_1xBLTS_subsequence(...
          samplesTm                = Vv.bltsSamplesTm(:, :, iBlts), ...
          zvNValidSamplesPerRecord = Vv.zvNValidSamplesPerRecord, ...
          ssid                     = Cv.bltsSsidArray(iBlts), ...
          iBlts                    = iBlts, ...
          hasSwfFormat             = Cv.hasSwfFormat, ...
          isAchg                   = Cv.isAchgFpa.logical2doubleNan(), ...
          iCalibL                  = Cv.iCalibL, ...
          iCalibH                  = Cv.iCalibH, ...
          iLsf                     = Cv.iLsf, ...
          dtSec                    = dtSec, ...
          isLfr                    = Cv.isLfr, ...
          isTdsCwf                 = Cv.isTdsCwf, ...
          zvcti                    = Cv.zvcti, ...
          ufv                      = Cv.ufv, ...
          Cal                      = A.Cal);
      end
    end    % calibrate_voltages_5xBLTS_subsequence



    % Calibrate one BLTS channel.
    function samplesAVolt = calibrate_1xBLTS_subsequence(A)
      arguments
        A.ssid
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
      % but the original code made calibrate_voltages_5xBLTS() to large and
      % unwieldy. Having many arguments also highlights the exact dependencies.

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
      irf.assert.sizes(A.samplesTm, [-1, -2])   % One BLTS channel. CWF/SWF.

      if isequaln(A.ssid, bicas.proc.L1L2.const.C.SSID_DICT("UNKNOWN"))
        % ==> Calibrated data set to NaN.
        samplesAVolt = nan(size(A.samplesTm));

      elseif isequaln(A.ssid, bicas.proc.L1L2.const.C.SSID_DICT("GND")) || ...
          isequaln(A.ssid, bicas.proc.L1L2.const.C.SSID_DICT("REF25V"))
        % ==> No calibration.
        % NOTE: samplesTm stores TM units using float!
        samplesAVolt = A.samplesTm;

      else
        assert(bicas.proc.L1L2.const.SSID_is_ASR(A.ssid))
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
          A.iBlts, A.ssid, A.isAchg, A.iCalibL, A.iCalibH, A.iLsf);
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
    end    % calibrate_1xBLTS_subsequence



    function AsrSamplesAVoltSrm = relabel_reconstruct_samples_5xBLTS_to_9xASR(...
        bltsSamplesAvolt, bltsSdidArray, L)
      % PROPOSAL: Automated tests.

      Tmk = bicas.utils.Timekeeper('bicas.proc.L1L2.dc.relabel_reconstruct_samples_5xBLTS_to_9xASR', L);

      [nRecTot, nSamplesPerRecordChannel] = irf.assert.sizes(...
        bltsSamplesAvolt, [-1, -2, bicas.const.N_BLTS], ...
        bltsSdidArray,    [-1,     bicas.const.N_BLTS]);

      % -----------------------------------------------------------------
      % Pre-allocate AsrSamplesAVoltSrm: All (ASID) channels, all records
      % -----------------------------------------------------------------
      % IMPLEMENTATION NOTE: Preallocation is very important for speeding
      % up LFR-SWF which tends to be broken into subsequences of 1 record.
      AsrSamplesAVoltSrm = bicas.utils.SameRowsMap(...
        "uint8", nRecTot, 'CONSTANT', ...
        nan(nRecTot, nSamplesPerRecordChannel), ...
        bicas.proc.L1L2.const.C.ASID_DICT.values);



      [iRec1Ar, iRec2Ar, nSs] = irf.utils.split_by_change(bltsSdidArray);
      for iSs = 1:nSs
        iRec1 = iRec1Ar(iSs);
        iRec2 = iRec2Ar(iSs);

        SsAsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.relabel_reconstruct_samples_5xBLTS_to_9xASR_subsequence(...
          bltsSdidArray(   iRec1,          :), ...
          bltsSamplesAvolt(iRec1:iRec2, :, :));

        % Set demuxed subsequence signals (some records) to the global arrays
        % (all records).
        AsrSamplesAVoltSrm.set_rows(SsAsrSamplesAVoltSrm, [iRec1:iRec2]');
      end

      Tmk.stop_log(nRecTot, 'record', nSs, 'subsequence')
    end



    % EXPERIMENTAL
    %
    % Derive VSIB from BLTS samples. Vectorized.
    %
    % RETURN VALUE
    % ============
    % bltsVsibAr
    %       N x 5. SWF: Set if at least one bit is set for any sample within
    %       a snapshot.
    function bltsVsibAr = get_VSIB_5xBLTS_NEW(...
        Bso, bltsSamplesAVoltAr, bltsSsidAr, isAchgFpa, L)
      % PROPOSAL: Test code.
      % PROPOSAL: SatSettings as argument.
      %   PRO: Simpler test code.

      Tmk = bicas.utils.Timekeeper('get_VSIB_5xBLTS_NEW', L);

      [nSpr, nRec] = irf.assert.sizes(...
        bltsSamplesAVoltAr, [-2, -1, bicas.const.N_BLTS], ...
        bltsSsidAr,         [-2,     bicas.const.N_BLTS], ...
        isAchgFpa,          [-2]);
      assert(bicas.proc.L1L2.const.is_SSID(bltsSsidAr))

      % Expand variables to be of the same size as bltsSamplesAVoltAr
      % -------------------------------------------------------------
      % NOTE: This could possibly lead to memory problems, which could be
      % mitigated by e.g. calling bicas.proc.L1L2.sat.get_threshold_VSIB_NEW
      % once per BLTS.
      isAchgFpa   = repmat(        isAchgFpa,              [1, nSpr, bicas.const.N_BLTS]);
      bltsSsidAr  = repmat(permute(bltsSsidAr, [1, 3, 2]), [1, nSpr, 1                 ]);

      SatSettings = bicas.proc.L1L2.sat.from_BSO_extract_saturation_settings(Bso);

      bltsVsibAr  = bicas.proc.L1L2.sat.get_threshold_VSIB_NEW(...
        SatSettings, bltsSamplesAVoltAr, bltsSsidAr, isAchgFpa);

      % Normalize CWF/SWF data to one array format.
      % N x M x 5 --> N x 1 x 5 --> N x 5
      % Logical OR over all VSIBs within snapshot.
      bltsVsibAr = any(    bltsVsibAr, 2);
      bltsVsibAr = permute(bltsVsibAr, [1, 3, 2]);

      Tmk.stop_log(nRec, 'record')
    end



    % EXPERIMENTAL.
    %
    % Convert samples stored as 5x BLTSs to 9x ASRs (without reconstructing
    % missing data).
    %
    % Intended as future partial conceptual replacement for
    % bicas.proc.L1L2.dc.relabel_reconstruct_samples_5xBLTS_to_9xASR().
    %
    function Cdac = convert_samples_5xBLTS_to_9xASR_NEW( ...
        bltsSamplesAVoltAr, bltsVsibAr, bltsSdidAr, L)

      % Tmk = bicas.utils.Timekeeper(...
      %   'bicas.proc.L1L2.dc.convert_samples_5xBLTS_to_9xASR_NEW', L);

      [nRec, nSamplesPerRecordChannel] = irf.assert.sizes(...
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

      %==================================================
      % Convert bltsSchdCa --> Cdac
      % ---------------------------
      % Copy values from 5x BLTSs into 9x SCHD (1x Cdac)
      % (but no reconstruction of missing values)
      %==================================================
      Cdac = bicas.proc.L1L2.ChannelDataAsrCollection(nRec);
      for asrSdid = bicas.proc.L1L2.const.C.SDID_ASR_AR'

        % ~Preallocate empty SCHD for the current ASR/SDID.
        AsrSchd = bicas.proc.L1L2.SingleChannelData(...
          nan(  nRec, nSamplesPerRecordChannel), ...
          false(nRec, 1));

        for iBlts = 1:bicas.const.N_BLTS
          bRecCopy          = ( bltsSdidAr(:, iBlts) == asrSdid );

          % Copy BLTS samples into selected records/rows in previously created
          % SCHD for ASR.
          BltsSchd          = bltsSchdCa{iBlts};
          AsrSchd(bRecCopy) = BltsSchd(bRecCopy);
        end
        Cdac = Cdac.set_channel(asrSdid, AsrSchd);
      end

      % Tmk.stop_log(nRec, 'record')
    end



  end    % methods(Static, Access=private)



end
