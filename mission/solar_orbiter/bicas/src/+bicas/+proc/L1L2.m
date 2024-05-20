%
% Class that collects miscellaneous functions for processing L1/L1R-->L2.
%
% This class is not meant to be instantiated.
%
%
% CODE CONVENTIONS
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like"
%   data, use the first MATLAB array index to represent CDF records.
%
%
% SOME INTERMEDIATE PROCESSING DATA FORMATS
% =========================================
% - PreDc = Pre-(Demuxing & Calibration) Data
%       Generic data format that can represent all forms of input datasets
%       before demuxing and calibration. Can use an arbitrary number of samples
%       per record. Some variables are therefore not used in CWF output
%       datasets.
% - PostDc = Post-(Demuxing & Calibration) Data
%       Data format that includes calibrated currents & calibrated & demuxed
%       voltages.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-10, with source code from data_manager_old.m.
%
classdef L1L2
  %#######################################################################################################################
  %
  % PROPOSAL: Move normalize_CALIBRATION_TABLE_INDEX() to some collection of utils.
  %   PROPOSAL: bicas.proc.utils
  %       CON: Function is too specific. Has inputDsi as argument.
  %           CON: Could be less bad than this file.
  %
  % PROPOSAL: Submit ZV attributes.
  %   PRO: Can interpret fill values.
  %       Ex: Can doublecheck TDS RSWF snapshot length using fill values and compare with zVar SAMPS_PER_CH (which seems
  %           to be bad).
  %
  % PROPOSAL: Return (to execute_SWM), global attributes.
  %   PRO: Needed for output datasets: CALIBRATION_TABLE, CALIBRATION_VERSION
  %       ~CON: CALIBRATION_VERSION refers to algorithm and should maybe be a SETTING.
  %
  % PROPOSAL: Class for HkSciTime.
  %
  %#######################################################################################################################



  %#############################
  %#############################
  methods(Static, Access=public)
    %#############################
    %#############################



    % Only converts relevant HK ZVs to be on SCI Epoch. Other code later
    % decides whether to actually use it (BDM).
    function HkSciTime = process_HK_CDF_to_HK_on_SCI_TIME(InSci, InHk, Bso, L)
      % PROPOSAL: Separate function for the actual interpolation of data
      %           (changing time array HK-->SCI).
      %
      % NOTE: Function only uses InSci fields InSci.Zv.Epoch and InSci.Zv.ACQUISITION_TIME
      %   PROPOSAL: Replace argument InSci --> SciEpoch and Sci_ACQUISITION_TIME



      % ASSERTIONS
      assert(isa(InSci, 'bicas.InputDataset'))
      assert(isa(InHk,  'bicas.InputDataset'))

      HkSciTime = [];



      %===================================================================
      % Select whether HK should use
      %   (1) (HK) Epoch, or
      %   (2) (HK) ACQUISITION_TIME (not always available).
      % ---------------------------------------------------
      % IMPLEMENTATION NOTE: Historically, there have been datasets where
      % Epoch contains errors, but ACQUISITION_TIME seems OK. This should
      % be phased out eventually.
      %===================================================================
      ACQUISITION_TIME_EPOCH_UTC = Bso.get_fv('INPUT_CDF.ACQUISITION_TIME_EPOCH_UTC');
      USE_ZV_ACQUISITION_TIME_HK = Bso.get_fv('PROCESSING.HK.USE_ZV_ACQUISITION_TIME');
      if USE_ZV_ACQUISITION_TIME_HK
        hkEpoch = bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
          InHk.Zv.ACQUISITION_TIME, ...
          ACQUISITION_TIME_EPOCH_UTC);

        L.logf('warning', 'Using HK zVar ACQUISITION_TIME instead of Epoch.')
      else
        hkEpoch = InHk.Zv.Epoch;
      end



      %==================================================================
      % Log time intervals to enable comparing available SCI and HK data
      %==================================================================
      TimeZvs = [];    % Temporary struct only used for logging.
      TimeZvs.HK_Epoch  = InHk.Zv.Epoch;
      TimeZvs.SCI_Epoch = InSci.Zv.Epoch;
      if isfield(InHk.Zv, 'ACQUISITION_TIME')
        TimeZvs.HK_ACQUISITION_TIME_tt2000 = ...
          bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
          InHk.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
      end
      if isfield(InSci.Zv, 'ACQUISITION_TIME') && ~isempty(InSci.Zv.ACQUISITION_TIME)
        TimeZvs.SCI_ACQUISITION_TIME_tt2000 = ...
          bicas.proc.utils.ACQUISITION_TIME_to_TT2000(...
          InSci.Zv.ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC);
      end
      bicas.utils.log_ZVs(TimeZvs, Bso, L);



      %===================
      % WARNINGS / ERRORS
      %===================
      if ~issorted(hkEpoch, 'strictascend')
        % Ex: zVar ACQUISITION_TIME in test file
        % TDS___TESTDATA_RGTS_TDS_CALBA_V0.8.6/
        % solo_HK_rpw-bia_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
        % is not monotonically increasing (in fact, it is completely
        % strange).
        error(...
          ['HK timestamps do not increase monotonically', ...
          ' (USE_ZV_ACQUISITION_TIME_HK=%g).'], ...
          USE_ZV_ACQUISITION_TIME_HK)

      end
      if ~irf.utils.ranges_intersect(InSci.Zv.Epoch, hkEpoch)
        %---------------------------------------
        % CASE: SCI does not overlap HK in time
        %---------------------------------------

        % NOTE: "WARNING" (rather than error) only makes sense if it is
        % possible to later meaningfully permit non-intersection.
        [settingValue, settingKey] = Bso.get_fv(...
          'PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY');
        bicas.default_anomaly_handling(L, ...
          settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
          'SCI and HK time ranges do not overlap in time.', ...
          'BICAS:SWMProcessing')

      elseif ~irf.utils.is_range_subset(InSci.Zv.Epoch, hkEpoch)
        %-------------------------------------------------
        % CASE: SCI does not cover a subset of HK in time
        %-------------------------------------------------
        % NOTE: This anomaly is obviously implied by the anomaly above
        % (SCI, HK do not overlap). It is therefore only meaningful to
        % detect it if the above anomaly is not detected.
        hk1RelativeSec = 1e-9 * (min(hkEpoch) - min(InSci.Zv.Epoch));
        hk2RelativeSec = 1e-9 * (max(hkEpoch) - max(InSci.Zv.Epoch));

        anomalyDescrMsg = sprintf(...
          ['HK time range is not a superset of SCI time range.', ...
          ' Can not reliably interpolate HK data for all of SCI.', ...
          ' HK begins %g s AFTER SCI begins. HK ends %g s BEFORE SCI ends.'], ...
          hk1RelativeSec, ...
          -hk2RelativeSec);

        [settingValue, settingKey] = Bso.get_fv(...
          'PROCESSING.HK.TIME_NOT_SUPERSET_OF_SCI_POLICY');
        bicas.default_anomaly_handling(L, ...
          settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
          anomalyDescrMsg, 'BICAS:DatasetFormat:SWMProcessing')

      end



      % Derive time margin within which the nearest HK value will be used.
      % NOTE: Requires >=2 records. 0 or 1 records ==> NaN (and MATLAB
      %       warning).
      %
      % BUG?: This looks as if it could yield bad results. If the margin
      %       is small, then interpolation would yield too many NaN.
      %   NOTE: HK nominally has a much lower "sampling rate" than science
      %         data. ==> The derived margin tends to be higher.
      %   NOTE: Time difference between HK samples varies, which means
      %         that mode() will not necessarily identify the de facto
      %         most common time difference.
      hkEpochExtrapMargin = mode(diff(hkEpoch)) / 2;



      %=============================================================
      % Derive BDM
      % ----------
      % NOTE: Only obtains one BDM per record
      %       ==> Can not change BDM in the middle of a record.
      % NOTE: Can potentially also obtain BDM from LFR SCI, but that
      %       decision should not be made here.
      %       See bicas.proc.L1L2.lfr.process_CDF_to_PreDc().
      %=============================================================
      bdmDoubleNan = bicas.utils.interpolate_nearest(...
        hkEpochExtrapMargin, ...
        hkEpoch, ...
        InHk.ZvFpa.HK_BIA_MODE_MUX_SET.int2doubleNan(), ...
        InSci.Zv.Epoch);
      HkSciTime.bdmFpa = bicas.utils.FPArray(bdmDoubleNan, 'FILL_VALUE', NaN).cast('uint8');



      %==================================================================
      % Derive DIFF_GAIN / isAchgFpa
      % ----------------------------
      % NOTE: Not perfect handling of time when 1 snapshot/record, since
      % one should ideally use time stamps for every LFR _sample_.
      %==================================================================
      HkSciTime.isAchgFpa = bicas.utils.FPArray.floatNan2logical(...
        bicas.utils.interpolate_nearest(...
        hkEpochExtrapMargin, ...
        hkEpoch, ...
        InHk.ZvFpa.HK_BIA_DIFF_GAIN.int2doubleNan(), ...
        InSci.Zv.Epoch));



      %=====================================
      % Derive HK_BIA_MODE_DIFF_PROBE / DLR
      %=====================================
      % NOTE: FPA uint8 --> float-NaN --> FPA logical
      HkSciTime.dlrFpa = bicas.utils.FPArray.floatNan2logical(...
        bicas.utils.interpolate_nearest(...
        hkEpochExtrapMargin, ...
        hkEpoch, ...
        InHk.ZvFpa.HK_BIA_MODE_DIFF_PROBE.int2doubleNan(), ...
        InSci.Zv.Epoch));



      %======================
      % Derive isSweepingFpa
      %======================
      isSweepingFpa = bicas.proc.L1L2.autodetect_sweeps(...
        hkEpoch, ...
        InHk.ZvFpa.HK_BIA_MODE_MUX_SET, ...
        [...
        InHk.ZvFpa.HK_BIA_BIAS1, ...
        InHk.ZvFpa.HK_BIA_BIAS2, ...
        InHk.ZvFpa.HK_BIA_BIAS3...
        ], ...
        Bso);
      HkSciTime.isSweepingFpa = bicas.utils.FPArray.floatNan2logical(...
        bicas.utils.interpolate_nearest(...
        hkEpochExtrapMargin, ...
        hkEpoch, ...
        isSweepingFpa.logical2doubleNan, ...
        InSci.Zv.Epoch));



      % ASSERTIONS
      irf.assert.struct(HkSciTime, {'bdmFpa', 'isAchgFpa', 'dlrFpa', 'isSweepingFpa'}, {})
    end



    % Utility function to shorten code.
    %
    % NOTE: Operates on entire ZvStruct since CALIBRATION_TABLE_INDEX exists
    % for L1R, but not L1, and the corresponding field may thus be or not be
    % present.
    function CALIBRATION_TABLE_INDEX = normalize_CALIBRATION_TABLE_INDEX(...
        ZvStruct, nRecords, inputDsi)

      C = bicas.classify_BICAS_L1_L1R_to_L2_DSI(inputDsi);

      if C.isL1r
        CALIBRATION_TABLE_INDEX = ZvStruct.CALIBRATION_TABLE_INDEX;
      elseif C.isL1
        CALIBRATION_TABLE_INDEX = nan(nRecords, 2);
      else
        error(...
          ['Can not normalize CALIBRATION_TABLE_INDEX', ...
          ' for this DSI classification.'])
      end

      irf.assert.sizes(CALIBRATION_TABLE_INDEX, [nRecords, 2])
    end



    % Convert PreDc+PostDc to something that
    % (1) represents a TDS dataset (hence the name), and
    % (2) ALMOST REPRESENTS an LFR dataset (the rest is done in a wrapper).
    %
    % This function only changes the data format (and selects data to send
    % to CDF).
    %
    % IMPLEMENTATION NOTE: This method is used by both LFR and TDS since the
    % L2 output datasets are very similar, despite that the input L1/L1R LFR
    % & TDS datasets are very dissimilar.
    %
    function [OutSci] = process_PostDc_to_CDF(SciPreDc, SciPostDc, outputDsi)
      % PROPOSAL: Rename to something shared between LFR and TDS, then use
      %           two wrappers.
      %   PROPOSAL: process_PostDc_to_LFR_TDS_CDF_core
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

    end    % process_PostDc_to_CDF



  end    % methods(Static, Access=public)



end
