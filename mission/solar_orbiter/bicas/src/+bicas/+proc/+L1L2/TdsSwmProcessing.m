%
% SWMP for processing TDS L1/L1R --> L2.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef TdsSwmProcessing < bicas.proc.SwmProcessing
  % PROPOSAL: Automatic test code.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=private)
    inputSciDsi
    inputSci    % Classification of type of processing (based on input dataset).
    outputDsi
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % inputSciDsi
    %       The science input dataset will be interpreted as having this
    %       DSI.
    %       RATIONALE: InputDatasetsMap should contain the same as a CDF
    %       global attribute but
    %       (1) it could be missing, or
    %       (2) sometimes one may want to read an ROC-SGSE dataset as if it
    %           was an RODP dataset or the other way around.
    %
    function obj = TdsSwmProcessing(inputSciDsi, outputDsi)
      obj.inputSciDsi = inputSciDsi;
      obj.inputSci    = bicas.classify_BICAS_L1_L1R_to_L2_DSI(inputSciDsi);

      obj.outputDsi   = outputDsi;
    end



    % OVERRIDE
    function OutputDatasetsMap = production_function(obj, ...
        InputDatasetsMap, rctDir, NsoTable, Bso, L)

      InputHkCdf  = InputDatasetsMap('HK_cdf');
      InputCurCdf = InputDatasetsMap('CUR_cdf');
      InputSciCdf = InputDatasetsMap('SCI_cdf');



      %==========================================
      % Configure bicas.proc.L1L2.cal.Cal object
      %==========================================
      % NOTE: TDS L1R never uses ZVCTI2.
      if obj.inputSci.isTdsCwf
        settingUseGactRct = 'PROCESSING.L1R.TDS.CWF.USE_GA_CALIBRATION_TABLE_RCTS';
        tdsRcttid         = 'TDS-CWF';
      else
        settingUseGactRct = 'PROCESSING.L1R.TDS.RSWF.USE_GA_CALIBRATION_TABLE_RCTS';
        tdsRcttid         = 'TDS-RSWF';
      end
      useGactRct = obj.inputSci.isL1r && Bso.get_fv(settingUseGactRct);
      useZvcti2  = false;    % Always false for TDS.

      % Create a synthetic zv_BW since it does not exist for TDS (only LFR).
      % --
      % NOTE: This should not be regarded as a hack but as ~normalization to
      % avoid later special cases.
      zv_BW = uint8(ones(...
        size(InputSciCdf.Zv.CALIBRATION_TABLE_INDEX, 1), ...
        1));

      Rctdc = bicas.proc.L1L2.cal.rct.findread.get_nominal_RCTDC(...
        useGactRct, tdsRcttid, rctDir, ...
        InputSciCdf.Ga.CALIBRATION_TABLE, ...
        InputSciCdf.Zv.CALIBRATION_TABLE_INDEX, ...
        zv_BW, ...
        min(InputSciCdf.Zv.Epoch), ...
        max(InputSciCdf.Zv.Epoch), ...
        L);

      Cal = bicas.proc.L1L2.cal.Cal(Rctdc, useGactRct, useZvcti2, Bso);



      %==============
      % Process data
      %==============
      HkSciTimePd  = bicas.proc.L1L2.process_HK_CDF_to_HK_on_SCI_TIME(InputSciCdf, InputHkCdf,  Bso, L);
      InputSciCdf  = obj.process_normalize_CDF(                       InputSciCdf,              Bso, L);
      SciDcip      = obj.process_CDF_to_DCIP(                         InputSciCdf, HkSciTimePd);
      SciDcop      = bicas.proc.L1L2.dc.process_calibrate_demux(      SciDcip, InputCurCdf, Cal, NsoTable, Bso, L);
      OutputSciCdf = bicas.proc.L1L2.process_DCOP_to_CDF(             SciDcip, SciDcop, obj.outputDsi);



      OutputDatasetsMap = containers.Map();
      RctdCa = Rctdc.get_global_RCTD_CA();
      OutputDatasetsMap('SCI_cdf') = bicas.OutputDataset(OutputSciCdf.Zv, OutputSciCdf.Ga, RctdCa);
    end



  end    % methods(Access=public)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Only "normalizes" data to account for technically
    % illegal input TDS datasets. It should try to:
    % ** modify L1 to look like L1R
    % ** mitigate historical bugs in the input datasets
    % ** mitigate for not yet implemented features in input datasets
    %
    function InSciNorm = process_normalize_CDF(obj, InSci, Bso, L)

      % Default behaviour: Copy values, except for values which are
      % modified later
      InSciNorm = InSci;

      nRecords = irf.assert.sizes(InSci.Zv.Epoch, [-1]);



      %===================================
      % Normalize CALIBRATION_TABLE_INDEX
      %===================================
      InSciNorm.Zv.CALIBRATION_TABLE_INDEX = bicas.proc.L1L2.normalize_ZVCTI(...
        InSci.Zv, nRecords, obj.inputSciDsi);



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
        fnChangeList, obj.inputSciDsi, Bso, L, ...
        'SYNCHRO_FLAG', 'INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY')



      %=========================
      % Normalize SAMPLING_RATE
      %=========================
      if any(InSci.Zv.SAMPLING_RATE == 255)
        [settingValue, settingKey] = Bso.get_fv(...
          'PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY');
        anomalyDescrMsg = ...
          ['Finds illegal, stated sampling frequency', ...
          ' 255 in TDS L1/L1R LFM-RSWF dataset.'];

        if obj.inputSci.isTdsRswf
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
                settingValue, settingKey, 'OTHER', anomalyDescrMsg)

            otherwise
              bicas.default_anomaly_handling(L, ...
                settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
                anomalyDescrMsg, ...
                'BICAS:DatasetFormat')
          end
        else
          error('BICAS:DatasetFormat', anomalyDescrMsg)
        end
      end



      if obj.inputSci.isTdsRswf
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

          [settingValue, settingKey] = Bso.get_fv(...
            'PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY');
          switch(settingValue)

            case 'ROUND'
              bicas.default_anomaly_handling(...
                L, settingValue, settingKey, 'OTHER', ...
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
                settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
                anomalyDescrMsg, ...
                'BICAS:Assertion:DatasetFormat')

          end    % switch
        end    % if
      end    % if

    end    % process_normalize_CDF



    % Convert TDS CDF data (PDs) to DCIP.
    function Dcip = process_CDF_to_DCIP(obj, InSci, HkSciTime)
      %
      % BUG?: Does not use CHANNEL_STATUS_INFO.
      % NOTE: BIAS output datasets do not have a variable for the length of
      % snapshots. Need to use NaN/fill value.

      % ASSERTIONS: VARIABLES
      assert(isa(InSci, 'bicas.InputDataset'))
      irf.assert.struct(HkSciTime, {'bdmFpa', 'isAchgFpa', 'dlrFpa', 'isSweepingFpa'}, {})



      % ASSERTIONS: CDF
      bicas.proc.utils.assert_increasing(...
        InSci.Zv.Epoch, true, 'BICAS:DatasetFormat', ...
        ['Voltage (science) dataset timestamps Epoch do not', ...
        ' increase monotonously.']...
        )
      [nRecords, WAVEFORM_DATA_nChannels, nCdfSamplesPerRecord] = irf.assert.sizes(...
        InSci.Zv.Epoch,         [-1], ...
        InSci.Zv.WAVEFORM_DATA, [-1, -2, -3]);
      if     obj.inputSci.isL1r,   WAVEFORM_DATA_nChannels_expected = 3;
      elseif obj.inputSci.isL1,    WAVEFORM_DATA_nChannels_expected = 8;
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
      if obj.inputSci.isTdsRswf
        assert(...
          nCdfSamplesPerRecord == solo.hwzv.const.TDS_RSWF_L1R_SAMPLES_PER_RECORD, ...
          'Unexpected number of samples per CDF record (%i). Expected %i.', ...
          nCdfSamplesPerRecord, solo.hwzv.const.TDS_RSWF_L1R_SAMPLES_PER_RECORD)
      else
        assert(nCdfSamplesPerRecord == 1)
      end



      % TODO-NI: Why convert to double? To avoid precision problems when
      % doing math with other variables?
      zvFreqHz = double(InSci.Zv.SAMPLING_RATE);



      Zv    = [];

      Zv.Epoch                   = InSci.Zv.Epoch;
      % NOTE: DELTA_PLUS_MINUS is only applies to Epoch, and must therefore have
      % consistent number of dimensions, regardless of CWF/SWF.
      Zv.DELTA_PLUS_MINUS        = bicas.proc.utils.derive_DELTA_PLUS_MINUS(...
        zvFreqHz, 1);
      Zv.freqHz                  = zvFreqHz;
      Zv.QUALITY_BITMASK         = InSci.ZvFpa.QUALITY_BITMASK;
      Zv.QUALITY_FLAG            = InSci.ZvFpa.QUALITY_FLAG;
      Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
      Zv.bdmFpa                  = HkSciTime.bdmFpa;
      Zv.isAchgFpa               = HkSciTime.isAchgFpa;
      Zv.dlrFpa                  = HkSciTime.dlrFpa;
      Zv.ufv                     = HkSciTime.isSweepingFpa.array(false);
      Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;



      %=====================================
      % Set Zv.nValidSamplesPerRecord
      %=====================================
      if obj.inputSci.isTdsRswf
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

      Zv.bltsSamplesTm(:, :, 1) = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,1), Zv.nValidSamplesPerRecord );
      Zv.bltsSamplesTm(:, :, 2) = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,2), Zv.nValidSamplesPerRecord );
      Zv.bltsSamplesTm(:, :, 3) = bicas.proc.utils.set_NaN_end_of_rows( zv_WAVEFORM_DATA_modif(:,:,3), Zv.nValidSamplesPerRecord );
      Zv.bltsSamplesTm(:, :, 4) = nan(nRecords, nCdfSamplesPerRecord);
      Zv.bltsSamplesTm(:, :, 5) = nan(nRecords, nCdfSamplesPerRecord);



      Ga = [];
      Ga.OBS_ID    = InSci.Ga.OBS_ID;
      Ga.SOOP_TYPE = InSci.Ga.SOOP_TYPE;

      % Only set because the code shared with LFR requires it.
      Zv.iLsf      = nan( nRecords, 1);
      Zv.lrx       = ones(nRecords, 1);
      Zv.BW        = true(nRecords, 1);

      Dcip = bicas.proc.L1L2.DemultiplexingCalibrationInput(Zv, Ga, obj.inputSci.isTdsRswf, false, obj.inputSci.isTdsCwf);
    end    % process_CDF_to_DCIP



  end    % methods(Access=private)



end
