%
% SWMP for processing LFR L1/L1R --> L2.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef LfrSwmProcessing < bicas.proc.SwmProcessing
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
    function obj = LfrSwmProcessing(inputSciDsi, outputDsi)
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
      useGactRct = obj.inputSci.isL1r && Bso.get_fv('PROCESSING.L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS');
      useZvcti2  = obj.inputSci.isL1r && Bso.get_fv('PROCESSING.L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2');

      Rctdc = bicas.proc.L1L2.cal.rct.findread.get_nominal_RCTDC(...
        useGactRct, 'LFR', rctDir, ...
        InputSciCdf.Ga.CALIBRATION_TABLE, ...
        InputSciCdf.Zv.CALIBRATION_TABLE_INDEX, ...
        InputSciCdf.Zv.BW, ...
        min(InputSciCdf.Zv.Epoch), ...
        max(InputSciCdf.Zv.Epoch), ...
        L);

      Cal = bicas.proc.L1L2.cal.Cal(Rctdc, useGactRct, useZvcti2, Bso);



      %==============
      % Process data
      %==============
      HkSciTimePd  = bicas.proc.L1L2.process_HK_CDF_to_HK_on_SCI_TIME(InputSciCdf, InputHkCdf,  Bso, L);
      InputSciCdf  = obj.process_normalize_CDF(                       InputSciCdf, Bso, L);
      SciDcip      = obj.process_CDF_to_DCIP(                         InputSciCdf, HkSciTimePd, Bso, L);
      SciDcop      = bicas.proc.L1L2.dc.process_calibrate_demux(      SciDcip, InputCurCdf, Cal, NsoTable, Bso, L);
      OutputSciCdf = obj.process_DCOP_to_CDF(                         SciDcip, SciDcop);



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



    % Only "normalizes" data to account for technically illegal input LFR
    % datasets. It should try to:
    % ** modify L1 data to look like L1R
    % ** mitigate historical bugs in input datasets
    % ** mitigate for not yet implemented features in input datasets
    %
    function InSciNorm = process_normalize_CDF(obj, InSci, Bso, L)
      % ASSERTIONS: VARIABLES
      assert(isa(InSci, 'bicas.InputDataset'))

      % Default behaviour: Copy values, except for values which are
      % modified later
      InSciNorm = InSci;

      nRecords = irf.assert.sizes(InSci.Zv.Epoch, [-1]);



      % NOTE: Very old L1R LFR-SBM1/2 (test) data datasets have been
      % observed to have QUALITY_BITMASK with illegal data type
      % CDF_UINT1/uint8 while newer ones do not. Could normalize for this
      % but it should be better to simply not support (and thus not use)
      % such datasets.
      % Ex: QUALITY_BITMASK uses CDF_UINT1/uint8 in
      %
      if ~strcmp(InSciNorm.ZvFpa.QUALITY_BITMASK.mc, 'uint16')
        error('BICAS:DatasetFormat', 'zVariable QUALITY_BITMASK is not uint16 (MATLAB class).')
      end



      %===================================
      % Normalize CALIBRATION_TABLE_INDEX
      %===================================
      InSciNorm.Zv.CALIBRATION_TABLE_INDEX = ...
        bicas.proc.L1L2.normalize_ZVCTI(...
        InSci.Zv, nRecords, obj.inputSciDsi);



      %========================
      % Normalize SYNCHRO_FLAG
      %========================
      has_SYNCHRO_FLAG      = isfield(InSci.Zv, 'SYNCHRO_FLAG');
      has_TIME_SYNCHRO_FLAG = isfield(InSci.Zv, 'TIME_SYNCHRO_FLAG');
      if      has_SYNCHRO_FLAG && ~has_TIME_SYNCHRO_FLAG

        % CASE: Everything nominal.
        InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.SYNCHRO_FLAG;

      elseif ~has_SYNCHRO_FLAG && has_TIME_SYNCHRO_FLAG

        % CASE: Input CDF uses wrong zVar name.
        [settingValue, settingKey] = ...
          Bso.get_fv('INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY');
        bicas.default_anomaly_handling(L, ...
          settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
          'Found zVar TIME_SYNCHRO_FLAG instead of SYNCHRO_FLAG.')
        L.log('warning', ...
          'Using illegally named zVar TIME_SYNCHRO_FLAG as SYNCHRO_FLAG.')
        InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.TIME_SYNCHRO_FLAG;

      elseif has_SYNCHRO_FLAG && has_TIME_SYNCHRO_FLAG

        % CASE: Input CDF has two ZVs: one with correct name, one with
        % incorrect name

        %------------------------
        % "Normal" normalization
        %------------------------
        % 2020-01-21: Based on skeletons (.skt; L1R, L2), SYNCHRO_FLAG
        % seems to be the correct zVar.
        if Bso.get_fv(...
            'INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED') ...
            && isempty(InSci.Zv.SYNCHRO_FLAG)
          %----------------------------------------------------------
          % Workaround: Normalize LFR data to handle variations that
          % should not exist
          %----------------------------------------------------------
          % Handle that SYNCHRO_FLAG (empty) and TIME_SYNCHRO_FLAG
          % (non-empty) may BOTH be present. "DEFINITION BUG" in
          % definition of datasets/skeleton?
          % Ex: LFR___TESTDATA_RGTS_LFR_CALBUT_V0.7.0/ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_4129f0b_CNE_V02.cdf /2020-03-17

          InSciNorm.Zv.SYNCHRO_FLAG = InSci.Zv.TIME_SYNCHRO_FLAG;
        else
          error('BICAS:DatasetFormat', ...
            ['Input dataset has both zVar SYNCHRO_FLAG and', ...
            ' TIME_SYNCHRO_FLAG.'])
        end
      else
        error('BICAS:DatasetFormat', ...
          'Input dataset does not have zVar SYNCHRO_FLAG as expected.')
      end



      %=======================================================================================================
      % Set QUALITY_BITMASK, QUALITY_FLAG:
      % Replace illegally empty data with fill values/NaN
      % ------------------------------------------------------------------
      % IMPLEMENTATION NOTE: QUALITY_BITMASK, QUALITY_FLAG have been found
      % empty in test data, but should have attribute DEPEND_0 = "Epoch"
      % ==> Should have same number of records as Epoch.
      %
      % Can not save CDF with zVar with zero records (crashes when reading
      % CDF). ==> Better create empty records.
      %
      % Examples of QUALITY_FLAG = empty:
      %  MYSTERIOUS_SIGNAL_1_2016-04-15_Run2__7729147__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-SWF_7729147_CNE_V01.cdf
      %  ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_4129f0b_CNE_V02.cdf (TESTDATA_RGTS_LFR_CALBUT_V1.1.0)
      %  ROC-SGSE_L1R_RPW-LFR-SBM2-CWF-E_6b05822_CNE_V02.cdf (TESTDATA_RGTS_LFR_CALBUT_V1.1.0)
      %=======================================================================================================
      % PROPOSAL: Move to the code that reads CDF datasets instead. Generalize to many zVariables.
      % PROPOSAL: Regard as "normalization" code. ==> Group together with other normalization code.
      %=======================================================================================================
      [settingValue, settingKey] = Bso.get_fv(...
        'PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY');

      InSciNorm.ZvFpa.QUALITY_BITMASK = bicas.proc.L1L2.LfrSwmProcessing.normalize_ZV_empty(...
        L, settingValue, settingKey, nRecords, ...
        InSci.ZvFpa.QUALITY_BITMASK, 'QUALITY_BITMASK');

      InSciNorm.ZvFpa.QUALITY_FLAG    = bicas.proc.L1L2.LfrSwmProcessing.normalize_ZV_empty(...
        L, settingValue, settingKey, nRecords, ...
        InSci.ZvFpa.QUALITY_FLAG,    'QUALITY_FLAG');

      % ASSERTIONS
      irf.assert.sizes(...
        InSciNorm.ZvFpa.QUALITY_BITMASK, [nRecords, 1], ...
        InSciNorm.ZvFpa.QUALITY_FLAG,    [nRecords, 1])

    end    % process_normalize_CDF



    % Convert LFR CDF data to DCIP.
    %
    % IMPLEMENTATION NOTE: Does not modify InSci in an attempt to save RAM
    % (should help MATLAB's optimization). Unclear if actually works.
    %
    function Dcip = process_CDF_to_DCIP(obj, InSci, HkSciTime, Bso, L)
      %
      % PROBLEM: Hard-coded CDF data types (MATLAB classes).
      % MINOR PROBLEM: Still does not handle LFR zVar TYPE for determining
      % "virtual snapshot" length. Should only be relevant for
      % V01_ROC-SGSE_L2R_RPW-LFR-SURV-CWF (not V02) which should expire.

      % ASSERTIONS: VARIABLES
      assert(isa(InSci, 'bicas.InputDataset'))
      irf.assert.struct(HkSciTime, {'bdmFpa', 'isAchgFpa', 'dlrFpa', 'isSweepingFpa'}, {})

      % ASSERTIONS: CDF
      bicas.proc.utils.assert_increasing(...
        InSci.Zv.Epoch, true, 'BICAS:DatasetFormat', ...
        ['Voltage (science) dataset timestamps Epoch do not', ...
        ' increase monotonously.']...
        )
      nRecords = irf.assert.sizes(InSci.Zv.Epoch, [-1]);



      %============
      % Set iLsfZv
      %============
      if     obj.inputSci.isLfrSbm1
        zvILsf = ones(nRecords, 1) * 2;   % Always value "2" (F1, "FREQ = 1").
      elseif obj.inputSci.isLfrSbm2
        zvILsf = ones(nRecords, 1) * 3;   % Always value "3" (F2, "FREQ = 2").
      else
        zvILsf = InSci.Zv.FREQ + 1;
        % NOTE: Translates from LFR's FREQ values (0=F0 etc) to LSF
        % index values (1=F0) used in loaded RCT data structs.
      end
      irf.assert.sizes(zvILsf, [nRecords])



      % NOTE: Needed also for 1 SPR.
      zvFreqHz = solo.hwzv.get_LSF( zvILsf );

      % Obtain the relevant values (one per record) from zVariables R0,
      % R1, R2, and the virtual "R3".
      zvLrx = solo.hwzv.get_LRX(...
        InSci.Zv.R0, ...
        InSci.Zv.R1, ...
        InSci.Zv.R2, ...
        zvILsf);



      %===================================================================
      % IMPLEMENTATION NOTE: E & V must be floating-point so that values
      % can be set to NaN.
      %
      % Switch last two indices of E.
      % ==> index 2 = "snapshot" sample index, including for CWF
      %               (sample/record, "snapshots" consisting of 1 sample).
      %     index 3 = E1/E2 component
      %               NOTE: 1/2=index into array; these are diffs but not
      %               equivalent to any particular diffs).
      %===================================================================
      E = single(permute(InSci.Zv.E, [1,3,2]));

      % ASSERTIONS
      nCdfSamplesPerRecord = irf.assert.sizes(...
        InSci.Zv.V, [nRecords, -1], ...
        E,          [nRecords, -1, 2]);
      if obj.inputSci.isLfrSurvSwf,   assert(nCdfSamplesPerRecord == solo.hwzv.const.LFR_SWF_SNAPSHOT_LENGTH)
      else,                           assert(nCdfSamplesPerRecord == 1)
      end



      Zv = [];

      Zv.bltsSamplesTm(:, :, 1) = single(InSci.Zv.V);
      % Copy values when there is actual data for that BLTS as determined
      % by zvLrx. Otherwise NaN.
      % zvLrx == 0: BLTS 4/5 contain data.
      % zvLrx == 1: BLTS 2/3 contain data.
      Zv.bltsSamplesTm(:, :, 2) = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvLrx==0 );
      Zv.bltsSamplesTm(:, :, 3) = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvLrx==0 );
      Zv.bltsSamplesTm(:, :, 4) = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvLrx==1 );
      Zv.bltsSamplesTm(:, :, 5) = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvLrx==1 );

      Zv.Epoch                   = InSci.Zv.Epoch;
      % NOTE: DELTA_PLUS_MINUS is only applies to Epoch, and must therefore have
      % consistent number of dimensions, regardless of CWF/SWF.
      Zv.DELTA_PLUS_MINUS        = bicas.proc.utils.derive_DELTA_PLUS_MINUS(...
        zvFreqHz, 1);
      Zv.freqHz                  = zvFreqHz;
      Zv.nValidSamplesPerRecord  = ones(nRecords, 1) * nCdfSamplesPerRecord;
      Zv.BW                      = InSci.Zv.BW;
      Zv.ufv                     = ~logical(InSci.Zv.BW) | HkSciTime.isSweepingFpa.array(false);
      Zv.isAchgFpa               = HkSciTime.isAchgFpa;
      Zv.dlrFpa                  = HkSciTime.dlrFpa;
      Zv.iLsf                    = zvILsf;

      Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
      Zv.CALIBRATION_TABLE_INDEX = InSci.Zv.CALIBRATION_TABLE_INDEX;

      Zv.QUALITY_BITMASK         = InSci.ZvFpa.QUALITY_BITMASK;
      Zv.QUALITY_FLAG            = InSci.ZvFpa.QUALITY_FLAG;

      Zv.lrx                     = zvLrx;



      %=====================================
      % Set BDM
      % -------
      % Select which source of BDM is used.
      %=====================================
      [bdmSrcSettingValue, bdmSrcSettingKey] = Bso.get_fv('PROCESSING.LFR.MUX_MODE_SOURCE');
      switch(bdmSrcSettingValue)
        case 'BIAS_HK'
          L.log('debug', 'Using BIAS HK mux mode (BDM).')
          bdmFpa = HkSciTime.bdmFpa;

        case 'LFR_SCI'
          L.log('debug', 'Using LFR SCI mux mode (BDM).')
          bdmFpa = InSci.ZvFpa.BIAS_MODE_MUX_SET;

        case 'BIAS_HK_LFR_SCI'
          L.log('debug', ...
            ['Using mux mode (BDM) from BIAS HK when available, and', ...
            ' from LFR SCI when the former is not available.'])

          bdmFpa = HkSciTime.bdmFpa.complement(InSci.ZvFpa.BIAS_MODE_MUX_SET);

        otherwise
          error('BICAS:ConfigurationBug', ...
            'Illegal settings value %s="%s"', bdmSrcSettingKey, bdmSrcSettingValue)
      end
      Zv.bdmFpa = bdmFpa;



      Ga = [];
      Ga.OBS_ID    = InSci.Ga.OBS_ID;
      Ga.SOOP_TYPE = InSci.Ga.SOOP_TYPE;

      Dcip = bicas.proc.L1L2.DemultiplexingCalibrationInput(Zv, Ga, obj.inputSci.isLfrSurvSwf, true, false);

    end    % process_CDF_to_DCIP



    function [OutSci] = process_DCOP_to_CDF(obj, SciDcip, SciDcop)
      % NOTE: Most processing is done in function shared between LFR and
      %       TDS.
      OutSci = bicas.proc.L1L2.process_DCOP_to_CDF(...
        SciDcip, SciDcop, obj.outputDsi);

      OutSci.Zv.BW = SciDcip.Zv.BW;
    end



  end    % methods(Access=private)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Local utility function to shorten & clarify code.
    %
    % ARGUMENTS
    % =========
    % ZvFpa1
    %       ZV-like FPA. Column vector (Nx1) or empty.
    %
    function ZvFpa2 = normalize_ZV_empty(...
        L, settingValue, settingKey, nRecords, ZvFpa1, zvName)

      if ~isempty(ZvFpa1)
        % Do nothing (except assertion later).
        ZvFpa2 = ZvFpa1;
      else
        anomalyDescrMsg = sprintf(...
          'zVar "%s" from the LFR SCI source dataset is empty.', ...
          zvName);

        switch(settingValue)
          case 'USE_FILL_VALUE'
            bicas.default_anomaly_handling(L, ...
              settingValue, settingKey, 'OTHER', ...
              anomalyDescrMsg, ...
              'BICAS:DatasetFormat:SWMProcessing')

            L.logf('warning', 'Using fill values for %s.', zvName)
            ZvFpa2 = bicas.utils.FPArray(...
              zeros(nRecords, 1, ZvFpa1.mc), 'ONLY_FILL_POSITIONS');

          otherwise
            bicas.default_anomaly_handling(L, ...
              settingValue, settingKey, 'ERROR_ILLEGAL_SETTING', ...
              anomalyDescrMsg, ...
              'BICAS:DatasetFormat:SWMProcessing')
        end
      end

      irf.assert.sizes(ZvFpa2, [NaN])
    end



  end    % methods(Static, Access=private)



end
