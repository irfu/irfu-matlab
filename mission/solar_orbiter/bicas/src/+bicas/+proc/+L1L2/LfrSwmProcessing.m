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
    inputSciDsid
    inputSci    % Classification of type of processing (based on input dataset).
    outputDsid
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % inputSciDsid
    %       The science input dataset will be interpreted as having this
    %       DSID.
    %       RATIONALE: InputDatasetsMap should contain the same as a CDF
    %       global attribute but
    %       (1) it could be missing, or
    %       (2) sometimes one may want to read an ROC-SGSE dataset as if it
    %           was an RODP dataset or the other way around.
    %
    function obj = LfrSwmProcessing(inputSciDsid, outputDsid)
      obj.inputSciDsid = inputSciDsid;
      obj.inputSci    = bicas.classify_BICAS_L1_L1R_to_L2_DSID(inputSciDsid);

      obj.outputDsid   = outputDsid;
    end



    % OVERRIDE
    function OutputDatasetsMap = production_function(obj, ...
        InputDatasetsMap, rctDir, NsoTable, Bso, L)

      InputHkCdf  = InputDatasetsMap('HK_cdf');
      InputCurCdf = InputDatasetsMap('CUR_cdf');
      InputSciCdf = InputDatasetsMap('SCI_cdf');



      %======================
      % Create VCAL and CCAL
      %======================
      useGactRct = obj.inputSci.isL1r;
      useNbci    = obj.inputSci.isL1r;

      Rctdc = bicas.proc.L1L2.cal.rct.findread.get_nominal_RCTDC(...
        useGactRct, 'LFR', rctDir, ...
        InputSciCdf.Ga.CALIBRATION_TABLE, ...
        InputSciCdf.Zv.CALIBRATION_TABLE_INDEX(:, 1) + 1, ...
        InputSciCdf.Zv.BW, ...
        min(InputSciCdf.Zv.Epoch), ...
        max(InputSciCdf.Zv.Epoch), ...
        L);
      BiasRctdCa = Rctdc.get_RCTD_CA('BIAS');

      Vcds = bicas.proc.L1L2.cal.VoltageCalibrationDataSupplierImpl(Rctdc, useNbci, Bso);
      Vcal = bicas.proc.L1L2.cal.VoltageCalibration(Vcds, useGactRct, Bso);
      Ccal = bicas.proc.L1L2.cal.CurrentCalibrationImpl(BiasRctdCa{1}, Bso);



      %==============
      % Process data
      %==============
      HkSciTimePd  = bicas.proc.L1L2.process_HK_CDF_to_HK_on_SCI_TIME(InputSciCdf, InputHkCdf,  Bso, L);
      InputSciCdf  = obj.process_normalize_CDF(                       InputSciCdf, Bso, L);
      SciDcip      = obj.process_CDF_to_DCIP(                         InputSciCdf, HkSciTimePd, Bso, L);
      SciDcop      = bicas.proc.L1L2.dc.process_calibrate_demux(      SciDcip, InputCurCdf, obj.outputDsid, Vcal, Ccal, NsoTable, Bso, L);
      OutputSciCdf = obj.process_DCOP_to_CDF(                         SciDcip, SciDcop);



      OutputDatasetsMap            = containers.Map();
      RctdCa                       = Rctdc.get_global_RCTD_CA();
      OutputDatasetsMap('SCI_cdf') = bicas.OutputDataset(...
        OutputSciCdf.Zv, OutputSciCdf.Ga, RctdCa);
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
        error('BICAS:DatasetFormat', ...
          'zVariable QUALITY_BITMASK is not uint16 (MATLAB class).')
      end



      %===================================
      % Normalize CALIBRATION_TABLE_INDEX
      %===================================
      InSciNorm.Zv.CALIBRATION_TABLE_INDEX = ...
        bicas.proc.L1L2.normalize_CALIBRATION_TABLE_INDEX(...
        InSci.Zv, nRecords, obj.inputSciDsid);



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
        % NOTE: Translates from LFR's FREQ values (0=F0 etc) to LSF index
        % values (1=F0) used in loaded RCT data structs.
      end
      irf.assert.sizes(zvILsf, [nRecords])



      % NOTE: Needed also for 1 SPR.
      samplRateHz = solo.hwzv.get_LSF( zvILsf );

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
      aspr = irf.assert.sizes(...
        InSci.Zv.V, [nRecords, -1], ...
        E,          [nRecords, -1, 2]);
      if obj.inputSci.isLfrSurvSwf,   assert(aspr == solo.hwzv.const.LFR_SWF_SNAPSHOT_LENGTH)
      else,                           assert(aspr == 1)
      end



      Zv = [];

      %======================
      % Set Zv.bltsVoltageTm
      %======================
      Zv.bltsVoltageTm(:, :, 1) = single(InSci.Zv.V);
      % Copy values when there is actual data for that BLTS as determined by
      % zvLrx. Otherwise NaN.
      % zvLrx == 0: BLTS 4/5 contain data.
      % zvLrx == 1: BLTS 2/3 contain data.
      Zv.bltsVoltageTm(:, :, 2) = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvLrx==0 );
      Zv.bltsVoltageTm(:, :, 3) = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvLrx==0 );
      Zv.bltsVoltageTm(:, :, 4) = bicas.proc.utils.set_NaN_rows( E(:,:,1), zvLrx==1 );
      Zv.bltsVoltageTm(:, :, 5) = bicas.proc.utils.set_NaN_rows( E(:,:,2), zvLrx==1 );

      Zv.tt2000                  = InSci.Zv.Epoch;
      % NOTE: DELTA_PLUS_MINUS only applies to Epoch, and must therefore have
      % consistent number of dimensions, regardless of CWF/SWF.
      Zv.DELTA_PLUS_MINUS        = bicas.proc.utils.derive_DELTA_PLUS_MINUS(...
        samplRateHz, 1);
      Zv.samplRateHz             = samplRateHz;
      Zv.uspr                    = ones(nRecords, 1) * aspr;
      % BW needed because it is copied to the output CDF.
      Zv.BW                      = InSci.Zv.BW;

      Zv.isAchgFpa               = HkSciTime.isAchgFpa;
      Zv.dlrFpa                  = HkSciTime.dlrFpa;
      Zv.iLsf                    = zvILsf;
      Zv.lrx                     = zvLrx;
      Zv.SYNCHRO_FLAG            = InSci.Zv.SYNCHRO_FLAG;
      Zv.L1qbmFpa                = InSci.ZvFpa.QUALITY_BITMASK;
      Zv.QflFpa                  = InSci.ZvFpa.QUALITY_FLAG;

      % QRCB arrayss for
      % (1) when LFR ZV BW says BIAS is OFF, and
      % (2) BIAS is sweeping.
      % so that "quality actions" can be taken later based on these.
      Zv.biasOffQrcb             = ~logical(InSci.Zv.BW);
      Zv.sweepQrcb               = HkSciTime.isSweepingFpa.array(false);



      % Replace CALIBRATION_TABLE_INDEX-->NbriFpa + NbciFpa
      bBlank                     = ~logical(InSci.Zv.BW);
      Zv.NbriFpa                 = bicas.utils.FPArray(...
        InSci.Zv.CALIBRATION_TABLE_INDEX(:, 1), 'FILL_POSITIONS', bBlank) + uint8(1);
      Zv.NbciFpa                 = bicas.utils.FPArray(...
        InSci.Zv.CALIBRATION_TABLE_INDEX(:, 2), 'FILL_POSITIONS', bBlank);



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

      Dcip = bicas.proc.L1L2.DemultiplexingCalibrationInput(...
        Zv, Ga, obj.inputSci.isLfrSurvSwf, true, false);
    end    % process_CDF_to_DCIP



    function [OutSci] = process_DCOP_to_CDF(obj, SciDcip, SciDcop)
      % NOTE: Most processing is done in a function which is shared between LFR
      %       and TDS. This wrapper is needed to handle the small difference
      %       between LFR and TDS.
      OutSci = bicas.proc.L1L2.process_DCOP_to_CDF(...
        SciDcip, SciDcop, obj.outputDsid);

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
        L, policySettingValue, policySettingKey, nRecords, ZvFpa1, zvName)

      if ~isempty(ZvFpa1)
        % Do nothing (except assertion later).
        ZvFpa2 = ZvFpa1;
      else
        anomalyDescrMsg = sprintf(...
          'zVar "%s" from the LFR SCI source dataset is empty.', ...
          zvName);

        switch(policySettingValue)
          case 'USE_FILL_VALUE'
            bicas.default_anomaly_handling(L, ...
              policySettingValue, policySettingKey, 'OTHER', ...
              anomalyDescrMsg, ...
              'BICAS:DatasetFormat:SWMProcessing')

            L.logf('warning', 'Using fill values for %s.', zvName)
            ZvFpa2 = bicas.utils.FPArray(...
              zeros(nRecords, 1, ZvFpa1.mc), 'ONLY_FILL_POSITIONS');

          otherwise
            bicas.default_anomaly_handling(L, ...
              policySettingValue, policySettingKey, 'ERROR_ILLEGAL_SETTING', ...
              anomalyDescrMsg, ...
              'BICAS:DatasetFormat:SWMProcessing')
        end
      end

      irf.assert.sizes(ZvFpa2, [NaN])
    end



  end    % methods(Static, Access=private)



end
