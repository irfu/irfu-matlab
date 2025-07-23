%
% SWMP for processing L2 LFR CWF to L3 OSR+DSR density, E field, and ScPot.
%
%
% CODE CONVENTIONS
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like"
%   data, use the first MATLAB array index to represent CDF records.
%
%
% PERFORMANCE
% ===========
% Simple tests shows that the duplicated downsampling of duplicated quality
% variables (bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template()) slows down L2->L3
% processing (takes ~20% more time). This was implemented for a bugfix. One can
% eliminate this (except for QUALITY_FLAG for density which is different), but
% it makes the implementation sensitive to any EFIELD/DENSITY/SCPOT-specific
% future modifications of their particular OSR quality variables (i.e. one has
% to change the implementation "a lot").
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef L3OsrDsrSwmProcessing < bicas.proc.SwmProcessing
  % PROPOSAL: Automatic test code.
  %   NOTE: There are limited tests.
  %
  % PROPOSAL: Split up processing between (a) density, and (b) E-field & SCPOT
  %           into separate SWMs.
  %   PRO: Faster processing when only processing subset of L3 DSIs.
  %       CON: Not very heavy operation.
  %   PRO: Leads to better organization of code.
  %       PRO: process_L2_to_L3() is too large and should be split up anyway.
  %   CON: DENSITY is a function EFIELD+SCPOT, and thus has to be processed
  %        after the latter.



  %###################
  %###################
  % PRIVATE CONSTANTS
  %###################
  %###################
  properties(Constant, Access=private)



    % Define length of bins, and relative position of corresponding
    % bin timestamps.
    BIN_LENGTH_WOLS_NS        = int64(10e9);
    BIN_TIMESTAMP_POS_WOLS_NS = int64(bicas.proc.L2L3.L3OsrDsrSwmProcessing.BIN_LENGTH_WOLS_NS / 2);



  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % OVERRIDE
    %
    % NOTE: Does not use NsoTable.
    %
    function OutputDatasetsMap = production_function(obj, ...
        InputDatasetsMap, rctDir, NsoTable, Bso, L)

      InputLfrCwfCdf = InputDatasetsMap('LFR-SURV-CWF-E_cdf');

      Excd = bicas.proc.L2L3.ExternalCodeImplementation();

      %==============
      % Process data
      %==============
      [...
        EfieldOsrCdf,  EfieldDsrCdf, ...
        ScpotOsrCdf,   ScpotDsrCdf, ...
        DensityOsrCdf, DensityDsrCdf ...
      ] = bicas.proc.L2L3.L3OsrDsrSwmProcessing.process_L2_to_L3(...
        InputLfrCwfCdf, Excd, Bso, L);

      OutputDatasetsMap = containers.Map();
      OutputDatasetsMap('EFIELD_OSR_cdf')  = EfieldOsrCdf;
      OutputDatasetsMap('EFIELD_DSR_cdf')  = EfieldDsrCdf;
      OutputDatasetsMap('SCPOT_OSR_cdf')   = ScpotOsrCdf;
      OutputDatasetsMap('SCPOT_DSR_cdf')   = ScpotDsrCdf;
      OutputDatasetsMap('DENSITY_OSR_cdf') = DensityOsrCdf;
      OutputDatasetsMap('DENSITY_DSR_cdf') = DensityDsrCdf;
    end



  end    % methods(Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % ~Process L2-->L3 (not VHT).
    %
    % NOTE: Function assumes that (some) fill values (the dedicated constant
    % values) for integer-valued zVariables are identical in input and output
    % datasets.
    %
    % NOTE: Function does not discard data with QUALITY_FLAG==fill value, as
    % opposed to QUALITY_FLAG < threshold.
    %
    % NOTE: Keeps quality zVariables also when all data in record is FV.
    % Both OSR and DSR. /EJ 2025-06-12
    %
    % IMPLEMENTATION NOTE: This function is separate from
    % bicas.proc.L2L3.L3OsrDsrSwmProcessing.production_function() to facilitate
    % automated tests, in particular by adding an explicit dependence on
    % bicas.proc.L2L3.ExternalCodeImplementation.
    %
    function [...
        OutEfieldOsr,  OutEfieldDsr, ...
        OutScpotOsr,   OutScpotDsr, ...
        OutDensityOsr, OutDensityDsr ...
        ] ...
        = process_L2_to_L3(InLfrCwf, Excd, Bso, L)

      if 0
        % TEST
        % PRETEND that input QUALITY_FLAG, L2_QUALITY_BITMASK have other values
        % (overwrite).
        InLfrCwf.ZvFpa.QUALITY_FLAG       = bicas.utils.FPArray(...
          uint8(ones(size(InLfrCwf.ZvFpa.QUALITY_FLAG))) * 2, ...
          'NO_FILL_POSITIONS');  % TEST
        InLfrCwf.ZvFpa.L2_QUALITY_BITMASK = bicas.utils.FPArray(...
          uint16(ones(size(InLfrCwf.ZvFpa.L2_QUALITY_BITMASK))) * 0, ...
          'NO_FILL_POSITIONS');  % TEST
      end



      % PROPOSAL: Split up into different parts for EFIELD, SCPOT, DENSITY
      %           (still combine non-downsampled and downsampled).
      %   CON: Slows down overall processing.
      %       PRO: Must read same L2 dataset multiple times.
      %       PRO: Must read L3 SCPOT dataset to produce L3 DENSITY dataset.
      %   CON: There is much shared functionality for 3 quality ZVs.
      %       PRO: Same ~constants
      %           Ex: INPUT_DSI, BIN_LENGTH_WOLS_NS, BIN_TIMESTAMP_POS_WOLS_NS
      %       PRO: Read setting QUALITY_FLAG_MIN
      %       PRO: Normalizing CWF zVar names.
      %       PRO: Preparations for downsampled.
      %           Bin locations, bundling of records,
      %           Downsampling of quality variables
      %               (QUALITY_FLAG, QUALITY_BITMASK, L2_QUALITY_BITMASK).
      %           DELTA_PLUS_MINUS_dsr
      %
      % NOTE: ROC BUG: https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/48
      %         L1 QUALITY_BITMASK seems to use the wrong value (255) as
      %         fill value (FILLVAL=65535). ==> A BICAS bug fix would not
      %         fix the entire issue.
      %   PROPOSAL: Use double also for CDF integer variables so NaN can
      %             represent fill value also for these.
      %   PROPOSAL: Implement MATLAB equivalent of the JUICE pipeline's
      %             FPA class.
      %
      % NOTE: L2 LFR-CWF-E skt previously had zVar
      %   QUALITY_BITMASK=CDF_UINT1, fill value=255 (wrong)
      % until skt V12 when it was changed to
      %   QUALITY_BITMASK=CDF_UINT2, fill value 65535 (correct).

      Tmk = bicas.utils.Timekeeper('bicas.proc.L2L3.L3OsrDsrSwmProcessing.process_L2_to_L3', L);
      assert(isa(Excd, 'bicas.proc.L2L3.ExternalCodeAbstract'))



      % EXPERIMENTAL
      % %#########################################
      % % NSO table, L2_QUALITY_BITMASK --> QRCBs
      % %#########################################
      % % Read NSO table into QRCBs ONCE, so that it does not need to be done
      % % later.
      % AllQrcbMap = bicas.proc.qual.NSO_table_to_QRCB_map(...
      %   bicas.const.Q.ALL_QRCID_AR, NsoTable, InLfrCwf.Zv.Epoch, L);
      % clear NsoTable
      % %
      % saturationQualitySchemeId = Bso.get_fv('PROCESSING.SATURATION.QUALITY_SCHEME');
      % l2qbmAr = InLfrCwf.ZvFpa.L2_QUALITY_BITMASK.array(uint16(0));
      % % TEST
      % if 0
      %   l2qbmAr = uint16(zeros(size(l2qbmAr)));
      %   l2qbmAr(1:100) = 1;
      % end
      % SaturationQrcbMap = bicas.proc.L2L3.qual.get_saturation_QRCBs_from_L2QBM(...
      %   l2qbmAr, saturationQualitySchemeId);
      % AllQrcbMap.union(SaturationQrcbMap)



      %=======================================
      % Call BICAS-external code to calculate
      % (1) EFIELD, SCPOT, and from that
      % (2) DENSITY.
      %=======================================
      LfrCwfZv = [];
      LfrCwfZv.Epoch            = InLfrCwf.Zv.Epoch;
      LfrCwfZv.VDC_Fpa          = InLfrCwf.ZvFpa.VDC;
      LfrCwfZv.EDC_Fpa          = InLfrCwf.ZvFpa.EDC;
      LfrCwfZv.QUALITY_FLAG_Fpa = InLfrCwf.ZvFpa.QUALITY_FLAG;
      R = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT_DENSITY(LfrCwfZv, Excd, Bso);



      % EXPERIMENTAL
      % % Update QRCB.
      % AllQrcbMap.set("BAD_DENSITY", R.NeScpQualityBitFpa.array(false));



      %====================================================================
      % Derive values for CDF global attribute "Misc_calibration_versions"
      %====================================================================
      % NOTE: Should not add BICAS version to GA
      % "Misc_calibration_versions" since this is already encoded in GA
      % "Software_version" (together with "Software_name").
      %
      % NOTE: Density "Misc_calibration_versions" contains all three
      % versions, since density is derived from PSP.
      vdccalStr    = ['solo.vdccal() code version ',     R.vdccalCodeVerStr];
      vdccalMatStr = ['solo.vdccal() calibration file ', R.vdccalMatVerStr];
      psp2neStr    = ['solo.psp2ne() code version ',     R.psp2neCodeVerStr];
      gaEfieldScpot_Misc_calibration_versions = {vdccalStr, vdccalMatStr};
      gaDensity_Misc_calibration_versions     = ...
        [gaEfieldScpot_Misc_calibration_versions, {psp2neStr}];



      %================================================================
      % Misc. variables shared between datasets and later modified for
      % specific datasets
      %================================================================
      TemplateOsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.get_OSR_template(InLfrCwf);

      %=======================================
      % Generate data structures for datasets
      %=======================================
      OutEfieldOsr  = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_efield( TemplateOsr, R.EdcSrfMvpmFpa,                       gaEfieldScpot_Misc_calibration_versions);
      OutScpotOsr   = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_scpot(  TemplateOsr, R.ScpotVoltFpa,  R.PspVoltFpa,         gaEfieldScpot_Misc_calibration_versions);
      OutDensityOsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_density(TemplateOsr, R.NeScpCm3Fpa,   R.NeScpQualityBitFpa, gaDensity_Misc_calibration_versions);

      OutEfieldDsr  = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_efield( OutEfieldOsr,  R.EdcSrfMvpmFpa,              L);
      %OutEfieldDsr  = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_efield( TemplateOsr,  R.EdcSrfMvpmFpa,              L);    % TEST
      OutScpotDsr   = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_scpot(  OutScpotOsr,   R.ScpotVoltFpa, R.PspVoltFpa, L);
      OutDensityDsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_density(OutDensityOsr, R.NeScpCm3Fpa,                L);



      nRecordsOsr = size(OutDensityOsr.Zv.Epoch, 1);
      nRecordsDsr = size(OutDensityDsr.Zv.Epoch, 1);
      Tmk.stop_log(nRecordsOsr, 'OSR record', nRecordsDsr, 'DSR record')
    end    % process_L2_to_L3



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Starting template for OSR datasets. Return value is modified.
    function TemplateOsr = get_OSR_template(InLfrCwf)
      Ga = struct();
      Ga.OBS_ID             = InLfrCwf.Ga.OBS_ID;
      Ga.SOOP_TYPE          = InLfrCwf.Ga.SOOP_TYPE;

      Zv = struct();
      Zv.Epoch              = InLfrCwf.Zv.Epoch;
      Zv.DELTA_PLUS_MINUS   = InLfrCwf.ZvFpa.DELTA_PLUS_MINUS;

      %-------------------------------------
      % Copy quality ZVs from input dataset
      %-------------------------------------
      Zv.QUALITY_FLAG       = InLfrCwf.ZvFpa.QUALITY_FLAG;
      % Zv.QUALITY_FLAG       = bicas.utils.FPArray(...
      %   uint8(ones(size(InLfrCwf.ZvFpa.QUALITY_FLAG))) * 2, ...
      %   'NO_FILL_POSITIONS');  % TEST
      Zv.QUALITY_BITMASK    = InLfrCwf.ZvFpa.QUALITY_BITMASK;
      Zv.L2_QUALITY_BITMASK = InLfrCwf.ZvFpa.L2_QUALITY_BITMASK;
      % Zv.L2_QUALITY_BITMASK       = bicas.utils.FPArray(...
      %   uint16(ones(size(InLfrCwf.ZvFpa.L2_QUALITY_BITMASK))) * 0, ...
      %   'NO_FILL_POSITIONS');  % TEST

      TemplateOsr = struct('Ga', Ga, 'Zv', Zv);
    end



    function Out = OSR_efield(TemplateOsr, EdcSrfMvpmFpa, gaMisc_calibration_versions)
      Out = TemplateOsr;
      Out.Ga.Misc_calibration_versions = gaMisc_calibration_versions;

      Out.Zv.EDC_SRF                   = EdcSrfMvpmFpa.cast('single');

      Out = bicas.OutputDataset(Out.Zv, Out.Ga, cell(0,1));
    end



    function Out = OSR_scpot(TemplateOsr, ScpotVoltFpa, PspVoltFpa, gaMisc_calibration_versions)
      Out = TemplateOsr;
      Out.Ga.Misc_calibration_versions = gaMisc_calibration_versions;

      Out.Zv.SCPOT                     = ScpotVoltFpa.cast('single');
      Out.Zv.PSP                       = PspVoltFpa.  cast('single');

      Out = bicas.OutputDataset(Out.Zv, Out.Ga, cell(0,1));
    end



    function Out = OSR_density(TemplateOsr, NeScpCm3Fpa, NeScpQualityBitFpa, gaMisc_calibration_versions)
      Out = TemplateOsr;
      Out.Ga.Misc_calibration_versions = gaMisc_calibration_versions;

      Out.Zv.DENSITY                   = NeScpCm3Fpa.cast('single');

      % NOTE: Behaviour w.r.t. FPs:
      %   Density bit FP ==> L3_QUALITY_BITMASK density bit=false
      %                      (since there is no FP for individual quality bits).
      [QUALITY_FLAG, L3_QUALITY_FLAG] = bicas.proc.L2L3.qual.get_quality_ZVs_density(...
        NeScpQualityBitFpa.array(false));
      Out.Zv.QUALITY_FLAG             = Out.Zv.QUALITY_FLAG.min(QUALITY_FLAG);
      Out.Zv.L3_QUALITY_BITMASK       = bicas.utils.FPArray(L3_QUALITY_FLAG);

      Out = bicas.OutputDataset(Out.Zv, Out.Ga, cell(0,1));
    end



    function Out = DSR_efield(OutEfieldOsr, EdcSrfMvpmOsrFpa, L)
      [TemplateDsrZv, iRecordsInBinCa] = bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template(...
        OutEfieldOsr.Zv.Epoch, ...
        OutEfieldOsr.Zv.QUALITY_FLAG, ...
        OutEfieldOsr.Zv.QUALITY_BITMASK, ...
        OutEfieldOsr.Zv.L2_QUALITY_BITMASK, ...
        bicas.proc.L2L3.L3OsrDsrSwmProcessing.BIN_LENGTH_WOLS_NS, ...
        bicas.proc.L2L3.L3OsrDsrSwmProcessing.BIN_TIMESTAMP_POS_WOLS_NS, L);
      Out = struct('Zv', TemplateDsrZv, 'Ga', OutEfieldOsr.Ga);

      [EdcSrfDsrFpa, EdcstdSrfDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
        EdcSrfMvpmOsrFpa, ...
        bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
        iRecordsInBinCa, L);
      Out.Zv.EDC_SRF    = EdcSrfDsrFpa.   cast('single');
      Out.Zv.EDCSTD_SRF = EdcstdSrfDsrFpa.cast('single');

      Out = bicas.OutputDataset(Out.Zv, Out.Ga, cell(0,1));
    end



    function Out = DSR_scpot(OutScpotOsr, ScpotVoltOsrFpa, PspVoltOsrFpa, L)
      [TemplateDsrZv, iRecordsInBinCa] = bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template(...
        OutScpotOsr.Zv.Epoch, ...
        OutScpotOsr.Zv.QUALITY_FLAG, ...
        OutScpotOsr.Zv.QUALITY_BITMASK, ...
        OutScpotOsr.Zv.L2_QUALITY_BITMASK, ...
        bicas.proc.L2L3.L3OsrDsrSwmProcessing.BIN_LENGTH_WOLS_NS, ...
        bicas.proc.L2L3.L3OsrDsrSwmProcessing.BIN_TIMESTAMP_POS_WOLS_NS, L);
      Out = struct('Zv', TemplateDsrZv, 'Ga', OutScpotOsr.Ga);

      % Downsample SCPOT
      [ScpotDsrFpa, ScpotstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
        ScpotVoltOsrFpa, ...
        bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
        iRecordsInBinCa, ...
        L);
      Out.Zv.SCPOT    = ScpotDsrFpa.   cast('single');
      Out.Zv.SCPOTSTD = ScpotstdDsrFpa.cast('single');

      % Downsample PSP
      [PspDsrFpa, PspstdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
        PspVoltOsrFpa, ...
        bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
        iRecordsInBinCa, L);
      Out.Zv.PSP    = PspDsrFpa.   cast('single');
      Out.Zv.PSPSTD = PspstdDsrFpa.cast('single');

      Out = bicas.OutputDataset(Out.Zv, Out.Ga, cell(0,1));
    end



    function Out = DSR_density(OutDensityOsr, NeScpCm3OsrFpa, L)
      [TemplateDsrZv, iRecordsInBinCa] = bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template(...
        OutDensityOsr.Zv.Epoch, ...
        OutDensityOsr.Zv.QUALITY_FLAG, ...
        OutDensityOsr.Zv.QUALITY_BITMASK, ...
        OutDensityOsr.Zv.L2_QUALITY_BITMASK, ...
        bicas.proc.L2L3.L3OsrDsrSwmProcessing.BIN_LENGTH_WOLS_NS, ...
        bicas.proc.L2L3.L3OsrDsrSwmProcessing.BIN_TIMESTAMP_POS_WOLS_NS, L);
      Out = struct('Zv', TemplateDsrZv, 'Ga', OutDensityOsr.Ga);

      [DensityDsrFpa, DensitystdDsrFpa] = bicas.proc.dsr.downsample_sci_ZV(...
        NeScpCm3OsrFpa, ...
        bicas.const.N_MIN_OSR_SAMPLES_PER_BIN, ...
        iRecordsInBinCa, L);
      Out.Zv.DENSITY    = DensityDsrFpa.   cast('single');
      Out.Zv.DENSITYSTD = DensitystdDsrFpa.cast('single');

      Out.Zv.L3_QUALITY_BITMASK = bicas.proc.dsr.downsample_ZV_bitmask(...
        OutDensityOsr.Zv.L3_QUALITY_BITMASK, iRecordsInBinCa);

      Out = bicas.OutputDataset(Out.Zv, Out.Ga, cell(0,1));
    end



  end    % methods(Static, Access=private)



end
