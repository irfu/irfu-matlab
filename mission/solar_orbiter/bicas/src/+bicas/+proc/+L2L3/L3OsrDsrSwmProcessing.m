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
  % PROPOSAL: Class name without "L3".
  %   PRO: "L3" is implicit from the package "L2L3".
  %
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



    % Define length of bins, and relative position of corresponding bin
    % timestamps.
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
        InputLfrCwfCdf, NsoTable, Excd, Bso, L);

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
        = process_L2_to_L3(InLfrCwf, NsoTable, Excd, Bso, L)

      % PROPOSAL: Abolish struct "Ga".
      % PROPOSAL: Abolish struct "Zv".
      %   PRO: Not passed to any function etc.
      %   PRO: Shorter.
      %   CON: Makes it less clear which variables are ZV-like and not. Prefix
      %        "Zv." effectively constitutes an alternative prefix stating
      %        this.
      %     CON: This convention is not followed anywhere else.

      Ga = struct();
      Ga.SOOP_TYPE              = InLfrCwf.Ga.SOOP_TYPE;
      Ga.OBS_ID                 = InLfrCwf.Ga.OBS_ID;
      Zv = struct();            % ----------------------
      Zv.Epoch                  = InLfrCwf.Zv.Epoch;
      Zv.VDC_Fpa                = InLfrCwf.ZvFpa.VDC;
      Zv.EDC_Fpa                = InLfrCwf.ZvFpa.EDC;
      Zv.QUALITY_FLAG_Fpa       = InLfrCwf.ZvFpa.QUALITY_FLAG;
      Zv.QUALITY_BITMASK_Fpa    = InLfrCwf.ZvFpa.QUALITY_BITMASK;
      Zv.L2_QUALITY_BITMASK_Fpa = InLfrCwf.ZvFpa.L2_QUALITY_BITMASK;
      Zv.DELTA_PLUS_MINUS_Fpa   = InLfrCwf.ZvFpa.DELTA_PLUS_MINUS;
      clear InLfrCwf

      assert(isa(NsoTable, "bicas.NsoTable"))



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

      Tmk = bicas.utils.Timekeeper(...
        'bicas.proc.L2L3.L3OsrDsrSwmProcessing.process_L2_to_L3', L);
      assert(isa(Excd, 'bicas.proc.L2L3.ExternalCodeAbstract'))



      %=========================================
      % NSO table, L2_QUALITY_BITMASK --> QRCBs
      %=========================================
      %------------------------------
      % NSO table-->L3 QRCBs
      %             L3 DENSITY QRCBs
      %------------------------------
      L3Qrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        bicas.const.qrc.Q.L3_QRCSM.qrcidAr, ...
        NsoTable, Zv.Epoch, L);
      L3DensityQrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        bicas.const.qrc.Q.L3_DENSITY_QRCSM.qrcidAr, ...
        NsoTable, Zv.Epoch, L);

      %-------------------------------
      % L2_QUALITY_BITMASK-->L3 QRCBs
      %-------------------------------
      % Obtain channel saturation QRCBs by reading L2_QUALITY_BITMASK *AND*
      % converting it to an QRCBM.
      ChannelSaturationQrcbm = ...
        bicas.proc.L2L3.qrc.L2QBM_to_channel_saturation_QRCBs(...
        Zv.L2_QUALITY_BITMASK_Fpa.array(uint16(0)), ...
        Bso.get_fv('PROCESSING.SATURATION.QUALITY_SCHEME'));
      L3Qrcbm.union(ChannelSaturationQrcbm)
      if 0
      % DEBUG
        L3Qrcbm.set("SATURATION_ZV_V2", true(size(Zv.Epoch)))
        bicas.debug.plot_QRCBM(L3Qrcbm, Zv.Epoch, "L3Qrcbm")
      end



      %-----------------------------------------------------------------------
      % Calculate what QUALITY_FLAG for L2 LFR CWF *SHOULD HAVE BEEN*, had it
      % not been for saturation or sweeps (sic!).
      %-----------------------------------------------------------------------
      % NOTE: Sweeps are blanked (in L2), so their QUALITY_FLAG values do not
      % matter.
      % IMPLEMENTATION NOTE: Must distinguish between
      % (1) the derived L2 non-saturation QUALITY_FLAG, and
      % (2) the (true) L2 input QUALITY FLAG.
      Zv.QUALITY_FLAG_nonsatFpa = ...
        bicas.proc.L2L3.qrc.get_L2_nonsaturation_nonsweep_QUALITY_FLAG( ...
        Zv.Epoch, NsoTable, Zv.QUALITY_FLAG_Fpa.fpAr, L);



      %==========================
      % Call BICAS-external code
      %==========================

      %--------------------------------
      % Blank input data based on QRCs
      %--------------------------------
      % bicas.debug.plot_VDC_EDC_FPA(Zv.VDC_Fpa, Zv.EDC_Fpa, Zv.Epoch, "Before QRC blanking")   % DEBUG
      VDC_Fpa = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        Zv.VDC_Fpa, L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "vdcFvIndexAr");
      EDC_Fpa = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        Zv.EDC_Fpa, L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "edcFvIndexAr");
      % bicas.debug.plot_VDC_EDC_FPA(VDC_Fpa, EDC_Fpa, Zv.Epoch, "After QRC blanking")    % DEBUG

      %-------------------------------------------------------
      % Blank input data when QUALITY_FLAG is below threshold
      %-------------------------------------------------------
      % IMPORTANT NOTE: Uses derived QUALITY_FLAG_nonsatFpa (not QUALITY_FLAG
      % from L2 input file)!
      % NOTE: It is unclear what is the best way to treat QUALITY_FLAG=FV.
      % NOTE: Treatment of this special case is documented in readme.txt.
      QUALITY_FLAG_minForUse = uint8(Bso.get_fv(...
        'PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN'));
      bDoNotUseFpa = Zv.QUALITY_FLAG_nonsatFpa < QUALITY_FLAG_minForUse;
      bDoNotUse    = bDoNotUseFpa.array(false);   % Is [FP==>false] wise?
      VDC_Fpa(bDoNotUse, :) = bicas.utils.FPArray.FP_SINGLE;
      EDC_Fpa(bDoNotUse, :) = bicas.utils.FPArray.FP_SINGLE;

      %---------------------------------------
      % Call BICAS-external code to calculate
      % (1) EFIELD, SCPOT, and from that
      % (2) DENSITY.
      %---------------------------------------
      R = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT_DENSITY(Excd, ...
        Epoch  =Zv.Epoch, ...
        VDC_Fpa=VDC_Fpa, ...
        EDC_Fpa=EDC_Fpa);

      %------------------------------------
      % Blank output EFIELD based on QRCBs
      %------------------------------------
      % IMPLEMENTATION NOTE: Can add blanking of other L3 variables too, if
      % needed.
      R.EdcSrfMvpmFpa = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        R.EdcSrfMvpmFpa, L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "efieldFvIndexAr");

      %----------------------------------------------------------------------
      % Set L3 density BAD_DENSITY QRCB using information from solo.psp2ne()
      %----------------------------------------------------------------------
      % NOTE: Behaviour w.r.t. FPs:
      %   Density bit FP ==> L3_QUALITY_BITMASK density bit=false
      %                      (since there is no FP for individual quality bits).
      L3DensityQrcbm.set("BAD_DENSITY", R.NeScpQualityBitFpa.array(false));



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
      TemplateOsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.get_OSR_template( ...
        OBS_ID                =Ga.OBS_ID, ...
        SOOP_TYPE             =Ga.SOOP_TYPE, ...
        Epoch                 =Zv.Epoch, ...
        DELTA_PLUS_MINUS_Fpa  =Zv.DELTA_PLUS_MINUS_Fpa, ...
        QUALITY_FLAG_Fpa      =Zv.QUALITY_FLAG_nonsatFpa, ...
        QUALITY_BITMASK_Fpa   =Zv.QUALITY_BITMASK_Fpa, ...
        L2_QUALITY_BITMASK_Fpa=Zv.L2_QUALITY_BITMASK_Fpa);

      %=======================================
      % Generate data structures for datasets
      %=======================================
      OutEfieldOsr  = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_efield( TemplateOsr, R.EdcSrfMvpmFpa,                 gaEfieldScpot_Misc_calibration_versions);
      OutScpotOsr   = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_scpot(  TemplateOsr, R.ScpotVoltFpa,  R.PspVoltFpa,   gaEfieldScpot_Misc_calibration_versions);
      OutDensityOsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.OSR_density(TemplateOsr, R.NeScpCm3Fpa,   L3DensityQrcbm, gaDensity_Misc_calibration_versions);

      OutEfieldDsr  = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_efield( OutEfieldOsr,  R.EdcSrfMvpmFpa,              L);
      OutScpotDsr   = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_scpot(  OutScpotOsr,   R.ScpotVoltFpa, R.PspVoltFpa, L);
      OutDensityDsr = bicas.proc.L2L3.L3OsrDsrSwmProcessing.DSR_density(OutDensityOsr, R.NeScpCm3Fpa,                L);



      nRecordsOsr = size(OutDensityOsr.Zv.Epoch, 1);
      nRecordsDsr = size(OutDensityDsr.Zv.Epoch, 1);
      Tmk.stop_log(nRecordsOsr, 'OSR record', nRecordsDsr, 'DSR record')
    end    % process_L2_to_L3()



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Starting template for OSR datasets. Return value is modified.
    function TemplateOsr = get_OSR_template(Ga, Zv)
      % PROPOSAL: Abolish function?!
      %   PRO: The function does very little.
      %     PRO: Call to function is as long as replacing the function call with
      %          code.
      %     CON: The function renames variables.
      %     CON: The function isolates a conceptual unit of code.
      arguments
        Ga.OBS_ID
        Ga.SOOP_TYPE
        Zv.Epoch
        Zv.DELTA_PLUS_MINUS_Fpa
        Zv.QUALITY_FLAG_Fpa
        Zv.QUALITY_BITMASK_Fpa
        Zv.L2_QUALITY_BITMASK_Fpa
      end
      % NOTE: Copies "Ga" directly.

      Zv2 = struct();
      Zv2.Epoch              = Zv.Epoch;
      Zv2.DELTA_PLUS_MINUS   = Zv.DELTA_PLUS_MINUS_Fpa;

      Zv2.QUALITY_FLAG       = Zv.QUALITY_FLAG_Fpa;
      Zv2.QUALITY_BITMASK    = Zv.QUALITY_BITMASK_Fpa;
      Zv2.L2_QUALITY_BITMASK = Zv.L2_QUALITY_BITMASK_Fpa;

      TemplateOsr = struct('Ga', Ga, 'Zv', Zv2);
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



    function Out = OSR_density(TemplateOsr, NeScpCm3Fpa, L3DensityQrcbm, gaMisc_calibration_versions)
      Out = TemplateOsr;
      Out.Ga.Misc_calibration_versions = gaMisc_calibration_versions;

      Out.Zv.DENSITY                   = NeScpCm3Fpa.cast('single');

      [QUALITY_FLAG, L3_QUALITY_FLAG] = bicas.proc.qrc.QRCB_arrays_to_quality_ZVs(...
        L3DensityQrcbm, bicas.const.qrc.Q.L3_DENSITY_QRCSM, "L3_QUALITY_BITMASK");

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
