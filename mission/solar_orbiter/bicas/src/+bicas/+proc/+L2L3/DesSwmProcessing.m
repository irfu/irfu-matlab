%
% SWMP for processing L2 LFR CWF to L3 OSR+DSR density, E field, and ScPot
% ("DES").
%
% NOTE: Is currently also used for experimental L2-->L3 SBMx DENSITY processing.
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
classdef DesSwmProcessing < bicas.proc.SwmProcessing
  % PROPOSAL: Better class name
  %   PROPOSAL: Class name without "L3".
  %     PRO: "L3" is implicit from the package "L2L3".
  %   NOTE: Should cover both SURV and SBMx.
  %   ~Density+Efield+Scpot
  %   ~CWF
  %   ~OSR, DSR
  %   DensityEfieldScpotSwmProcessing
  %   DensEfieldScpotSwmProcessing
  %     CON: Long
  %     CON: Not in order of processing
  %   DesSwmProcessing -- DONE
  %     CON: There is no abbreviation "DES".
  %   CwfDesProcessing
  %
  % PROPOSAL: Automatic test code.
  %   NOTE: There are limited tests.
  %
  % PROPOSAL: Split up processing between (a) density, and (b) E-field & SCPOT
  %           into separate SWMs.
  %   PRO: Faster processing when only processing subset of L3 DSIDs.
  %       CON: Not very heavy operation.
  %   PRO: Leads to better organization of code.
  %       PRO: process_L2_to_L3() is too large and should be split up anyway.
  %   CON: DENSITY is a function EFIELD+SCPOT, and thus has to be processed
  %        after the latter.
  %
  % PROPOSAL: Abolish iSbm.
  %   PRO: Seems unused and not needed.



  %###################
  %###################
  % PRIVATE CONSTANTS
  %###################
  %###################
  properties(Constant, Access=private)



    % Define length of bins, and relative position of corresponding bin
    % timestamps.
    BIN_LENGTH_WOLS_NS        = int64(10e9);
    BIN_TIMESTAMP_POS_WOLS_NS = int64(bicas.proc.L2L3.DesSwmProcessing.BIN_LENGTH_WOLS_NS / 2);



  end



  %############################
  %############################
  % PRIVATE INSTANCE VARIABLES
  %############################
  %############################
  properties(GetAccess=private)
    % Whether to process
    % (true)  L2 SBMx CWF --> L3 SBMx DENSITY, or
    % (false) L2 SURF CWF --> L3 DENSITY+EFIELD+SCPOT
    bSbmx
    % 1/2 = Whether processing SBM1 or SBM2.
    iSbm
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = DesSwmProcessing(bSbmx, iSbm)
      % ASSERTIONS
      assert(islogical(bSbmx) & isscalar(bSbmx))
      if bSbmx
        assert(ismember(iSbm, [1, 2]))
      else
        assert(isempty(iSbm))
      end

      obj.bSbmx = bSbmx;
      obj.iSbm  = iSbm;
    end



    % OVERRIDE
    %
    function OutputDatasetsMap = production_function(obj, ...
        InputDatasetsMap, rctDir, NsoTable, Bso, L)

      % IMPLEMENTATION NOTE: Currently executes all DENSITY+EFIELD+SCPOT OSR+DSR
      % processing both for SURV-CWF (old) and SBMx (experimental) DENSITY
      % processing, and rejects the parts which are not needed. This is
      % technically inefficient, but simplifies the implementation and should be
      % accurate and should probably be kept until there is a good reason to
      % change it(?).
      %   PROPOSAL: Make
      %             bicas.proc.L2L3.DesSwmProcessing.process_L2_to_L3()
      %             skip DSR.

      % PROPOSAL: Create instance variable obj.pfiid in constructor instead.
      if ~obj.bSbmx
        pfiid = 'LFR-SURV-CWF-E_cdf';
      else
        pfiid = 'LFR-SBMx-CWF-E_cdf';
      end
      InputLfrCwfCdf = InputDatasetsMap(pfiid);

      Excd = bicas.proc.L2L3.ExternalCodeImplementation();

      %==============
      % Process data
      %==============
      [...
        EfieldOsrCdf,  EfieldDsrCdf, ...
        ScpotOsrCdf,   ScpotDsrCdf, ...
        DensityOsrCdf, DensityDsrCdf ...
        ] = obj.process_L2_to_L3(InputLfrCwfCdf, NsoTable, Excd, Bso, L);

      OutputDatasetsMap = containers.Map();
      if ~obj.bSbmx
        OutputDatasetsMap('EFIELD_OSR_cdf')   = EfieldOsrCdf;
        OutputDatasetsMap('EFIELD_DSR_cdf')   = EfieldDsrCdf;
        OutputDatasetsMap('SCPOT_OSR_cdf')    = ScpotOsrCdf;
        OutputDatasetsMap('SCPOT_DSR_cdf')    = ScpotDsrCdf;
        OutputDatasetsMap('DENSITY_OSR_cdf')  = DensityOsrCdf;
        OutputDatasetsMap('DENSITY_DSR_cdf')  = DensityDsrCdf;
      else
        % NOTE: Ignoring all other return data.
        OutputDatasetsMap('SBMx_EFIELD_cdf')  = EfieldOsrCdf;
        OutputDatasetsMap('SBMx_SCPOT_cdf')   = ScpotOsrCdf;
        OutputDatasetsMap('SBMx_DENSITY_cdf') = DensityOsrCdf;
      end
    end



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
    % bicas.proc.L2L3.DesSwmProcessing.production_function() to facilitate
    % automated tests, in particular by adding an explicit dependence on
    % bicas.proc.L2L3.ExternalCodeImplementation.
    %
    function [...
        OutEfieldOsr,  OutEfieldDsr, ...
        OutScpotOsr,   OutScpotDsr, ...
        OutDensityOsr, OutDensityDsr ...
        ] ...
        = process_L2_to_L3(obj, InLfrCwf, NsoTable, Excd, Bso, L)



      saturationSchemeId  = Bso.get_fv('PROCESSING.SATURATION.QUALITY_SCHEME');
      inputSkeletonVerNbr = str2double(InLfrCwf.Ga.Skeleton_version{1});
      % --------------------------------------------------------------------
      % ASSERTION: Check that the input datasets use the correct quality bit
      %            saturation scheme.
      % --------------------------------------------------------------------
      % Official SOLO_L2_RPW-LFR-SURV/SBM1/SBM2-CWF-E datasets use
      % GLOBAL_SATURATION starting at Skeleton_version = 18.
      %
      % IMPLEMENTATION NOTE: ROC has not yet reprocessed SBMx to have
      % inputSkeletonVerNbr>=18. This is due to a mistake in the EJ-ROC
      % communication. Not yet checking for SBMx so that code can be tested at
      % least. /2026-01-14
      if ~obj.bSbmx
        switch saturationSchemeId
          case "CHANNEL_SATURATION"
            if inputSkeletonVerNbr < 18
              error( ...
                "The input CDF has Skeleton_version=%d indicating that it" + ...
                " does not use CHANNEL_SATURATION and is therefore" + ...
                " incompatible.", ...
                inputSkeletonVerNbr)
            end

          case "GLOBAL_SATURATION"
            if inputSkeletonVerNbr >= 18
              error( ...
                "The input CDF has Skeleton_version=%d indicating that it" + ...
                " does not use GLOBAL_SATURATION and is therefore" + ...
                " incompatible.", ...
                inputSkeletonVerNbr)
            end

          otherwise
            error("Illegal value of saturationSchemeId=%s.", saturationSchemeId)
        end    % switch
      end    % if



      % PROPOSAL: Abolish struct "Ga".
      % PROPOSAL: Abolish struct "Zv".
      %   PRO: Not passed to any function etc.
      %   PRO: Shorter.
      %   CON: Makes it less clear which variables are ZV-like and not. Prefix
      %        "Zv." effectively constitutes an alternative prefix stating
      %        this.
      %     CON: This convention is not followed anywhere else.
      %   CON-PROPOSAL: Rename argument InLfrCwf to something short and never
      %                 rename it in the code.

      Ga = struct();
      Ga.SOOP_TYPE            = InLfrCwf.Ga.SOOP_TYPE;
      Ga.OBS_ID               = InLfrCwf.Ga.OBS_ID;
      %
      Zv = struct();
      Zv.tt2000               = InLfrCwf.Zv.Epoch;
      Zv.VDC_Fpa              = InLfrCwf.ZvFpa.VDC;
      Zv.EDC_Fpa              = InLfrCwf.ZvFpa.EDC;
      Zv.QflFpa               = InLfrCwf.ZvFpa.QUALITY_FLAG;
      Zv.L1qbmFpa             = InLfrCwf.ZvFpa.QUALITY_BITMASK;
      Zv.L2qbmFpa             = InLfrCwf.ZvFpa.L2_QUALITY_BITMASK;
      Zv.DELTA_PLUS_MINUS_Fpa = InLfrCwf.ZvFpa.DELTA_PLUS_MINUS;
      clear InLfrCwf

      assert(isa(NsoTable, "bicas.NsoTable"))



      % PROPOSAL: Split up into different parts for EFIELD, SCPOT, DENSITY
      %           (still combine non-downsampled and downsampled).
      %   CON: Slows down overall processing.
      %       PRO: Must read same L2 dataset multiple times.
      %       PRO: Must read L3 SCPOT dataset to produce L3 DENSITY dataset.
      %   CON: There is much shared functionality for 3 quality ZVs.
      %       PRO: Same ~constants
      %           Ex: INPUT_DSID, BIN_LENGTH_WOLS_NS, BIN_TIMESTAMP_POS_WOLS_NS
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
        'bicas.proc.L2L3.DesSwmProcessing.process_L2_to_L3', L);
      assert(isa(Excd, 'bicas.proc.L2L3.ExternalCodeAbstract'))



      %------------------------------------
      % Derive QRCBMs and synthetic L2 QFL
      %------------------------------------
      [L3Qrcbm, L3DensityQrcbm, SyntheticL2QflFpa] = ...
        bicas.proc.L2L3.DesSwmProcessing.get_QRCBMs_synthetic_L2_QFL(...
        ...
        Zv.L2qbmFpa, Zv.QflFpa, Zv.tt2000, NsoTable, ...
        Bso.get_fv('PROCESSING.SATURATION.QUALITY_SCHEME'), L);

      %---------------------------------
      % Blank selected VDC, EDC samples
      %---------------------------------
      qflMinForUse = uint8(Bso.get_fv('PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN'));
      [VDC_Fpa, EDC_Fpa] = ...
        bicas.proc.L2L3.DesSwmProcessing.set_VDC_EDC_FPs_before_processing( ...
        ...
        Zv.VDC_Fpa, Zv.EDC_Fpa, SyntheticL2QflFpa, ...
        L3Qrcbm, qflMinForUse);



      %---------------------------------------
      % Call BICAS-external code to calculate
      % (1) EFIELD, SCPOT, and from that
      % (2) DENSITY.
      %---------------------------------------
      R = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT_DENSITY(Excd, ...
        tt2000  = Zv.tt2000, ...
        VDC_Fpa = VDC_Fpa, ...
        EDC_Fpa = EDC_Fpa, ...
        L       = L);

      %-----------------------------
      % Blank output based on QRCBs
      %-----------------------------
      % Blank EFIELD.
      R.EdcSrfMvpmFpa = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        R.EdcSrfMvpmFpa, L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "efieldFvIndexAr");
      % Blank DENSITY (but not density quality bit).
      R.NeScpCm3Fpa   = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        R.NeScpCm3Fpa,   L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "densityFvIndexAr");
      % Blank PSP and SCPOT.
      R.PspVoltFpa    = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        R.PspVoltFpa,    L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "scpotFvIndexAr");
      R.ScpotVoltFpa  = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        R.ScpotVoltFpa,  L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "scpotFvIndexAr");



      %----------------------------------------------------------------------
      % Set L3 density BAD_DENSITY QRCB using information from solo.psp2ne()
      %----------------------------------------------------------------------
      % NOTE: Behaviour w.r.t. FPs:
      %   Density bit FP ==> L3_QUALITY_BITMASK density bit=false
      %                      (since there is no FP for individual quality bits).
      L3DensityQrcbm.set("BAD_DENSITY", R.NeScpQualityBitFpa.array(false));



      %================================================
      % Derive GA values which apply to all 6x L3 CDFs
      %================================================
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
      %
      gaL3Caveats = bicas.proc.qrc.QRCB_arrays_to_GA_CAVEATS( ...
        L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM);



      %================================================================
      % Misc. variables shared between datasets and later modified for
      % specific datasets
      %================================================================
      TemplateOsr = bicas.proc.L2L3.DesSwmProcessing.get_OSR_template( ...
        OBS_ID               = Ga.OBS_ID, ...
        SOOP_TYPE            = Ga.SOOP_TYPE, ...
        CAVEATS              = gaL3Caveats, ...
        tt2000               = Zv.tt2000, ...
        DELTA_PLUS_MINUS_Fpa = Zv.DELTA_PLUS_MINUS_Fpa, ...
        QflFpa               = SyntheticL2QflFpa, ...
        L1qbmFpa             = Zv.L1qbmFpa, ...
        L2qbmFpa             = Zv.L2qbmFpa);

      %=======================================
      % Generate data structures for datasets
      %=======================================
      OutEfieldOsr  = bicas.proc.L2L3.DesSwmProcessing.OSR_efield( TemplateOsr, R.EdcSrfMvpmFpa,                 gaEfieldScpot_Misc_calibration_versions);
      OutScpotOsr   = bicas.proc.L2L3.DesSwmProcessing.OSR_scpot(  TemplateOsr, R.ScpotVoltFpa,  R.PspVoltFpa,   gaEfieldScpot_Misc_calibration_versions);
      OutDensityOsr = bicas.proc.L2L3.DesSwmProcessing.OSR_density(TemplateOsr, R.NeScpCm3Fpa,   L3DensityQrcbm, gaDensity_Misc_calibration_versions);

      OutEfieldDsr  = bicas.proc.L2L3.DesSwmProcessing.DSR_efield( OutEfieldOsr,  R.EdcSrfMvpmFpa,              L);
      OutScpotDsr   = bicas.proc.L2L3.DesSwmProcessing.DSR_scpot(  OutScpotOsr,   R.ScpotVoltFpa, R.PspVoltFpa, L);
      OutDensityDsr = bicas.proc.L2L3.DesSwmProcessing.DSR_density(OutDensityOsr, R.NeScpCm3Fpa,                L);



      nRecordsOsr = size(OutDensityOsr.Zv.Epoch, 1);
      nRecordsDsr = size(OutDensityDsr.Zv.Epoch, 1);
      Tmk.stop_log(nRecordsOsr, 'OSR record', nRecordsDsr, 'DSR record')
    end    % process_L2_to_L3()



  end    % methods(Access=public)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Function to group together related code in process_L2_to_L3().
    %
    % IMPLEMENTATION NOTE: In principle, the call to
    % bicas.proc.L2L3.qrc.get_synthetic_L2_QFL() is
    % independent of the function, and does in a sense does not belong there,
    % but the need to group together code is greater and takes precedence.
    %
    function [L3Qrcbm, L3DensityQrcbm, SyntheticL2QflFpa] = ...
        get_QRCBMs_synthetic_L2_QFL(...
        ...
        L2qbmFpa, L2QflFpa, tt2000Ar, NsoTable, ...
        saturationQualityScheme, L)
      % PROPOSAL: Move to bicas.proc.L2L3.qrc.

      %=========================================
      % NSO table, L2_QUALITY_BITMASK --> QRCBs
      %=========================================
      %------------------------------
      % NSO table-->L3 QRCBs
      %             L3 DENSITY QRCBs
      %------------------------------
      L3Qrcbm        = bicas.proc.qrc.NSO_table_to_QRCBM(...
        bicas.const.qrc.Q.L3_QRCSM.qrcidAr,         NsoTable, tt2000Ar, L);
      L3DensityQrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        bicas.const.qrc.Q.L3_DENSITY_QRCSM.qrcidAr, NsoTable, tt2000Ar, L);

      %-------------------------------
      % L2_QUALITY_BITMASK-->L3 QRCBs
      %-------------------------------
      % Reconstruct L2 channel saturation QRCBs by reading L2_QUALITY_BITMASK
      % *AND* converting it to an QRCBM.
      ChannelSaturationQrcbm = ...
        bicas.proc.L2L3.qrc.L2QBM_to_channel_saturation_QRCBs(...
        ...
        L2qbmFpa.array(uint16(0)), ...
        saturationQualityScheme);
      L3Qrcbm.union(ChannelSaturationQrcbm)



      %=========================================
      % Calculate a "synthetic L2 QUALITY_FLAG"
      %=========================================
      SyntheticL2QflFpa = bicas.proc.L2L3.qrc.get_synthetic_L2_QFL(...
        tt2000Ar, NsoTable, ChannelSaturationQrcbm, L2QflFpa, L);
    end



    % Function to group together related code in process_L2_to_L3().
    % Remove (blank) L2 data before sending it to the external L2-to-L3
    % calibration.
    %
    function [VDC_Fpa, EDC_Fpa] = set_VDC_EDC_FPs_before_processing( ...
        VDC_Fpa, EDC_Fpa, SyntheticL2QflFpa, ...
        L3Qrcbm, qflMinForUse)
      % PROPOSAL: Test code.

      % irf.assert.sizes(...
      %   VDC_Fpa, [-1, 3], ...
      %   EDC_Fpa, [-1, 3])
      % bicas.utils.validate_QFL(QUALITY_FLAG_minForUse)
      % assert(isscalar(QUALITY_FLAG_minForUse))
      % assert(isa(L3Qrcbm, "bicas.proc.QrcbMap"))

      %--------------------------------
      % Blank input data based on QRCs
      %--------------------------------
      % bicas.debug.plot_VDC_EDC_FPA(Zv.VDC_Fpa, Zv.EDC_Fpa, Zv.tt2000, "Before QRC blanking")   % DEBUG
      VDC_Fpa = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        VDC_Fpa, L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "vdcFvIndexAr");
      EDC_Fpa = bicas.proc.L2L3.qrc.set_FPA_samples_FP(...
        EDC_Fpa, L3Qrcbm, bicas.const.qrc.Q.L3_QRCSM, "edcFvIndexAr");
      % bicas.debug.plot_VDC_EDC_FPA(VDC_Fpa, EDC_Fpa, Zv.tt2000, "After QRC blanking")    % DEBUG

      %-------------------------------------------------------
      % Blank input data when QUALITY_FLAG is below threshold
      %-------------------------------------------------------
      % IMPORTANT NOTE: Uses derived QUALITY_FLAG_nonsatFpa (not QUALITY_FLAG
      % from L2 input file)!
      % NOTE: It is unclear what is the best way to treat QUALITY_FLAG=FV.
      % NOTE: Treatment of this special case is documented in readme.txt.
      bDoNotUseFpa = SyntheticL2QflFpa < qflMinForUse;
      bDoNotUse    = bDoNotUseFpa.array(false);   % Is [FP==>false] wise?
      VDC_Fpa(bDoNotUse, :) = bicas.utils.FPArray.FP_SINGLE;
      EDC_Fpa(bDoNotUse, :) = bicas.utils.FPArray.FP_SINGLE;
    end



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
        Ga.CAVEATS
        Zv.tt2000
        Zv.DELTA_PLUS_MINUS_Fpa
        Zv.QflFpa
        Zv.L1qbmFpa
        Zv.L2qbmFpa
      end
      % NOTE: Copies "Ga" directly.

      Zv2 = struct();
      Zv2.Epoch              = Zv.tt2000;
      Zv2.DELTA_PLUS_MINUS   = Zv.DELTA_PLUS_MINUS_Fpa;

      Zv2.QUALITY_FLAG       = Zv.QflFpa;
      Zv2.QUALITY_BITMASK    = Zv.L1qbmFpa;
      Zv2.L2_QUALITY_BITMASK = Zv.L2qbmFpa;

      %bicas.debug.plot_FPA(Zv.tt2000, Zv.QflFpa, uint8(-1), "QflFpa")

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

      [qfl, l3qbm] = bicas.proc.qrc.QRCB_arrays_to_quality_ZVs(...
        L3DensityQrcbm, bicas.const.qrc.Q.L3_DENSITY_QRCSM, "L3_QUALITY_BITMASK");
      gaCaveats = bicas.proc.qrc.QRCB_arrays_to_GA_CAVEATS( ...
        L3DensityQrcbm, bicas.const.qrc.Q.L3_DENSITY_QRCSM);

      Out.Zv.QUALITY_FLAG             = Out.Zv.QUALITY_FLAG.min(qfl);
      Out.Zv.L3_QUALITY_BITMASK       = bicas.utils.FPArray(l3qbm);
      Out.Ga.CAVEATS                  = [Out.Ga.CAVEATS; gaCaveats];

      Out = bicas.OutputDataset(Out.Zv, Out.Ga, cell(0,1));
    end



    function Out = DSR_efield(OutEfieldOsr, EdcSrfMvpmOsrFpa, L)
      [TemplateDsrZv, iRecordsInBinCa] = bicas.proc.dsr.get_LFR_CWF_DSR_ZVs_template(...
        OutEfieldOsr.Zv.Epoch, ...
        OutEfieldOsr.Zv.QUALITY_FLAG, ...
        OutEfieldOsr.Zv.QUALITY_BITMASK, ...
        OutEfieldOsr.Zv.L2_QUALITY_BITMASK, ...
        bicas.proc.L2L3.DesSwmProcessing.BIN_LENGTH_WOLS_NS, ...
        bicas.proc.L2L3.DesSwmProcessing.BIN_TIMESTAMP_POS_WOLS_NS, L);
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
        bicas.proc.L2L3.DesSwmProcessing.BIN_LENGTH_WOLS_NS, ...
        bicas.proc.L2L3.DesSwmProcessing.BIN_TIMESTAMP_POS_WOLS_NS, L);
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
        bicas.proc.L2L3.DesSwmProcessing.BIN_LENGTH_WOLS_NS, ...
        bicas.proc.L2L3.DesSwmProcessing.BIN_TIMESTAMP_POS_WOLS_NS, L);
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
