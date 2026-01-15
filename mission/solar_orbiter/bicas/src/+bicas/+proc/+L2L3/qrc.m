%
% Collection of code relating to QRCs for L2 to L3 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qrc



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Given some L2 QRCSs which set unique L2QBM bits, set (L3) QRCBs (with the
    % same QRCIDs) from those quality bits. Asserts that this is easily doable
    % by requiring all QRCSs to set exactly one quality bit which is unique.
    %
    %
    % ARGUMENTS
    % =========
    % L2QualityBitSettingQrcsm
    %       Must ONLY contain those QRCSs which quality bits should be read
    %       from L2QBM.
    %
    function Qrcbm = L2QBM_to_QRCBs(l2qbmAr, L2QualityBitSettingQrcsm)
      % IMPLEMENTATION NOTE: Function is designed for
      % SOLO_L2_RPW-LFR-SURV-CWF-E's L2_QUALITY_BITMASK but is generalized to
      % arbitrary QRCSs, partly to simplify testing and not mix algorithms with
      % global constants.

      assert(isa(L2QualityBitSettingQrcsm, "bicas.proc.QrcSettingsMap"))
      assert(iscolumn(l2qbmAr) & isa(l2qbmAr, "uint16"))

      Qrcbm       = bicas.proc.QrcbMap(numel(l2qbmAr));
      % Collection of quality bit positions for all used QRCDSs. Only used for
      % asserting against collisions.
      allBitPosAr = zeros(0, 1);

      for qrcid = L2QualityBitSettingQrcsm.qrcidAr'
        Qrcs = L2QualityBitSettingQrcsm.get(qrcid);
        assert(isa(Qrcs, "bicas.proc.QrcSettingL2"))

        qrcLxqbm = Qrcs.("l2qbm");
        bitPosAr = bicas.proc.qrc.LxQBM_to_bit_positions(qrcLxqbm);
        assert(isscalar(bitPosAr), "QRC does not set exactly one bit.")

        allBitPosAr(end+1, 1) = bitPosAr;

        qrcbAr = logical(bitand(l2qbmAr, qrcLxqbm));
        Qrcbm.add(qrcid, qrcbAr)
      end

      % ASSERTION: No overlap in (specified) QRCS quality bits.
      % NOTE: Can not verify against quality bits in QRCS not submitted to this
      %       function.
      irf.assert.number_set(allBitPosAr)
    end



    % Given L2QBM, derive the corresponding CHANNEL_SATURATION QRCBs.
    %
    %
    % RATIONALE
    % =========
    % Thisis necessary for obtaining QRCBs used for blanking science data sent
    % to EXCD (solo.vdccal(), solo.psp2ne()).
    %
    %
    % NOTE: Always creates QRCBs for CHANNEL_SATURATION QRCIDs but never for
    % GLOBAL_SATURATION (since no action is taken for them in L2-->L3
    % processing no matter what). Sets the QRCB arrays differently depending on
    % "saturationQualitySchemeId".
    %
    function SaturationQrcbm = L2QBM_to_channel_saturation_QRCBs(...
        l2qbmAr, saturationQualitySchemeId)

      assert(iscolumn(l2qbmAr) & isa(l2qbmAr, "uint16"))

      SaturationQrcbm = bicas.proc.QrcbMap(numel(l2qbmAr));
      SaturationQrcbm.add_false(bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR)

      switch saturationQualitySchemeId
        case "CHANNEL_SATURATION"

          % Update CHANNEL_SATURATION QRCBs.
          ChannelSaturationQrcbm = bicas.proc.L2L3.qrc.L2QBM_to_QRCBs(...
            l2qbmAr, bicas.const.qrc.Q.L2_CHANNEL_SATURATION_QRCSM);
          SaturationQrcbm.union(ChannelSaturationQrcbm)

        case "GLOBAL_SATURATION"

          ;   % Do nothing. All saturation QRCBs are false.

          % IMPLEMENTATION NOTE: The functionality is technically incorrect
          % since this option will not set sensible channel saturation QRCBs.
          % This is not overly important though, since
          % (1) One can not use the global saturation for anything sensible
          % (can not blank selected channels), and (2) the GLOBAL_SATURATION
          % functionality is supposed to be phased out anyway.
          %
          % PROPOSAL: Have FULL_SATURATION ==> Set all CHANNEL_SATURATION for
          %           all channels.

        otherwise

          error("BICAS:ConfigurationBug", ...
            "Illegal argument saturationQualitySchemeId=""%s"".", ...
            saturationQualitySchemeId)
      end

      assert(isequal( ...
        sort(SaturationQrcbm.qrcidAr), ...
        sort(bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR)))
    end



    % Derive a (partly) "synthetic" L2 LFR CWF QUALITY_FLAG which is similar to
    % the true L2 QUALITY_FLAG but which tries to emulate what the true L2
    % QUALITY_FLAG would have been if it had not been lowered due to (1)
    % SATURATION_ZV_V*, and (2) ANT3_FAILING.
    %
    %
    % RATIONALE
    % =========
    % The value is intended to be used as input to processing L2-->L3, both
    % (a) for filtering data (treshold; settting
    %     "PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN")), and
    % (b) as starting point for deriving the L3 QFL which may be higher than the
    %     corresponding actual L2 QUALITY FLAG when L3 contains valid values
    %     derived from a subset of good data, but ignores a subset of bad data
    %     which was the reason for lowering L2 QFL, but which is therefore not a
    %     reason for lowering L3 QFL.
    %     Ex: Using non-saturated L2 channels in the presence of other saturated
    %         L2 channels.
    %
    % NOTE: This function is a "better hack" since it tries to solve a currently
    % impossible task. It is to be used only while waiting for a better
    % solution to become available.
    %
    %
    % STRENGTHS, LIMITATIONS, RISKS
    % =============================
    % The actual ideal synthetic QFL value is impossible to derive, since not
    % all the information used when deriving L1R-->L2 is available when
    % processing L2-->L3.
    %
    % QRC events which *CAN* be accounted for are:
    % (1) QRC events in the NSO table.
    % (2) QRC events which are stored in the L2 CDF
    %     Ex: SATURATION_ZV_*
    %
    % QRC events which can *NOT* be accounted for are:
    % (1) QRC events which are autodetected during processing L1R-->L2 but which
    %     are not stored in the L2 CDF.
    %     Ex: SWEEP:
    %       This is OK since the samples are blanked anyway, though the QFL
    %       could be off.
    % (2) QRC(B)s which are read from data sources which are available during
    %     processing L1R-->L2, but not L2-->L3
    %     Ex: BIAS_HW_OFF (QRCB read from LFR L1R CDF):
    %           This is OK since the samples are blank anyway, though the QFL
    %           could be off.
    %
    % NOTE/BUG: This function can potentially set a maximum QUALITY_FLAG value
    % which is higher than what L1R could possibly have (not the QUALITY_FLAG
    % max=4, but the highest QUALITY_FLAG which ROC wants to produce).
    % QUALITY_FLAG is eventually capped by setting
    % "PROCESSING.ZV_QUALITY_FLAG_MAX" when BICAS writes the CDF anyway, but
    % this setting then has to use the appropriate value.
    %
    % NOTE: The main failure mode happens when L1R caps QFL at the same time as
    % BICAS caps L2 QFL (to the same value). Then this function will still
    % "uncap" the synthetic L2 QFL and raise the value.
    %
    %
    % RATIONALE
    % =========
    % This synthetic L2 QUALITY FLAG is needed for
    % (1) filtering data L2-->L3 (settting
    % "PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN"), and
    % (2) an input value that can be used for deriving the final L3 QUALITY_FLAG
    % which may be higher than the corresponding actual L2 QUALITY FLAG when L3
    % contains valid values derived from e.g. non-saturated L2 channels in the
    % presence of other saturated L2 channels.
    %
    function SyntheticL2QflFpa = get_synthetic_L2_QFL( ...
        tt2000Ar, NsoTable, ChannelSaturationQrcbm, L2QflFpa, L)

      % NOTE: One could consider also removing/excluding ANT3_FAILING, since
      % the unaffected channels should be OK. Has no instructions to do so yet
      % though. /2025-08-27

      % QRCs which, when they apply, imply that the (completely) synthetic QFL
      % should be used.
      EXCEPTIONS_QRCID_AR = [...
        bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR; ...
        "ANT3_FAILING"
        ];
      NONEXCEPTIONS_QRCID_AR = setdiff( ...
        bicas.const.qrc.Q.L2_QRCSM.qrcidAr, EXCEPTIONS_QRCID_AR);



      % ASSERTIONS
      assert(isa(L2QflFpa, "bicas.utils.FPArray"))
      assert(isa(ChannelSaturationQrcbm, "bicas.proc.QrcbMap"))
      assert(isequal( ...
        ChannelSaturationQrcbm.qrcidAr, ...
        bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR))



      %--------------------------------------------------------------------
      % Create synthetic L2 QRCBs, which includes as many QRCs as possible
      %--------------------------------------------------------------------
      % NOTE: QRCB values are not set for L2 QRCs (events) which can not be
      % recreated because they are autodetected or read from unavailable CDFs
      % (L1R). If they are specified in the NSO table, then they (QRC events)
      % are included though.
      AllSynthL2Qrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        bicas.const.qrc.Q.L2_QRCSM.qrcidAr, NsoTable, tt2000Ar, L);
      AllSynthL2Qrcbm.union(ChannelSaturationQrcbm)



      %-------------------------------------------------------------------------
      % Determine when L2 QRCs apply for which L2 QFL should (ideally) never be
      % lowered
      %-------------------------------------------------------------------------
      bQrcExceptionAr = false(AllSynthL2Qrcbm.nRecords, 1);
      for qrcid = EXCEPTIONS_QRCID_AR'
        bQrcExceptionAr = bQrcExceptionAr | AllSynthL2Qrcbm.get(qrcid);
      end



      % Local utility function to remove repetition.
      %
      % NOTE: SynthQrcbm2 is only returned for debugging purposes.
      function [SynthL2Qrcbm2, SynthL2QflFpa2] = get_synthetic_QFL(qrcidRemoveAr)
        SynthQrcsm2 = copy(bicas.const.qrc.Q.L2_QRCSM);
        SynthQrcsm2.remove_many(qrcidRemoveAr)

        SynthL2Qrcbm2 = copy(AllSynthL2Qrcbm);
        SynthL2Qrcbm2.remove_many(qrcidRemoveAr)

        [SynthL2Qfl2, ~] = bicas.proc.qrc.QRCB_arrays_to_quality_ZVs(...
          SynthL2Qrcbm2, SynthQrcsm2, "L2_QUALITY_BITMASK");
        SynthL2QflFpa2   = bicas.utils.FPArray(...
          SynthL2Qfl2, 'FILL_POSITIONS', L2QflFpa.fpAr);
      end



      %--------------------------------------
      % Create 2x "always synthetic" L2 QFLs
      %--------------------------------------
      % NOTE: These have synthetic values for all timestamps.

      % Version which only includes the QRC exceptions.
      [ExceptionsSynthL2Qrcbm, ExceptionsSynthL2QflFpa] = ...
        get_synthetic_QFL(NONEXCEPTIONS_QRCID_AR);
      % --
      % Version which includes as many QRCs as possible, except for the QRC
      % exceptions.
      [NonexceptionsSynthL2Qrcbm, NonexceptionsSynthL2QflFpa] = ...
        get_synthetic_QFL(EXCEPTIONS_QRCID_AR);
      %bicas.debug.plot_QRCBM(tt2000Ar, ExceptionsSynthL2Qrcbm,                "ExceptionsSynthL2Qrcbm");
      %bicas.debug.plot_FPA(  tt2000Ar, ExceptionsSynthL2QflFpa,    uint8(-1), "ExceptionsSynthL2QflFpa");
      %bicas.debug.plot_QRCBM(tt2000Ar, NonexceptionsSynthL2Qrcbm,             "NonexceptionsSynthL2Qrcbm");
      %bicas.debug.plot_FPA(  tt2000Ar, NonexceptionsSynthL2QflFpa, uint8(-1), "NonexceptionsSynthL2QflFpa");



      %======================================================================
      % Determine when to use synthetic L2 QFL values instead of true L2 QFL
      % values
      % --------------------------------------------------------------------
      % Tries to minimize the number of synthetic L2 QFL values used.
      %======================================================================

      % Determine when the L2 QFL is lower than L2 QFL should be if only
      % considering exception QRCs. When this is true, the L2 QFL must have been
      % lowered due to some reason other than the exception QRCs. ==> Must never
      % use a synthetic QFL value.
      % IMPLEMENTATION NOTE: This is a trick to further minimize the number of
      % timestamps which use synthetic L2 QFL values (in the return value).
      % Ex: ANT3_FAILING set L2 QFL=1, but L2 QFL=0 for unknown reason. ==> Use
      %     the (true) L2 QFL=0 value.
      BExceptionsSynthQflHigherFpa = ExceptionsSynthL2QflFpa > L2QflFpa;
      bExceptionsSynthQflHigherAr  = BExceptionsSynthQflHigherFpa.array(false);

      bUseSynthQfl = ~bExceptionsSynthQflHigherAr & bQrcExceptionAr;

      %bicas.debug.plot(    tt2000Ar, bQrcExceptionAr,     "bQrcQflExceptionAr");
      %bicas.debug.plot_FPA(tt2000Ar, L2QflFpa, uint8(-1), "L2QflFpa");

      %========================
      % Construct return value
      %========================
      SyntheticL2QflFpa               = L2QflFpa;
      SyntheticL2QflFpa(bUseSynthQfl) = NonexceptionsSynthL2QflFpa(bUseSynthQfl);

      %bicas.debug.plot_FPA(tt2000Ar, SyntheticL2QflFpa, uint8(-1), "SyntheticL2QflFpa");
    end



    % Selectively blank one 2D FPA based on QRCBs. Intended for blanking e.g. VDC
    % and EDC before they are passed to solo.vdccal() and solo.psp2ne(), or
    % their return values.
    %
    %
    % ARGUMENTS
    % =========
    % Fpa
    % Qrcbm
    % Qrcsm
    %       Must contain the same keys as Qrcbm.
    % qrcsFieldName
    %       String. Must refer to QRCS field containing column array of channel
    %       indices (corresponding to the second dimension in Fpa).
    %
    function Fpa = set_FPA_samples_FP(Fpa, Qrcbm, Qrcsm, qrcsFieldName)
      % PROPOSAL: Somehow support QRCS fields which are scalar logical.
      %   PRO: Natural for describing removal of density, scpot.
      %   CON-PROPOSAL: Can describe scalar values as having one channel. -- IMPLEMENTED
      %     CON: Ugly.

      assert(isa(Qrcbm, "bicas.proc.QrcbMap"))
      assert(isa(Qrcsm, "bicas.proc.QrcSettingsMap"))
      assert(isequal(Qrcbm.qrcidAr, Qrcsm.qrcidAr))
      assert(isa(Fpa,   "bicas.utils.FPArray"))
      assert(isstring(qrcsFieldName))
      nRecords  = Qrcbm.nRecords;
      nChannels = irf.assert.sizes(Fpa, [nRecords, -1]);

      bFv            = false(nRecords, nChannels);
      channelIndexAr = repmat(1:nChannels, [nRecords, 1]);

      for qrcid = Qrcbm.qrcidAr'
        qrcbAr = Qrcbm.get(qrcid);    % (nRecords, 1)
        Qrcs   = Qrcsm.get(qrcid);
        assert(isa(Qrcs, "bicas.proc.QrcSettingL3"))

        fvChannelIndexAr = Qrcs.(qrcsFieldName);

        % ASSERTIONS
        assert(iscolumn(fvChannelIndexAr))

        % Arrays of the same size as Fpa.
        bQrcbAr         = repmat(qrcbAr, [1, nChannels]);
        bChannelMatchAr = ismember(channelIndexAr, fvChannelIndexAr);

        bFv = bFv | (bQrcbAr & bChannelMatchAr);
      end

      Fpa(bFv) = bicas.utils.FPArray.get_scalar_FP(Fpa.mc);
    end



  end    % methods(Static)



end
