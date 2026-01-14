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



    % Derive a "synthetic" L2 LFR CWF QUALITY_FLAG (the input to L3) equivalent
    % to the L2 QUALITY_FLAG
    % (1) in the absence of
    %         QRC SATURATION_ZV_V*
    %     (desired effect)
    % (2) in the absence of QRCs which are not specified in NSO table
    %         BIAS_HW_OFF  (from LFR)
    %         SWEEP        (autodetected)
    %     (undesired effect, but does not matter since data is blanked),
    % (3) if L1R QUALITY_FLAG had not been lowered below max (undesired effect).
    %
    % NOTE: In practice, this should be the L2 QUALITY_FLAG derived from the
    % NSO table (minus saturation), i.e. also without autodetected sweeps.
    %
    % NOTE/BUG: This function can set a maximum QUALITY_FLAG value which is
    % higher than what L1R could possibly have (not the QUALITY_FLAG max=4, but
    % the highest QUALITY_FLAG which ROC wants to produce). QUALITY_FLAG is
    % eventually capped by setting "PROCESSING.ZV_QUALITY_FLAG_MAX" when writing
    % the CDF anyway, but this setting then has to use the appropriate value.
    %
    % BUG: This function does not consider the L2 QUALITY_FLAG when there is no
    % saturation. Therefore, when L2 QUALITY_FLAG is lowered only due to L1R
    % QUALITY_FLAG being lowered (in the absence of saturation or any other
    % BICAS mechanism), this function will still generate a higher QUALITY_FLAG
    % value.
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
        tt2000Ar, NsoTable, qflFpAr, L)

      assert(islogical(qflFpAr))

      % NOTE: One could consider also removing/excluding ANT3_FAILING, since
      % the unaffected channels should be OK. Has no instructions to do so yet
      % though. /2025-08-27
      SyntheticL2Qrcsm = copy(bicas.const.qrc.Q.L2_QRCSM);
      SyntheticL2Qrcsm.remove_many(bicas.const.qrc.Q.SATURATION_QRCID_AR);

      SyntheticL2Qrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        SyntheticL2Qrcsm.qrcidAr, NsoTable, tt2000Ar, L);

      [SyntheticL2Qfl, ~] = bicas.proc.qrc.QRCB_arrays_to_quality_ZVs(...
        SyntheticL2Qrcbm, SyntheticL2Qrcsm, "L2_QUALITY_BITMASK");

      SyntheticL2QflFpa = bicas.utils.FPArray(...
        SyntheticL2Qfl, 'FILL_POSITIONS', qflFpAr);
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
