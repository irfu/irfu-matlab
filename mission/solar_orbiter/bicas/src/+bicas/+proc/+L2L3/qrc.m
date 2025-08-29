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

      Qrcbm     = bicas.proc.QrcbMap(numel(l2qbmAr));
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
        case 'CHANNEL_SATURATION'

          % Update CHANNEL_SATURATION QRCBs.
          ChannelSaturationQrcbm = bicas.proc.L2L3.qrc.L2QBM_to_QRCBs(...
            l2qbmAr, bicas.const.qrc.Q.L2_CHANNEL_SATURATION_QRCSM);
          SaturationQrcbm.union(ChannelSaturationQrcbm)

        case 'GLOBAL_SATURATION'

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



    % Better "hack" for obtaining the QUALITY_FLAG for L2 LFR CWF (the input to
    % L3) in the absence of saturation.
    %
    % NOTE: In practice, this should be the L2 QUALITY_FLAG derived from the
    % NSO table (minus saturation), i.e. also without autodetected sweeps.
    %
    % RATIONALE
    % =========
    % This is needed for deriving the L3 QUALITY_FLAG which may be higher than
    % the corresponding L2 QUALITY FLAG when L3 contains valid values derived
    % from non-saturated L2 channels in the presence of other saturated L2
    % channels.
    %
    function SyntheticL2QflFpa = get_synthetic_L2_QFL( ...
        tt2000Ar, NsoTable, qflFpAr, L)
      % PROPOSAL: Separate function for deriving QUALITY_FLAG.
      %   PRO: Also "needed" for EFIELD+SCPOT which do not use QUALITY_BITMASK.
      %   CON: Can ignore return value.
      %     CON: bicas.proc.qrc.QRCB_arrays_to_quality_ZVs() still requires
      %          lxqbmName and QRCSs which contain some LxQBM value.
      %       CON-PROPOSAL: Special value to ignore retrieving a QRCS LxQBM value.

      assert(islogical(qflFpAr))

      % NOTE: One could consider also removing/excluding ANT3_FAILING, since
      % the unaffected channels should be OK. Has no instructions to do so yet
      % though. /2025-08-27
      NonsaturationL2Qrcsm = copy(bicas.const.qrc.Q.L2_QRCSM);
      NonsaturationL2Qrcsm.remove_many(bicas.const.qrc.Q.SATURATION_QRCID_AR);

      L2Qrcbm = bicas.proc.qrc.NSO_table_to_QRCBM(...
        NonsaturationL2Qrcsm.qrcidAr, NsoTable, tt2000Ar, L);

      [NonsaturationL2Qfl, ~] = bicas.proc.qrc.QRCB_arrays_to_quality_ZVs(...
        L2Qrcbm, NonsaturationL2Qrcsm, "L2_QUALITY_BITMASK");

      SyntheticL2QflFpa = bicas.utils.FPArray(...
        NonsaturationL2Qfl, 'FILL_POSITIONS', qflFpAr);
    end



    % Blank 2D FPA based on QRCBs. Intended for blanking e.g. VDC and EDC
    % before they are passed to solo.vdccal() and solo.psp2ne(), or their
    % return values.
    %
    %
    % ARGUMENTS
    % =========
    % Qrcbm
    % Qrcsm
    %       Must contain the same keys as Qrcbm.
    % qrcsFieldName
    %       String. Must refer to QRCS field containing column array of channel
    %       indices (second dimension in Fpa).
    %
    function Fpa = set_FPA_samples_FP(Fpa, Qrcbm, Qrcsm, qrcsFieldName)
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
