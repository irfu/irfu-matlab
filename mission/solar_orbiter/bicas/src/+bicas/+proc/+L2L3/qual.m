%
% Collection of code relating to quality ZVs for L2 to L3 processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qual
  % PROPOSAL: Automatic test code.



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % PLANNED TO BE PHASED OUT
    %
    % NOTE: Sets quality ZVs based on BAD_DENSITY QRCB *ONLY*.
    % NOTE: ONLY FOR DENSITY.
    function [QUALITY_FLAG, L3_QUALITY_BITMASK] = ...
        get_quality_ZVs_density(badDensityQrcbAr)
      % NOTE: Can not get NSO events from NSO table.

      % IMPLEMENTATION NOTE: Function is (as of 2023-12-18) in principle more
      % complicated than necessary w.r.t. L3_QUALITY_BITMASK but the
      % architecture is chosen to be analogue with
      % bicas.proc.L1L2.qual.get_quality_ZVs() so that it can be expanded in a
      % similar way if needed.

      assert(islogical(badDensityQrcbAr))

      nRecords = size(badDensityQrcbAr, 1);

      QrcbMap = bicas.proc.QrcbMap(nRecords);
      QrcbMap.add("BAD_DENSITY", badDensityQrcbAr);

      [QUALITY_FLAG, L3_QUALITY_BITMASK] = ...
        bicas.proc.qual.QRCB_arrays_to_quality_ZVs(...
        QrcbMap, bicas.const.Q.L3DENSITY_QRCSM, "L3_QUALITY_BITMASK");
    end



    % UNUSED. EXPERIMENTAL
    %
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
    function QrcbMap = L2QBM_to_QRCBs(l2qbmAr, L2QualityBitSettingQrcsm)
      % IMPLEMENTATION NOTE: Function is designed for
      % SOLO_L2_RPW-LFR-SURV-CWF-E's L2_QUALITY_BITMASK but is generalized to
      % arbitrary QRCSs, partly to simplify testing and not mix algorithms with
      % global constants.

      assert(isa(L2QualityBitSettingQrcsm, "bicas.proc.QrcSettingsMap"))
      assert(iscolumn(l2qbmAr) & isa(l2qbmAr, "uint16"))

      QrcbMap     = bicas.proc.QrcbMap(numel(l2qbmAr));
      % Collection of quality bit positions for all used QRCDSs. Only used for
      % asserting against collisions.
      allBitPosAr = zeros(0, 1);

      for qrcid = L2QualityBitSettingQrcsm.get_QRCIDs()'
        Qrcs = L2QualityBitSettingQrcsm.get(qrcid);
        assert(isa(Qrcs, "bicas.proc.QrcSettingL2"))

        qrcLxqbm = Qrcs.("L2_QUALITY_BITMASK");
        bitPosAr = bicas.proc.qual.LxQBM_to_bit_positions(qrcLxqbm);
        assert(isscalar(bitPosAr), "QRC does not set exactly one bit.")

        allBitPosAr(end+1, 1) = bitPosAr;

        qrcbAr = logical(bitand(l2qbmAr, qrcLxqbm));
        QrcbMap.add(qrcid, qrcbAr)
      end

      % ASSERTION: No overlap in (specified) QRCS quality bits.
      % NOTE: Can not verify against quality bits in QRCS not submitted to this
      %       function.
      irf.assert.number_set(allBitPosAr)
    end



    % EXPERIMENTAL
    %
    % NOTE: Always sets the GLOBAL_SATURATION and CHANNEL_SATURATION QRCIDs but
    % sets the QRCB arrays differently depending on "saturationQualitySchemeId".
    %
    function SaturationQrcbMap = L2QBM_to_saturation_QRCBs(...
        l2qbmAr, saturationQualitySchemeId)

      assert(iscolumn(l2qbmAr) & isa(l2qbmAr, "uint16"))

      SaturationQrcbMap = bicas.proc.QrcbMap(numel(l2qbmAr));
      SaturationQrcbMap.add_false(bicas.const.Q.SATURATION_QRCID_AR)

      switch saturationQualitySchemeId
        case 'CHANNEL_SATURATION'

          ChannelSaturationQrcbMap = bicas.proc.L2L3.qual.L2QBM_to_QRCBs(...
            l2qbmAr, bicas.const.Q.L2_CHANNEL_SATURATION_QRCSM);
          SaturationQrcbMap.add_map(ChannelSaturationQrcbMap)

        case 'GLOBAL_SATURATION'

          ;   % Do nothing. All saturation QRCBs are false.
          % IMPLEMENTATION NOTE: (1) Can not use the global saturation for
          % anything sensible (can not blank selected channels), and (2) the
          % GLOBAL_SATURATION functionality is supposed to be phased out
          % anyway.

        otherwise

          error("BICAS:ConfigurationBug", ...
            "Illegal argument saturationQualitySchemeId=""%s"".", ...
            saturationQualitySchemeId)
      end

      assert(isequal( ...
        sort(SaturationQrcbMap.qrcidAr), ...
        sort(bicas.const.Q.SATURATION_QRCID_AR)))
    end



    % UNUSED SO FAR.
    %
    % For blanking VDC and EDC before they are passed to EXCD
    %
    function [VDC_Fpa, EDC_Fpa] = set_VDC_EDC_samples_FV(...
        VDC_Fpa, EDC_Fpa, QrcbMap, Qrcsm)

      assert(isa(QrcbMap, "bicas.proc.QrcbMap"))
      assert(isa(Qrcsm,   "bicas.proc.QrcSettingsMap"))
      assert(isa(VDC_Fpa, "bicas.utils.FPArray"))
      assert(isa(EDC_Fpa, "bicas.utils.FPArray"))
      nRecords = QrcbMap.nRecords;
      irf.assert.sizes(...
        VDC_Fpa, [nRecords, 3], ...
        EDC_Fpa, [nRecords, 3])

      bVdcFv         = false(nRecords, 3);
      bEdcFv         = false(nRecords, 3);
      channelIndexAr = repmat([1:3], [nRecords, 1]);

      for qrcid = QrcbMap.qrcidAr'
        qrcbAr = QrcbMap.get(qrcid);    % (nRecords, 1)
        Qrcs   = Qrcsm.get(  qrcid);
        assert(isa(Qrcs, "bicas.proc.QrcSettingL3"))

        % Arrays of the same size as VDC_Fpa & EDC_Fpa.
        bQrcbAr          = repmat(qrcbAr, [1, 3]);
        bVdcChannelMatch = ismember(channelIndexAr, Qrcs.vdcFvIndexAr);
        bEdcChannelMatch = ismember(channelIndexAr, Qrcs.edcFvIndexAr);

        bVdcFv           = bVdcFv | (bQrcbAr & bVdcChannelMatch);
        bEdcFv           = bEdcFv | (bQrcbAr & bEdcChannelMatch);
      end

      VDC_Fpa(bVdcFv) = bicas.utils.FPArray.FP_SINGLE;
      EDC_Fpa(bEdcFv) = bicas.utils.FPArray.FP_SINGLE;
    end



  end    % methods(Static)



end
