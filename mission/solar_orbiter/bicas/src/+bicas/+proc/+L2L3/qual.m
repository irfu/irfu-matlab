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
    % function QrcbMap = get_saturation_QRCBs_from_L2QBM(l2qbmAr, saturationQualitySchemeId)
    %   assert(iscolumn(l2qbmAr) & isa(l2qbmAr, "uint16"))
    %
    %   switch saturationQualitySchemeId
    %     case 'CHANNEL_SATURATION'
    %
    %       QrcbMap = bicas.proc.qual.LxQBM_to_QRCB_maps(...
    %         l2qbmAr, 'SOLO_L2_RPW-LFR-SURV-CWF-E', "L2_QUALITY_BITMASK", ...
    %         bicas.const.Q.CHANNEL_SATURATION_QRCS_MAP);
    %
    %     case 'GLOBAL_SATURATION'
    %
    %       QrcbMap = bicas.proc.QrcbMap(numel(l2qbmAr));
    %       for qrcid = bicas.const.Q.CHANNEL_SATURATION_QRCID_AR
    %         QrcbMap.add(qrcid, false(size(l2qbmAr)))
    %       end
    %
    %     otherwise
    %
    %       error("BICAS:ConfigurationBug", ...
    %         "Illegal argument saturationQualitySchemeId=""%s"".", ...
    %         saturationQualitySchemeId)
    %   end
    % end



    % UNUSED SO FAR. UNTESTED.
    % function [VDC_Fpa, EDC_Fpa] = set_VDC_EDC_samples_FV(VDC_Fpa, EDC_Fpa, QrcbMap, QrcsMap)
    %   assert(isa(QrcbMap, "bicas.proc.QrcbMap"))
    %   assert(isa(Qrcsmap, "containers.Map"))
    %   assert(isa(VDC_Fpa, "bicas.utils.FPArray"))
    %   assert(isa(EDC_Fpa, "bicas.utils.FPArray"))
    %   nRecords = QrcbMap.nRecords;
    %   irf.assert.size(...
    %     VDC_Fpa, [nRecords, 3], ...
    %     EDC_Fpa, [nRecords, 3])
    %
    %
    %
    %   bVdcFv         = false(nRecords, 3);
    %   bEdcFv         = false(nRecords, 3);
    %   channelIndexAr = repmat([1:3], [nRecords, 1]);
    %
    %   for qrcid = QrcbMap.qrcidAr'
    %     qrcbAr     = QrcbMap.get(qrcid);    % (nRecords, 1)
    %     Qrcs       = QrcsMap(    qrcid);
    %     assert(isa(Qrcs, "bicas.proc.QrcSetting"))
    %     Qrcds      = Qrcs.get(dsi);
    %     if isempty(Qrcds)
    %       continue
    %     end
    %
    %     % Arrays of the same size as voltageAr.
    %     bQrcbAr          = repmat(qrcbAr, [1, 3]);
    %     bVdcChannelMatch = ismember(channelIndexAr, Qrcds.vdcIndexAr);
    %     bEdcChannelMatch = ismember(channelIndexAr, Qrcds.edcIndexAr);
    %
    %     bVdcFv           = bVdcFv | (bQrcbAr & bVdcChannelMatch);
    %     bEdcFv           = bEdcFv | (bQrcbAr & bEdcChannelMatch);
    %   end
    %
    %   VDC_Fpa(bVdcFv) = bicas.utils.FPArray.FP_SINGLE;
    %   EDC_Fpa(bEdcFv) = bicas.utils.FPArray.FP_SINGLE;
    % end



  end    % methods(Static)



end
