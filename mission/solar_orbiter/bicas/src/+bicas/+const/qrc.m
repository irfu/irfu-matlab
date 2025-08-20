%
% Class for defining QRC-related constants, in particular QRCSs (stored inside
% QRCSMs).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qrc
  % PROPOSAL: Rename QRCS or QRCSM.
  %   PRO: QRCS and QRCSM are the only publically visible products of this class.
  %   CON: The subject of the class is QRCs.
  %
  % PROPOSAL: Replace bicas.const.qrc.Q with separate constants.
  %   PRO: Shorter variables paths.
  %   CON: Can not shorten/copy in code: Q = bicas.const.qrc.Q.
  %     CON: Has not been needed so far.



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant)



    % Constants for common values
    % ---------------------------
    % NOTE: Useful for e.g. testing.
    LxQBM_NONE       = uint16(0);

    % Absolute min & max for ZV QUALITY_FLAG, according to the definition in
    % external metadata standards.
    QUALITY_FLAG_MIN = uint8(0);
    QUALITY_FLAG_MAX = uint8(4);



    Q = bicas.const.qrc.init_QRC_constants();



  end
  properties(Constant, Access=private)



    % QUALITY_FLAG value for both
    % (1) global saturation's "full saturation" QRC, and
    % (2) channel saturation QRCs.
    QUALITY_FLAG_SATURATION = uint8(0);



  end



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Function for the initializing QRCSM with all QRCSs for channel saturation
    % in L2 CDFs and L3 CDFs.
    %
    function [L2Qrcsm, L3Qrcsm] = init_L2_L3_CHANNEL_SATURATION_QRCSM()
      L2Qrcsm = bicas.proc.QrcSettingsMap();
      L3Qrcsm = bicas.proc.QrcSettingsMap();

      %====================
      % CHANNEL SATURATION
      %====================
      % NOTE: DC/AC diffs use the same QRCIDs and settings.
      % NOTE: Channel saturation QRCIDs technically represent saturation on the
      %       corresponding ZVs, not the specified ASR SSIDs. This is only
      %       important if adding automatic saturation detection for non-ASRs
      %       (which is highly unlikely).
      function add_L2_channel_saturation_QRCS(qrcid, l2qbm)
        L2Qrcs = bicas.proc.QrcSettingL2(...
          QUALITY_FLAG       = bicas.const.qrc.QUALITY_FLAG_SATURATION, ...
          L2_QUALITY_BITMASK = uint16(l2qbm));
        L2Qrcsm.add(qrcid, L2Qrcs);
      end

      add_L2_channel_saturation_QRCS("SATURATION_ZV_V1",   1)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V2",   2)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V3",   4)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V12",  8)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V13", 16)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V23", 32)

      L3Qrcsm.add("SATURATION_ZV_V1",  bicas.proc.QrcSettingL3(vdcFvIndexAr=1))
      L3Qrcsm.add("SATURATION_ZV_V2",  bicas.proc.QrcSettingL3(vdcFvIndexAr=2))
      L3Qrcsm.add("SATURATION_ZV_V3",  bicas.proc.QrcSettingL3(vdcFvIndexAr=3))
      L3Qrcsm.add("SATURATION_ZV_V12", bicas.proc.QrcSettingL3(edcFvIndexAr=1))
      L3Qrcsm.add("SATURATION_ZV_V13", bicas.proc.QrcSettingL3(edcFvIndexAr=2))
      L3Qrcsm.add("SATURATION_ZV_V23", bicas.proc.QrcSettingL3(edcFvIndexAr=3))

      assert(isequal(L2Qrcsm.qrcidAr, L3Qrcsm.qrcidAr))
    end



    function [L2Qrcsm] = init_L2_GLOBAL_SATURATION_QRCSM()
      L2Qrcsm = bicas.proc.QrcSettingsMap();

      % Global saturation quality variable scheme:
      % NOTE: L2QBM_BIT_PARTIAL_SATURATION is used for two different QRCIDs.
      L2QBM_BIT_PARTIAL_SATURATION = uint16(1);
      L2QBM_BIT_FULL_SATURATION    = uint16(2);

      %====================
      % PARTIAL_SATURATION
      %====================
      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG       = uint8(1), ...
        L2_QUALITY_BITMASK = L2QBM_BIT_PARTIAL_SATURATION);
      L2Qrcsm.add("PARTIAL_SATURATION", Qrcs)

      %=================
      % FULL_SATURATION
      %=================
      % NOTE: Also set PARTIAL saturation bit when FULL
      % saturation. /YK 2020-10-02.
      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG       = bicas.const.qrc.QUALITY_FLAG_SATURATION, ...
        L2_QUALITY_BITMASK = L2QBM_BIT_FULL_SATURATION + L2QBM_BIT_PARTIAL_SATURATION);
      L2Qrcsm.add("FULL_SATURATION", Qrcs);
    end



    % Function for initializing QRCSMs containing all QRCSs for producing
    % (1) all (official) L2 CDFs, and
    % (3) all L3 density CDFs.
    %
    % NOTE: The quality bit definitions here must be consistent with the
    % definitions in the corresponding CDF skeletons. Definitions in general
    % must be consistent with documentation.
    %
    function Q = init_L2_L3_L3Density_QRCSM()

      [Q.L2_CHANNEL_SATURATION_QRCSM, ...
        Q.L3_CHANNEL_SATURATION_QRCSM] = ...
        bicas.const.qrc.init_L2_L3_CHANNEL_SATURATION_QRCSM();

      Q.L2_GLOBAL_SATURATION_QRCSM = ...
        bicas.const.qrc.init_L2_GLOBAL_SATURATION_QRCSM();

      L2Qrcsm        = bicas.proc.QrcSettingsMap();
      L3Qrcsm        = bicas.proc.QrcSettingsMap();
      L3DensityQrcsm = bicas.proc.QrcSettingsMap();

      L2Qrcsm.add_QRCSM(Q.L2_CHANNEL_SATURATION_QRCSM);
      L3Qrcsm.add_QRCSM(Q.L3_CHANNEL_SATURATION_QRCSM);
      L2Qrcsm.add_QRCSM(Q.L2_GLOBAL_SATURATION_QRCSM);

      %=================
      % Local constants
      %=================
      % No circular dependence?!! bicas.proc.L1L2.const <-> bicas.const.qrc
      S = bicas.proc.L1L2.const.C.SSID_DICT;



      %=================
      % THRUSTER_FIRING
      %=================
      % NOTE: There will be an L1 QUALITY_BITMASK bit for thruster firings in
      % the future according to
      % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
      % Therefore(?) not setting any bit in L2_QUALITY_BITMASK.
      % (YK 2020-11-03 did not ask for any to be set.)
      Qrcs = bicas.proc.QrcSettingL2(QUALITY_FLAG = uint8(1));
      L2Qrcsm.add("THRUSTER_FIRING", Qrcs);

      %=============
      % BIAS_HW_OFF
      %=============
      % This condition corresponds when LFR ZV "BW" says that BIAS h/w is off.
      % --
      % IMPLEMENTATION NOTE: It could be argued that this condition is not
      % really an "error", or "quality problem" at all and that the
      % corresponding functionality should not be implemented using QRCs. It is
      % an edge case, conceptually. However, it has been implemented as a QRC
      % anyway to make the hardcoding/documentation/configuration of BICAS's
      % behaviour when this condition applies clearer. It is "obvious" (?) that
      % the voltages should be blanked, but not the currents, and that should
      % (preferably) be configured clearly somewhere.
      % --
      % NOTE: Deletes all values, both voltages & currents.
      Qrcs = bicas.proc.QrcSettingL2(...
        voltageFvSsidAr = S.values, ...
        currentFvIantAr = [1:3]');
      L2Qrcsm.add("BIAS_HW_OFF", Qrcs);

      %=======
      % SWEEP
      %=======
      Qrcs = bicas.proc.QrcSettingL2(...
        voltageFvSsidAr = S.values, ...   % Blank all SSIDs.
        currentFvIantAr = [1:3]');
      L2Qrcsm.add("SWEEP", Qrcs);

      %==============
      % ANTx_FAILING
      %==============
      % NOTE: There is currently (2025-08-26) only a need for ANT3_FAILING, but
      % in case other antennas fail, one can just create the corresponding
      % QRCSs for ANT1 & ANT2 and fill the NSO table with the corresponding
      % entries.

      % TEMPORARY VALUES. Good values not decided.
      % NOTE: QRC only affects some channels.
      QUALITY_FLAG_ANTx_FAILING = bicas.const.qrc.QUALITY_FLAG_MAX;
      L2QBM_ANTx_FAILING        = bicas.const.qrc.LxQBM_NONE;

      % Qrcs = bicas.proc.QrcSettingL2(...
      %   QUALITY_FLAG       = QUALITY_FLAG_ANTx_FAILING, ...
      %   voltageFvSsidAr    = S(["DC_V1" "DC_V12" "DC_V13" "AC_V12" "AC_V13"]'));
      % L2Qrcsm.add("ANT1_FAILING", Qrcs);
      %
      % Qrcs = bicas.proc.QrcSettingL2(...
      %   QUALITY_FLAG       = QUALITY_FLAG_ANTx_FAILING, ...
      %   voltageFvSsidAr    = S(["DC_V2" "DC_V12" "DC_V23" "AC_V12" "AC_V23"]'));
      % L2Qrcsm.add("ANT2_FAILING", Qrcs);

      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG       = QUALITY_FLAG_ANTx_FAILING, ...
        voltageFvSsidAr    = S(["DC_V3" "DC_V13" "DC_V23" "AC_V13" "AC_V23"]'));
      L2Qrcsm.add("ANT3_FAILING", Qrcs);



      %=============
      % BAD_DENSITY
      %=============
      % NOTE: L3_QUALITY_BITMASK currently only exists for L3 DENSITY datasets
      % (i.e. not for other L3 datasets). If L3_QUALITY_BITMASK is used for other
      % L3 datasets in the future, then the bits may have different meanings for
      % those datasets.
      Qrcs = bicas.proc.QrcSettingL3Density(...
        QUALITY_FLAG       = uint8(1), ...
        L3_QUALITY_BITMASK = uint16(1));
      L3DensityQrcsm.add("BAD_DENSITY", Qrcs);



      %=============================
      % V3_UNINTENTIONALLY_FLOATING
      %=============================
      % V3 is unintentionally floating after sweeps due to bad commanding.
      % https://github.com/irfu/irfu-matlab/issues/156
      % NOTE: Removes EFIELD output from solo.vdccal() (i.e. not by removing
      % input to solo.vdccal()).
      Qrcs = bicas.proc.QrcSettingL3(efieldFvIndexAr=[1 2 3]');
      L3Qrcsm.add("V3_UNINTENTIONALLY_FLOATING", Qrcs);



      %===============
      % REMOVE_EFIELD
      %===============
      % Intended for when there are arbitrary reasons to remove EFIELD data.
      % NOTE: Removes EFIELD output from solo.vdccal() (i.e. not by removing
      % input to solo.vdccal()).
      Qrcs = bicas.proc.QrcSettingL3(efieldFvIndexAr=[1 2 3]');
      L3Qrcsm.add("REMOVE_EFIELD", Qrcs);



      Q.L2_QRCSM         = L2Qrcsm;
      Q.L3_QRCSM         = L3Qrcsm;
      Q.L3_DENSITY_QRCSM = L3DensityQrcsm;
    end



    function Q = init_QRC_constants()
      Q = bicas.const.qrc.init_L2_L3_L3Density_QRCSM();

      Q.CHANNEL_SATURATION_QRCID_AR = string(Q.L2_CHANNEL_SATURATION_QRCSM.qrcidAr);
      Q.GLOBAL_SATURATION_QRCID_AR  = string(Q.L2_GLOBAL_SATURATION_QRCSM.qrcidAr);
      Q.SATURATION_QRCID_AR = [...
        Q.CHANNEL_SATURATION_QRCID_AR; ...
        Q.GLOBAL_SATURATION_QRCID_AR];

      % All legal QRCIDs, or all kinds of processing. This defines the set of
      % legal QRCIDs, including ones that can be used in the NSO table file.
      % IMPLEMENTATION NOTE: This is required for asserting QRCIDs in the NSO
      % table file.
      Q.ALL_QRCID_AR = [...
        Q.L2_QRCSM.qrcidAr; ...
        Q.L3_QRCSM.qrcidAr; ...
        Q.L3_DENSITY_QRCSM.qrcidAr];
    end



  end    % methods(Static, Access=private)



end
