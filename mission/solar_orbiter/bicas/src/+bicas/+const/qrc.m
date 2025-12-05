%
% Class for defining QRC-related constants, in particular QRCSs (stored inside
% QRCSMs).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef qrc
  % PROPOSAL: Replace bicas.const.qrc.Q with separate constants.
  %   PRO: Shorter variables paths.
  %   CON: Can not shorten/copy in code: Q = bicas.const.qrc.Q.
  %     CON: Has not been needed so far.
  % PROPOSAL: Constants for QFL values (0..4).



  % QUALITY_FLAG values according to
  % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
  % ===========================================================================
  % 0	Bad data
  % 1	Known problems, use at your own risk
  % 2	Survey data, possibly not publication-quality
  % 3	Good for publication, subject to PI approval
  % 4	Excellent data which has received special treatment
  % """"""""
  % The QUALITY_FLAG values of the RPW L1 CDF are all set to "3" by default.
  %
  % From this initial status, the RCS must then decide in an autonomous way to
  % change or not the value of QUALITY_FLAG records, when transferring the
  % quality information to the output RPW CDF. The RCS decision must rely on
  % relevant data (i.e, QUALITY_BITMASK, other data, parent data, etc.) as
  % explained in the next section.
  % """"""""



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant)



    % Constants for common values
    % ---------------------------
    % NOTE: Useful for e.g. testing.
    LxQBM_NONE = uint16(0);

    % Absolute min & max for ZV QUALITY_FLAG, according to the definition in
    % external metadata standards.
    QFL_MIN = uint8(0);
    QFL_MAX = uint8(4);



    Q = bicas.const.qrc.init_QRC_constants();



  end
  properties(Constant, Access=private)



    % QUALITY_FLAG value for both
    % (1) global saturation's "full saturation" QRC, and
    % (2) channel saturation QRCs.
    QFL_SATURATION = uint8(0);



  end



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Function for the initializing a QRCSM with all QRCSs for CHANNEL
    % saturation in L2 CDFs and L3 CDFs.
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
          qfl   = bicas.const.qrc.QFL_SATURATION, ...
          l2qbm = uint16(l2qbm));

        L2Qrcsm.add(qrcid, L2Qrcs);
      end

      % NOTE: Setting QFL and L2QBM. Not blanking anything.
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V1",   1)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V2",   2)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V3",   4)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V12",  8)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V13", 16)
      add_L2_channel_saturation_QRCS("SATURATION_ZV_V23", 32)

      % NOTE: *NOT* setting QFL and L3QBM. Blanking L2 input data before using
      % it for deriving L3* (i.e. not alter the actual content of L2 CDFs).
      L3Qrcsm.add("SATURATION_ZV_V1",  bicas.proc.QrcSettingL3(vdcFvIndexAr=1))
      L3Qrcsm.add("SATURATION_ZV_V2",  bicas.proc.QrcSettingL3(vdcFvIndexAr=2))
      L3Qrcsm.add("SATURATION_ZV_V3",  bicas.proc.QrcSettingL3(vdcFvIndexAr=3))
      L3Qrcsm.add("SATURATION_ZV_V12", bicas.proc.QrcSettingL3(edcFvIndexAr=1))
      L3Qrcsm.add("SATURATION_ZV_V13", bicas.proc.QrcSettingL3(edcFvIndexAr=2))
      L3Qrcsm.add("SATURATION_ZV_V23", bicas.proc.QrcSettingL3(edcFvIndexAr=3))

      assert(isequal(L2Qrcsm.qrcidAr, L3Qrcsm.qrcidAr))
    end



    % Function for the initializing QRCSM with all QRCSs for GLOBAL saturation
    % in L2 CDFs.
    %
    function [L2Qrcsm] = init_L2_GLOBAL_SATURATION_QRCSM()
      L2Qrcsm = bicas.proc.QrcSettingsMap();

      % NOTE: L2QBM_BIT_PARTIAL_SATURATION is used for two different QRCIDs.
      L2QBM_BIT_PARTIAL_SATURATION = uint16(1);
      L2QBM_BIT_FULL_SATURATION    = uint16(2);

      %====================
      % PARTIAL_SATURATION
      %====================
      Qrcs = bicas.proc.QrcSettingL2(...
        qfl   = uint8(1), ...
        l2qbm = L2QBM_BIT_PARTIAL_SATURATION);
      L2Qrcsm.add("PARTIAL_SATURATION", Qrcs)

      %=================
      % FULL_SATURATION
      %=================
      % NOTE: Also set PARTIAL saturation bit when FULL
      % saturation. /YK 2020-10-02.
      Qrcs = bicas.proc.QrcSettingL2(...
        qfl   = bicas.const.qrc.QFL_SATURATION, ...
        l2qbm = L2QBM_BIT_FULL_SATURATION + L2QBM_BIT_PARTIAL_SATURATION);
      L2Qrcsm.add("FULL_SATURATION", Qrcs);
    end



    % Function for initializing QRCSMs in three groups representing different
    % output CDFs/DSIDs. The three groups are:
    % (1) all (official) L2 CDFs,
    % (2) all L3 CDFs, and
    % (3) all L3 density CDFs.
    % NOTE: Groups (2) and (3) overlap.
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
      % NOTE: Empirically, ROC's thruster firing quality bit covers much longer
      % time (several hours) than the actual thruster firing according to the
      % actual effect on data. /Erik P G Johansson 2025-08-29
      Qrcs = bicas.proc.QrcSettingL2(qfl = uint8(1));
      L2Qrcsm.add("THRUSTER_FIRING", Qrcs);

      %=============
      % BIAS_HW_OFF
      %=============
      % This condition corresponds when LFR ZV "BW" states that BIAS h/w is off.
      % --
      % IMPLEMENTATION NOTE: It could be argued that this condition is not
      % really an "error", or "quality problem" at all and that the
      % corresponding functionality should not be implemented using QRCs. It is
      % an edge case, conceptually. However, it has been implemented as a QRC
      % anyway to make the hardcoding/documentation/configuration of BICAS's
      % behaviour when this condition applies clearer. It is "obvious" (?) that
      % the voltages should be blanked, but not the currents, and that should
      % (preferably) be configured clearly somewhere and doing it using QRCSs
      % makes sense. Rather, this use case is an argument for that QRC* should
      % be renamed.
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
      % in case other antennas fail, one can just create the corresponding QRCSs
      % for ANT1 & ANT2 and fill the NSO table with the corresponding entries.
      %
      % NOTE: This QRC only affects some channels.
      %
      % NOTE: As of 2025-09-02, ROC modifies CDFs created by BICAS (at ROC;
      % i.e. L2) to cap QUALITY_FLAG<=1="Known problems, use at your own risk"
      % for ANT3 failing as defined by ROC. Erik P G Johansson (2025-09-02) is
      % not aware of any formal document specifying this value (though maybe it
      % is obvious considering the QUALITY_FLAG value definitions). It was
      % mentioned in an e-mail thread (below) and is consistent with
      % LIRA-generated CDFs.
      % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
      %
      % """"""""
      % The RPW data pipeline at LESIA includes the capability to update a
      % given set of CDF files with specific information for given data
      % products and time ranges. This function will be for instance used to
      % notify RPW data users about degraded science data impact during the
      % periods when the ANT3 failed.
      %
      % In the ANT3 case, the idea to is to set the value of the two following
      % items:
      %
      % - CAVEATS global attribute —> Exact message has to be defined, but
      %   should be something like "RPW electrical antenna 3 [MY] failure
      %   reported during the current day (see where QUALITY_FLAG=1)."
      % - QUALITY_FLAG zVariable —> Will be set to 1 when ANT3 failure is
      %   reported
      %
      % These changes will be performed for the L1 CDF by the ROC, but same
      % information shall be also found in the child L2 (and L3 when needed).
      % """"""""
      % /Xavier Bonnin e-mail 2024-07-12, 11:24
      %
      % NOTE: It is not conceptually obvious whether QUALITY_FLAG should be
      % capped if data affected by ANT3_FAILING is entirely removed. However,
      % since ROC caps QUALITY_FLAG<=1 for ANT3_FAILING in L1/L1R (see above),
      % in principle, capping QFL<=1 should have no influence here, unless
      % artificially raising it (is allowed for L2-->L3).

      % NOTE: Using the *EXACT SAME* string which ROC uses for setting GA
      % CAVEATS in CDFs output by RCSs (when agreed upon). /2025-09-02
      ANT3_FAILING_GA_CAVEATS = "RPW electrical antenna 3 [MY] failure" + ...
        " reported during the current day (see QUALITY_FLAG=1).";

      % IMPLEMENTATION NOTE: Creating QRCSs for both 6xL2 and 6xL3 datasets.
      % This is needed since GA CAVEATS values are *NOT* inherited (or modified)
      % from parent CDFs as opposed to for QFL and LxQBM.
      % --
      % Blanking data in accordance with
      % https://github.com/irfu/irfu-matlab/issues/142 .
      Qrcs = bicas.proc.QrcSettingL2(...
        qfl             = uint8(1), ...
        voltageFvSsidAr = S(["DC_V3" "DC_V13" "DC_V23" "AC_V13" "AC_V23"]'), ...
        gaCaveats       = ANT3_FAILING_GA_CAVEATS);
      L2Qrcsm.add("ANT3_FAILING", Qrcs);
      %
      Qrcs = bicas.proc.QrcSettingL3(...
        gaCaveats       = ANT3_FAILING_GA_CAVEATS);
      L3Qrcsm.add("ANT3_FAILING", Qrcs);



      %=============
      % BAD_DENSITY
      %=============
      % NOTE: L3_QUALITY_BITMASK currently only exists for L3 DENSITY datasets
      % (i.e. not for other L3 datasets). If L3_QUALITY_BITMASK is used for
      % other L3 datasets in the future, then the bits may have different
      % meanings for those datasets.
      Qrcs = bicas.proc.QrcSettingL3Density(...
        qfl   = uint8(1), ...
        l3qbm = uint16(1));
      L3DensityQrcsm.add("BAD_DENSITY", Qrcs);



      %===============================
      % ANT3_UNINTENTIONALLY_FLOATING
      %===============================
      % ANT3 is unintentionally floating after sweeps due to bad commanding.
      % https://github.com/irfu/irfu-matlab/issues/156
      % NOTE: Separate QRCSs for L1/L1R-->L2 and L2-->L3 processing.
      % NOTE: 2025-12-03: ROC's SOLO_L1_RPW-BIA-CURRENT datasets contain
      %       non-zero ANT3 bias currents for ANT3_UNINTENTIONALLY_FLOATING,
      %       which is then copied to the BIAS L2 datasets. Therefore, the L2
      %       bias currents are also non-zero.
      % YK 2025-12-03: Keep samples, but blank bias current.
      ANT3_UNINTENTIONALLY_FLOATING_CAVEATS = ...
        "BIAS unintentionally sets zero bias current on ANT3 (see QUALITY_FLAG=1) during the current day.";
      Qrcs = bicas.proc.QrcSettingL2(...
        qfl             = uint8(1), ...
        gaCaveats       = ANT3_UNINTENTIONALLY_FLOATING_CAVEATS, ...
        currentFvIantAr = [3]);
      L2Qrcsm.add("ANT3_UNINTENTIONALLY_FLOATING", Qrcs);
      % --
      % NOTE: Removes EFIELD output from solo.vdccal() (i.e. not by removing
      %       *INPUT* to solo.vdccal()).
      % NOTE: Removing EFIELD should be unnecessary due to (synthetic)
      %       L2 QUALITY_FLAG < PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN.
      % NOTE: Unclear what should be the correct QFL value when data is removed.
      %       Therefore not setting it.
      %   PROPOSAL: QFL=0
      Qrcs = bicas.proc.QrcSettingL3(...
        efieldFvIndexAr = [1 2 3]', ...
        gaCaveats       = ANT3_UNINTENTIONALLY_FLOATING_CAVEATS);
      L3Qrcsm.add("ANT3_UNINTENTIONALLY_FLOATING", Qrcs);



      %=============================
      % REMOVE_EFIELD/DENSITY/SCPOT
      %=============================
      % Intended for when there are arbitrary reasons to remove L3 data.
      %   Ex: Ensure that L3 data is never removed outside the time interval for
      %   which there is data. (Should ideally be handled by solo.vdccal(),
      %   solo.psp2ne() themselves, but might not.)
      %   Ex: Some part of calibration has been found to be wrong. Has happened
      %       for EFIELD and this was the quickfix.
      % NOTE: Removes output from solo.vdccal()/solo.psp2ne() (i.e. not
      % by removing input to those functions).
      % NOTE: Does not cap QFL (yet).
      Qrcs = bicas.proc.QrcSettingL3(efieldFvIndexAr  = [1 2 3]');
      L3Qrcsm.add("REMOVE_EFIELD", Qrcs);
      %
      Qrcs = bicas.proc.QrcSettingL3(densityFvIndexAr = [1]);
      L3Qrcsm.add("REMOVE_DENSITY", Qrcs);
      %
      Qrcs = bicas.proc.QrcSettingL3(scpotFvIndexAr   = [1]');
      L3Qrcsm.add("REMOVE_SCPOT", Qrcs);



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
