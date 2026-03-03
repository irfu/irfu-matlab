%
% Collection of code relating to QRCs for L1/L1R to L2 processing.
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



    % NOTE: Always returns both the GLOBAL_SATURATION and CHANNEL_SATURATION
    % QRCIDs but sets the QRCB arrays differently depending on
    % "saturationQualitySchemeId".
    %
    % ARGUMENTS
    % =========
    % VsibZvm
    %       ZVM with VSIBs (not VSTBs, not VSQBs, not QRCBs) for the respective
    %       ASR channels.
    %
    function SaturationQrcbm = VSIBs_to_saturation_QRCBs(...
        tt2000Ar, ...
        VsibZvm, isSwf, vstbFractionThreshold, cwfSlidingWindowLengthSec)

      % PROPOSAL: Change order of arguments.

      SaturationQrcbm = bicas.proc.QrcbMap(numel(tt2000Ar));

      %---------------------------------
      % Obtain CHANNEL_SATURATION QRCBs
      %---------------------------------
      ChannelSaturationQrcbm = ...
        bicas.proc.L1L2.qrc.VSIBs_to_channel_saturation_QRCBs(...
        VsibZvm, tt2000Ar, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec);
      SaturationQrcbm.union(ChannelSaturationQrcbm)

      %--------------------------------
      % Obtain GLOBAL_SATURATION QRCBs
      %--------------------------------
      GlobalSaturationQrcbm = ...
        bicas.proc.L1L2.qrc.channel_saturation_to_global_saturation_QRCBs(...
        ChannelSaturationQrcbm);
      SaturationQrcbm.union(GlobalSaturationQrcbm)

      assert(isequal( ...
        sort(SaturationQrcbm.qrcidAr), ...
        sort(bicas.const.qrc.Q.SATURATION_QRCID_AR)))
    end



    % Given VSIBs for 9x ASRs (incl. reconstruction), derive channel saturation
    % QRCBs.
    %
    % NOTE: One can only obtain channel saturation on the 9x ASRs *AFTER*
    %       reconstruction, using the already set VSIBs contained within
    %       the SCHDs.
    %       Ex: V1 and V2 are saturated by being high ==> V12 = 0 (approx.)
    %           ==> V12 does not appear saturated, but the reconstructed V12
    %               value is still affected by the saturation on V1 and V2.
    %
    %
    % IMPLEMENTATION NOTE: Using moving window algo. (CWF) on ASRs instead of
    %                      on BLTSs
    % =======================================================================
    % This function applies the moving window algorithm
    % (bicas.utils.sliding_window_over_fraction()) on the ASR (CWF)
    % channels rather than on the BLTSs, since one can then avoid possible bad
    % VSQB behaviour.
    %
    % Ex 1: Consider a reconstructed channel when a VSQB (after moving window)
    % ends on one underlying source channel and begins on another, at roughly
    % the same time: If there is a non-saturated hole (VSQB=0) between the end
    % and beginning, then the reconstructed channel will have a non-saturated
    % hole (VSQB=0) that might never have existed if moving window was applied
    % on the reconstructed channel's threshold saturation bits instead of the
    % two source channels.
    %
    % Ex 2: If the two source channels for a reconstructed channel separately
    % contain too few threshold saturated samples for the moving window algo. to
    % set VSQB=1 (e.g. 30% for a window fraction 50%), but enough when combined
    % (e.g. 30%+30% > 50%), then the reconstructed channel's saturation bits
    % would be zero, despite being very much affected by saturation.
    %
    function Qrcbm = VSIBs_to_channel_saturation_QRCBs(...
        VsibZvm, tt2000Ar, isSwf, ...
        vstbFractionThreshold, cwfSlidingWindowLengthSec)

      assert(isa(VsibZvm, "bicas.utils.ZvMap"))
      assert(isscalar(isSwf) & islogical(isSwf))

      % IMPLEMENTATION NOTE: containers.Map does not support string-valued keys
      % (sic!)
      Qrcbm = bicas.proc.QrcbMap(numel(tt2000Ar));



      % Update Qrcbm wrt. the corresponding arguments.
      function handle_channel(sdidStr, channelSaturationQrcid)
        sdid       = bicas.proc.L1L2.const.C.SDID_DICT(sdidStr);
        vsibSdidAr = VsibZvm.get(sdid);

        if ~isSwf
          vsibSdidAr = bicas.utils.sliding_window_over_fraction(...
            tt2000Ar, vsibSdidAr, ...
            vstbFractionThreshold, cwfSlidingWindowLengthSec);
        end

        % IMPLEMENTATION NOTE: The (nested) function is called for the same
        % QRCID up to two times.
        if Qrcbm.has_QRCID(channelSaturationQrcid)
          Qrcbm.set(channelSaturationQrcid, ...
            Qrcbm.get(channelSaturationQrcid) | vsibSdidAr)
        else
          Qrcbm.add(channelSaturationQrcid, vsibSdidAr);
        end
      end



      handle_channel("DC_V1",  "SATURATION_ZV_V1")
      handle_channel("DC_V2",  "SATURATION_ZV_V2")
      handle_channel("DC_V3",  "SATURATION_ZV_V3")
      handle_channel("DC_V12", "SATURATION_ZV_V12")
      handle_channel("DC_V13", "SATURATION_ZV_V13")
      handle_channel("DC_V23", "SATURATION_ZV_V23")
      handle_channel("AC_V12", "SATURATION_ZV_V12")
      handle_channel("AC_V13", "SATURATION_ZV_V13")
      handle_channel("AC_V23", "SATURATION_ZV_V23")
    end



    % Convert channel saturation QRCBs to global saturation QRCBs.
    %
    % NOTE: Always sets QRCID="PARTIAL_SATURATION" QRCB=false since it can not
    %       be autodetected.
    % NOTE: Produces QRCB values using both the old scheme ("GLOBAL_SATURATION")
    %       and new scheme ("CHANNEL_SATURATION") for backward-compatibility
    %       during development and testing. "GLOBAL_SATURATION" should be phased
    %       out eventually. /Erik P G Johansson, 2025-07-08
    %
    % ARGUMENTS
    % =========
    % ChannelSaturationQrcbm
    %       Must contain at least all the CHANNEL_SATURATION QRCIDs, but may
    %       contain more which are then ignored.
    %
    function GlobalSaturationQrcbm = ...
        channel_saturation_to_global_saturation_QRCBs( ...
        ChannelSaturationQrcbm)

      assert(isa(ChannelSaturationQrcbm, 'bicas.proc.QrcbMap'))

      nRecords = ChannelSaturationQrcbm.nRecords;

      fullSaturationQrcbAr = false(nRecords, 1);
      for qrcid = bicas.const.qrc.Q.CHANNEL_SATURATION_QRCID_AR'
        fullSaturationQrcbAr = ...
          fullSaturationQrcbAr | ChannelSaturationQrcbm.get(qrcid);
      end

      GlobalSaturationQrcbm = bicas.proc.QrcbMap(nRecords);

      GlobalSaturationQrcbm.add("FULL_SATURATION",    fullSaturationQrcbAr);
      % NOTE: Always false since partial saturation (QRC) can not be
      % autodetected.
      GlobalSaturationQrcbm.add_false("PARTIAL_SATURATION");
    end



    % Set samples based on NSO table and SSIDs.
    %
    % ARGUMENTS
    % =========
    % voltageAr
    %       Float. 1D or more dimensions. First dimension is CDF records.
    %       NOTE: Does not have to have any particular unit.
    % ssidAr
    %       Same size as samplesAr. Same number of rows as Qrcbm.
    % Qrcbm
    % Qrcsm
    %       Must contain the same keys as Qrcbm.
    %
    function voltageAr = set_5xBLTS_voltage_samples_FV(...
        voltageAr, ssidAr, Qrcbm, Qrcsm)

      % IMPLEMENTATION NOTE: Input arrays samplesAr & ssidAr must have same
      % arbitrary size (not arbitrary for first dimension).
      %   PRO: It simplifies testing.
      %   PRO: It simplifies the implementation.
      %   PRO: It makes the function more generic.
      % IMPLEMENTATION NOTE: Argument Qrcsm is there instead of global
      % constant to make testing simpler & more robust.

      % ASSERTIONS
      assert(isfloat(voltageAr))
      assert(isequal(size(voltageAr), size(ssidAr)))
      assert(isa(Qrcbm, "bicas.proc.QrcbMap"))
      assert(isa(Qrcsm, "bicas.proc.QrcSettingsMap"))
      assert(isequal(Qrcbm.qrcidAr, Qrcsm.qrcidAr))
      assert(Qrcbm.nRecords == size(voltageAr, 1))    % Nbr. of records.

      sizeAr = size(voltageAr);
      bFv    = false(sizeAr);
      for qrcid = Qrcbm.qrcidAr'
        qrcbAr = Qrcbm.get(qrcid);    % (nRecords, 1)
        Qrcs   = Qrcsm.get(qrcid);

        % Arrays of the same size as voltageAr.
        bQrcbAr    = repmat(qrcbAr, [1, sizeAr(2:end)]);
        bSsidMatch = ismember(ssidAr, Qrcs.voltageFvSsidAr);

        bFv        = bFv | (bQrcbAr & bSsidMatch);
      end
      voltageAr(bFv) = NaN;
    end



    % Overwrite records of current with FVs as specified in QRCBs.
    % Cf. bicas.proc.L1L2.qrc.set_5xBLTS_voltage_samples_FV().
    %
    % ARGUMENTS
    % =========
    % currentAr
    %       Float. Size (nRecords, 3).
    % Qrcbm
    % Qrcsm
    %       Must contain the same keys as Qrcbm.
    %
    function currentAr = set_current_samples_FV(currentAr, Qrcbm, Qrcsm)
      % IMPLEMENTATION NOTE: Argument Qrcsm is there (instead of using a
      % global constant) to make test code simpler & more robust.

      assert(isfloat(currentAr))
      assert(isa(Qrcbm, 'bicas.proc.QrcbMap'))
      assert(isa(Qrcsm, 'bicas.proc.QrcSettingsMap'))
      assert(isequal(Qrcbm.qrcidAr, Qrcsm.qrcidAr))
      irf.assert.sizes(currentAr, [Qrcbm.nRecords, 3])

      % PROPOSAL: Create and use bAntennas = logical size (1, 3) + repmat.
      %   PRO: Faster?

      iAntennaAr = repmat([1:3], [Qrcbm.nRecords, 1]);
      bFv        = false(size(currentAr));
      for qrcid = Qrcbm.qrcidAr'
        qrcbAr = Qrcbm.get(qrcid);    % (nRecords, 1)
        Qrcs   = Qrcsm.get(qrcid);

        % Arrays of the same size as currentAr.
        bQrcbAr       = repmat(qrcbAr, [1, 3]);
        bAntennaMatch = ismember(iAntennaAr, Qrcs.currentFvIantAr);

        bFv           = bFv | (bQrcbAr & bAntennaMatch);
      end
      currentAr(bFv) = NaN;
    end






  end    % methods(Static)



end
