%
% "Encode" the demultiplexer part of the BIAS subsystem.
% See
%   bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_all_ASRs()
%   bicas.proc.L1L2.AntennaSignalId
%   bicas.proc.L1L2.SignalSourceId
%   bicas.proc.L1L2.SignalDestinationId
%
%
% NOTE
% ====
% It is in principle arbitrary (probably) how the GND and "2.5V Ref" signals,
% which are generated by the instrument, should be represented in the datasets,
% since the datasets assume that only assumes signals from the antennas. The
% implementation classifies them as antennas, including for diffs, but the
% signalTypeCategory specifies that they should be calibrated differently. In
% principle, one could represent unknown signals (unknown mux mode) as antenna
% signals too.
% --
% Demultiplexer is designed to not be aware of that TDS only digitizes BLTS 1-3
% (not 4-5) and does not need to be.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-18
%
classdef demuxer



  methods(Static, Access=public)



    % Function that "encodes" the demultiplexer part of the BIAS subsystem.
    % For a specified BDM and DLR setting, it determines, for every BLTS, it
    % associates
    % (1) which (physical) input signal (Antennas, GND, "2.5V Ref",
    %     unknown), and
    % (2) as what ASR (if any) should the BLTS be represented in the
    %     datasets.
    %
    %
    % RATIONALE
    % =========
    % Meant to collect all hard-coded information about the demultiplexer
    % ROUTING of signals in the BIAS specification, Table 4.
    %
    %
    % EDGE CASES
    % ==========
    % Function must be able to handle:
    % ** bdm = NaN
    %    Unknown BDM, e.g. due to insufficient HK time coverage.
    % ** BLTS 1-3 signals labelled as "GND" or "2.5V Ref" in BDMs 5-7.
    %
    %
    % ARGUMENTS
    % =========
    % bdm
    %       Scalar value. Demultiplexer mode.
    %       NOTE: Can be NaN to represent unknown BDM.
    %       Implies that AsrSamplesVolt fields are correctly
    %       sized with NaN values.
    % dlrFpa
    %       Scalar value.
    %
    %
    % RETURN VALUES
    % =============
    % RoutingArray
    %       Array of bicas.proc.L1L2.Routing objects, one per BLTS.
    %       (iBlts).
    function RoutingArray = get_routings(bdmFpa, dlrFpa)
      % switch-case checks bdmFpa values.
      assert(isscalar(bdmFpa) && isa(bdmFpa, 'bicas.utils.FPArray') && strcmp(bdmFpa.mc, 'uint8'))
      assert(isscalar(dlrFpa) && isa(dlrFpa, 'bicas.utils.FPArray') && strcmp(dlrFpa.mc, 'logical'))

      R = bicas.proc.L1L2.Routing.C;

      dlrFloat = dlrFpa.logical2doubleNan();
      if isnan(dlrFloat)
        R.DC_V1x = R.UNKNOWN_TO_NOWHERE;
        R.AC_V1x = R.UNKNOWN_TO_NOWHERE;
      elseif dlrFloat
        R.DC_V1x = R.DC_V13;
        R.AC_V1x = R.AC_V13;
      else
        R.DC_V1x = R.DC_V12;
        R.AC_V1x = R.AC_V12;
      end

      % IMPLEMENTATION NOTE: switch-case statement does not work
      % for NaN. Therefore not using NaN for
      bdmInt = bdmFpa.array(uint8(255));
      switch(bdmInt)

        case 0   % "Standard operation" : We have all information.

          % Summarize the routing.
          RoutingArray(1) = R.DC_V1;
          RoutingArray(2) = R.DC_V1x;
          RoutingArray(3) = R.DC_V23;

        case 1   % Probe 1 fails

          RoutingArray(1) = R.DC_V2;
          RoutingArray(2) = R.DC_V3;
          RoutingArray(3) = R.DC_V23;

          % NOTE: Can not derive anything extra for DC. BLTS 1-3
          % contain redundant data (regardless of DLR).

        case 2   % Probe 2 fails

          RoutingArray(1) = R.DC_V1;
          RoutingArray(2) = R.DC_V3;
          RoutingArray(3) = R.DC_V1x;

        case 3   % Probe 3 fails

          RoutingArray(1) = R.DC_V1;
          RoutingArray(2) = R.DC_V2;
          RoutingArray(3) = R.DC_V1x;

        case 4   % Calibration mode 0

          RoutingArray(1) = R.DC_V1;
          RoutingArray(2) = R.DC_V2;
          RoutingArray(3) = R.DC_V3;

        case {5,6,7}   % Calibration mode 1/2/3

          switch(bdmInt)
            case 5
              RoutingArray(1) = R.REF25V_TO_DC_V1;
              RoutingArray(2) = R.REF25V_TO_DC_V2;
              RoutingArray(3) = R.REF25V_TO_DC_V3;
            case {6,7}
              RoutingArray(1) = R.GND_TO_DC_V1;
              RoutingArray(2) = R.GND_TO_DC_V2;
              RoutingArray(3) = R.GND_TO_DC_V3;
          end

        case 255
          % NOTE: Could route unknown DC signals to V1-V3, but
          % since this behaviour is probably not very obvious to
          % the user, the code effectively deletes the information
          % instead.
          RoutingArray(1) = R.UNKNOWN_TO_NOWHERE;
          RoutingArray(2) = R.UNKNOWN_TO_NOWHERE;
          RoutingArray(3) = R.UNKNOWN_TO_NOWHERE;

          % NOTE: The routing of BLTS 4 & 5 is identical for all BDMs
          % (but does depend on the DLR). Can therefore route them
          % also when the BDM is unknown.

        otherwise
          error('BICAS:Assertion:IllegalArgument:DatasetFormat', ...
            'Illegal argument value bdm=%g.', bdm)
      end    % switch

      RoutingArray(4) = R.AC_V1x;
      RoutingArray(5) = R.AC_V23;
    end



    % (1) Given demultiplexer routings, convert the (already calibrated)
    %     BLTSs to (subset of) ASRs.
    % (2) Derive the remaining ASRs (samples) from those ASRs which have
    %     already been set.
    %       NOTE: This derivation from fully calibrated ASR samples only
    %       requires addition/subtraction of ASRs. It does not require any
    %       sophisticated/non-trivial calibration since the relationships
    %       between the ASRs are so simple. The only consideration is that
    %       DC diffs have higher accurracy than DC singles, and should have
    %       precedence when deriving ASRs in the event of redundant
    %       information.
    %
    % NOTE: This code does NOT handle the equivalent of demultiplexer
    % multiplication of the BLTS signal (alpha, beta, gamma in the BIAS
    % specification). It is assumed that the supplied BLTS samples have been
    % calibrated to account for this already.
    %
    %
    % ARGUMENTS
    % =========
    % SdidArray
    %       Length-5 array of SSIDs. One SDID per BLTS.
    % bltsSamplesAVolt
    %       Cell array of vectors/matrices, length 5.
    %       {iBlts} = Vector/matrix with sample values
    %                 for that BLTS channel, calibrated as for ASR.
    %
    %
    % RETURN VALUES
    % =============
    % AsrSamplesAVoltSrm
    %       Samples for all ASRs (singles, diffs) which can
    %       possibly be derived from the BLTS (BIAS_i). Those
    %       which can not be derived are correctly sized
    %       containing only NaN. Struct with fields.
    % --
    % NOTE: Separate names bltsSamplesAVolt & AsrSamplesAVoltSrm to denote
    % that they are organized by BLTS and ASRs respectively.
    %
    function AsrSamplesAVoltSrm = calibrated_BLTSs_to_all_ASRs(SdidArray, bltsSamplesAVolt)
      % PROPOSAL: Log message for BDM=NaN.

      % ASSERTIONS
      assert(numel(SdidArray) == bicas.const.N_BLTS)
      assert(isa(SdidArray, 'bicas.proc.L1L2.SignalDestinationId'))
      assert(isnumeric(bltsSamplesAVolt))
      irf.assert.sizes(bltsSamplesAVolt, [-1, -2, bicas.const.N_BLTS])

      % Set only those ASRs for which there is data.
      AsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.assign_ASR_samples_from_BLTS(...
        bltsSamplesAVolt, SdidArray);

      % Set those ASRs for which there is NO data.
      % NOTE: Argument (handle object) is modified.
      bicas.proc.L1L2.demuxer.complement_ASR(AsrSamplesAVoltSrm);
    end



    % Given an incomplete SRM of ASID-labelled samples, derive the missing
    % ASRs.
    %
    % NOTE: In the event of redundant FIELDS, but not redundant DATA
    % (non-fill value), the code can NOT make intelligent choice of only
    % using available data to replace fill values.
    %   Ex: BDM=1: Fields (V1, V2, V12) but V1 does not contain any data
    %   (fill values). Could in principle derive V1=V2-V12 but code does
    %   not know this.
    %   NOTE: Unlikely that this will ever happen, or that the
    %   instrument RPW is even able to return data for this situation.
    %
    % NOTE: Only public for the purpose of automatic testing.
    %
    %
    % ARGUMENTS
    % =========
    % AsrSamplesAVoltSrm
    %       NOTE: Modifies argument.
    %
    function complement_ASR(AsrSamplesAVoltSrm)
      assert(isa(AsrSamplesAVoltSrm, 'bicas.utils.SameRowsMap'))

      % Shorten variable names.
      C     = bicas.proc.L1L2.AntennaSignalId.C;
      AsSrm = AsrSamplesAVoltSrm;

      %================
      % Derive AC ASRs
      %================
      % AC ASRs are separate from DC. Does not have to be in loop.
      % IMPLEMENTATION NOTE: Must be executed before DC loop. Otherwise
      % nFnAfter == 9 condition does not work.
      AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, C.AC_V13, C.AC_V12, C.AC_V23);

      %================
      % Derive DC ASRs
      %================
      nAsidBefore = AsSrm.length;
      while true
        % NOTE: Relation DC_V13 = DC_V12 + DC_V23 has precedence for
        % deriving diffs since it is better to derive a diff from
        % (initially available) diffs rather than singles, directly or
        % indirectly, if possible.
        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, C.DC_V13, C.DC_V12, C.DC_V23);

        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, C.DC_V1,  C.DC_V12, C.DC_V2);
        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, C.DC_V1,  C.DC_V13, C.DC_V3);
        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, C.DC_V2,  C.DC_V23, C.DC_V3);
        nAsidAfter = AsSrm.length;

        if (nAsidBefore == nAsidAfter) || (nAsidAfter == 9)
          break
        end
        nAsidBefore = nAsidAfter;
      end

      %===================================================================
      % Add all ASIDs which have not yet been assigned
      % ----------------------------------------------
      % IMPLEMENTATION NOTE: This is needed to handle for situations when
      % the supplied fields can not be used to determine all nine fields.
      %   Ex: bdm=1,2,3
      %===================================================================

      keysCa = AsSrm.keys;

      % IMPLEMENTATION NOTE: Can not use bicas.utils.SameRowsMap methods
      % for deriving the entire size (samples per record), until possibly
      % using a future bicas.utils.SameSizeTypeMap instead.
      %tempNaN = nan(AsMap.nRows(), 1);
      tempNaN = nan(size(AsrSamplesAVoltSrm(keysCa{1})));

      for asidNameCa = bicas.proc.L1L2.AntennaSignalId.C.ALL_ASID_NAMES_CA'
        asidName = asidNameCa{1};
        if ~AsSrm.isKey(asidName)
          AsSrm.add(asidName, tempNaN);
        end
      end

    end



  end    % methods(Static, Access=public)



  %###########################################################################



  methods(Static, Access=private)



    % Given FIVE BLTS sample arrays, copy those which correspond to ASRs
    % (five or fewer!) into a bicas.utils.SameRowsMap.
    function AsrSamplesSrm = assign_ASR_samples_from_BLTS(...
        bltsSamplesAVolt, SdidArray)

      % ASSERTIONS
      assert(numel(SdidArray) == bicas.const.N_BLTS)
      assert(isnumeric(bltsSamplesAVolt))
      irf.assert.sizes(bltsSamplesAVolt, [-1, -2, bicas.const.N_BLTS])

      nRows = size(bltsSamplesAVolt, 1);
      AsrSamplesSrm = bicas.utils.SameRowsMap('char', nRows, 'EMPTY');
      for iBlts = 1:bicas.const.N_BLTS
        if ~SdidArray(iBlts).isNowhere
          AsrSamplesSrm.add(...
            SdidArray(iBlts).Asid.s, ...
            bltsSamplesAVolt(:, :, iBlts));
        end
      end
    end



    % Utility function. Derive missing ASR fields from other fields. If
    % exactly two of the Map keys exist in S, then derive the third using
    % the relationship AsMap(Asid1.s) == AsMap(Asid2.s) + AsMap(Asid3.s).
    %
    % ARGUMENTS
    % =========
    % Asid1, Asid2, Asid3
    %       ASIDs whose ID strings may or may not be keys in AsMap. If
    %       exactly one of them is missing in "As", then the key+value is
    %       created with values assuming that the field contents are related
    %       through the relationship value1 = value2 + value3. In other
    %       cases, "AsMap" is returned unmodified.
    %
    function AsSrm = complete_relation(AsSrm, Asid1, Asid2, Asid3)
      assert(isa(AsSrm, 'bicas.utils.SameRowsMap'))

      e1 = AsSrm.isKey(Asid1.s);
      e2 = AsSrm.isKey(Asid2.s);
      e3 = AsSrm.isKey(Asid3.s);

      if     ~e1 &&  e2 &&  e3   AsSrm.add(Asid1.s, AsSrm(Asid2.s) + AsSrm(Asid3.s));
      elseif  e1 && ~e2 &&  e3   AsSrm.add(Asid2.s, AsSrm(Asid1.s) - AsSrm(Asid3.s));
      elseif  e1 &&  e2 && ~e3   AsSrm.add(Asid3.s, AsSrm(Asid1.s) - AsSrm(Asid2.s));
      end
    end



  end    % methods(Static, Access=public)



end    % classdef
