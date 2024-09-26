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
    % For a specified BDM and DLR setting, for every BLTS, it determines
    % (1) the (physical) input signal (Antennas, GND, "2.5V Ref",
    %     unknown), and
    % (2) as what ASR (SDID) (if any) should the BLTS be represented in
    %     the datasets.
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
    % ** bdm = fill position
    %    Ex: Unknown BDM, e.g. due to insufficient HK time coverage.
    % ** BLTS 1-3 labelled as "GND" or "2.5V Ref" in BDMs 5-7.
    %
    %
    % ARGUMENTS
    % =========
    % bdmFpa
    %       Scalar value.
    %       NOTE: Can be fill position to represent unknown BDM.
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
      assert(isscalar(bdmFpa) && isa(bdmFpa, 'bicas.utils.FPArray') && strcmp(bdmFpa.mc, 'uint8'))
      assert(isscalar(dlrFpa) && isa(dlrFpa, 'bicas.utils.FPArray') && strcmp(dlrFpa.mc, 'logical'))

      R = bicas.sconst.C.S_ROUTING_DICT;

      dlrFloat = dlrFpa.logical2doubleNan();
      if isnan(dlrFloat)
        DC_V1x = R("UNKNOWN_TO_NOWHERE");
        AC_V1x = R("UNKNOWN_TO_NOWHERE");
      elseif dlrFloat
        DC_V1x = R("DC_V13");
        AC_V1x = R("AC_V13");
      else
        DC_V1x = R("DC_V12");
        AC_V1x = R("AC_V12");
      end

      % IMPLEMENTATION NOTE: switch-case statement does not work for NaN.
      % Therefore using local non-NaN fill value.
      BDM_INT_FV = uint8(255);
      bdmInt = bdmFpa.array(BDM_INT_FV);
      switch(bdmInt)

        case 0   % "Standard operation" : We have all information.

          % Summarize the routing.
          RoutingArray(1) = R("DC_V1");
          RoutingArray(2) =    DC_V1x;
          RoutingArray(3) = R("DC_V23");

        case 1   % Probe 1 fails

          RoutingArray(1) = R("DC_V2");
          RoutingArray(2) = R("DC_V3");
          RoutingArray(3) = R("DC_V23");

          % NOTE: Can not derive anything extra for DC. BLTS 1-3
          % contain redundant data (regardless of DLR).

        case 2   % Probe 2 fails

          RoutingArray(1) = R("DC_V1");
          RoutingArray(2) = R("DC_V3");
          RoutingArray(3) =   DC_V1x;

        case 3   % Probe 3 fails

          RoutingArray(1) = R("DC_V1");
          RoutingArray(2) = R("DC_V2");
          RoutingArray(3) =    DC_V1x;

        case 4   % Calibration mode 0

          RoutingArray(1) = R("DC_V1");
          RoutingArray(2) = R("DC_V2");
          RoutingArray(3) = R("DC_V3");

        case {5,6,7}   % Calibration mode 1/2/3

          switch(bdmInt)
            case 5
              RoutingArray(1) = R("REF25V_TO_DC_V1");
              RoutingArray(2) = R("REF25V_TO_DC_V2");
              RoutingArray(3) = R("REF25V_TO_DC_V3");
            case {6,7}
              RoutingArray(1) = R("GND_TO_DC_V1");
              RoutingArray(2) = R("GND_TO_DC_V2");
              RoutingArray(3) = R("GND_TO_DC_V3");
          end

        case BDM_INT_FV
          % CASE: Fill position
          % -------------------
          % NOTE: Could route unknown DC signals to V1-V3, but
          % since this behaviour is probably not very obvious to
          % the user, the code effectively deletes the information
          % instead.
          RoutingArray(1) = R("UNKNOWN_TO_NOWHERE");
          RoutingArray(2) = R("UNKNOWN_TO_NOWHERE");
          RoutingArray(3) = R("UNKNOWN_TO_NOWHERE");

          % NOTE: The routing of BLTS 4 & 5 is identical for all BDMs
          % (but does depend on the DLR). Can therefore route them
          % also when the BDM is unknown.

        otherwise
          error('BICAS:Assertion:IllegalArgument:DatasetFormat', ...
            'Illegal argument value bdm=%g.', bdm)
      end    % switch

      RoutingArray(4) = AC_V1x;
      RoutingArray(5) = R("AC_V23");
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
      assert(isa(SdidArray, 'bicas.proc.L1L2.SignalDestinationId'))
      assert(isnumeric(bltsSamplesAVolt))
      irf.assert.sizes(...
        bltsSamplesAVolt, [-1, -2, bicas.const.N_BLTS], ...
        SdidArray,        [ 1,     bicas.const.N_BLTS])

      % Assign arrays only for those ASIDs for which there is data.
      AsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.assign_ASR_samples_from_BLTS(...
        bltsSamplesAVolt, SdidArray);

      % Assign arrays for the remaining ASIDs. Reconstruct data when possible.
      % NOTE: The function modifies the ARGUMENT (handle object).
      bicas.proc.L1L2.demuxer.reconstruct_missing_ASR_samples(AsrSamplesAVoltSrm);
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
    %   PROBLEM: Can happen if blanking data in future implementations due to
    %   NSOs.
    %
    % NOTE: Only public for the purpose of automatic testing.
    %
    %
    % ARGUMENTS
    % =========
    % AsrSamplesAVoltSrm
    %       NOTE: Function modifies the argument (handle class)!
    %
    function reconstruct_missing_ASR_samples(AsrSamplesAVoltSrm)
      % PROPOSAL: Better name
      %   NOTE: Can not always reconstruct all signals.
      %   --
      %   ASR
      %   ASID
      %   complete
      %   complement
      %   reconstruct
      %   signals
      %   missing (antenna) signals
      %   singles, diffs
      %   --
      %   reconstruct_missing_signals
      %   reconstruct_missing_antenna_signals
      %   reconstruct_missing_ASRs
      %   add_missing_ASRs
      %   add_missing_ASIDs
      %   add_reconstruct_missing_ASRs
      %   --
      %   If using SDID keys (in the future), then
      %     add_missing_SDIDs

      assert(isa(AsrSamplesAVoltSrm, 'bicas.utils.SameRowsMap'))

      % Shorten variable names.
      A     = bicas.sconst.C.S_ASID_DICT;
      AsSrm = AsrSamplesAVoltSrm;

      %================
      % Derive AC ASRs
      %================
      % AC ASRs are separate from DC. Does not have to be in loop.
      % IMPLEMENTATION NOTE: Must be executed before DC loop. Otherwise
      % nFnAfter == 9 condition does not work.
      AsSrm = bicas.proc.L1L2.demuxer.complete_relation(...
        AsSrm, A("AC_V13"), A("AC_V12"), A("AC_V23"));

      %================
      % Derive DC ASRs
      %================
      nAsidBefore = AsSrm.numEntries;
      while true
        % NOTE: Relation DC_V13 = DC_V12 + DC_V23 has precedence for
        % deriving diffs since it is better to derive a diff from
        % (initially available) diffs rather than singles, directly or
        % indirectly, if possible.
        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, A("DC_V13"), A("DC_V12"), A("DC_V23"));

        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, A("DC_V1"),  A("DC_V12"), A("DC_V2"));
        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, A("DC_V1"),  A("DC_V13"), A("DC_V3"));
        AsSrm = bicas.proc.L1L2.demuxer.complete_relation(AsSrm, A("DC_V2"),  A("DC_V23"), A("DC_V3"));
        nAsidAfter = AsSrm.numEntries;

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

      keyArray = AsSrm.keys;

      % IMPLEMENTATION NOTE: Can not use bicas.utils.SameRowsMap methods
      % for deriving the entire size (samples per record), until possibly
      % using a future bicas.utils.SameSizeTypeMap instead.
      tempNaN = nan(size(AsrSamplesAVoltSrm(keyArray)));

      for Asid = bicas.sconst.C.S_ASID_DICT.values'
        if ~AsSrm.isKey(Asid)
          AsSrm.add(Asid, tempNaN);
        end
      end

    end



  end    % methods(Static, Access=public)



  %###########################################################################



  methods(Static, Access=private)



    % Given FIVE BLTS sample arrays, copy those which correspond to ASRs (five
    % or fewer!) into an SRM with correct ASID keys for the corresponding
    % arrays.
    function AsrSamplesSrm = assign_ASR_samples_from_BLTS(...
        bltsSamplesAVolt, SdidArray)

      % ASSERTIONS
      assert(isnumeric(bltsSamplesAVolt))
      nRows = irf.assert.sizes( ...
        bltsSamplesAVolt, [-1, -2, bicas.const.N_BLTS], ...
        SdidArray,        [ 1,     bicas.const.N_BLTS]);

      AsrSamplesSrm = bicas.utils.SameRowsMap( ...
        "bicas.proc.L1L2.AntennaSignalId", nRows, 'EMPTY');
      for iBlts = 1:bicas.const.N_BLTS
        if ~SdidArray(iBlts).isNowhere
          % NOTE: Converting from SDID to ASID and using ASID as key. Not sure
          % if conceptually sensible.
          Asid = SdidArray(iBlts).Asid;
          AsrSamplesSrm.add(Asid, bltsSamplesAVolt(:, :, iBlts));
        end
      end
    end



    % Utility function. Derive missing ASR fields/channels from other
    % fields/channels. If exactly two of the SRM keys exist in AsrSrm, then
    % derive the third one using the relationship
    % AsSrm(Asid1) == AsSrm(Asid2) + AsSrm(Asid3).
    %
    % ARGUMENTS
    % =========
    % Asid1, Asid2, Asid3
    %       ASIDs whose ID strings may or may not be keys in AsSrm. If
    %       exactly one of them is missing in "As", then the key+value is
    %       created with values assuming that the field contents are related
    %       through the relationship value1 = value2 + value3. In other
    %       cases, "AsSrm" is returned unmodified.
    %
    function AsSrm = complete_relation(AsSrm, Asid1, Asid2, Asid3)
      % PROPOSAL: Handle propagating quality bits derived from BLTSs.
      %   NOTE: BLTS saturation bits should propagate.
      %   This is thus a qualitatively different behaviour from samples which
      %   are reconstructed.
      % PROPOSAL: More generic implementation independent of SRMs and ASIDs.
      %   PROBLEM: Needs way to detect missing data (channels) and a way to add
      %            the reconstructed data.
      %     PROPOSAL: Input is FPAs for all ASRs/SDIDs.
      % PROPOSAL: Handle both samples and quality bits.
      %
      % NOTE: Truly generic vectorized algorithm. FP = Fill Position as in FPAs.
      %   NOTE: Can handle both
      %     (1) reconstructing data (e.g. diffs from singles), and
      %     (2) "spread" of data (e.g. saturation bits).
      %   bSet1 =  bFp1 & ~bFp2 & ~bFp3
      %   bSet2 = ~bFp1 &  bFp2 & ~bFp3
      %   bSet3 = ~bFp1 & ~bFp2 &  bFp3
      %   A1(bSet1) = func_1_from_23(A2(bSet1), A3(bSet1))
      %   A2(bSet2) = func_2_from_13(A1(bSet2), A3(bSet2))
      %   A3(bSet3) = func_3_from_12(A1(bSet3), A2(bSet3))
      assert(isa(AsSrm, 'bicas.utils.SameRowsMap'))

      e1 = AsSrm.isKey(Asid1);
      e2 = AsSrm.isKey(Asid2);
      e3 = AsSrm.isKey(Asid3);

      if     ~e1 &&  e2 &&  e3   AsSrm.add(Asid1, AsSrm(Asid2) + AsSrm(Asid3));
      elseif  e1 && ~e2 &&  e3   AsSrm.add(Asid2, AsSrm(Asid1) - AsSrm(Asid3));
      elseif  e1 &&  e2 && ~e3   AsSrm.add(Asid3, AsSrm(Asid1) - AsSrm(Asid2));
      end
    end



  end    % methods(Static, Access=public)



end    % classdef
