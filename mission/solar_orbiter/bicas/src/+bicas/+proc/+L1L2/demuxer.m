%
% "Encode" the demultiplexer part of the BIAS subsystem.
% See
%   bicas.proc.L1L2.demuxer.relabel_reconstruct_samples_BLTS_to_ASR_subsequence()
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
% The demultiplexer code is designed to not be aware of that TDS only digitizes
% BLTS 1-3 (not 4-5) and does not need to be.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-18
%
classdef demuxer
% PROPOSAL: Better name.
%   PRO: Code really implements "relabel BLTS to ASR" and "reconstructing
%        ASRs".



  methods(Static, Access=public)



    % For a specified BDM and DLR setting, for every BLTS, determine
    % (1) the (type of physical) input signal (Antennas, GND, "2.5V Ref",
    %     unknown), and
    % (2) as what ASR (SDID) (if any) should the corresponding samples be stored
    %     as in the datasets.
    %     NOTE: Output datasets only contain dedicated zVariables for ASRs, not
    %     e.g. "2.5V REF". Such non-ASR samples must still be routed to one of
    %     those same zVariabes.
    %
    % This function therefore "encodes" the demultiplexer part of the BIAS
    % subsystem (and more).
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
    %    Ex: Unknown BDM, e.g. due to insufficient HK time coverage (must use HK
    %    for TDS).
    % ** dlr = fill position
    %    Ex: Unknown DLR, e.g. due to insufficient HK time coverage.
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
    %
    function RoutingArray = get_routings(bdmFpa, dlrFpa)
      assert(isscalar(bdmFpa) && isa(bdmFpa, 'bicas.utils.FPArray') && strcmp(bdmFpa.mc, 'uint8'))
      assert(isscalar(dlrFpa) && isa(dlrFpa, 'bicas.utils.FPArray') && strcmp(dlrFpa.mc, 'logical'))

      R = bicas.proc.L1L2.const.C.ROUTING_DICT;

      dlrFloat = dlrFpa.logical2doubleNan();
      if isnan(dlrFloat)
        DC_V1x_DLR = R("UNKNOWN_TO_NOWHERE");
        AC_V1x_DLR = R("UNKNOWN_TO_NOWHERE");
      elseif dlrFloat
        DC_V1x_DLR = R("DC_V13");
        AC_V1x_DLR = R("AC_V13");
      else
        DC_V1x_DLR = R("DC_V12");
        AC_V1x_DLR = R("AC_V12");
      end

      % IMPLEMENTATION NOTE: switch-case statement does not work for NaN.
      % Therefore using local non-NaN fill value.
      BDM_INT_FV = uint8(255);
      bdmInt = bdmFpa.array(BDM_INT_FV);
      switch(bdmInt)

        case 0   % "Standard operation" : We have all information.

          % Summarize the routing.
          RoutingArray(1) = R("DC_V1");
          RoutingArray(2) =    DC_V1x_DLR;
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
          RoutingArray(3) =    DC_V1x_DLR;

        case 3   % Probe 3 fails

          RoutingArray(1) = R("DC_V1");
          RoutingArray(2) = R("DC_V2");
          RoutingArray(3) =    DC_V1x_DLR;

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

      RoutingArray(4) =    AC_V1x_DLR;
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
    % sdidArray
    %       1x5 array of SSIDs. One SDID per BLTS.
    % bltsSamplesAVolt
    %       1x5 cell array.
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
    function AsrSamplesAVoltSrm = relabel_reconstruct_samples_BLTS_to_ASR_subsequence(...
        sdidArray, bltsSamplesAVolt)
      % PROPOSAL: Log message for BDM=NaN.

      % ASSERTIONS
      assert(isa(sdidArray, 'uint8'))
      assert(isnumeric(bltsSamplesAVolt))
      irf.assert.sizes(...
        bltsSamplesAVolt, [-1, -2, bicas.const.N_BLTS], ...
        sdidArray,        [ 1,     bicas.const.N_BLTS])

      % Assign arrays only for those ASIDs for which there is data.
      AsrSamplesAVoltSrm = bicas.proc.L1L2.demuxer.relabel_samples_BLTS_to_ASR_subsequence(...
        bltsSamplesAVolt, sdidArray);

      % Assign arrays for the remaining ASIDs. Reconstruct data when possible.
      % NOTE: The function modifies the ARGUMENT (handle object).
      bicas.proc.L1L2.demuxer.reconstruct_ASR_samples_subsequence(AsrSamplesAVoltSrm);
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
    function reconstruct_ASR_samples_subsequence(AsrSamplesAVoltSrm)
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
      A     = bicas.proc.L1L2.const.C.ASID_DICT;
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

      useAsidArray = AsSrm.keys;

      % IMPLEMENTATION NOTE: Can not use bicas.utils.SameRowsMap methods
      % for deriving the entire size (samples per record), until possibly
      % using a future bicas.utils.SameSizeTypeMap instead.
      tempNaN = nan(size(AsrSamplesAVoltSrm(useAsidArray(1))));

      for asid = bicas.proc.L1L2.const.C.ASID_DICT.values'
        if ~AsSrm.isKey(asid)
          AsSrm.add(asid, tempNaN);
        end
      end

    end



    % EXPERIMENTAL, UNUSED FUNCTION
    %
    % Intended as future conceptual replacement for
    % bicas.proc.L1L2.demuxer.complete_relation().
    %
    % Complement/derive redundant data in three same-sized arrays. If one
    % element is missing (labelled as a fill position) while the corresponding
    % elements in the two other arrays are not, then the missing element is
    % derived from the other two.
    %
    % The function does not verify that pre-existing data is consistent with
    % the specified functions which define the relationship between the
    % functions.
    %
    %
    % ARGUMENTS
    % =========
    % A1, A2, A3
    %       Data arrays of the same MATLAB class and size.
    % bFp1, bFp2, bFp3
    %       Logical arrays of the same size as A1.
    %       True=fill position in corresponding A* array.
    % fh23to1, fh13to2, fh12to3
    %       Function handles z=f(x,y) for how to derive missing elements. Must
    %       be vectorized.
    %
    function [A1,A2,A3] = reconstruct_missing_data(A1,A2,A3, bFp1,bFp2,bFp3, fh23to1, fh13to2, fh12to3)
      % TODO-NI: Will lead to memory problems for large arrays?
      %   TODO: Test.
      %   NOTE: Using cell array for return value may cause more in-memory data
      %         copying.

      assert(strcmp(class(A1), class(A2)))
      assert(strcmp(class(A1), class(A3)))

      assert(islogical(bFp1))
      assert(islogical(bFp2))
      assert(islogical(bFp3))

      assert(isequal(size(A1), size(A2)))
      assert(isequal(size(A1), size(A3)))
      assert(isequal(size(A1), size(bFp1)))
      assert(isequal(size(A1), size(bFp2)))
      assert(isequal(size(A1), size(bFp3)))

      bDerive1 =  bFp1 & ~bFp2 & ~bFp3;
      bDerive2 = ~bFp1 &  bFp2 & ~bFp3;
      bDerive3 = ~bFp1 & ~bFp2 &  bFp3;

      % NOTE: If using actual arrays for samples, then need at least one more
      % dimension for SWF (1 record=1 snapshot): Ax(b, :) = func(Ay(b, :),
      % Az(b, :).
      A3(bDerive3) = fh12to3(A1(bDerive3), A2(bDerive3));
      A2(bDerive2) = fh13to2(A1(bDerive2), A3(bDerive2));
      A1(bDerive1) = fh23to1(A2(bDerive1), A3(bDerive1));
    end



    % EXPERIMENTAL, UNUSED FUNCTION
    %
    % Intended as future conceptual replacement for
    % bicas.proc.L1L2.demuxer.reconstruct_ASR_samples_subsequence() (though
    % "global", not for a subsequence).
    %
    function SdcdDict = reconstruct_ASR_samples2(SdcdDict)

      assert(isa(SdcdDict, 'bicas.proc.L1L2.SdChannelDataDict'))

      % Shorten variable names.
      D = bicas.proc.L1L2.const.C.SDID_DICT;



      % Reconstructs values using relationship SDID_1 = SDID_2 + SDID_3.
      function reconstruct_missing_data_helper(sumSdidStr1, termSdidStr2, termSdidStr3)
        % NOTE: Printout is very useful for being able to follow how values are
        % being reconstructed, e.g. when debugging and verifying automated
        % tests.
        % fprintf("%-6s = %-6s + %-6s\n", sumSdidStr1, termSdidStr2, termSdidStr3)

        Sdcd1 = SdcdDict.get(D( sumSdidStr1));
        Sdcd2 = SdcdDict.get(D(termSdidStr2));
        Sdcd3 = SdcdDict.get(D(termSdidStr3));

        [...
          Sdcd1, ...
          Sdcd2, ...
          Sdcd3 ...
        ] = ...
          bicas.proc.L1L2.demuxer.reconstruct_missing_data(...
          Sdcd1, ...
          Sdcd2, ...
          Sdcd3, ...
          Sdcd1.bFp, ...
          Sdcd2.bFp, ...
          Sdcd3.bFp, ...
          @(x,y) (x+y), ...
          @(x,y) (x-y), ...
          @(x,y) (x-y) ...
          );

        SdcdDict = SdcdDict.set(D( sumSdidStr1), Sdcd1);
        SdcdDict = SdcdDict.set(D(termSdidStr2), Sdcd2);
        SdcdDict = SdcdDict.set(D(termSdidStr3), Sdcd3);
      end



      %================
      % Derive AC ASRs
      %================
      % AC ASRs are separate from DC ASRs and only satisfy one relationship
      % since there are only three of them. Therefore does not have to be in
      % loop.
      reconstruct_missing_data_helper("AC_V13",   "AC_V12", "AC_V23")

      %================
      % Derive DC ASRs
      %================
      nFp0 = SdcdDict.nFp;
      while true
        % NOTE: Relation DC_V13 = DC_V12 + DC_V23 has precedence for deriving
        % diffs (i.e. it should come first) since it is better to derive a diff
        % from (initially available) diffs rather than singles, directly or
        % indirectly, if possible.
        % Ex: The is only information on V1, V12, V23.
        % ==> Derive V13 first, then use V1 to derive V2, V3.
        % ==> If V1 is lower-accuracy, then it will not affect V13.
        %     If V1 is saturated,      then it will not affect V13.
        % Note that this is due to that sets of three diffs are always
        % redundant (contain redundant information), as opposed to sets of
        % three singles which are not.

        reconstruct_missing_data_helper("DC_V13",  "DC_V12", "DC_V23");

        reconstruct_missing_data_helper("DC_V1",   "DC_V12", "DC_V2");
        reconstruct_missing_data_helper("DC_V1",   "DC_V13", "DC_V3");
        reconstruct_missing_data_helper("DC_V2",   "DC_V23", "DC_V3");

        % NOTE: Impossible to get Dcd.nFp == 0...
        if (SdcdDict.nFp == nFp0) || (SdcdDict.nFp == 0)
          break
        end
        nFp0 = SdcdDict.nFp;
      end
    end



  end    % methods(Static, Access=public)



  %###########################################################################



  methods(Static, Access=private)



    % Given FIVE BLTS sample arrays, copy those which correspond to ASRs (five
    % or fewer!) into an SRM with correct ASID keys for the corresponding
    % arrays.
    function AsrSamplesSrm = relabel_samples_BLTS_to_ASR_subsequence(...
        bltsSamplesAVolt, sdidArray)

      % ASSERTIONS
      assert(isnumeric(bltsSamplesAVolt))
      nRows = irf.assert.sizes( ...
        bltsSamplesAVolt, [-1, -2, bicas.const.N_BLTS], ...
        sdidArray,        [ 1,     bicas.const.N_BLTS]);

      AsrSamplesSrm = bicas.utils.SameRowsMap("uint8", nRows, 'EMPTY');
      for iBlts = 1:bicas.const.N_BLTS
        if ~bicas.proc.L1L2.const.SDID_is_nowhere(sdidArray(iBlts))
          % NOTE: Converting from SDID to ASID and using ASID as key. Not sure
          % if conceptually sensible.
          asid = bicas.proc.L1L2.const.SDID_ASR_to_ASID(sdidArray(iBlts));

          AsrSamplesSrm.add(asid, bltsSamplesAVolt(:, :, iBlts));
        end
      end
    end



    % Utility function. Derive missing ASR fields/channels from other
    % fields/channels. If exactly two of the SRM keys exist in AsrSrm, then
    % derive the third one using the relationship
    % AsSrm(asid1) == AsSrm(asid2) + AsSrm(asid3).
    %
    % ARGUMENTS
    % =========
    % asid1, asid2, asid3
    %       ASIDs whose ID strings may or may not be keys in AsSrm. If
    %       exactly one of them is missing in "As", then the key+value is
    %       created with values assuming that the field contents are related
    %       through the relationship value1 = value2 + value3. In other
    %       cases, "AsSrm" is returned unmodified.
    %
    function AsSrm = complete_relation(AsSrm, asid1, asid2, asid3)
      % PROPOSAL: Handle propagating quality bits derived from BLTSs.
      %   NOTE: BLTS saturation bits should propagate.
      %   This is thus a qualitatively different behaviour from samples which
      %   are reconstructed.
      %
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

      e1 = AsSrm.isKey(asid1);
      e2 = AsSrm.isKey(asid2);
      e3 = AsSrm.isKey(asid3);

      if     ~e1 &&  e2 &&  e3   AsSrm.add(asid1, AsSrm(asid2) + AsSrm(asid3));
      elseif  e1 && ~e2 &&  e3   AsSrm.add(asid2, AsSrm(asid1) - AsSrm(asid3));
      elseif  e1 &&  e2 && ~e3   AsSrm.add(asid3, AsSrm(asid1) - AsSrm(asid2));
      end
    end



  end    % methods(Static, Access=public)



end    % classdef
