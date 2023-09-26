%
% "Encode" the demultiplexer part of the BIAS subsystem.
% See
%   bicas.proc.L1L2.demuxer.calibrated_BLTSs_to_ASRs()
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
% DEFINITIONS
% ===========
% See bicas.proc.L1L2.cal.Cal, bicas.proc.L1L2.demuxer_latching_relay.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-18
%
classdef demuxer
    
    
    
    methods(Static, Access=public)



        % Function that "encodes" the demultiplexer part of the BIAS subsystem.
        % For a specified mux mode and demuxer latching relay setting, it
        % determines, for every BLTS, it associates
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
        % ** demuxMode = NaN                                   
        %    Unknown demux mode, e.g. due to insufficient HK time coverage.
        % ** BLTS 1-3 signals labelled as "GND" or "2.5V Ref" in demux modes 5-7.
        % NOTE: Can not hande unknown dlrUsing12.
        %
        %
        % ARGUMENTS
        % =========
        % demuxMode
        %       Scalar value. Demultiplexer mode.
        %       NOTE: Can be NaN to represent unknown demux mode.
        %       Implies that AsrSamplesVolt fields are correctly
        %       sized with NaN values.
        % dlrUsing12
        %       See bicas.proc.L1L2.demuxer_latching_relay().
        %
        %
        % RETURN VALUES
        % =============
        % RoutingArray
        %       Array of bicas.proc.L1L2.Routing objects, one per BLTS.
        %       (iBlts).
        function RoutingArray = get_routings(demuxMode, dlrUsing12)
            assert(isscalar(demuxMode))   % switch-case checks values.
            assert(isscalar(dlrUsing12))

            R = bicas.proc.L1L2.Routing.C;



            if dlrUsing12
                R.DC_V1x = R.DC_V12;
                R.AC_V1x = R.AC_V12;
            else
                R.DC_V1x = R.DC_V13;
                R.AC_V1x = R.AC_V13;
            end
            
            switch(demuxMode)
                
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
                    % contain redundant data (regardless of latching relay
                    % setting).
                    
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

                    switch(demuxMode)
                        case 5
                            RoutingArray(1) = R.REF25V_1_TO_DC_V1;
                            RoutingArray(2) = R.REF25V_2_TO_DC_V2;
                            RoutingArray(3) = R.REF25V_3_TO_DC_V3;
                        case {6,7}
                            RoutingArray(1) = R.GND_TO_DC_V1;
                            RoutingArray(2) = R.GND_TO_DC_V2;
                            RoutingArray(3) = R.GND_TO_DC_V3;
                    end
                    
                otherwise
                    % IMPLEMENTATION NOTE: switch-case statement does not work
                    % for NaN. Therefore using "otherwise".
                    if isnan(demuxMode)
                        
                        % NOTE: Could route unknown DC signals to V1-V3, but
                        % since this behaviour is probably not very obvious to
                        % the user, the code effectively deletes the information
                        % instead.
                        RoutingArray(1) = R.UNKNOWN_TO_NOWHERE;
                        RoutingArray(2) = R.UNKNOWN_TO_NOWHERE;
                        RoutingArray(3) = R.UNKNOWN_TO_NOWHERE;
                        
                        % NOTE: The routing of BLTS 4 & 5 is identical for all
                        % mux modes (but does depend on the latching relay). Can
                        % therefore route them also when the mux mode is
                        % unknown.
                        
                    else
                        error('BICAS:Assertion:IllegalArgument:DatasetFormat', ...
                            'Illegal argument value demuxMode=%g.', demuxMode)
                    end
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
        % bltsSamplesAVolt
        %       Cell array of vectors/matrices, length 5.
        %       {iBlts} = Vector/matrix with sample values
        %                 for that BLTS channel, calibrated as for ASR.
        %
        %
        % RETURN VALUES
        % =============
        % SsidArray
        %       Length-5 array of SSIDs. One SSID per BLTS.
        % AsrSamplesVolt
        %       Samples for all ASRs (singles, diffs) which can
        %       possibly be derived from the BLTS (BIAS_i). Those
        %       which can not be derived are correctly sized
        %       containing only NaN. Struct with fields.
        % --
        % NOTE: Separate names bltsSamplesAVolt & AsrSamplesAVolt to denote that
        % they are organized by BLTS and ASRs respectively.
        %
        function AsrSamplesAVoltMap = calibrated_BLTSs_to_ASRs(SdidArray, bltsSamplesAVoltCa)
            % PROPOSAL: Log message for mux=NaN.
            
            % ASSERTIONS
            assert(numel(SdidArray) == 5)
            assert(isa(SdidArray, 'bicas.proc.L1L2.SignalDestinationId'))
            assert(iscell(bltsSamplesAVoltCa))            
            assert(numel(bltsSamplesAVoltCa)==5)
            % Should ideally check for all indices, but one helps.
            assert(isnumeric(bltsSamplesAVoltCa{1}))
            
            AsrSamplesAVoltMap = bicas.proc.L1L2.demuxer.assign_ASR_samples_from_BLTS(...
                bltsSamplesAVoltCa, SdidArray);
            AsrSamplesAVoltMap = bicas.proc.L1L2.demuxer.complement_ASR(AsrSamplesAVoltMap);
        end
        
        
        
        % Given an incomplete containers.Map of ASR-labelled samples, derive the
        % missing ASRs.
        %
        % NOTE: In the event of redundant FIELDS, but not redundant DATA
        % (non-fill value), the code can NOT make intelligent choice of only
        % using available data to replace fill values.
        %   Ex: mux=1: Fields (V1, V2, V12) but V1 does not contain any data
        %   (fill values). Could in principle derive V1=V2-V12 but code does
        %   not know this.
        %   NOTE: Unlikely that this will ever happen, or that the
        %   instrument RPW is even able to return data for this situation.
        %
        % NOTE: Only public for the purpose of automatic testing.
        %
        function AsrSamplesAVoltMap = complement_ASR(AsrSamplesAVoltMap)
            assert(isa(AsrSamplesAVoltMap, 'bicas.utils.SameRowsMap'))
            
            % Shorten variable names.
            C  = bicas.proc.L1L2.AntennaSignalId.C;
            AsMap = AsrSamplesAVoltMap;
            
            %================
            % Derive AC ASRs
            %================
            % AC ASRs are separate from DC. Does not have to be in loop.
            % IMPLEMENTATION NOTE: Must be executed before DC loop. Otherwise
            % nFnAfter == 9 condition does not work.
            AsMap = bicas.proc.L1L2.demuxer.complete_relation(AsMap, C.AC_V13, C.AC_V12, C.AC_V23);

            %================
            % Derive DC ASRs
            %================
            nAsidBefore = AsMap.length;
            while true
                % NOTE: Relation DC_V13 = DC_V12 + DC_V23 has precedence for
                % deriving diffs since it is better to derive a diff from
                % (initially available) diffs rather than singles, directly or
                % indirectly, if possible.
                AsMap = bicas.proc.L1L2.demuxer.complete_relation(AsMap, C.DC_V13, C.DC_V12, C.DC_V23);
                
                AsMap = bicas.proc.L1L2.demuxer.complete_relation(AsMap, C.DC_V1,  C.DC_V12, C.DC_V2);
                AsMap = bicas.proc.L1L2.demuxer.complete_relation(AsMap, C.DC_V1,  C.DC_V13, C.DC_V3);
                AsMap = bicas.proc.L1L2.demuxer.complete_relation(AsMap, C.DC_V2,  C.DC_V23, C.DC_V3);
                nAsidAfter = AsMap.length;
                
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
            %   Ex: mux=1,2,3
            %===================================================================            
            tempNaN = nan(AsMap.nRows(), 1);

            for asidNameCa = bicas.proc.L1L2.AntennaSignalId.C.ALL_ASID_NAMES_CA'
                asidName = asidNameCa{1};
                if ~AsMap.isKey(asidName)
                    AsMap.add(asidName, tempNaN);
                end
            end
            
            AsrSamplesAVoltMap = AsMap;
        end
        
        
        
    end    % methods(Static, Access=public)
    
    
    
    %###########################################################################
    
    
    
    methods(Static, Access=private)
    

    
        % Given FIVE BLTS sample arrays, copy those which correspond to ASRs
        % (five or fewer!) into a bicas.utils.SameRowsMap.
        function AsrSamplesMap = assign_ASR_samples_from_BLTS(...
                BltsSamplesCa, SdidArray)

            % ASSERTIONS
            assert(numel(BltsSamplesCa) == 5 && iscell(BltsSamplesCa))
            assert(numel(SdidArray) == 5)

            nRows = size(BltsSamplesCa{1}, 1);
            AsrSamplesMap = bicas.utils.SameRowsMap('char', nRows, 'empty');
            for iBlts = 1:5
                if ~isequal(SdidArray(iBlts).value, 'Nowhere')
                    AsrSamplesMap.add(...
                        SdidArray(iBlts).value.s, ...
                        BltsSamplesCa{iBlts});
                end
            end
        end
        

        
        % Utility function. Derive missing ASR fields from other fields. If
        % exactly two of the Map keys exist in S, then derive the third using
        % the relationship AsMap(asid1) == AsMap(asid2) + AsMap(asid3)
        %
        % ARGUMENTS
        % =========
        % AsMap
        %       containers.Map
        % asid1, asid2, asid3
        %       ASID names key strings which may or may not be keys in AsMap. If
        %       exactly one of them is missing in "As", then the key+value is
        %       created with values assuming that the field contents are related
        %       through the relationship value1 = value2 + value3.
        %       In other cases, "AsMap" is returned unmodified.
        %
        function AsMap = complete_relation(AsMap, asid1, asid2, asid3)
            e1 = AsMap.isKey(asid1.s);
            e2 = AsMap.isKey(asid2.s);
            e3 = AsMap.isKey(asid3.s);

            if     ~e1 &&  e2 &&  e3   AsMap.add(asid1.s, AsMap.get(asid2.s) + AsMap.get(asid3.s));
            elseif  e1 && ~e2 &&  e3   AsMap.add(asid2.s, AsMap.get(asid1.s) - AsMap.get(asid3.s));
            elseif  e1 &&  e2 && ~e3   AsMap.add(asid3.s, AsMap.get(asid1.s) - AsMap.get(asid2.s));
            end
        end



    end    % methods(Static, Access=public)



end    % classdef
