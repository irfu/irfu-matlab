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
    
    
    
    properties(Access=private, Constant)
        % IMPLEMENTATION NOTE: Defining one constant struct, which contains
        % multiple constants as fields. Advantages:
        % (1) Many constants need to be defined using other constants (in this
        %     class) and MATLAB then requires that one uses the full qualifiers,
        %     i.e. bicas.proc.L1L2.demuxer.* which makes the code very long. One
        %     can also not use "import".
        % (2) Makes it possible to access constants through a variable copy of
        %     this constant rather than using the long qualifiers.
        C = bicas.proc.L1L2.demuxer.init_const();
    end



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
        %       See bicas.proc.L1L2.demuxer_latching_relay.
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
        %     BLTSs to ASRs.
        % (2) Derive the ASRs (samples) from those ASRs which have not already
        %     been set.
        %       NOTE: This derivation from fully calibrated ASR samples only
        %       requires addition/subtraction of ASRs. It does not require any
        %       sophisticated/non-trivial calibration since the relationships
        %       between the ASRs are so simple. The only consideration is that
        %       DC diffs have higher accurracy than DC singles, and should have
        %       precedence when deriving ASRs.
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
        % AsrSamplesVolt
        %       Samples for all ASRs (singles, diffs) which can
        %       possibly be derived from the BLTS (BIAS_i). Those
        %       which can not be derived are correctly sized
        %       containing only NaN. Struct with fields.
        % --
        % NOTE: Separate names bltsSamplesAVolt & AsrSamplesAVolt to denote that
        % they are organized by BLTS and ASRs respectively.
        %
        function AsrSamplesAVolt = calibrated_BLTSs_to_ASRs(RoutingArray, bltsSamplesAVoltCa)
            % PROPOSAL: Log message for mux=NaN.
            
            % ASSERTIONS
            assert(numel(RoutingArray) == 5)
            assert(isa(RoutingArray, 'bicas.proc.L1L2.Routing'))
            assert(iscell(bltsSamplesAVoltCa))            
            irf.assert.vector(bltsSamplesAVoltCa)
            assert(numel(bltsSamplesAVoltCa)==5)
            % Should ideally check for all indices, but one helps.
            assert(isnumeric(bltsSamplesAVoltCa{1}))
            
            % AS = "ASR Samples" (avolt)
            As = struct();
            As = bicas.proc.L1L2.demuxer.assign_ASR_samples_from_BLTS(...
                As, bltsSamplesAVoltCa, RoutingArray);
            As = bicas.proc.L1L2.demuxer.complement_ASR(As);

            % ASSERTIONS
            irf.assert.struct(As, bicas.proc.L1L2.demuxer.C.ASR_FIELDNAMES_CA, {})
            
            % Assign return values.
            AsrSamplesAVolt = As;
        end
        
        
        
        % Given an incomplete ASR struct, derive the missing ASRs using the
        % existing ones.
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
        function AsrSamplesAVolt = complement_ASR(AsrSamplesAVolt)
            % Shorten variable names.
            C  = bicas.proc.L1L2.demuxer.C;
            As = AsrSamplesAVolt;   
            
            % NOTE: A relation can be complemented at most once, but the
            % algorithm will try to complement relations multiple times.
            %
            % PROPOSAL: Add empty/NaN same-sized fields when can not derive any more fields. Assertion on fieldnames
            %           at the end.
            %
            % PROPOSAL: Implement on a record-by-record basis, with b indices to
            %           select which records should be derived.
            %   PRO: Can apply over multiple demultiplexer modes.
            %       CON: Not necessarily for non-antenna data. Might want to
            %       keep V12=NaN, but V1,V2=non-NaN for e.g. 2.5V ref data.
            %       Maybe...
            %   Ex: b1 = ~isnan(As.dcV1), ...
            %       b = b1 & b2 & ~b12;
            %       As.dcV12(b, :) = As.dcV1(b, :) - As.dcV2(b, :);
            
            N_ASR_FIELDNAMES = numel(C.ASR_FIELDNAMES_CA);
            
            %================
            % Derive AC ASRs
            %================
            % AC ASRs are separate from DC. Does not have to be in loop.
            % IMPLEMENTATION NOTE: Must be executed before DC loop. Otherwise
            % nFnAfter == 9 condition does not work.
            As = bicas.proc.L1L2.demuxer.complete_relation(As, 'acV13', 'acV12', 'acV23');

            %================
            % Derive DC ASRs
            %================
            nFnBefore = numel(fieldnames(As));
            while true
                % NOTE: Relation dcV13 = dcV12 + dcV23 has precedence for
                % deriving diffs since it is better to derive a diff from
                % (initially available) diffs rather than singles, directly or
                % indirectly, if possible.
                As = bicas.proc.L1L2.demuxer.complete_relation(As, 'dcV13', 'dcV12', 'dcV23');
                
                As = bicas.proc.L1L2.demuxer.complete_relation(As, 'dcV1',  'dcV12', 'dcV2');
                As = bicas.proc.L1L2.demuxer.complete_relation(As, 'dcV1',  'dcV13', 'dcV3');
                As = bicas.proc.L1L2.demuxer.complete_relation(As, 'dcV2',  'dcV23', 'dcV3');
                nFnAfter = numel(fieldnames(As));
                
                if (nFnBefore == nFnAfter) || (nFnAfter == 9)
                    break
                end
                nFnBefore = nFnAfter;
            end
                        
            %===================================================================
            % Assign all ASRs/fields which have not yet been assigned
            % -------------------------------------------------------
            % IMPLEMENTATION NOTE: This is needed to handle for situations when
            % the supplied fields can not be used to determine all fields (five
            % fields are supplied but the system of equations is
            % overdetermined/redundant).
            %   Ex: mux=1,2,3
            %===================================================================            
            fnCa    = fieldnames(As);
            tempNaN = nan(size(As.(fnCa{1})));
            for iFn = 1:N_ASR_FIELDNAMES
                fn = C.ASR_FIELDNAMES_CA{iFn};
                if ~isfield(As, fn)
                    As.(fn) = tempNaN;
                end
            end
            
            irf.assert.struct(As, C.ASR_FIELDNAMES_CA, {})
            
            AsrSamplesAVolt = As;
        end
        
        
        
    end    % methods(Static, Access=public)
    
    
    
    %###########################################################################
    
    
    
    methods(Static, Access=private)
    %methods(Static, Access=public)
    
    
        function R = init_const()
            SDID = bicas.proc.L1L2.SignalDestinationId.C;
            % Table that associates SDIDs with ASR struct fieldnamnes (FN).
            R.SDID_ASR_FN_TABLE = {...
                SDID.DC_V1,  'dcV1'; ...
                SDID.DC_V2,  'dcV2'; ...
                SDID.DC_V3,  'dcV3'; ...
                ...
                SDID.DC_V12, 'dcV12'; ...
                SDID.DC_V13, 'dcV13'; ...
                SDID.DC_V23, 'dcV23'; ...
                ...
                SDID.AC_V12, 'acV12'; ...
                SDID.AC_V13, 'acV13'; ...
                SDID.AC_V23, 'acV23'; ...
            };

            R.ASR_FIELDNAMES_CA = R.SDID_ASR_FN_TABLE(:, 2);
        end
    
    

        % Overwrite ASR struct fields from all FIVE BLTS, given specified
        % routings. Does not touch other struct fields.
        function AsrSamples = assign_ASR_samples_from_BLTS(...
                AsrSamples, BltsSamplesAVolt, RoutingArray)
            
            % ASSERTIONS
            irf.assert.all_equal(...
                [numel(BltsSamplesAVolt), numel(RoutingArray), 5])
            
            for iBlts = 1:5
                if ~isequal(RoutingArray(iBlts).dest.value, 'Nowhere')
                    asrFn = bicas.proc.L1L2.demuxer.get_ASR_fieldname(...
                        RoutingArray(iBlts).dest);
                    AsrSamples.(asrFn) = BltsSamplesAVolt{iBlts};
                end
            end
        end
        
        
        
        % Convert a bicas.proc.L1L2.SignalDestinationId object (nominally
        % representing an ASR) to a corresponding struct fieldname.
        function fn = get_ASR_fieldname(Sdid)
            % PROPOSAL: New name implying "destination".
            
            SDID_ASR_FN_TABLE = bicas.proc.L1L2.demuxer.C.SDID_ASR_FN_TABLE;
            
            for i =1:size(SDID_ASR_FN_TABLE, 1)
                if isequal(Sdid, SDID_ASR_FN_TABLE{i, 1})
                    fn = SDID_ASR_FN_TABLE{i, 2};
                    return
                end
            end
            error('BICAS:Assertion:IllegalArgument', ...
                'Illegal argument BltsDest.')
        end
        
        
        
        % Utility function. Derive missing ASR fields from other fields. If
        % exactly two of the fieldnames exist in S, then derive the third using
        % the relationship As.(fn1) == As.(fn2) + As.(fn3)
        %
        % ARGUMENTS
        % =========
        % As            
        %       Struct.
        % fn1, fn2, fn3
        %       Existent or non-existent fieldnames in s. If exactly one of them
        %       is missing in "As", then the field is created with values
        %       assuming that the field contents are related through the
        %       relationship fn1 = fn2 + fn3.
        %       In other cases, "As" is returned unmodified.
        %
        function As = complete_relation(As, fn1, fn2, fn3)
            e1 = isfield(As, fn1);
            e2 = isfield(As, fn2);
            e3 = isfield(As, fn3);
            if     ~e1 &&  e2 &&  e3     As.(fn1) = As.(fn2) + As.(fn3);
            elseif  e1 && ~e2 &&  e3     As.(fn2) = As.(fn1) - As.(fn3);
            elseif  e1 &&  e2 && ~e3     As.(fn3) = As.(fn1) - As.(fn2);
            end
        end



    end    % methods(Static, Access=public)



end    % classdef
