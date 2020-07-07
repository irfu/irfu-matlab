%
% "Encode" the demultiplexer part of the BIAS subsystem.
% See
%   bias.demultiplexer.main.
%   bicas.BLTS_src_dest
%
%
% NOTE
% ====
% It is in principle arbitrary (probably) how the GND and "2.5V Ref" signals, which are generated by the instrument,
% should be represented in the datasets, since the datasets assume that only assumes signals from the antennas. The
% implementation classifies them as antennas, including for diffs, but the signalTypeCategory specifies that they should
% be calibrated differently. In principle, one could represent unknown signals (unknown mux mode) as antenna signals
% too.
% --
% Demultiplexer is not aware of that TDS only digitizes BLTS 1-3 (not 4-5) and does not need to be. 
%
%
% DEFINITIONS
% ===========
% See bicas.calib, bicas.demultiplexer_latching_relay.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-18
%
classdef demultiplexer
    % PROPOSAL: Change fieldname <routing>.src
    %   NOTE: Can not use ASR/AS ID. Does not cover all sources.
    %   PROPOSAL: .signalSource
    % PROPOSAL: Change fieldname <routing>.dest
    %   PROOPSAL: .datasetAsId
    
    
    
    properties(Access=private, Constant)
        
        % IMPLEMENTATION NOTE: Reasons to pre-define
        %   (1) bicas.BLTS_src_dest objects, and
        %   (2) structs created by
        % bicas.demultiplexer.routing (consisting of bicas.BLTS_src_dest objects):
        % ** Only makes constructor calls when bicas.demultiplexer is (statically) initialized, not every time the
        %    bicas.demultiplexer.main function is called. ==> Faster
        % ** Makes all necessary calls bicas.BLTS_src_dest constructor immediately, regardless of for which arguments
        %    the bicas.demultiplexer.main function is called (in particular mux mode). ==> Tests more code immediately.
        %    ==> Finds "routing bugs" faster.
        
        SRC_DC_V1   = bicas.BLTS_src_dest('DC single', [1]);
        SRC_DC_V2   = bicas.BLTS_src_dest('DC single', [2]);
        SRC_DC_V3   = bicas.BLTS_src_dest('DC single', [3]);
        SRC_DC_V12  = bicas.BLTS_src_dest('DC diff',   [1,2]);
        SRC_DC_V13  = bicas.BLTS_src_dest('DC diff',   [1,3]);
        SRC_DC_V23  = bicas.BLTS_src_dest('DC diff',   [2,3]);
        SRC_AC_V12  = bicas.BLTS_src_dest('AC diff',   [1,2]);
        SRC_AC_V13  = bicas.BLTS_src_dest('AC diff',   [1,3]);
        SRC_AC_V23  = bicas.BLTS_src_dest('AC diff',   [2,3]);

        % Define constants each one representing different routings.
        ROUTING_DC_V1  = bicas.demultiplexer.routing('DC single', [1]);
        ROUTING_DC_V2  = bicas.demultiplexer.routing('DC single', [2]);
        ROUTING_DC_V3  = bicas.demultiplexer.routing('DC single', [3]);
        ROUTING_DC_V12 = bicas.demultiplexer.routing('DC diff',   [1,2]);
        ROUTING_DC_V13 = bicas.demultiplexer.routing('DC diff',   [1,3]);
        ROUTING_DC_V23 = bicas.demultiplexer.routing('DC diff',   [2,3]);
        ROUTING_AC_V12 = bicas.demultiplexer.routing('AC diff',   [1,2]);
        ROUTING_AC_V13 = bicas.demultiplexer.routing('AC diff',   [1,3]);
        ROUTING_AC_V23 = bicas.demultiplexer.routing('AC diff',   [2,3]);

        ROUTING_25VREF_1_TO_DC_V1     = bicas.demultiplexer.routing('2.5V Ref',  [], 'DC single', [1]);
        ROUTING_25VREF_2_TO_DC_V2     = bicas.demultiplexer.routing('2.5V Ref',  [], 'DC single', [2]);
        ROUTING_25VREF_3_TO_DC_V3     = bicas.demultiplexer.routing('2.5V Ref',  [], 'DC single', [3]);
        ROUTING_GND_1_TO_V1_LF        = bicas.demultiplexer.routing('GND',       [], 'DC single', [1]);
        ROUTING_GND_2_TO_V2_LF        = bicas.demultiplexer.routing('GND',       [], 'DC single', [2]);
        ROUTING_GND_3_TO_V3_LF        = bicas.demultiplexer.routing('GND',       [], 'DC single', [3]);
        ROUTING_UNKNOWN_1_TO_NOTHING  = bicas.demultiplexer.routing('Unknown',   [], 'Nowhere', []);
        ROUTING_UNKNOWN_2_TO_NOTHING  = bicas.demultiplexer.routing('Unknown',   [], 'Nowhere', []);
        ROUTING_UNKNOWN_3_TO_NOTHING  = bicas.demultiplexer.routing('Unknown',   [], 'Nowhere', []);
%         ROUTING_UNKNOWN_1_TO_V1_LF    = bicas.demultiplexer.routing('Unknown',   [], 'DC single', [1]);
%         ROUTING_UNKNOWN_2_TO_V2_LF    = bicas.demultiplexer.routing('Unknown',   [], 'DC single', [2]);
%         ROUTING_UNKNOWN_3_TO_V3_LF    = bicas.demultiplexer.routing('Unknown',   [], 'DC single', [3]);
        % NOTE: BLTS 4 & 5 are never unknown.
    end



    methods(Static, Access=public)
        
        
        
        % NEW FUNCTION. NOT USED YET BUT MEANT TO REPLACE OLD FUNCTION "simple_demultiplex_subsequence_OLD".
        %
        %
        % Function that "encodes" the demultiplexer part of the BIAS subsystem.For a specified mux mode and demuxer
        % latching relay setting, it determines/encodes
        % (1) which (physical) input signal (Antennas, GND, "2.5V Ref", unknown) is routed to which physical output signal (BLTS)
        %       NOTE: This is needed for calibration.
        % (2) as what ASR (if any) should the BLTS be represented in the datasets,
        % (3) derive the ASRs (samples) from the ASRs which have not been set.
        %       NOTE: This derivation from fully calibrated ASR samples only requires addition/subtraction of ASRs.
        %             It does not require any sophisticated/non-trivial calibration since the relationships between the ASRs are
        %             so simple. The only consideration is that DC diffs have higher accurracy than DC singles, and should have
        %             precedence when deriving ASRs.
        %
        % NOTE: This code does NOT handle the equivalent of demultiplexer multiplication of the BLTS signal (alpha,
        % beta, gamma in the BIAS specification). It is assumed that the supplied BLTS samples have been calibrated to
        % account for this already.
        % 
        % 
        % USAGE
        % =====
        % The function is meant to be called in two different ways, typically twice for any single time period with
        % samples:
        % (1) To obtain signal type info needed for how to calibrate every BIAS-LFR/TDS signal (BIAS_i) signal given any
        %     demux mode.
        % (2) To derive the complete set of ASR samples from the given BLTS samples.
        %
        %
        % RATIONALE
        % =========
        % Meant to collect all hardcoded information about the demultiplexer routing of signals in the BIAS
        % specification, Table 4.
        %
        %
        % EDGE CASES
        % ==========
        % Function must be able to handle:
        % demuxMode = NaN                                  : Unknown demux mode, e.g. due to insufficient HK time coverage.
        % BLTS 1-3 signals labelled as "GND" or "2.5V Ref". : Demux modes 5-7.
        %
        %
        % ARGUMENTS
        % =========
        % demuxMode          : Scalar value. Demultiplexer mode.
        %                      NOTE: Can be NaN to represent unknown demux mode. Implies that AsrSamplesVolt fields are
        %                      correctly sized with NaN values.
        % dlrUsing12         : See bicas.demultiplexer_latching_relay.
        % bltsSamplesAVolt   : Cell array of vectors/matrices, length 5. {iBlts} = Vector/matrix with sample values
        %                      for that BLTS channel, calibrated as for ASR.
        % --
        % NOTE: There is no argument for diff gain since this function does not calibrate/multiply signals by factor.
        %
        %
        % RETURN VALUES
        % =============
        % BltsSrcArray   : Array of bicas.BLTS_src_dest objects. (iBlts) = Represents the origin of the corresponding
        %                  BLTS.
        % AsrSamplesVolt : Samples for all ASRs (singles, diffs) which can possibly be derived from the BLTS (BIAS_i). Those which can not be
        %                  derived are correctly sized containing only NaN. Struct with fields.
        %                  NOTE: See "EDGE CASES".
        %
        %
        % NOTE: Separate names bltsSamplesAVolt & AsrSamplesAVolt to denote that they are organized by BLTS and ASRs.
        function [BltsSrcArray, AsrSamplesAVolt] = main(demuxMode, dlrUsing12, bltsSamplesAVolt)
            % PROPOSAL: Function name that implies constant settings (MUX_SET at least; DIFF_GAIN?!).
            %   invert
            %   main
            %   demux, demultiplex
            %
            % PROPOSAL: Split into two functions.
            %   (1) Function that returns routing and calibration information
            %   (2) Function that derives the missing ASRs (complements struct)
            %       NOTE: DC and AC are separate groups of signals. ==> Simplifies.
            %   NOTE: Needs to use (e.g.) non-existance of field as indication that field has not been derived.
            %           ==> Must set only non-existent fields to NaN.
            %   CON: Might be harder than it seems.
            %       PRO: Unclear which assumptions to make without using knowledge of the mux table, in which case one
            %       does not want to split up the function in two (does not want to split up the knowledge of the mux
            %       table in two).
            %       CON: Only four relationships between DC ASRs. Only Vxy=f(Vyz,Vxz) (diff as a function of diffs) has
            %            preference over other relationships due to higher precision.
            
            % ASSERTIONS
            assert(isscalar(demuxMode))
            assert(isscalar(dlrUsing12))
            assert(iscell(bltsSamplesAVolt))            
            EJ_library.assert.vector(bltsSamplesAVolt)
            assert(numel(bltsSamplesAVolt)==5)
            assert(isnumeric(bltsSamplesAVolt{1}))   % Shuold ideally check for all indices, but one helps.
            
            % AS = "ASR Samples" (AVolt)
            NAN_VALUES = ones(size(bltsSamplesAVolt{1})) * NaN;
            As.dcV1  = NAN_VALUES;
            As.dcV2  = NAN_VALUES;
            As.dcV3  = NAN_VALUES;
            As.dcV12 = NAN_VALUES;
            As.dcV13 = NAN_VALUES;
            As.dcV23 = NAN_VALUES;
            As.acV12 = NAN_VALUES;
            As.acV13 = NAN_VALUES;
            As.acV23 = NAN_VALUES;
            
            
            
            if dlrUsing12
                ROUTING_DC_V1x = bicas.demultiplexer.ROUTING_DC_V12;
                ROUTING_AC_V1x = bicas.demultiplexer.ROUTING_AC_V12;
            else
                ROUTING_DC_V1x = bicas.demultiplexer.ROUTING_DC_V13;
                ROUTING_AC_V1x = bicas.demultiplexer.ROUTING_AC_V13;
            end

            import bicas.demultiplexer.*

            % IMPLEMENTATION NOTE: BLTS 4 & 5 are routed independently of mux mode, but the code hardcodes this
            % separately for every case (i.e. multiple times) for completeness.
            switch(demuxMode)
                
                case 0   % "Standard operation" : We have all information.
                    
                    % Summarize the routing.
                    RoutingArray(1) = bicas.demultiplexer.ROUTING_DC_V1;
                    RoutingArray(2) =                     ROUTING_DC_V1x;
                    RoutingArray(3) = bicas.demultiplexer.ROUTING_DC_V23;
                    RoutingArray(4) =                     ROUTING_AC_V1x;
                    RoutingArray(5) = bicas.demultiplexer.ROUTING_AC_V23;
                    As = assign_ASR_samples_from_BLTS(As, bltsSamplesAVolt, RoutingArray);
                    
                    % Derive the ASR:s not in the BLTS.
                    As.dcV2 = As.dcV1 - As.dcV12;
                    As.dcV3 = As.dcV2 - As.dcV23;
                    if dlrUsing12
                        As.dcV13 = As.dcV12 + As.dcV23;
                    else
                        As.dcV12 = As.dcV13 - As.dcV23;
                    end
                    
                case 1   % Probe 1 fails
                    
                    RoutingArray(1) = bicas.demultiplexer.ROUTING_DC_V2;
                    RoutingArray(2) = bicas.demultiplexer.ROUTING_DC_V3;
                    RoutingArray(3) = bicas.demultiplexer.ROUTING_DC_V23;
                    RoutingArray(4) =                     ROUTING_AC_V1x;
                    RoutingArray(5) = bicas.demultiplexer.ROUTING_AC_V23;
                    As = assign_ASR_samples_from_BLTS(As, bltsSamplesAVolt, RoutingArray);
                    
                    % NOTE: Can not derive anything extra for DC. BLTS 1-3 contain redundant data (regardless of
                    % latching relay setting).
                    
                case 2   % Probe 2 fails
                    
                    RoutingArray(1) = bicas.demultiplexer.ROUTING_DC_V1;
                    RoutingArray(2) = bicas.demultiplexer.ROUTING_DC_V3;
                    RoutingArray(3) =                     ROUTING_DC_V1x;
                    RoutingArray(4) = bicas.demultiplexer.ROUTING_AC_V13;
                    RoutingArray(5) = bicas.demultiplexer.ROUTING_AC_V23;
                    As = assign_ASR_samples_from_BLTS(As, bltsSamplesAVolt, RoutingArray);
                    
                    if dlrUsing12
                        % CASE: Know V1, V3, V12
                        As.dcV2  = As.dcV1 - As.dcV12;
                        As.dcV13 = As.dcV1 - As.dcV3;
                        As.dcV23 = As.dcV2 - As.dcV3;
                    else
                        % NOTE: Can not derive anything extra for DC. BLTS 1-3 contain redundant data.
                    end
                    
                case 3   % Probe 3 fails

                    RoutingArray(1) = bicas.demultiplexer.ROUTING_DC_V1;
                    RoutingArray(2) = bicas.demultiplexer.ROUTING_DC_V2;
                    RoutingArray(3) =                     ROUTING_DC_V1x;
                    RoutingArray(4) =                     ROUTING_AC_V1x;
                    RoutingArray(5) = bicas.demultiplexer.ROUTING_AC_V23;
                    As = assign_ASR_samples_from_BLTS(As, bltsSamplesAVolt, RoutingArray);

                    if dlrUsing12
                        % NOTE: Can not derive anything extra for DC. BLTS 1-3 contain redundant data.
                    else
                        % CASE: Know V1, V2, V13
                        As.dcV3  = As.dcV1 - As.dcV13;
                        As.dcV12 = As.dcV1 - As.dcV2;
                        As.dcV23 = As.dcV2 - As.dcV3;
                    end
                    
                case 4   % Calibration mode 0
                    
                    RoutingArray(1) = bicas.demultiplexer.ROUTING_DC_V1;
                    RoutingArray(2) = bicas.demultiplexer.ROUTING_DC_V2;
                    RoutingArray(3) = bicas.demultiplexer.ROUTING_DC_V3;
                    RoutingArray(4) =                     ROUTING_AC_V1x;
                    RoutingArray(5) = bicas.demultiplexer.ROUTING_AC_V23;
                    As = assign_ASR_samples_from_BLTS(As, bltsSamplesAVolt, RoutingArray);
                    
                    As.dcV12 = As.dcV1 - As.dcV2;
                    As.dcV13 = As.dcV1 - As.dcV3;
                    As.dcV23 = As.dcV2 - As.dcV3;

                case {5,6,7}   % Calibration mode 1/2/3

                    switch(demuxMode)
                        case 5
                            RoutingArray(1) = bicas.demultiplexer.ROUTING_25VREF_1_TO_DC_V1;
                            RoutingArray(2) = bicas.demultiplexer.ROUTING_25VREF_2_TO_DC_V2;
                            RoutingArray(3) = bicas.demultiplexer.ROUTING_25VREF_3_TO_DC_V3;
                        case {6,7}
                            RoutingArray(1) = bicas.demultiplexer.ROUTING_GND_1_TO_V1_LF;
                            RoutingArray(2) = bicas.demultiplexer.ROUTING_GND_2_TO_V2_LF;
                            RoutingArray(3) = bicas.demultiplexer.ROUTING_GND_3_TO_V3_LF;
                    end
                    RoutingArray(4) =                     ROUTING_AC_V1x;
                    RoutingArray(5) = bicas.demultiplexer.ROUTING_AC_V23;
                    As = assign_ASR_samples_from_BLTS(As, bltsSamplesAVolt, RoutingArray);
                    
                    As.dcV12 = As.dcV1 - As.dcV2;
                    As.dcV13 = As.dcV1 - As.dcV3;
                    As.dcV23 = As.dcV2 - As.dcV3;

                otherwise
                    % IMPLEMENTATION NOTE: switch-case statement does not work for NaN. Therefore using "otherwise".
                    if isnan(demuxMode)
                        
                        % NOTE: The routing of BLTS 4 & 5 is identical for all mux modes (but does depend on the
                        % latching relay). Can therefore route them also when the mux mode is unknown.
                        RoutingArray(1) = bicas.demultiplexer.ROUTING_UNKNOWN_1_TO_NOTHING;
                        RoutingArray(2) = bicas.demultiplexer.ROUTING_UNKNOWN_2_TO_NOTHING;
                        RoutingArray(3) = bicas.demultiplexer.ROUTING_UNKNOWN_3_TO_NOTHING;
                        RoutingArray(4) =                     ROUTING_AC_V1x;
                        RoutingArray(5) = bicas.demultiplexer.ROUTING_AC_V23;
                        As = assign_ASR_samples_from_BLTS(As, bltsSamplesAVolt, RoutingArray);
                        
                        % BUG/NOTE: Does not set any DC samples. Therefore just keeping the defaults. Can not derive any
                        % DC signals from other DC signals.
                        
                    else
                        error('BICAS:demultiplexer:main:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value demuxMode=%g.', demuxMode)
                    end
            end    % switch
            
            % IMPLEMENTATION NOTE: Can be placed outside switch-case since BLTS 4 & 5 are routed identically for all mux
            % modes.
            if dlrUsing12
                As.acV13 = As.acV12 + As.acV23;
            else
                As.acV12 = As.acV13 - As.acV23;
            end
            
            AsrSamplesAVolt = As;
            
            % ASSERTIONS
            assert(numel(RoutingArray) == 5)
            assert(isstruct(RoutingArray))
            
            BltsSrcArray = [RoutingArray.src];
        end
        
        
        
    end   % methods(Static, Access=public)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    methods(Static, Access=private)
    %methods(Static, Access=public)
        
        
        
        % Create struct for a given BLTS consisting of two BLTS_src_dest objects.
        %   .src  : Where the BLTS comes from. This is used to determine how the signal should be calibrated.
        %   .dest : How the BLTS should be stored in the datasets.
        %
        function R = routing(srcCategory, srcAntennas, varargin)
            % PROPOSAL: Make into class.
            
            R     = struct();
            R.src = bicas.BLTS_src_dest(srcCategory, srcAntennas);
            
            if numel(varargin) == 0
                R.dest = R.src;
            elseif numel(varargin) == 2
                R.dest = bicas.BLTS_src_dest(varargin{1}, varargin{2});
            else
                error('BICAS:demultiplexer:Assertion:IllegalArgument', '')
            end
        end
        
        
        
        % Assign ASR variables from all five BLTS, given specified routings.
        function AsrSamples = assign_ASR_samples_from_BLTS(AsrSamples, BltsSamples, RoutingArray)
            % ASSERTIONS
            EJ_library.assert.all_equal([numel(BltsSamples), numel(RoutingArray), 5])
            
            for iBlts = 1:5
                if ~strcmp(RoutingArray(iBlts).dest.category, 'Nowhere')
                    AsrFn = bicas.demultiplexer.get_ASR_fieldname( RoutingArray(iBlts).dest );
                    AsrSamples.(AsrFn) = BltsSamples{iBlts};
                end
            end
        end
        
        
        
        % Convert a BLTS_src_dest (representing an ASR) to a corresponding fieldname.
        function fn = get_ASR_fieldname(BltsDest)
            % PROPOSAL: New name implying "destination".
            
            if     isequal(BltsDest, bicas.demultiplexer.SRC_DC_V1)    fn = 'dcV1';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_DC_V2)    fn = 'dcV2';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_DC_V3)    fn = 'dcV3';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_DC_V12)   fn = 'dcV12';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_DC_V13)   fn = 'dcV13';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_DC_V23)   fn = 'dcV23';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_AC_V12)   fn = 'acV12';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_AC_V13)   fn = 'acV13';
            elseif isequal(BltsDest, bicas.demultiplexer.SRC_AC_V23)   fn = 'acV23';
            else
                error('BICAS:demultiplexer:Assertion:IllegalArgument', 'Illegal argument BltsDest.')
            end
        end
        
        
        
%         % EXPERIMENTAL. NOT USED.
%         function AsrSamplesVolt = complement_DC_signals(AsrSamplesVolt)
%             As = AsrSamplesVolt;
%             
%             % BUG??!!! Algorithm should ideally start over with highest-precedence relation after any field has been
%             % derived?
%             % NOTE: A relation can be complemented at most once, but algorithm will try to complement relations
%             % multiple times.            
%             % PROPOSAL: complete_relation returns successFlag. Use list of function pointers. Remove functions that have
%             % executed.
%             %   CON: Overkill.
%             % PROPOSAL: Add empty/NaN same-sized fields when can not derive any more fields. Assertion on fieldnames
%             %           at the end.
%             
%             % AC ASRs are separate from DC. Does not have to be in loop.
%             % IMPLEMENTATION NOTE: Must be executed before DC loop. Otherwise nFnAfter == 9 condition does not work.
%             As = complete_relation(As, 'acV13', 'acV12', 'acV23');
% 
%             nFnBefore = numel(fieldnames(As));
%             while true
%                 % NOTE: This relation has precedence since it is better to derive a diff from (initially available)
%                 % diffs rather than singles, directly or indirectly, if possible.
%                 %
%                 % Example of bad derivation:
%                 %   Mux mode=0, dlr12=true
%                 %   Initially known values (BLTS): V1,V12,V23
%                 %   V2=f(V1,V12)
%                 %   V3=f(V2,V23)
%                 %   V13=f(V1,V3)   ## Bad since could use V13=f(V12,V23), only using initial values.
%                 As = complete_relation(As, 'dcV13',    'dcV12',    'dcV23');
%                 
%                 As = complete_relation(As, 'dcV1',     'dcV12',    'dcV2');
%                 As = complete_relation(As, 'dcV2',     'dcV23',    'dcV3');
%                 As = complete_relation(As, 'dcV1',     'dcV13',    'dcV3');
%                 nFnAfter = numel(fieldnames(As));
%                 
%                 if (nFnBefore == nFnAfter) || (nFnAfter == 9)
%                     break
%                 end
%                 nFnBefore = nFnAfter;
%             end
%             
%             AsrSamplesVolt = As;
%         end
%         
%         
%         
%         % EXPERIMENTAL. NOT USED.
%         %
%         % s             : Struct
%         % fn1, fn2, fn3 : Existent or non-existent fieldnames in s.
%         %
%         % If exactly two of the fieldnames exist in s, then derive the third using the relationship s.(fn1) == s.(fn2) + s.(fn3)
%         function S = complete_relation(S, fn1, fn2, fn3)
%             e1 = isfield(S, fn1);
%             e2 = isfield(S, fn2);
%             d3 = isfield(S, fn3);
%             if     ~e1 &&  e2 &&  e3     S.(fn1) = S.(fn2) + S.(fn3);
%             elseif  e1 && ~e2 &&  e3     S.(fn2) = S.(fn1) - S.(fn3);
%             elseif  e1 &&  e2 && ~e3     S.(fn3) = S.(fn1) - S.(fn2);
%             end
%         end
                
        
        
    end    % methods(Static, Access=public)
    
    
    
end    % classdef
