classdef demultiplexer

        methods(Static, Access=public)

        
        
        % NEW FUNCTION. NOT USED YET BUT MEANT TO REPLACE OLD FUNCTION "simple_demultiplex_subsequence_OLD".
        %
        % (1) Returns the demultiplexer mode-dependent information needed for calibrating a BLTS interface signal (BIAS_i),
        % (2) Derives as much as possible of all the antenna singles and diffs from the available BIAS-LFR/TDS signals
        % (BIAS_i), except the calibration (i.e. only addition and subtraction).
        %
        % Meant to be called in two different ways, typically twice for any time period with samples.
        % (1) To obtain signal type info needed for how to calibrate every BIAS-LFR/TDS signal (BIAS_i) signal given any demux mode. 
        % (2) To derive the complete set of ASR samples from the given BLTS samples.
        %
        % RATIONALE: Meant to collect all hardcoded information about the demultiplexer routing of signals.
        % NOTE: Does not perform any calibration. The closest is to calculate diffs and singles from diffs and singles.
        % 
        % 
        % ARGUMENTS
        % =========
        % demuxMode            : Scalar value. Demultiplexer mode.
        % dlrUsing12           : 0/1, true/false. DLR = Demultiplexer Latching Relay.
        %                         False=0 = Using diffs V13_DC, V13_AC
        %                         True =1 = Using diffs V12_DC, V12_AC
        % BltsSamplesCalibVolt : Cell array of matrices, length 5. {iBlts} = Vector with sample values for that channel.
        %                        BIAS calibrated volts.
        % --
        % NOTE: No argument for diff gain since this function does not calibrate.
        %
        %
        % RETURN VALUES
        % =============
        % BltsAsrType : Struct array. (iBlts) = Number representing ASR type of the BLTS data, which depends on the mux mode.
        %               Has fields
        %                   .antennas = Numeric vector of length 0, 1 or 2.
        %                           Either [] (no signal, e.g. BIAS_4/5 for TDS), [iAnt] (single), or [iAnt1, iAnt2] (diff).
        %                           NOTE: iAnt1 < iAnt2. iAnt/iAnt1/iAnt2 = {1,2,3}.
        %                           Represents the current routing of signals.
        %                   .category = String constant representing the category/type of signal on the channel.
        %                           DC single, DC diff, AC low-gain, AC high-gain, no signal
        % AsrSamplesVolt
        %             : All representations of antenna signals which can possibly be derived from the BLTS (BIAS_i).
        %               Struct with fields named as in the BIAS specification: .Vi_LF, .Vij_LF, .Vij_LF_AC
        %               NOTE: Calibration signals GND and 2.5V Ref which are generated internally by BIAS are also
        %               stored in these variables although they are technically not ASRs. See implementation.
        %
        %
        % DEFINITIONS, NAMING CONVENTIONS
        % ===============================
        % See bicas.calib.
        %
        function [BltsAsrType, AsrSamplesVolt] = main(demuxMode, dlrUsing12, BltsSamplesCalibVolt)
            % PROPOSAL: Function name that implies constant settings (MUX_SET at least; DIFF_GAIN?!).
            %   invert
            %   main
            %   demux, demultiplex
            %   
            % PROPOSAL/NOTE: BIAS calibrated volts = ASR volts (automatically for those ASR for which there is BLTS data)
            % TODO-DECISION: How handle calibration modes with fixed, constant BIAS-LFR/TDS signals ("GND", "2.5 V Ref")?
            %
            % PROBLEM: BltsAsrType.category for AC can not include low-gain/high-gain which leads to different set of
            % alternatives than used for selecting transfer functions.
            % PROPOSAL: "Assertion" for using good combination of mux mode and latching relay. Log warning if assertion
            %           fails.
            % PROPOSAL: Assertions for returned string constants (one of listed legal alternatives).
            
            % ASSERTIONS
            assert(isscalar(demuxMode))
            assert(isscalar(dlrUsing12))
            assert(iscell(BltsSamplesCalibVolt))
            EJ_library.utils.assert.vector(BltsSamplesCalibVolt)
            assert(numel(BltsSamplesCalibVolt)==5)
            
            % CV = (BIAS) Calibrated (BLTS) volt
            BIAS_1_Cv = BltsSamplesCalibVolt{1};
            BIAS_2_Cv = BltsSamplesCalibVolt{2};
            BIAS_3_Cv = BltsSamplesCalibVolt{3};
            BIAS_4_Cv = BltsSamplesCalibVolt{4};
            BIAS_5_Cv = BltsSamplesCalibVolt{5};
            
            % AS = "ASR Samples"
            NAN_VALUES = ones(size(BIAS_1_Cv)) * NaN;
            As.V1_LF     = NAN_VALUES;
            As.V2_LF     = NAN_VALUES;
            As.V3_LF     = NAN_VALUES;
            As.V12_LF    = NAN_VALUES;
            As.V13_LF    = NAN_VALUES;
            As.V23_LF    = NAN_VALUES;
            As.V12_LF_AC = NAN_VALUES;
            As.V13_LF_AC = NAN_VALUES;
            As.V23_LF_AC = NAN_VALUES;



            if dlrUsing12;   iAntB = 2;
            else         ;   iAntB = 3;
            end

            import bicas.demultiplexer.routing
            
            % NOTE: BLTS 5 = V23_LF_AC for all modes, but the code hardcodes this separately for every case for
            % completeness.
            switch(demuxMode)
                case 0   % "Standard operation" : We have all information.

                    % Summarize the routing.
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [1,iAntB], 'DC diff',   BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [2,3],     'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % Derive the ASR:s not in the BLTS.
                    if dlrUsing12
                        As.V13_LF    = As.V12_LF    + As.V23_LF;
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF    = As.V13_LF    - As.V23_LF;
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end
                    As.V2_LF     = As.V1_LF     - As.V12_LF;
                    As.V3_LF     = As.V2_LF     - As.V23_LF;
                    
                case 1   % Probe 1 fails

                    [BltsAsrType(1), As] = routing(As, [2],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [3],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [2,3],     'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % NOTE: Can not derive anything extra for DC. BLTS 1-3 contain redundant data.
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end
                    
                case 2   % Probe 2 fails
                    
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [3],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [1,iAntB], 'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % NOTE: Can not derive anything extra for DC. BLTS 1-3 contain redundant data.
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end
                    
                case 3   % Probe 3 fails
                    
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [2],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [1,iAntB], 'DC diff',   BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    % NOTE: Can not derive anything extra for DC. BLTS 1-3 contain redundant data.
                    if dlrUsing12
                        As.V13_LF_AC = V12_LF_AC + V23_LF_AC;
                    else
                        As.V12_LF_AC = V13_LF_AC - V23_LF_AC;
                    end
                    
                case 4   % Calibration mode 0
                    
                    [BltsAsrType(1), As] = routing(As, [1],       'DC single', BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [2],       'DC single', BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [3],       'DC single', BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',        BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',        BIAS_5_Cv);
                    
                    As.V12_LF    = As.V1_LF    - As.V2_LF;
                    As.V13_LF    = As.V1_LF    - As.V3_LF;
                    As.V23_LF    = As.V2_LF    - As.V3_LF;
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end

                case {5,6,7}   % Calibration mode 1/2/3
                    
                    switch(demuxMode)
                        case 5
                            signalTypeCategory = '2.5V Ref';
                        case {6,7}
                            signalTypeCategory = 'GND';
                    end
                    
                    % NOTE: It is in principle arbitrary (probably) how the GND and 2.5V Ref signals, which are
                    % generated by the instrument, should be represented in the datasets, since the datasets assume that
                    % only assumes signals from the antennas. The implementation classifies them as antennas, including
                    % for diffs, but the signalTypeCategory specifies that they should be calibrated differently.
                    [BltsAsrType(1), As] = routing(As, [1],       signalTypeCategory, BIAS_1_Cv);
                    [BltsAsrType(2), As] = routing(As, [2],       signalTypeCategory, BIAS_2_Cv);
                    [BltsAsrType(3), As] = routing(As, [3],       signalTypeCategory, BIAS_3_Cv);
                    [BltsAsrType(4), As] = routing(As, [1,iAntB], 'AC',               BIAS_4_Cv);
                    [BltsAsrType(5), As] = routing(As, [2,3],     'AC',               BIAS_5_Cv);

                    As.V12_LF    = As.V1_LF    - As.V2_LF;
                    As.V13_LF    = As.V1_LF    - As.V3_LF;
                    As.V23_LF    = As.V2_LF    - As.V3_LF;
                    if dlrUsing12
                        As.V13_LF_AC = As.V12_LF_AC + As.V23_LF_AC;
                    else
                        As.V12_LF_AC = As.V13_LF_AC - As.V23_LF_AC;
                    end

                otherwise
                    error('BICAS:proc_sub:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value for mux_set.')
            end    % switch
            
            AsrSamplesVolt = As;
            
            assert(numel(BltsAsrType) == 5)
        end
        
        
        
    end   % methods(Static, Access=public)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    methods(Static, Access=private)



        % Utility function for "demultiplexer".
        % (1) Assign BLTS samples to the specified ASR variable, but only if they should be assigned (not/ambiguous "2.5V Ref", "GND").
        % (2) Return BLTS ASR type on standard format for the signal.
        %
        % NOTE: AsrSamplesSet (struct) = All ASRs, of which one (or none) is assigned.
        %
        function [BltsAsrType, AsrSamplesSet] = routing(AsrSamplesSet, antennas, category, BltsSamples)
            % PROPOSAL: Use "translate" function.
            %   CON: Not designed to replace this type of code: Multiple variables to compare & non-string values to compare.
            
            % Normalize vector to row vector since "isequal" is sensitive to row/column vectors.
            antennas = antennas(:)';
            
            % Assign BltsType.
            BltsAsrType.antennas = antennas;
            BltsAsrType.category = category;
            
            % Modify AsrSamples (and assertion on arguments).
            if     isequal(antennas, [1])   && strcmp(category, 'DC single')   AsrSamplesSet.V1_LF     = BltsSamples;
            elseif isequal(antennas, [2])   && strcmp(category, 'DC single')   AsrSamplesSet.V2_LF     = BltsSamples;
            elseif isequal(antennas, [3])   && strcmp(category, 'DC single')   AsrSamplesSet.V3_LF     = BltsSamples;
            elseif isequal(antennas, [1,2]) && strcmp(category, 'DC diff')     AsrSamplesSet.V12_LF    = BltsSamples;
            elseif isequal(antennas, [1,3]) && strcmp(category, 'DC diff')     AsrSamplesSet.V13_LF    = BltsSamples;
            elseif isequal(antennas, [2,3]) && strcmp(category, 'DC diff')     AsrSamplesSet.V23_LF    = BltsSamples;
            elseif isequal(antennas, [1,2]) && strcmp(category, 'AC')          AsrSamplesSet.V12_LF_AC = BltsSamples;
            elseif isequal(antennas, [1,3]) && strcmp(category, 'AC')          AsrSamplesSet.V13_LF_AC = BltsSamples;
            elseif isequal(antennas, [2,3]) && strcmp(category, 'AC')          AsrSamplesSet.V23_LF_AC = BltsSamples;
            elseif isequal(antennas, [1])   && strcmp(category, '2.5V Ref')    % Do nothing. AsrSamples.V1_LF = BltsSamples;  ?
            elseif isequal(antennas, [2])   && strcmp(category, '2.5V Ref')    % Do nothing. AsrSamples.V2_LF = BltsSamples;  ?
            elseif isequal(antennas, [3])   && strcmp(category, '2.5V Ref')    % Do nothing. AsrSamples.V3_LF = BltsSamples;  ?
            elseif isequal(antennas, [1])   && strcmp(category, 'GND')         % Do nothing. AsrSamples.V1_LF = BltsSamples;  ?
            elseif isequal(antennas, [2])   && strcmp(category, 'GND')         % Do nothing. AsrSamples.V2_LF = BltsSamples;  ?
            elseif isequal(antennas, [3])   && strcmp(category, 'GND')         % Do nothing. AsrSamples.V3_LF = BltsSamples;  ?
            else
                error('BICAS:proc_SUB:Assertion:IllegalArgument', 'Illegal combination of arguments "antennas" and "category".')
            end
        end
        
        
        
    end    % methods(Static, Access=private)
        
end