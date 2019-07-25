classdef calibration
% Class for (1) loading calibration data from file, and (2) library/utility functions that calibrate data.
%
% UNFINISHED. NOT USED BY MAIN PROGRAM YET. NOT ENTIRELY CLEAR WHAT IT SHOULD INCLUDE.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-15



% BOGIQ:
% ------
% PROPOSAL: Need multiple functions
%   Apply TF in frequency domain, in time domain(?)   (in both directions = inverse & non-inverse)
%   Convert TF from frequency domain to time domain (and maybe reverse for debugging)
%   Invert TF function.
%   Multiply/combine transfer functions.
%   Convert files data to time-domain transfer functions?
%   (Select one BIAS TF of many? Interpolate between transfer functions?!!!)
%   Select a TDS/LFR TF
%   Combine TFs: TDS/LFR, BIAS, possibly capacitance-TF.
%   Apply TF (inverse)
%   Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
%   Special treatment of 0 Hz?!
%
% QUESTION: How implement the parasitic capacitance TF? Use analytical TF directly, or convert to TF table?
% TODO: Read RCT
%
% PROPOSAL: apply_transfer_function_in_freq as separate utility function?

    properties(Access=private)
        
        % Different TFs. Initialized in constructor.
        % IMPLEMENTATION NOTE: Variable names include the prefix "BIAS" to indicate that these are the TFs for BIAS
        % as opposed to TDS and LFR which will have to be added later.
        %biasDcAbsGain        % BIAS gain for DC absolute (non-diff) data (instead of TF).
        %biasDcDiffGain       % BIAS gain for DC diff data (instead of TF).
        %biasDcAbsOffsets     % Unclear definition so far.
        
        biasAcLoGainTf
        biasAcHiGainTf
        
        BiasDcRecordList = [];
    end

    %###################################################################################################################

    methods(Access=public)

        function obj = calibration()
        
            % QUESTION: Is it wise to specify the paths in the constructor? Read the filenames (and relative directory) from the constants instead?
            
            %obj.biasDcAbsGain  = 
            %obj.biasDcDiffGain = 
            
            obj.BiasDcRecordList    = init_biasDcRecordList();

        end


        
%         % TEST
%         function inputSignalVolt = calibrate_DC(obj, temperatureCelsius, stimuli, antennaCh, lfrCh, outputSignalVolt)
%             %if numel(antennaCh) == 2
%             %    antennaCh = sort(antennaCh);
%             %end
%             
%             if     isequal(antennaCh, [1])
%             elseif isequal(antennaCh, [2])
%             elseif isequal(antennaCh, [3])
%             elseif isequal(antennaCh, [1,2])
%             elseif isequal(antennaCh, [1,3])
%             elseif isequal(antennaCh, [1,3])
%             else
%                 error('BICAS:Assertion', 'Can not interpret antenna channels.')
%             end
%             
%             inputSignalVolt = offset + slope .* outputSignalVolt;
%         end
%         
% 
    end    % methods(Access=public)
    
    %###################################################################################################################

    %methods(Static, Access=private)
    methods(Static, Access=public)   % Temporarily public for testing purposes.
        
%         function tfZ = parasitic_capacitance_TF(tfOmega)
%             % Calucalte Z(omega) values for TF representing parasitic capacitances (based on analytic function).
% 
%             % Function name? "Input capacitance"?
%             % Not read R & C from constants here? Submit as arguments?
%             capacitanceFarad =
%             impedanceOhm     =
%             
%             % Correct for a TF?!
%             tfZ = 1 / (1 + 1j*tfOmega*capacitanceFarad*impedanceOhm);
%             
%             error('BICAS:calibration:OperationNotImplemented', 'Function not implemented Yet.')
%         end



%         function biasCurrentAmpere = HK_TM_bias_to_bias_current(hkBiasTm)
%         % Convert HK diagnostic TM bias current values to physical units.
%         %
%         % NOTE: The HK bias current values are the measured onboard analog bias current values but are only meant as
%         % diagnostic values, not the proper bias current values. Therefore the calibrated values should only be seen as
%         % approximate.
%         % NOTE: The conversion function can be found in the BIAS spec, under BIAS1/2/3 under "Telemetry". (Not to be confused with the corresponding telecommands.) 
%         % NOTE: The conversion functions are identical for all three probes.
%
%         % NOTE: 
%         end

    end

    %###################################################################################################################

    methods(Static, Access=public)




        
    end    % methods(Static, Access=public)

end
