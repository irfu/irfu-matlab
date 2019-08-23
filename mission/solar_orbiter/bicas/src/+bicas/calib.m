classdef calib
% Class for (1) loading calibration data from file, and (2) library/utility functions that calibrate data.
%
% NOTE: RCT reading functions assume that the same type of RCT (BIAS, LFR, TDS-CWF or TDS-RSWF) is identical (in all
% relevant parts) for RODP and ROC-SGSE.
%
%
% DEFINITIONS
% ===========
% CPV = counts/volt
% VPC = volt/count
%
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
% PROPOSAL: Change name.
%   PRO: "calibration" is too long.
%   PRO: "calibration" overlaps with function names.
%   PROPOSAL: "calib", "cal".



%     properties(Access=private)
%         
%         % Different TFs. Initialized in constructor.
%         % IMPLEMENTATION NOTE: Variable names include the prefix "BIAS" to indicate that these are the TFs for BIAS
%         % as opposed to TDS and LFR which will have to be added later.
%         %biasDcAbsGain        % BIAS gain for DC absolute (non-diff) data (instead of TF).
%         %biasDcDiffGain       % BIAS gain for DC diff data (instead of TF).
%         %biasDcAbsOffsets     % Unclear definition so far.
%         
%         biasAcLoGainTf
%         biasAcHiGainTf
%         
%         BiasDcRecordList = [];
%     end

    %###################################################################################################################

%     methods(Access=public)
% 
%         function obj = calib()
%         
%             % QUESTION: Is it wise to specify the paths in the constructor? Read the filenames (and relative directory) from the constants instead?
%             
%             %obj.biasDcAbsGain  = 
%             %obj.biasDcDiffGain = 
%             
%             obj.BiasDcRecordList    = init_biasDcRecordList();
% 
%         end
% 
% 
%         
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
%     end    % methods(Access=public)
    
    %###################################################################################################################

    methods(Static, Access=public)
        

        
        % Find the path to the RCT to use, using the filenames specified in the documentation. If there are multiple matching candidates, choose the latest one as indicated
        % by the filename.
        %
        % RCT filenaming convention is described in:	ROC-PRO-DAT-NTT-00006-LES, 1/1 draft, Sect 4.3.2-3.
        %
        % IMPLEMENTATION NOTE: BICAS is not following the documented filenaming convention here since (1) TDS and
        % LFR teams do not seem to use it and (2) it is uncertain how it can be applied to BIAS RCTs (which receiver
        % should the BIAS RCT specify?). Therefore using a regular expression that covers more files.
        %     NOTE: LFR RCTs use 2+12 digits instead of 10 in the timestamps (they add seconds=2 digits).
        %     NOTE: TDS RCTs use 2+6  digits instead of 10 in the timestamps (the have no time of day, only date)
        %
        % IMPLEMENTATION NOTE: Useful to have this as separate functionality so that the chosen RCT to use can be
        % explicitly overridden via e.g. settings.
        %
        function path = find_RCT(calibrationDir, pipelineId, rctId)            
            % PROPOSAL: Different regexp (not just the middle substring) for different types of RCT to handle different
            % naming conventions (e.g. different number of digits, whether to include substrings RCT, RPW).
            % PROPOSAL: Put type-dependent regexp. naming conventions in settings.
            %   PROPOSAL: Exclude pipelineId prefix and CDF suffix.
            
            % Examples of RCT filenames
            % -------------------------
            % BIAS:
            %       ROC-SGSE_CAL_RCT-BIAS_V201803211625.cdf   (old implemented convention)
            %       ROC-SGSE_CAL_RPW_BIAS_V201908231028.cdf   (new implemented convention, closer to documentation)
            %           SOLO_CAL_RCT-BIAS_V201901141146.cdf   (old implemented convention)
            % LFR:
            %       ROC-SGSE_CAL_RCT-LFR-BIAS_V20180724165443.cdf
            %           SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf
            % TDS:
            %           SOLO_CAL_RCT-TDS-LFM-CWF-E_V20190128.cdf
            %           SOLO_CAL_RCT-TDS-LFM-RSWF-E_V20190128.cdf
            %         (Two types of calibration files, but only RODP versions)

            %============================
            % Create regexp for filename
            %============================
            switch(pipelineId)
                case {'ROC-SGSE', 'RGTS'}
                    filenamePrefix = 'ROC-SGSE';
                case 'RODP'
                    filenamePrefix = 'SOLO';
                otherwise
                    error('calib:Assertion:IllegalArgument', 'Illegal pipelineId="%s"', pipelineId)
            end
            clear pipelineId
            
            % IMPLEMENTATION NOTE: This switch statement
            % (1) verifies the input/output, AND
            % (2) concentrates BICAS' knowledge of RCT filenaming to this function (the caller does not need to know
            %     anything about filenames).
            switch(rctId)
                case 'BIAS'
                    %rctDependentSubRegexp = 'RCT-BIAS_V20[0-9]{10}';   % Old convention.
                    rctDependentSubRegexp = 'RPW_BIAS_V20[0-9]{10}';
                case 'LFR'
                    rctDependentSubRegexp = 'RCT-LFR-BIAS_V20[0-9]{12}';
                case 'TDS-CWF'
                    rctDependentSubRegexp = 'RCT-TDS-LFM-CWF-E_V20[0-9]{6}';
                case 'TDS-RSWF'
                    rctDependentSubRegexp = 'RCT-TDS-LFM-RSWF-E_V20[0-9]{6}';
                otherwise
                    error('calib:Assertion:IllegalArgument', 'Illegal rctId="%s"', rctId)
            end
            clear rctId
            
            %filenameRegexp = sprintf('%s_CAL_[A-Z_-]*%s_V20[0-9]{6,12}.(cdf|CDF)', filenamePrefix, rctDependentSubRegexp);
            filenameRegexp = sprintf('%s_CAL_%s.(cdf|CDF)', filenamePrefix, rctDependentSubRegexp);



            %=================================================
            % Find candidate files and select the correct one
            %=================================================
            dirObjectList = dir(calibrationDir);
            dirObjectList([dirObjectList.isdir]) = [];    % Eliminate directories.
            filenameList = {dirObjectList.name};
            filenameList(~EJ_library.utils.regexpf(filenameList, filenameRegexp)) = [];    % Eliminate non-matching filenames.
            
            % ASSERTION / WARNING
            if numel(filenameList) == 0
                % ERROR
                error('Can not find any calibration file that matches regular expression "%s" in directory "%s".', filenameRegexp, calibrationDir);
            end
            % CASE: There is at least on candidate file.
            
            filenameList = sort(filenameList);
            filename     = filenameList{end};
            path = fullfile(calibrationDir, filename);
            
            if numel(filenameList) > 1
                % WARNING/INFO/NOTICE
                msg = sprintf('Found multiple calibration files matching regular expression "%s"\nin directory "%s".\nSelecting the latest one as indicated by the filename: "%s".\n', ...
                    filenameRegexp, calibrationDir, filename);
%                 for i = 1:numel(filenameList)
%                     msg = [msg, sprintf('    %s\n', filenameList{i})];
%                 end
                bicas.log('warning', msg)
            end
            
            % IMPLEMENTATION NOTE: Not logging which calibration file is selected, since this function is not supposed
            % to actually load the content.
        end



        % BOGIQ: RCT-reading functions
        % ============================
        % PROPOSAL: Use same code/function for reading calibration table, as for reading dataset (and master cdfs)?
        % PROPOSAL: Assert CDF skeleton/master version number.
        % PROPOSAL: Assert skeleton/master.
        % PROPOSAL: Convert from file units to "mathematical units"? (radians/s, radians (for phase shift))
        % PROPOSAL: Assert/warn (depending on setting?) file units.
        % PROPOSAL: Only use units in variable names.
        % PROPOSAL: Use utility function for reading every zVariable.
        %   PROPOSAL: Assert units from zVar attributes.
            
        function BiasCalibData = read_BIAS_RCT(filePath)
            
            Do = dataobj(filePath);
            
            Bcd.Epoch_L                  = Do.data.Epoch_L.data;
            Bcd.Epoch_H                  = Do.data.Epoch_H.data;
            Bcd.BIAS_CURRENT_OFFSET      = Do.data.BIAS_CURRENT_OFFSET.data;
            Bcd.BIAS_CURRENT_GAIN        = Do.data.BIAS_CURRENT_GAIN.data;
            Bcd.V_OFFSET                 = Do.data.V_OFFSET.data;
            Bcd.E_OFFSET                 = Do.data.E_OFFSET.data;
            Bcd.TRANSFER_FUNCTION_COEFFS = Do.data.TRANSFER_FUNCTION_COEFFS.data;
            
            BiasCalibData = Bcd;
        end
        
        
        
        function LfrCalibData = read_LFR_RCT(filePath)
            
            Do = dataobj(filePath);
            
            Lcd.Freqs_F0_Hz                 = Do.data.Freqs_F0.data;
            Lcd.Freqs_F1_Hz                 = Do.data.Freqs_F1.data;
            Lcd.Freqs_F2_Hz                 = Do.data.Freqs_F2.data;
            Lcd.Freqs_F3_Hz                 = Do.data.Freqs_F3.data;
            
            Lcd.BIAS_12345_amplitude_F0_Cpv = Do.data.TF_BIAS_12345_amplitude_F0.data;
            Lcd.BIAS_12345_amplitude_F1_Cpv = Do.data.TF_BIAS_12345_amplitude_F1.data;
            Lcd.BIAS_12345_amplitude_F2_Cpv = Do.data.TF_BIAS_12345_amplitude_F2.data;
            Lcd.BIAS_123_amplitude_F3_Cpv   = Do.data.TF_BIAS_123_amplitude_F3.data;
            
            Lcd.BIAS_12345_phase_F0_Deg     = Do.data.TF_BIAS_12345_phase_F0.data;
            Lcd.BIAS_12345_phase_F1_Deg     = Do.data.TF_BIAS_12345_phase_F1.data;
            Lcd.BIAS_12345_phase_F2_Deg     = Do.data.TF_BIAS_12345_phase_F2.data;
            Lcd.BIAS_123_phase_F3_Deg       = Do.data.TF_BIAS_123_phase_F3.data;
            
            LfrCalibData = Lcd;
        end
        
        
        
        function TdsCwfCalibData = read_TDS_CWF_RCT(filePath)
            
            Do = dataobj(filePath);
            
            % NOTE: Undocumented in CDF: zVar CALIBRATION_TABLE is volt/count for just multiplying the TDS signal.
            TdsCwfCalibData.calibrationFactorsVpc = Do.data.CALIBRATION_TABLE.data;
            %TdsCwfCalibData.CHANNEL_LABEL = Do.data.CHANNEL_LABEL.data;     % Unnecessary.
        end
        
        
        
        function TdsRswfCalibData = read_TDS_RSWF_RCT(filePath)
            Do = dataobj(filePath);
            
            Trcd.CALIBRATION_FREQUENCY_Hz  = Do.data.CALIBRATION_FREQUENCY.data;
            Trcd.CALIBRATION_AMPLITUDE_Vpc = Do.data.CALIBRATION_AMPLITUDE.data;
            Trcd.CALIBRATION_PHASE_Deg     = Do.data.CALIBRATION_PHASE.data;
            %Trcd.CHANNEL_ID_LABEL = Do.data.CHANNEL_ID_LABEL.data;    % Unnecessary
            
            TdsRswfCalibData = Trcd;
        end
        
        
        
%         function tfZ = parasitic_capacitance_TF(tfOmega)
%             % Calculate Z(omega) values for TF representing parasitic capacitances (based on analytic function).
% 
%             % Function name? "Input capacitance"?
%             % Not read R & C from constants here? Submit as arguments?
%             capacitanceFarad =
%             impedanceOhm     =
%             
%             % Correct for a TF?!
%             tfZ = 1 / (1 + 1j*tfOmega*capacitanceFarad*impedanceOhm);
%             
%             error('BICAS:calib:OperationNotImplemented', 'Function not implemented Yet.')
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

    end    % methods(Static, Access=public)

end
