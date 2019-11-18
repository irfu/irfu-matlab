classdef calib
%
% Class for (1) loading calibration data from file, and (2) library/utility functions that calibrate data.
% An instance contains all loaded RCTs.
%
% NOTE: RCT reading functions assume that the same type of RCT (BIAS, LFR, TDS-CWF or TDS-RSWF) is identical (in all
% relevant parts) for both the RODP and ROC-SGSE pipeline.
%
%
% NOTE: UNFINISHED. NOT USED BY MAIN PROGRAM YET.
%
%
% IMPLEMENTATION NOTES
% ====================
% Coded with extra many assertions since
% (1) undiscovered calibration bugs could be considered extra bad
% (2) it is expected to be hard to detect certain bugs
% (3) it could be hard to use automatic testing here.
% (4) to detect changing RCT formats, in particular from external teams.
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% Deg = Degrees (angle). 360 degrees=2*pi radians.
% CPV = counts/volt
% VPC = volt/count
% APC = ampere/count
% RPS = radians/second
% (count = TM/TC unit)
% LSF = LFR Sampling Frequency (F0...F3)
% TF = Transfer function (Z=Z(omega))
% FTF = Forward Transfer Function = TF that describes physical input-to-output (not the reverse)
% ITF = Inverse Transfer Function = TF that describes physical output-to-input (not the reverse)
% 
% BLTS = BIAS-LFR/TDS Signals
% ---------------------------
% Signals somewhere between the LFR/TDS ADCs and the antenna side of the BIAS demuxer
% including the BIAS transfer functions. Like BIAS_i, i=1..5, but includes various stages of calibration/non-calibration,
% including in particular
%   - TM units (inside LFR/TDS),
%   - BLTS interface volt (at the physical boundary BIAS-LFR/TDS (BIAS_i)), and
%   - calibrated values inside BIAS but without demuxer addition and subtraction inside BIAS (i.e. including
%     using BIAS offsets, BIAS transfer functions; volt).
% NOTE: Definition is partly created to avoid using term "BIAS_i" since it is easily confused with other things (the
% subsystem BIAS, bias currents), partly to include various stages of calibration.
%
% ASR = Antenna Signal Representation
% -----------------------------------
% Those measured signals which are ultimately derived/calibrated by BICAS,     i.e. Vi_LF, Vij_LF, Vij_LF_AC (i,j=1..3)
% in the BIAS specification.
% NOTE: This is different from the physical antenna signals (Vi_LF) which are essentially a subset of ASR, except for
% calibration errors and filtering.
% NOTE: This is different from the set Vi_DC, Vij_DC, Vij_AC of which a subset are equal to BIAS_i (which subset it is
% depends on the demux mode) and which is always in LFR/TDS calibrated volts.
%
% BIAS_i, i=1..5
% --------------
% Defined in BIAS specifications document. Equal to the physical signal at the physical boundary between BIAS and
% LFR/TDS. Unit: LFR/TDS calibrated volt. Mostly replaced by BLTS+specified unit in the code.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-15
%



% BOGIQ:
% ------
% TODO-NEED-INFO: Where does the parasitic capacitance TF fit into the calibration formulas?
% TODO-NEED-INFO: Parasitic capacitance values?
%
% PROPOSAL: Need multiple functions
%   Select a TDS/LFR TF
%   Combine TFs: TDS/LFR, BIAS, possibly capacitance-TF.
%   Apply TF (inverse)
%   Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
%
% PROPOSAL: Create general-purpose read_CDF function which handles indices correctly (1 vs many records).
% PROPOSAL: Function for permuting indices to handle dataobj's handling of 1 record-case.
%
% PROPOSAL: Derive the alpha, beta and gamma_lg/hg values from the BIAS transfer functions and log them.
%   NOTE: gamma values apply to AC (highpass filter) and their comparable amplitude is not at 0 Hz but at some other frequency.
%   NOTE: TFs may change over time.
%
% PROPOSAL: Calibration functions that do not work on a single sequence, but a list of sequences (with the same
% settings; not necessarily the same calibration time).
%
% TODO-DECISION: How distribute the "calibration formulas/algorithms between
%   (1) calibrate_* functions, 
%   (2) functions that select relevant calibration data (get_BIAS_calib_data)
%   (2) RCT reading functions,
%   (3) classes that store TFs?
%   Ex: Invert the (rat.func., tabulated) TFs
%   Ex: Extrapolate the tabulated TFs
%   Ex: If one modifies the TFs before applying the inverse (hypothetical; not implemented)
%   --
%   PROPOSAL: Separate modifications/choices/code that
%       (1) only need to be done/run once during the execution: modification of calibration data,
%       (2) are done every calibration (call to calibrate_*):   exact algorithms/formulas through which data is fed
%   PROPOSAL: read_*_RCT should not modify any calibration data, just store it: Not invert TFs, not extrapolate TFs to 0 Hz.
%       CON: Must separate unmodified and modified calibration data.
%       CON: "Must" label variables that they are modified/unmodified.
%       CON-PROBLEM: No clear line for what is modification or not.
%           Ex: TF frequency Hz-->omega rad/s
%           Ex: TF amplitude+phase-->Z
%       TODO-DECISION: Where is it natural to modify calibration data then?!
%   PROPOSAL: General philosophy should be that calibrate_* chooses as much as possible, and thus chooses different
%             functions to call.
%
% PROPOSAL: Move RCT functions to separate class.
%
% PROPOSAL: Simplify/clean-up calibrate_*
%   PROPOSAL: Function for dtSec.
%   PROPOSAL: Function for iCalibTimeL/H
%   PROPOSAL: Assertions for same number of Epoch as samples.
% PROPOSAL: SETTINGS for scalar calibrations.
%
%
% BOGIQ: RCT-reading functions
% ============================
% PROPOSAL: Use same code/function for reading calibration table, as for reading dataset (and master cdfs)?
% PROPOSAL: Assert CDF skeleton/master version number.
% PROPOSAL: Assert skeleton/master.
% PROPOSAL: Assert/warn (depending on setting?) file units.
% PROPOSAL: Only use units in variable names.
% PROPOSAL: Use utility function for reading every zVariable.
%   PROPOSAL: Assert units from zVar attributes.



    properties(Access=private)
        
        Bias
        LfrItfTable
        tdsCwfFactorsVpc
        TdsRswfItfList
        
        SETTINGS
    end

    %###################################################################################################################

    methods(Access=public)

        

        function obj = calib(calibrationDir, pipelineId, SETTINGS)
            % TODO-DECISION: Is it wise to specify the paths in the constructor? Read the filenames (and relative
            %                directory) from the constants instead?
            % TODO-DECISION: Appropriate to use SETTINGS this way? Submit calibration data directly?
            
            % IMPLEMENTATION NOTE: Must assign obj.SETTINGS before calling methods that rely on it having been set, e.g.
            % find_RCT.
            obj.SETTINGS = SETTINGS;
            
            obj.Bias             = obj.read_log_RCT(calibrationDir, pipelineId, 'BIAS');
            obj.LfrItfTable      = obj.read_log_RCT(calibrationDir, pipelineId, 'LFR');
            obj.tdsCwfFactorsVpc = obj.read_log_RCT(calibrationDir, pipelineId, 'TDS-CWF');
            obj.TdsRswfItfList   = obj.read_log_RCT(calibrationDir, pipelineId, 'TDS-RSWF');
        end



        % Convert/calibrate TC bias current: TM units --> physical units.
        %
        % NOTE: This is the normal way of obtaining bias current in physical units (as opposed to HK bias current).
        function biasCurrentAmpere = calibrate_TC_bias_TM_to_bias_current(obj, Epoch, tcBiasTm, iAntenna)
            % BUG: Can not handle time.
            
            % ASSERTION
            bicas.proc_utils.assert_Epoch(Epoch)   % NOTE: Asserts column vector ("zvar"-like).
            assert(isa(tcBiasTm, 'uint16'))
            
            %==============================
            % Obtain calibration constants
            %==============================
            iEpochListL = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochL);   iCalibTimeL = iEpochListL{1}(1);
            offsetAmpere = obj.Bias.Current.offsetsAmpere(iCalibTimeL, iAntenna);
            gainApc      = obj.Bias.Current.gainsApc(     iCalibTimeL, iAntenna);

            % CALIBRATE
            biasCurrentAmpere = offsetAmpere + gainApc .* tcBiasTm;    % LINEAR FUNCTION
        end



        % Convert/calibrate diagnostic HK TM bias current values to physical units.
        % Refers to BIAS HK zVars HK_BIA_BIAS1/2/3.
        %
        %
        % NOTES
        % =====
        % IMPORTANT NOTE: The HK bias current values are measured onboard but are only meant as DIAGNOSTIC values, NOT
        % AS THE PROPER BIAS CURRENT values for nominal use. Therefore the values should only be seen as approximate.
        %
        % NOTE: Walter Puccio, IRFU 2019-09-06: Values are measured on the order of once per second (and sent back as HK
        % even more rarely). Expect errors on the order of 5%.
        %
        % NOTE: The calibration data are NOT stored in the BIAS RCT.
        %
        % NOTE: The conversion function can be found in the BIAS specification, sections 3.4.4.{1-3} ("BIAS1" etc) under
        % "Telemetry". (Not to be confused with the corresponding telecommands.). The conversion functions are identical
        % for all three probes.
        %
        function biasCurrentAmpere = calibrate_HK_bias_TM_to_bias_current(obj, hkBiasCurrentTm, iAntenna)
            % ASSERTION: zVar HK_BIA_BIAS1/2/3's class in BIAS HK. Not strictly required, but the variable has to be
            % some integer which can contain the information.
            assert(isa(hkBiasCurrentTm, 'uint16'))
            
            offsetTm = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.HK_BIAS_CURRENT.OFFSET_TM');
            gainApc  = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.HK_BIAS_CURRENT.GAIN_APC');
            
            %===================================================================================================
            % CALIBRATE
            % ---------
            % Unsigned integer which represents ~signed integer.
            % ==> Intervals 0..0x7FFF and 0x8000...0xFFFF need to "change places".
            % ==> Need to flip bit representing sign to have one interval 0...0xFFFF with monotonic function
            %     TM-to-calibrated values.
            %===================================================================================================
            hkBiasCurrentTm   = bitxor(hkBiasCurrentTm, hex2dec('8000'));                       % FLIP BIT
            biasCurrentAmpere = gainApc(iAntenna) * (hkBiasCurrentTm + offsetTm(iAntenna));     % LINEAR FUNCTION
        end



        % NOTE: TEMPORARY IMPLEMENTATION(?) Only uses first Epoch value for determining calibration values.
        %
        function asrSamplesVolt = calibrate_LFR(obj, Epoch, lfrSamplesTm, iBltsChannel, BltsAsrType, iLfrFreq, biasHighGain)
            
            bicas.proc_utils.assert_Epoch(Epoch)
            assert(numel(lfrSamplesTm) == numel(Epoch))
            assert((1 <= iBltsChannel) && (iBltsChannel <= 5))
            assert((1 <= iLfrFreq)     && (iLfrFreq     <= 4), 'Illegal iLfrFreq value.')
            
            dtSec = double(Epoch(end) - Epoch(1)) / (numel(Epoch)-1) * 1e-9;   % Unit: s   (Epoch unit: ns)
            enableDetrending = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.DETRENDING_ENABLED');

            %==============================
            % Obtain calibration constants
            %==============================
            iEpochListL = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochL);   iCalibTimeL = iEpochListL{1}(1);
            iEpochListH = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochH);   iCalibTimeH = iEpochListH{1}(1);
            BiasCalibData = obj.get_BIAS_calib_data(BltsAsrType, biasHighGain, iCalibTimeL, iCalibTimeH);
            LfrTf = obj.LfrItfTable{iLfrFreq}{iBltsChannel};
            
            %============================================================
            % CALIBRATE: LFR TM --> LFR/BIAS interface volt --> ASR volt
            %============================================================
            % Create combined TF for LFR and BIAS.
            tf = @(omega) (...
                LfrTf.eval_linear(omega) ...
                .* ...
                BiasCalibData.Itf.eval(omega));

            % APPLY TRANSFER FUNCTION
            asrSamplesVolt = bicas.utils.apply_transfer_function(dtSec, lfrSamplesTm, tf, ...
                'enableDetrending', enableDetrending);
            
            asrSamplesVolt = asrSamplesVolt + BiasCalibData.offsetVolt;
        end



        % ARGUMENTS
        % =========
        % tdsCwfSamplesTm : Samples.
        % iBltsChannel    : 1..3.
        % BltsAsrType     : Struct with fields
        %       .antennas : [iAntenna] for single antenna, or [iAntenna, jAntenna] for diff antenna.
        %       .category : String constant
        %
        %
        % NOTE: TEMPORARY IMPLEMENTATION(?) Only uses first Epoch value for determining calibration values.
        %
        function asrSamplesVolt = calibrate_TDS_CWF(obj, Epoch, tdsCwfSamplesTm, iBltsChannel, BltsAsrType)
            % PROPOSAL: Some kind of assertion (assumption of) constant sampling frequency.
            
            % ASSERTIONS
            bicas.proc_utils.assert_Epoch(Epoch)
            assert(numel(tdsCwfSamplesTm) == numel(Epoch))
            assert((1 <= iBltsChannel) && (iBltsChannel <= 3))
            
            dtSec = double(Epoch(end) - Epoch(1)) / (numel(Epoch)-1) * 1e-9;   % Unit: s   ("Epoch" unit: ns)
            enableDetrending = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.DETRENDING_ENABLED');

            %==============================
            % Obtain calibration constants
            %==============================
            iEpochListL = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochL);   iCalibTimeL = iEpochListL{1}(1);
            iEpochListH = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochH);   iCalibTimeH = iEpochListH{1}(1);
            % NOTE: Low/high gain is irrelevant for TDS. Argument value arbitrary.
            BiasCalibData = obj.get_BIAS_calib_data(BltsAsrType, 0, iCalibTimeL, iCalibTimeH);

            %===============================================
            % CALIBRATE: TDS TM --> TDS/BIAS interface volt
            %===============================================
            bltsInterfVolt = obj.tdsCwfFactorsVpc(iBltsChannel) * tdsCwfSamplesTm;    % MULTIPLICATION
            
            %=====================================================
            % CALIBRATE: TDS/BIAS interface volt --> antenna volt
            %=====================================================
            % Create TF for BIAS.
            % APPLY TRANSFER FUNCTION
            asrSamplesVolt = bicas.utils.apply_transfer_function(...
                dtSec, ...
                bltsInterfVolt, ...
                @(omega) BiasCalibData.Itf.eval(omega), ...
                'enableDetrending', enableDetrending);                  

            asrSamplesVolt = asrSamplesVolt + BiasCalibData.offsetVolt;
        end



        function asrSamplesVolt = calibrate_TDS_RSWF(obj, Epoch, tdsRswfSamplesTm, iBltsChannel, BltsAsrType)
            % ASSERTIONS
            bicas.proc_utils.assert_Epoch(Epoch)
            assert(numel(tdsRswfSamplesTm) == numel(Epoch))
            assert((1 <= iBltsChannel) && (iBltsChannel <= 3))
            
            dtSec = double(Epoch(end) - Epoch(1)) / (numel(Epoch)-1) * 1e-9;   % Unit: s   ("Epoch" unit: ns)
            enableDetrending = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.DETRENDING_ENABLED');
            
            %==============================
            % Obtain calibration constants
            %==============================
            iEpochListL = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochL);   iCalibTimeL = iEpochListL{1}(1);
            iEpochListH = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochH);   iCalibTimeH = iEpochListH{1}(1);
            % NOTE: Low/high gain is irrelevant for TDS. Argument value arbitrary.
            BiasCalibData = obj.get_BIAS_calib_data(BltsAsrType, 0, iCalibTimeL, iCalibTimeH);
            
            itf = @(omega) (...
                obj.TdsRswfItfList{iBltsChannel}.eval_linear(omega) .* ...
                BiasCalibData.Itf.eval(omega));

            %====================================
            % CALIBRATE: TDS TM --> antenna volt
            %====================================
            asrSamplesVolt = bicas.utils.apply_transfer_function(...
                dtSec, ...
                tdsRswfSamplesTm, ...
                itf, ...
                'enableDetrending', enableDetrending);
            
            asrSamplesVolt = asrSamplesVolt + BiasCalibData.offsetVolt;
        end



    end    % methods(Access=public)

    
    
    %###################################################################################################################

    
    
    methods(Access=private)
        
        
        
        % Read any single RCT file, and log it. Effectively wraps the different RCT-reading functions.
        % 
        % IMPLEMENTATION NOTE: This method exists to
        % (1) run shared code that should be run when reading any RCT (logging, algorithm for finding file),
        % (2) separate logging from the RCT-reading code, so that one can read RCTs without BICAS.
        % IMPLEMENTATION NOTE: This method is an instance method only because of find_RCT.
        %
        %
        % ARGUMENTS
        % =========
        % pipelineId, rctId : String constants representing pipeline and RCT to be read.
        %
        function RctCalibData = read_log_RCT(obj, calibrationDir, pipelineId, rctId)
            
            filePath = bicas.RCT.find_RCT_by_SETTINGS_regexp(calibrationDir, pipelineId, rctId, obj.SETTINGS);
            bicas.logf('info', 'Reading %4s %-8s RCT: "%s"', pipelineId, rctId, filePath)

            switch(rctId)
                case 'BIAS'     ; Rcd = bicas.RCT.read_BIAS_RCT(filePath);
                case 'LFR'      ; Rcd = bicas.RCT.read_LFR_RCT(filePath);
                case 'TDS-CWF'  ; Rcd = bicas.RCT.read_TDS_CWF_RCT(filePath);
                case 'TDS-RSWF' ; Rcd = bicas.RCT.read_TDS_RSWF_RCT(filePath);
                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal rctId="%s"', rctId);
            end
            RctCalibData = Rcd;
        end



        % Return subset of already loaded BIAS calibration data, for specified settings.
        %
        function BiasCalibData = get_BIAS_calib_data(obj, BltsAsrType, biasHighGain, iCalibTimeL, iCalibTimeH)
            
            switch(BltsAsrType.category)
                case 'DC single'
                    assert(isscalar(BltsAsrType.antennas))
                    BiasItfSet    = obj.Bias.ItfSet.DcSingle;
                    offsetVolt = obj.Bias.dcSingleOffsetsVolt(iCalibTimeH, BltsAsrType.antennas);
                case 'DC diff'
                    BiasItfSet    = obj.Bias.ItfSet.DcDiff;
                    if     isequal(BltsAsrType.antennas(:)', [1,2]);   offsetVolt = obj.Bias.DcDiffOffsets.E12Volt(iCalibTimeH);
                    elseif isequal(BltsAsrType.antennas(:)', [2,3]);   offsetVolt = obj.Bias.DcDiffOffsets.E23Volt(iCalibTimeH);
                    elseif isequal(BltsAsrType.antennas(:)', [1,3]);   offsetVolt = obj.Bias.DcDiffOffsets.E13Volt(iCalibTimeH);
                    else
                        error('BICAS:calib:Assertion:IllegalArgument', 'Illegal BltsAsrType.antennas.');
                    end
                case 'AC'
                    if biasHighGain ; BiasItfSet = obj.Bias.ItfSet.AcHighGain;   offsetVolt = 0;
                    else            ; BiasItfSet = obj.Bias.ItfSet.AcLowGain;    offsetVolt = 0;
                    end
                otherwise
                    error('BICAS:calib:IllegalArgument:Assertion', 'Illegal argument BltsAsrType.category=%s', BltsAsrType.category)
            end
            BiasCalibData.Itf        = BiasItfSet{iCalibTimeL};
            BiasCalibData.offsetVolt = offsetVolt;
        end



    end    % methods(Access=private)

    %###################################################################################################################

    
    
    %methods(Static, Access=public)
    methods(Static, Access=private)



        % Internal utility function.
        %
        % Find which BIAS calibration data to use for specific Epoch values.
        % Assumes that the result is exactly one interval in Epoch so that transfer functions can be applied to it.
        % Therefore assumes that time values are increasing.
        %
        %
        % ARGUMENTS AND RETURN VALUES
        % ===========================
        % Epoch         : Column vector with CDF Epoch values. Must be monotonically increasing.
        % EpochEdgeList : Bias.epochL or Bias.epochH.          Must be monotonically increasing.
        % iEpochList    : Cell array. iEpochList{jCalib} = 1D vector of indices into Epoch for which to use jCalib as time
        %                 index into BIAS calibration data. Every vector should describe a continuous interval in the
        %                 original order. One can thus safely apply transfer functions on data selected this way.
        %                 NOTE: Intervals can be empty.
        %
        function [iEpochList] = get_calibration_time_interval(Epoch, EpochEdgeList)
            % PROPOSAL: Re-create as generic function that splits sorted vector into sub-arrays.
            
            % ASSERTIONS
            bicas.proc_utils.assert_Epoch(Epoch)
            bicas.proc_utils.assert_Epoch(EpochEdgeList)
            validateattributes(Epoch,         {'numeric'}, {'increasing'})
            validateattributes(EpochEdgeList, {'numeric'}, {'increasing'})

            % NOTE: int64(Inf) = int64 max value.
            % NOTE: Adds Inf to edges. "discretize" assigns NaN to values outside the list of edges, which is hereby avoided.
            jEpochCalib = discretize(Epoch, [EpochEdgeList; Inf], 'IncludedEdge', 'left');    % NaN for Epoch < epochL(1)

            % ASSERTION: All Epoch values were assigned a calibration time index.
            assert(~any(isnan(jEpochCalib)), 'BICAS:calib:get_calibration_time_interval:Assertion:IllegalArgument', ...
                'Found Epoch value(s) which could not be assigned to a time interval with calibration data.')

            iEpochList = {};
            for j = 1:numel(EpochEdgeList)
                iEpochList{j} = find(jEpochCalib == j);
            end
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



    end    %methods(Static, Access=private)

end
