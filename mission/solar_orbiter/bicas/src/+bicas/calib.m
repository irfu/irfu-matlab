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



    properties(Access=private, Constant)
        
        % Minimum number of numerator or denominator coefficients in the BIAS RCT.
        N_MIN_TF_NUMER_DENOM_COEFFS = 8;
        
        % Minimum number of entries in tabulated transfer functions in RCTs.
        TF_TABLE_MIN_LENGTH = 10;
        
    end


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



        % Find the path to the RCT to use, using the filenames specified in the documentation. If there are multiple
        % matching candidates, choose the latest one as indicated by the filename.
        %
        %
        % ARGUMENTS
        % =========
        % pipelineId, rctId : String constants representing pipeline and RCT to be read.
        %
        %
        % IMPLEMENTATION NOTES
        % ====================
        % Useful to have this as separate functionality so that the chosen RCT to use can be
        % explicitly overridden via e.g. settings.
        % This method is only an instance method so that it can use SETTINGS.
        %
        function path = find_RCT(obj, calibrationDir, pipelineId, rctId)

            %============================
            % Create regexp for filename
            %============================
            pipelineSettingsSegm = EJ_library.utils.translate({...
            {'ROC-SGSE', 'RGTS'}, 'RGTS';...
            {'RODP'},             'RODP'}, ...
            pipelineId, 'BICAS:calib:Assertion:IllegalArgument', sprintf('Illegal pipelineId="%s"', pipelineId));
            clear pipelineId
            
            % IMPLEMENTATION NOTE: This translation statement
            % (1) verifies the argument, AND
            % (2) separates the argument string constants from the SETTINGS naming convention.
            analyzerSettingsSegm = EJ_library.utils.translate({...
                {'BIAS'},     'BIAS'; ...
                {'LFR'},      'LFR'; ...
                {'TDS-CWF'},  'TDS-LFM-CWF'; ...
                {'TDS-RSWF'}, 'TDS-LFM-RSWF'}, ...
            rctId, 'BICAS:calib:Assertion:IllegalArgument', sprintf('Illegal rctId="%s"', rctId));
            clear rctId
            filenameRegexp = obj.SETTINGS.get_fv(sprintf('PROCESSING.RCT_REGEXP.%s.%s', pipelineSettingsSegm, analyzerSettingsSegm));



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
                error('BICAS:calib:CannotFindRegexMatchingRCT', ...
                    'Can not find any calibration file that matches regular expression "%s" in directory "%s".', ...
                    filenameRegexp, calibrationDir);
            end
            % CASE: There is at least one candidate file.
            
            filenameList = sort(filenameList);
            filename     = filenameList{end};
            path = fullfile(calibrationDir, filename);
            
            if numel(filenameList) > 1
                % WARNING/INFO/NOTICE
                msg = sprintf(...
                    ['Found multiple calibration files matching regular expression "%s"\n', ...
                     'in directory "%s".\n', ...
                     'Selecting the latest one as indicated by the filename: "%s".\n'], ...
                    filenameRegexp, calibrationDir, filename);
                for i = 1:numel(filenameList)
                    msg = [msg, sprintf('    %s\n', filenameList{i})];
                end
                bicas.log('debug', msg)
            end
            
            % IMPLEMENTATION NOTE: Not logging which calibration file is selected, since this function is not supposed
            % to actually load the content.
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
        function RctCalibData = read_log_RCT(obj, calibrationDir, pipelineId, rctId)
            
            filePath = obj.find_RCT(calibrationDir, pipelineId, rctId);
            bicas.logf('info', 'Reading %-4s %-8s RCT: "%s"', pipelineId, rctId, filePath)
            
            % ~EXPERIMENTAL. Correct code, but EXPERIMENTAL wrt. convenience/clarity etc.
%             funcHandle = EJ_library.utils.translate({...
%                 {'BIAS'    } @() (bicas.calib.read_BIAS_RCT(filePath));
%                 {'LFR'     } @() (bicas.calib.read_LFR_RCT(filePath));
%                 {'TDS-CWF' } @() (bicas.calib.read_TDS_CWF_RCT(filePath));
%                 {'TDS-RSWF'} @() (bicas.calib.read_TDS_RSWF_RCT(filePath))}, ...
%                 rctId, 'BICAS:calib:Assertion:IllegalArgument', sprintf('Illegal rctId="%s"', rctId));
%             
%             RctCalibData = funcHandle();

            switch(rctId)
                case 'BIAS'     ; Rcd = bicas.calib.read_BIAS_RCT(filePath);
                case 'LFR'      ; Rcd = bicas.calib.read_LFR_RCT(filePath);
                case 'TDS-CWF'  ; Rcd = bicas.calib.read_TDS_CWF_RCT(filePath);
                case 'TDS-RSWF' ; Rcd = bicas.calib.read_TDS_RSWF_RCT(filePath);
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
        
        
        
        function [Bias] = read_BIAS_RCT(filePath)
            % TODO-DECISION: How handle time?
            %   PROPOSAL: "Only" access the BIAS values (trans.func and other) through a function instead of selecting indices in a data struct.
            %       PROPOSAL: (private method) [omegaRps, zVpc] = get_transfer_func(epoch, signalType)
            %           signalType = 'DC single' etc
            
            Do = dataobj(filePath);
            
            % Constants for interpreting the array indices in the CDF.
            NUMERATOR   = 1;
            DENOMINATOR = 2;
            %
            DC_SINGLE = 1;
            DC_DIFF   = 2;
            AC_LG     = 3;
            AC_HG     = 4;
            
            try
                % NOTE: Assumes 1 CDF record or many (time-dependent values).
                % ==> Must handle that dataobj assigns differently for these two cases.
                epochL                   = bicas.calib.norm_do_zv(Do.data.Epoch_L);
                epochH                   = bicas.calib.norm_do_zv(Do.data.Epoch_H);
                biasCurrentOffsetsAmpere = bicas.calib.norm_do_zv(Do.data.BIAS_CURRENT_OFFSET);      % DEPEND_0 = Epoch_L
                biasCurrentGainsApc      = bicas.calib.norm_do_zv(Do.data.BIAS_CURRENT_GAIN);        % DEPEND_0 = Epoch_L
                dcSingleOffsetsVolt      = bicas.calib.norm_do_zv(Do.data.V_OFFSET);                 % DEPEND_0 = Epoch_H
                dcDiffOffsetsVolt        = bicas.calib.norm_do_zv(Do.data.E_OFFSET);                 % DEPEND_0 = Epoch_H
                tfCoeffs                 = bicas.calib.norm_do_zv(Do.data.TRANSFER_FUNCTION_COEFFS); % DEPEND_0 = Epoch_L

                nEpochL = size(epochL, 1);
                nEpochH = size(epochH, 1);

                % IMPLEMENTATION NOTE: Corrects for what seems to be a bug in dataobj. dataobj permutes/removes indices,
                % and permutes them differently depending on the number of CDF records (but wrong in all cases).
                %
                % 1 CDF record : cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       F/TTT"   # 3=number of dimensions/record
                % 2 CDF records: cdfdump: "TRANSFER_FUNCTION_COEFFS CDF_DOUBLE/1   3:[2,8,4]       T/TTT"
                % 1 CDF record:   size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [  4 2 8]
                % 2 CDF records:  size(Do.data.TRANSFER_FUNCTION_COEFFS.data) == [2 4 2 8]                
                tfCoeffs = permute(tfCoeffs, [1, 4,3,2]);



                %=======================================================
                % ASSERTIONS: Size of tfCoeffs/TRANSFER_FUNCTION_COEFFS
                %=======================================================
                assert(size(tfCoeffs, 1) == nEpochL)
                assert(size(tfCoeffs, 2) >= bicas.calib.N_MIN_TF_NUMER_DENOM_COEFFS)
                assert(size(tfCoeffs, 3) == 2)
                assert(size(tfCoeffs, 4) == 4)

                %================================
                % Assign struct that is returned
                %================================
                Bias.epochL = epochL;
                Bias.epochH = epochH;
                
                Bias.Current.offsetsAmpere   = biasCurrentOffsetsAmpere;
                Bias.Current.gainsApc        = biasCurrentGainsApc;
                Bias.dcSingleOffsetsVolt     = dcSingleOffsetsVolt;
                Bias.DcDiffOffsets.E12Volt   = dcDiffOffsetsVolt(:, 1);
                Bias.DcDiffOffsets.E13Volt   = dcDiffOffsetsVolt(:, 2);
                Bias.DcDiffOffsets.E23Volt   = dcDiffOffsetsVolt(:, 3);
                
                % NOTE: Using name "ItfSet" only to avoid "Itfs" (plural). (List, Table would be wrong? Use "ItfTable"?)
                Bias.ItfSet.DcSingle = bicas.calib.sequence_of_ITFs(...
                    tfCoeffs(:, :, NUMERATOR,   DC_SINGLE), ...
                    tfCoeffs(:, :, DENOMINATOR, DC_SINGLE));
                
                Bias.ItfSet.DcDiff = bicas.calib.sequence_of_ITFs(...
                    tfCoeffs(:, :, NUMERATOR,   DC_DIFF), ...
                    tfCoeffs(:, :, DENOMINATOR, DC_DIFF));
                
                Bias.ItfSet.AcLowGain = bicas.calib.sequence_of_ITFs(...
                    tfCoeffs(:, :, NUMERATOR,   AC_LG), ...
                    tfCoeffs(:, :, DENOMINATOR, AC_LG));
                
                Bias.ItfSet.AcHighGain = bicas.calib.sequence_of_ITFs(...
                    tfCoeffs(:, :, NUMERATOR,   AC_HG), ...
                    tfCoeffs(:, :, DENOMINATOR, AC_HG));
                
                % ASSERTION
                EJ_library.utils.assert.all_equal(...
                   [numel(Bias.ItfSet.DcSingle), ...
                    numel(Bias.ItfSet.DcDiff), ...
                    numel(Bias.ItfSet.AcLowGain), ...
                    numel(Bias.ItfSet.AcHighGain)])
                
                %==========================================================================
                % ASSERTIONS: All variables NOT based on tfCoeffs/TRANSFER_FUNCTION_COEFFS
                %==========================================================================
                bicas.proc_utils.assert_Epoch(Bias.epochL)
                bicas.proc_utils.assert_Epoch(Bias.epochH)
                validateattributes(Bias.epochL, {'numeric'}, {'increasing'})
                validateattributes(Bias.epochH, {'numeric'}, {'increasing'})
                
                assert(ndims(Bias.Current.offsetsAmpere)    == 2)
                assert(size( Bias.Current.offsetsAmpere, 1) == nEpochL)
                assert(size( Bias.Current.offsetsAmpere, 2) == 3)
                assert(ndims(Bias.Current.gainsApc)         == 2)
                assert(size( Bias.Current.gainsApc, 1)      == nEpochL)
                assert(size( Bias.Current.gainsApc, 2)      == 3)
                assert(ndims(Bias.dcSingleOffsetsVolt)      == 2)
                assert(size( Bias.dcSingleOffsetsVolt, 1)   == nEpochH)
                assert(size( Bias.dcSingleOffsetsVolt, 2)   == 3)
                for fn = fieldnames(Bias.DcDiffOffsets)'
                    assert(iscolumn(Bias.DcDiffOffsets.(fn{1}))           )
                    assert(length(  Bias.DcDiffOffsets.(fn{1})) == nEpochH)
                end
                
            catch Exc
                error('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath)
            end
        end



        % LfrItfTable : {iFreq}{iBiasChannel}, iFreq=1..4 representing LFR sampling frequencies F0...F3,
        %                   iFreq=1..3 : iBiasChannel=1..5 for BIAS_1..BIAS_5
        %                   iFreq=4    : iBiasChannel=1..3 for BIAS_1..BIAS_3
        %                  NOTE: This is different from LFR zVar FREQ.
        function LfrItfTable = read_LFR_RCT(filePath)
            Do = dataobj(filePath);
            
            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.

                % NOTE: There are separate TFs for each BLTS channel, not just separate LFR sampling frequencies, i.e.
                % there are 5+5+5+3 TFs.
                freqsHz{1}  = shiftdim(Do.data.Freqs_F0.data);    % NOTE: Index {iLfrFreq}.
                freqsHz{2}  = shiftdim(Do.data.Freqs_F1.data);
                freqsHz{3}  = shiftdim(Do.data.Freqs_F2.data);
                freqsHz{4}  = shiftdim(Do.data.Freqs_F3.data);

                amplCpv{1}  = shiftdim(Do.data.TF_BIAS_12345_amplitude_F0.data);
                amplCpv{2}  = shiftdim(Do.data.TF_BIAS_12345_amplitude_F1.data);
                amplCpv{3}  = shiftdim(Do.data.TF_BIAS_12345_amplitude_F2.data);
                amplCpv{4}  = shiftdim(Do.data.TF_BIAS_123_amplitude_F3.data);

                phaseDeg{1} = shiftdim(Do.data.TF_BIAS_12345_phase_F0.data);
                phaseDeg{2} = shiftdim(Do.data.TF_BIAS_12345_phase_F1.data);
                phaseDeg{3} = shiftdim(Do.data.TF_BIAS_12345_phase_F2.data);
                phaseDeg{4} = shiftdim(Do.data.TF_BIAS_123_phase_F3.data);

                for iLsf = 1:4
                    if iLsf ~= 4
                        nBltsChannels = 5;
                    else
                        nBltsChannels = 3;
                    end
                    
                    lsfFreqsHz  = freqsHz{iLsf};
                    lsfAmplCpv  = amplCpv{iLsf};
                    lsfPhaseDeg = phaseDeg{iLsf};

                    % ASSERTIONS: Check CDF array sizes, and implicitly that the CDF format is the expected one.
                    assert(iscolumn(freqsHz{iLsf}))
                    
                    assert(ndims(lsfAmplCpv)  == 2)
                    assert(ndims(lsfPhaseDeg) == 2)                    
                    assert(size( lsfAmplCpv,  1) >= bicas.calib.TF_TABLE_MIN_LENGTH)
                    assert(size( lsfPhaseDeg, 1) >= bicas.calib.TF_TABLE_MIN_LENGTH)
                    assert(size( lsfAmplCpv,  2) == nBltsChannels)
                    assert(size( lsfPhaseDeg, 2) == nBltsChannels)

                    for iBltsChannel = 1:nBltsChannels
                        
                        % NOTE: Inverting the tabulated TF.
                        Itf = EJ_library.utils.tabulated_transform(...
                            lsfFreqsHz * 2*pi, ...
                            1 ./ lsfAmplCpv(      :, iBltsChannel), ...
                            - deg2rad(lsfPhaseDeg(:, iBltsChannel)), ...
                            'extrapolatePositiveFreqZtoZero', 1);
                        
                        % ASSERTION: ITF
                        assert(~Itf.toward_zero_at_high_freq())
                        
                        LfrItfTable{iLsf}{iBltsChannel} = Itf;
                    end
                end
                
            catch Exc1
                Exc2 = MException('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        
        
        
        function tdsCwfFactorsVpc = read_TDS_CWF_RCT(filePath)
            
            Do = dataobj(filePath);
            
            try                
                % NOTE: Undocumented in CDF: zVar CALIBRATION_TABLE is volt/count for just multiplying the TDS signal (for
                % this kind of data). Is not a frequency-dependent transfer function.
                
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.
                
                tdsCwfFactorsVpc = shiftdim(Do.data.CALIBRATION_TABLE.data);
                
                % ASSERTIONS: Check CDF array sizes, no change in format.
                assert(iscolumn(tdsCwfFactorsVpc))
                assert(size(    tdsCwfFactorsVpc, 1) == 3)
                
            catch Exc1
                Exc2 = MException('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        
        
        
        function TdsRswfItfList = read_TDS_RSWF_RCT(filePath)
            
            Do = dataobj(filePath);
            
            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.
                freqsHz  = shiftdim(Do.data.CALIBRATION_FREQUENCY.data);
                amplVpc  = shiftdim(Do.data.CALIBRATION_AMPLITUDE.data);
                phaseDeg = shiftdim(Do.data.CALIBRATION_PHASE.data);
                
                % ASSERTIONS: Check CDF array sizes, no change in format.
                assert(iscolumn(freqsHz));                
                assert(ndims(amplVpc)     == 2)
                assert(ndims(phaseDeg)    == 2)
                assert(size( amplVpc,  1) == 3)
                assert(size( phaseDeg, 1) == 3)
                assert(size( amplVpc,  2) >= bicas.calib.TF_TABLE_MIN_LENGTH)
                assert(size( phaseDeg, 2) >= bicas.calib.TF_TABLE_MIN_LENGTH)
                
                EJ_library.utils.assert.all_equal([...
                    length(freqsHz), ...
                    size(amplVpc,  2), ...
                    size(phaseDeg, 2) ]);
                
                for iBltsChannel = 1:3
                    Itf = EJ_library.utils.tabulated_transform(...
                        freqsHz * 2*pi, ...
                        amplVpc(         iBltsChannel, :), ...
                        deg2rad(phaseDeg(iBltsChannel, :)), ...
                        'extrapolatePositiveFreqZtoZero', 1);
                    
                    % ASSERTION: INVERTED TF
                    assert(~Itf.toward_zero_at_high_freq(), ...
                        ['TDS RSWF transfer function appears to go toward zero at high frequencies. Has it not been', ...
                        ' inverted/backward in time, i.e. is it not physical output-to-input?'])
                    
                    TdsRswfItfList{iBltsChannel} = Itf;
                end
                
            catch Exc1
                Exc2 = MException('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        


        % Utility function
        %
        % Function for normalizing the indices of dataobj zVariables.
        % dataobj zVariable arrays have different meanings for their indices depending on whether there are one record
        % or many. If there is one record, then there is not record index. If there are multiple records, then the first
        % index represents the record number. This function inserts a size-one index as the first index.
        % 
        % DO   = dataobj(...)
        % data = Do.data.TRANSFER_FUNCTION_COEFFS.data
        %
        % NOTE: Not well tested on different types of zvar array sizes.
        % 
        function data = norm_do_zv(DataobjZVar)
            % PROPOSAL: Move to utils.
            % PROPOSAL: Shorter name:
            %   norm_dataobj_zvar
            %   norm_do_zv_data
            %   norm_do_zv
            
            data = DataobjZVar.data;            
            
            if DataobjZVar.nrec == 1
                %nDims = ndims(data);
                %order = [nDims + 1, 1:nDims];
                %data = permute(data, order);
                data = shiftdim(data, -1);
            end
        end
        
        
        
        % Utility function.
        %
        % NOTE: Arguments describe FTFs. Return value describes ITFs.
        function ItfArray = sequence_of_ITFs(ftfNumCoeffs, ftfDenomCoeffs)
            assert(size(ftfNumCoeffs, 1) == size(ftfDenomCoeffs, 1))
            ItfArray = {};
            
            for i = 1:size(ftfNumCoeffs, 1)
                % IMPORTANT NOTE: Invert TF: FTF --> ITF
                Itf = EJ_library.utils.rational_func_transform(...
                    ftfDenomCoeffs(i,:), ...
                    ftfNumCoeffs(i,:));
                
                % ASSERTIONS
                assert(Itf.has_real_impulse_response())
                % Assert ITF. Can not set proper error message.
                assert(~Itf.zero_in_high_freq_limit(), 'Transfer function is not inverted, i.e. not physical output-to-input.')
                
                ItfArray{end+1} = Itf;
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
