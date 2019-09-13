classdef calib
% Class for (1) loading calibration data from file, and (2) library/utility functions that calibrate data.
%
% NOTE: RCT reading functions assume that the same type of RCT (BIAS, LFR, TDS-CWF or TDS-RSWF) is identical (in all
% relevant parts) for RODP and ROC-SGSE.
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
% TF = Transfer function (Z=Z(omega))
% --
% BLTS : BIAS-LFR/TDS Signals. Signals somewhere between the LFR/TDS ADCs and the antenna side of the BIAS demuxer
%        including the BIAS transfer functions. Like BIAS_i, i=1..5, but includes various stages of calibration/non-calibration,
%        including in particular
%           - TM units (inside LFR/TDS),
%           - BLTS interface volt (at the physical boundary BIAS-LFR/TDS (BIAS_i)), and
%           - calibrated values inside BIAS but before addition and subtraction inside BIAS (i.e. after using BIAS offsets, BIAS
%        transfer functions; volt).
%        NOTE: Partly created to avoid using term "BIAS_i" since it is easily
%        confused with other things (the subsystem BIAS, bias currents), partly to include various stages of
%        calibration.
% ASR  : Antenna Signal Representations. Those measured signals which are ultimately derived/calibrated by BICAS,
%        i.e. Vi_LF, Vij_LF, Vij_LF_AC (i,j=1..3) in the BIAS specification.
%        NOTE: This is different from the physical antenna signals which are
%        essentially subset of ASR (Vi_LF), was it not for calibration errors and filtering.
%        NOTE: This is different from the set Vi_DC, Vij_DC, Vij_AC of which a subset are equal to BIAS_i
%        (which subset it is depends on the demux mode) and which is always in LFR/TDS calibrated volts.
% BIAS_i, i=1..5 : Defined in BIAS specifications document. Equal to the physical signal at the physical boundary
%        between BIAS and LFR/TDS. LFR/TDS calibrated volt. Mostly replaced by BLTS+unit in the code.
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
%
% PROPOSAL: Function for converting frequency list, amplitude list, phase list, into struct (omega, Z).
%
% PROPOSAL: Create general-purpose read_CDF function which handles indices correctly (1 vs many records).
% PROPOSAL: Function for permuting indices to handle dataobj's handling of 1 record-case.
%
% PROPOSAL: Somehow tie BIAS calibration variables closer to EpochL and EpochH. Should not be hardcoded in every
%           instance.
%
% PROPOSAL: Derive the alpha, beta and gamma_lg/hg values from the BIAS transfer functions and log them.
%   NOTE: gamma values apply to AC (highpass filter) and their comparable amplitude is not at 0 Hz but at some other frequency.
%   NOTE: TFs may change over time.
%
% PROPOSAL: Automatically detect in which "direction" a TF is: Delays, advances, or delays & advances the signal.
%   PRO: Can be used for detecting whether a TF is used correctly. Delaying TF should be inverted, advancing TF should
%        be used as is, delaying & advancing TF must be wrong.
%   NOTE: Can probably be read out from phase shifts alone.
%       TODO-DECISION/PROBLEM: How handle periodicity of angles? Must convert to absolute (non-periodic) angles representing rotation.
%           PROPOSAL: (For tabulated TFs) Use EJ's unwrap function.
%   PROBLEM: How do for analytical TFs?!
%       PROPOSAL: Define (numerical, anonymous) function for the phase shift. Numerically search for global min & max.
%           NOTE: Must limit min & max frequency.
%           NOTE: Can stop searching for min/max when found one example of negative/positive value. Does not really need
%                 global min & max.
%   PROPOSAL: Apply TF to numerical spike/delta/Dirac function.
%       PROBLEM: At what sampling rate?
% 
% PROPOSAL: Calibration functions that do not work on a single sequence, but a list of sequences (with the same
% settings; not necessarily the same calibration time).
%
%
%
% BOGIQ: RCT-reading functions
% ============================
% PROPOSAL: Use same code/function for reading calibration table, as for reading dataset (and master cdfs)?
% PROPOSAL: Assert CDF skeleton/master version number.
% PROPOSAL: Assert skeleton/master.
% PROPOSAL: Convert from file units to "mathematical units"? (radians/s, radians (for phase shift))
% PROPOSAL: Convert file transfer function to the right direction (invert it).
%   NOTE: Currently uncertain if such case exists, but it might.
% PROPOSAL: Assert/warn (depending on setting?) file units.
% PROPOSAL: Only use units in variable names.
% PROPOSAL: Use utility function for reading every zVariable.
%   PROPOSAL: Assert units from zVar attributes.



    properties(Access=private)
        
        % Different TFs. Initialized in constructor.
        % IMPLEMENTATION NOTE: Variable names include the prefix "BIAS" to indicate that these are the TFs for BIAS
        % as opposed to TDS and LFR which will have to be added later.
        %biasDcAbsGain        % BIAS gain for DC absolute (non-diff) data (instead of TF).
        %biasDcDiffGain       % BIAS gain for DC diff data (instead of TF).
        %biasDcAbsOffsets     % Unclear definition so far.
        
        %biasAcLoGainTf
        %biasAcHiGainTf
        
        %BiasDcRecordList = [];
        
        Bias
        LfrTfs;
        tdsCwfFactors;
        TdsRswfTfTableList;
        
        SETTINGS;
    end

    %###################################################################################################################

    methods(Access=public)

        function obj = calib(calibrationDir, pipelineId, SETTINGS)
            % TODO-DECISION: Is it wise to specify the paths in the constructor? Read the filenames (and relative directory) from the constants instead?
            % TODO-DECISION: Good to use SETTINGS this way? Submit calibration data directly?
            
            filePath = bicas.calib.find_RCT(calibrationDir, pipelineId, 'BIAS');
            bicas.logf('info', 'Reading BIAS     RCT "%s"', filePath)
            obj.Bias = bicas.calib.read_BIAS_RCT(filePath);
            
            filePath   = bicas.calib.find_RCT(calibrationDir, pipelineId, 'LFR');
            bicas.logf('info', 'Reading LFR      RCT "%s"', filePath)
            obj.LfrTfs = bicas.calib.read_LFR_RCT(filePath);
            
            filePath          = bicas.calib.find_RCT(calibrationDir, pipelineId, 'TDS-CWF');
            bicas.logf('info', 'Reading TDS-CWF  RCT "%s"', filePath)
            obj.tdsCwfFactors = bicas.calib.read_TDS_CWF_RCT(filePath);
            
            filePath               = bicas.calib.find_RCT(calibrationDir, pipelineId, 'TDS-RSWF');
            bicas.logf('info', 'Reading TDS-RSWF RCT "%s"', filePath)
            obj.TdsRswfTfTableList = bicas.calib.read_TDS_RSWF_RCT(filePath);
            
            obj.SETTINGS = SETTINGS;
        end



        % Convert/calibrate from TC bias current in TM units to physical units.
        % This is the normal way of obtaining bias current in physical units.
        function biasCurrentAmpere = TC_bias_TM_to_bias_current(obj, tcBiasTm, Epoch, iAntenna)
            % BUG: Can not handle time.
            
            bicas.proc_utils.assert_Epoch(Epoch)   % NOTE: Asserts column vector ("zvar"-like).
            
            biasCurrentAmpere = obj.Bias.Currents.offsetAmpere(iAntenna) + obj.Bias.Currents.gainApc(iAntenna) * tcBiasTm;
        end

        

        % Convert/calibrate diagnostic HK TM bias current values to physical units.
        % Refers to BIAS HK zVars HK_BIA_BIAS1/2/3.
        %
        % NOTES
        % =====
        % IMPORTANT NOTE: The HK bias current values are measured onboard but are only meant as diagnostic values, not
        % as the proper bias current values for nominal use. Therefore the values should only be seen as approximate.
        % NOTE: Walter Puccio, IRFU 2019-09-06: Values are measured on the order of once per second (and sent back as HK
        % even more rarely). Expect errors on the order of 5%.
        %
        % NOTE: The calibration data are NOT stored in the BIAS RCT.
        %
        % NOTE: The conversion function can be found in the BIAS specification, sections 3.4.4.{1-3} ("BIAS1" etc) under
        % "Telemetry". (Not to be confused with the corresponding telecommands.). The conversion functions are identical
        % for all three probes.
        %
        function biasCurrentAmpere = bias_HK_TM_to_bias_current(obj, hkBiasCurrentTm, iAntenna)
            % ASSERTION: zVar HK_BIA_BIAS1/2/3's class in BIAS HK. Not strictly required, but the variable has to be
            % some integer which can contain the information.
            assert(isa(hkBiasCurrentTm, 'uint16'))
            
            offsetTm = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.HK_BIAS_CURRENT.OFFSET_TM');
            gainApc  = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.HK_BIAS_CURRENT.GAIN_APC');
            
            % Unsigned integer which represents ~signed integer.
            % ==> Intervals 0..0x7FFF and 0x8000...0xFFFF need to "change places".
            % ==> Need to flip bit representing sign to have one interval 0...0xFFFF with monotonic function to calibrated values.
            hkBiasCurrentTm   = bitxor(hkBiasCurrentTm, hex2dec('8000'));
            biasCurrentAmpere = gainApc(iAntenna) * (hkBiasCurrentTm + offsetTm(iAntenna));
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
        function asrSamplesVolt = calibrate_TDS_CWF(obj, Epoch, tdsCwfSamplesTm, iBltsChannel, BltsAsrType)
            % PROPOSAL: Some kind of assertion (assumption of) constant sampling frequency.
            
            bicas.proc_utils.assert_Epoch(Epoch)
            assert(numel(tdsCwfSamplesTm) == numel(Epoch))
            assert((1 <= iBltsChannel) && (iBltsChannel <= 3))
            
            dt = double(Epoch(end) - Epoch(1)) / numel(Epoch) * 1e-9;   % Unit: s   (Epoch unit: ns)

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
            bltsInterfVolt = obj.tdsCwfFactors(iBltsChannel) * tdsCwfSamplesTm;
            
            %=====================================================
            % CALIBRATE: TDS/BIAS interface volt --> antenna volt
            %=====================================================
            tf = @(omega) bicas.calib.eval_analytical_transfer_func(...
                omega, BiasCalibData.tfNumerCoeffs, BiasCalibData.tfDenomCoeffs);
            asrSamplesVolt = bicas.utils.apply_transfer_function(dt, bltsInterfVolt, tf, ...
                'enableDetrending', obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.DETRENDING_ENABLED'));

            asrSamplesVolt = asrSamplesVolt + BiasCalibData.offsetVolt;
        end



        % NOTE: TEMPORARY IMPLEMENTATION(?) Only uses first Epoch value for determining calibration values.
        function asrSamplesVolt = calibrate_LFR(obj, Epoch, lfrSamplesTm, iBltsChannel, BltsAsrType, iLfrFreq, biasHighGain)
            
            bicas.proc_utils.assert_Epoch(Epoch)
            assert(numel(lfrSamplesTm) == numel(Epoch))
            assert((1 <= iBltsChannel) && (iBltsChannel <= 5))
            assert((1 <= iLfrFreq)     && (iLfrFreq     <= 4), 'Illegal iLfrFreq value.')
            
            dt = double(Epoch(end) - Epoch(1)) / numel(Epoch) * 1e-9;   % Unit: s   (Epoch unit: ns)

            %==============================
            % Obtain calibration constants
            %==============================
            iEpochListL = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochL);   iCalibTimeL = iEpochListL{1}(1);
            iEpochListH = bicas.calib.get_calibration_time_interval(Epoch, obj.Bias.epochH);   iCalibTimeH = iEpochListH{1}(1);
            BiasCalibData = obj.get_BIAS_calib_data(BltsAsrType, biasHighGain, iCalibTimeL, iCalibTimeH);
            LfrTf = obj.LfrTfs{iLfrFreq}{iBltsChannel};
            
            %============================================================
            % CALIBRATE: LFR TM --> LFR/BIAS interface volt --> ASR volt
            %============================================================
            tf = @(omega) (...
                bicas.calib.eval_tabulated_transfer_func(omega, LfrTf.omegaRps, LfrTf.zVpc) ...
                .* ...
                bicas.calib.eval_analytical_transfer_func(...
                omega, BiasCalibData.tfNumerCoeffs, BiasCalibData.tfDenomCoeffs));
            asrSamplesVolt = bicas.utils.apply_transfer_function(dt, lfrSamplesTm, tf, ...
                'enableDetrending', obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.DETRENDING_ENABLED'));
            
            asrSamplesVolt = asrSamplesVolt + BiasCalibData.offsetVolt;
        end
        
        
        
    end    % methods(Access=public)
    
    %###################################################################################################################

    methods(Access=private)



        % Internal utility function
        function BiasCalibData = get_BIAS_calib_data(obj, BltsAsrType, biasHighGain, iCalibTimeL, iCalibTimeH)
            
            switch(BltsAsrType.category)
                case 'DC single'
                    assert(isscalar(BltsAsrType.antennas))
                    BiasTfCoeffs = obj.Bias.Tfs.DcSingle;
                    offsetVolt   = obj.Bias.dcSingleOffsetsVolt(iCalibTimeH, BltsAsrType.antennas);
                case 'DC diff'
                    BiasTfCoeffs = obj.Bias.Tfs.DcDiff;
                    if     isequal(BltsAsrType.antennas(:)', [1,2]);   offsetVolt = obj.Bias.DcDiffOffsets.E12Volt(iCalibTimeH);
                    elseif isequal(BltsAsrType.antennas(:)', [2,3]);   offsetVolt = obj.Bias.DcDiffOffsets.E23Volt(iCalibTimeH);
                    elseif isequal(BltsAsrType.antennas(:)', [1,3]);   offsetVolt = obj.Bias.DcDiffOffsets.E13Volt(iCalibTimeH);
                    else
                        error('calib:calibrate_TDS_CWF:Assertion:IllegalArgument', 'Illegal BltsAsrType.antennas.');
                    end
                case 'AC'
                    if biasHighGain
                        BiasTfCoeffs = obj.Bias.Tfs.AcHighGain;
                        offsetVolt = 0;
                    else
                        BiasTfCoeffs = obj.Bias.Tfs.AcLowGain;
                        offsetVolt = 0;
                    end
                otherwise
                    error('BICAS:calib:IllegalArgument:Assertion', 'Illegal argument BltsAsrType.category=%s', BltsAsrType.category)
            end
            BiasCalibData.tfNumerCoeffs = BiasTfCoeffs.numerCoeffs(iCalibTimeL, :);
            BiasCalibData.tfDenomCoeffs = BiasTfCoeffs.denomCoeffs(iCalibTimeL, :);
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
            assert(~any(isnan(jEpochCalib)), 'Found Epoch value(s) which could not be assigned to a time interval with calibration data.')

            iEpochList = {};
            for j = 1:numel(EpochEdgeList)
                iEpochList{j} = find(jEpochCalib == j);
            end
        end
        
        
        
        function [Bias] = read_BIAS_RCT(filePath)
            % NOTE: UNFINISHED. DOES NOT SAVE/RETURN OFFSETS.
            
            % TODO-DECISION: How handle time?
            %   PROPOSAL: "Only" access the BIAS values (trans.func and other) through a function instead of selecting indices in a data struct.
            %       PROPOSAL: (private method) [omegaRps, zVpc] = get_transfer_func(epoch, signalType)
            %           signalType = 'DC single' etc
            
            Do = dataobj(filePath);
            
            % Constants for interpreting the array indices in the CDF.
            NUMERATOR   = 1;
            DENOMINATOR = 2;
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
                assert(size(tfCoeffs, 2) >= 8)
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

                Bias.Tfs.DcSingle.numerCoeffs   = tfCoeffs(:, :, NUMERATOR,   DC_SINGLE);
                Bias.Tfs.DcSingle.denomCoeffs   = tfCoeffs(:, :, DENOMINATOR, DC_SINGLE);
                
                Bias.Tfs.DcDiff.numerCoeffs     = tfCoeffs(:, :, NUMERATOR,   DC_DIFF);
                Bias.Tfs.DcDiff.denomCoeffs     = tfCoeffs(:, :, DENOMINATOR, DC_DIFF);
                
                Bias.Tfs.AcLowGain.numerCoeffs  = tfCoeffs(:, :, NUMERATOR,   AC_LG);
                Bias.Tfs.AcLowGain.denomCoeffs  = tfCoeffs(:, :, DENOMINATOR, AC_LG);
                
                Bias.Tfs.AcHighGain.numerCoeffs = tfCoeffs(:, :, NUMERATOR,   AC_HG);
                Bias.Tfs.AcHighGain.denomCoeffs = tfCoeffs(:, :, DENOMINATOR, AC_HG);
                
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



        % LfrTfs : LfrTfs{iFreq}{iBiasChannel}, iFreq=1..4 for F0..F3,
        %                   iFreq=1..3 : iBiasChannel=1..5 for BIAS_1..BIAS_5
        %                   iFreq=4    : iBiasChannel=1..3 for BIAS_1..BIAS_3
        %                  NOTE: This is different from LFR zVar FREQ.
        function LfrTfs = read_LFR_RCT(filePath)
            % BUG: 5+5+5+3 TFs, NOT 1+1+1+1.
            
            Do = dataobj(filePath);
            
            try
                % ASSUMPTION: Exactly 1 CDF record.
                % IMPLEMENTATION NOTE: Does not want to rely one dataobj special behaviour for 1 record case
                % ==> Remove leading singleton dimensions, much assertions.

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

                for iLfrFreq = 1:4
                    if iLfrFreq ~= 4
                        nBltsChannels = 5;
                    else
                        nBltsChannels = 3;
                    end

                    % ASSERTIONS: Check CDF array sizes, that no change in format.
                    assert(iscolumn(freqsHz{iLfrFreq}))
                    assert(ndims(amplCpv)  == 2)
                    assert(ndims(phaseDeg) == 2)
                    assert(size( amplCpv{iLfrFreq}, 2) == nBltsChannels)
                    assert(size(phaseDeg{iLfrFreq}, 2) == nBltsChannels)

                    for iBltsChannel = 1:nBltsChannels
                        omegaRps = freqsHz{iLfrFreq} * 2*pi;    % NOTE: Duplicates the same frequency list.
                        zVpc     = amplCpv{iLfrFreq}(:,iBltsChannel) .* exp(1i * deg2rad(phaseDeg{iLfrFreq}(:,iBltsChannel)));
                        
                        assert(isvector(omegaRps))
                        assert(isvector(zVpc))
                        assert(numel(omegaRps) == numel(zVpc))
                        
                        LfrTfs{iLfrFreq}{iBltsChannel}.omegaRps = omegaRps;
                        LfrTfs{iLfrFreq}{iBltsChannel}.zVpc     = zVpc;
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
        
        
        
        function TdsRswfTfList = read_TDS_RSWF_RCT(filePath)
            
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
                assert(size( amplVpc, 1)  == 3)
                assert(ndims(phaseDeg)    == 2)
                assert(size( phaseDeg, 1) == 3)
                
                EJ_library.utils.assert.all_equal([...
                    length(freqsHz), ...
                    size(amplVpc,  2), ...
                    size(phaseDeg, 2) ]);
                
                for iBltsChannel = 1:3
                    omegaRps = freqsHz;    % NOTE: Duplicates the same frequency list.
                    zVpc     = shiftdim(amplVpc(iBltsChannel, :) .* exp(1i * deg2rad(phaseDeg(iBltsChannel, :))));
                    
                    assert(iscolumn(omegaRps))
                    assert(iscolumn(zVpc))
                    
                    TdsRswfTfList(iBltsChannel).omegaRps = omegaRps;
                    TdsRswfTfList(iBltsChannel).zVpc     = zVpc;
                end
                
            catch Exc1
                Exc2 = MException('BICAS:calib:CannotInterpretRCT', 'Can not interpret calibration file (RCT) "%s"', filePath);
                Exc2 = Exc2.addCause(Exc1);
                throw(Exc2);
            end
        end
        

        
        % Find the path to the RCT to use, using the filenames specified in the documentation. If there are multiple
        % matching candidates, choose the latest one as indicated by the filename.
        %
        % RCT filenaming convention is described in:	ROC-PRO-DAT-NTT-00006-LES, 1/1 draft, Sect 4.3.2-3.
        %
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
        %
        % ARGUMENTS
        % =========
        % pipelineId, rctId : String constants representing pipeline and RCT to be read.
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
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal pipelineId="%s"', pipelineId)
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
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal rctId="%s"', rctId)
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
                error('BICAS:calib:CannotFindRegexMatchingRCT', ...
                    'Can not find any calibration file that matches regular expression "%s" in directory "%s".', ...
                    filenameRegexp, calibrationDir);
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
                bicas.log('debug', msg)
            end
            
            % IMPLEMENTATION NOTE: Not logging which calibration file is selected, since this function is not supposed
            % to actually load the content.
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
        
        
        
        % Evaluate (complex) analytical transfer function for arbitrary number of omega values.
        % 
        %     nc(1)*s^0 + ... + nc(nNc)*s^(nNc-1)
        % Z = -----------------------------------,   s = i*omega
        %     dc(1)*s^0 + ... + dc(nDc)*s^(nDc-1)
        %
        % nc = numerCoeffs
        % dc = denomCoeffs
        %
        % NOTE: MATLAB's "polyval" uses the opposite convention for the order/index of coefficients.
        %
        function Z = eval_analytical_transfer_func(omega, numerCoeffs, denomCoeffs)
            % PROPOSAL: Reclassify as generic function and move it.
            
            assert(isvector(numerCoeffs))
            assert(isvector(denomCoeffs))
            numerCoeffs = numerCoeffs(:);
            denomCoeffs = denomCoeffs(:);
            
            % ASSERTIONS.
            % NOTE: Column vectors are required for using "flipud".
            assert(iscolumn(numerCoeffs))
            assert(iscolumn(denomCoeffs))
            % ASSERTION: Try to detect accidently switching the coefficients.
            % Denominator polynomial should have at least as high a degree as the numerator polynomial. <==> Z should
            % not diverge as omega-->inf.
            nNc = find(numerCoeffs,   1, 'last' );   % Detect highest-order non-zero coefficient.
            nDc = find(denomCoeffs, 1, 'last' );
            assert(nDc >= nNc)
            
            % Calculate Z
            s = 1i * omega;
            Z = polyval(flipud(numerCoeffs), s) ./ polyval(flipud(denomCoeffs), s);
        end
        
        
        
        function Z = eval_tabulated_transfer_func(evalOmega, tfOmega, tfZ)
            % PROPOSAL: Extend TF to infinite frequency.
            %   PRO: Can play around with high sampling rates for testing.
            
            % Extend TF to zero radians/s.
            % NOTE: Can not just extrapolate/use the lowest-frequency Z to 0 Hz since it has to be real (not complex).
            tfOmega = [0;           tfOmega(:)];
            tfZ     = [abs(tfZ(1)); tfZ(:)];   
            
            Z = interp1(tfOmega, tfZ, evalOmega, 'linear');
        end



        % Function for normalizing the indices of dataobj zVariables.
        % dataobj zVariable arrays have different meanings for their indices depending on whether there are one record
        % or many. If there is one record, then there is not record index. If there are multiple records, then the first
        % index represents the record number. This function inserts a size-one index as the first index.
        % 
        % Do = dataobj(...)
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

            

    end    %methods(Static, Access=private)

end
