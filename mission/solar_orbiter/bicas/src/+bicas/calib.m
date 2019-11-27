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
% Class is implemented with extra many assertions since
% (1) undiscovered calibration bugs could be considered extra bad
% (2) it is expected to be hard to detect certain bugs
% (3) it could be hard to use automatic testing here.
% (4) to detect changing RCT formats, in particular from external teams.
% --
% All calibration functions of measured data are assumed to accept data from all BLTS (1-5), i.e. including TDS, in
% order to reduce the number assumptions that the calling code needs to make.
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
% TF  = Transfer function (Z=Z(omega))
% FTF = Forward Transfer Function = TF that describes physical input-to-output (not the reverse)
% ITF = Inverse Transfer Function = TF that describes physical output-to-input (not the reverse)
% ASR Volt = Calibrated to volt as for ASR, i.e. the final calibrated (measured) value.
%
% 
% BLTS = BIAS-LFR/TDS Signals
% ---------------------------
% Signals somewhere between the LFR/TDS ADCs and the non-antenna side of the BIAS demuxer
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
% Those "physical antenna signals" which BIAS-LFR/TDS is trying to measure, or a measurement thereof. In reality, the
% terminology is:
% ASR         : Pointer to a specific physical antenna signal, e.g. V12_LF (DC diff, antenna 1-2)
% ASR samples : Samples representing a specific ASR (as opposed to BLTS)
% NOTE: There are 9 ASRs, i.e. they can refer also to signals not represented in any single BLTS.
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
% PROPOSAL: Simplify/clean-up calibrate_*.
%   PROPOSAL: Function for iCalibTimeL/H.
%   PROPOSAL: Assertions for same number of Epoch as samples.
% PROPOSAL: SETTINGS for scalar calibrations.
%
% ~DOCUMENTATION BUG?!!:
%   In practice, BLTS refers to not only signals at the BIAS-LFT/TDS boundary, but also the BLTS number for
%   calibrated signals (at the antenna).
%       PROPOSAL: Change definition?
%       PROPOSAL: Change terminology.
%   In practice, ASR refers to both samples at antennas, and the type of signal.
%   PROPOSAL: Define
%       BLTS = BIAS-LFR/TDS Signal, and always refer to signal/samples. Variable with samples: BltsSamples
%       BLTS ID (BltsId) = Number that identifies one BLTS among other BLTS (BltsNbr, iBlts)
%   PROPOSAL: Abolish ASR. Define
%       AS ID (AsId) = Antenna signal ID (DC single/diff, AC diff; antennas)
%       AS (As)      = Antenna Signal
%       CON: AS to short?
%           PROPOSAL: A=Antenna --> RA=RPW Antenna
%               CON: RAS bad acronym.
%       NOTE: BLTS_src_dest is a superset of AS ID.
%
% PROPOSAL: Other name for bicas.BLTS_src_dest that does not reference BLTS.
%   PRO: Reference to BLTS is confusing.
%   PROPOSAL: Define acronym for all physical signal sources which is a superset of ASR/AS ID.
%       PROPOSAL: Be able to use for both physical signal sources, and for where to place in dataset.
%           PRO: Can use acronym for both, and for class bicas.BLTS_src_dest.
%       PROPOSAL: PSS  = Physical Signal Source
%       PROPOSAL: PS   = Physical Signal
%       PROPOSAL: PSSD = Physical Signal Source or Destination
%       PROPOSAL: PSSR = Physical Signal Source or Representation
%       PROPOSAL: PSSR = Physical Signal Source or Dataset Representation
%   PROPOSAL: Have different classes and acronyms for (1) physical signal sources and (2) dataset representation
%       ("BLTS src" and "BLTS dest") where (2) is in practice a subset of (1).
%   PROPOSAL: Use terms based on BLTS as formal terms (both samples and ID)
%       PROPOSAL: BLTS Source/Physical Signal Source
%           BLTSS, BltsSrc
%           BLTS-PSS, BltsPss
%       PROPOSAL: and BLTS Destination/Store/Dataset Variable
%           BLTSD, BltsDest
%           BLTSS, Bltss
%           BLTSDV, BltsDv
%
% NOTE: Should separately denote:
%   (1) The calibration/units of samples,
%   (2) The "identification" of signals (antenna or BLTS)
%   Ex: A variable may contain
%       (a) antenna-calibrated samples, but labelled as BLTS ID (before code demuxes), or
%       (b) non-calibrated/TM samples but labelled by antennas (disabled calibration for debugging).



    properties(Access=private)
        
        Bias
        LfrItfTable
        tdsCwfFactorsVpc
        TdsRswfItfList
        
        SETTINGS
        
        % Corresponds to SETTINGS key-value.
        % In principle unnecessary, but it shortens/clarifies the code.
        enableDetrending
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
            
            obj.enableDetrending = obj.SETTINGS.get_fv('PROCESSING.CALIBRATION.DETRENDING_ENABLED');
        end



        % Convert/calibrate TC bias current: TM units --> physical units.
        %
        % NOTE: This is the normal way of obtaining bias current in physical units (as opposed to HK bias current).
        %
        function biasCurrentAmpere = calibrate_TC_bias_TM_to_bias_current(obj, tcBiasCurrentTm, iAntenna, iCalibTimeL)
            % BUG: Can not handle time.
            
            % ASSERTION
            assert(isa(tcBiasCurrentTm, 'uint16'))
            
            %==============================
            % Obtain calibration constants
            %==============================
            offsetAmpere = obj.Bias.Current.offsetsAmpere(iCalibTimeL, iAntenna);
            gainApc      = obj.Bias.Current.gainsApc(     iCalibTimeL, iAntenna);

            % CALIBRATE
            biasCurrentAmpere = offsetAmpere + gainApc .* tcBiasCurrentTm;    % LINEAR FUNCTION
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
            
            % ASSERTION: zVar HK_BIA_BIAS1/2/3's class in BIAS HK.
            % Not strictly required, but the variable has to be some integer.
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
        % lfrSamplesTm   : 1D cell array of numeric 1D arrays.
        % asrSamplesVolt : 1D cell array of numeric 1D arrays.
        %
        function asrSamplesVolt = calibrate_LFR(obj, dtSec, lfrSamplesTm, iBltsChannel, BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH, iLfrFreq)
            
            % ASSERTIONS
            assert(iscell(lfrSamplesTm))
            EJ_library.utils.assert.vector(lfrSamplesTm)
            EJ_library.utils.assert.vector(dtSec)
            assert(numel(lfrSamplesTm) == numel(dtSec))
            assert((1 <= iBltsChannel) && (iBltsChannel <= 5))
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            assert((1 <= iLfrFreq)     && (iLfrFreq     <= 4), 'Illegal argument iLfrFreq=%g.', iLfrFreq)

            %==============================
            % Obtain calibration constants
            %==============================
            BiasCalibData = obj.get_BIAS_calib_data(BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            lfrItf        = obj.get_LFR_ITF(iBltsChannel, iLfrFreq);

            %=====================================
            % Create combined TF for LFR and BIAS
            %=====================================
            itf = @(omega) (...
                lfrItf(omega) ...
                .* ...
                BiasCalibData.itf(omega));

            %============================================================
            % CALIBRATE: LFR TM --> LFR/BIAS interface volt --> ASR volt
            %============================================================
            asrSamplesVolt = cell(size(lfrSamplesTm));
            for i = 1:numel(lfrSamplesTm)
                
                % APPLY TRANSFER FUNCTION
                tempAsrSamplesVolt = bicas.utils.apply_transfer_function(dtSec(i), lfrSamplesTm{i}(:), itf, ...
                    'enableDetrending', obj.enableDetrending);
            
                asrSamplesVolt{i} = tempAsrSamplesVolt + BiasCalibData.offsetVolt;
            end
        end



        % ARGUMENTS
        % =========
        % tdsCwfSamplesTm : Samples.
        % iBlts           : 1..5.
        %                   NOTE: TDS does not supply any signal for BLTS 4-5. Therefore the calibration procedure is
        %                   undefined and algorithm should return only NaN. BICAS should supply NaN samples for that
        %                   case anyway though.
        % BltsSrc         : bicas.BLTS_src_dest describing where the signal comes from.
        %
        function asrSamplesVolt = calibrate_TDS_CWF(obj, dtSec, tdsCwfSamplesTm, iBlts, BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH)

            % ASSERTIONS
            EJ_library.utils.assert.vector(dtSec)
            assert(numel(tdsCwfSamplesTm) == numel(dtSec))
            assert(iscell(tdsCwfSamplesTm))
            assert((1 <= iBlts) && (iBlts <= 5))
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            
            asrSamplesVolt = cell(size(tdsCwfSamplesTm));   % Initialize empty output variable.
                
            if ismember(iBlts, [1,2,3])
                
                %==============================
                % Obtain calibration constants
                %==============================
                % NOTE: Low/high gain is irrelevant for TDS. Argument value arbitrary.
                BiasCalibData = obj.get_BIAS_calib_data(BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
                
                for i = 1:numel(tdsCwfSamplesTm)
                    
                    %===============================================
                    % CALIBRATE: TDS TM --> TDS/BIAS interface volt
                    %===============================================
                    tempBltsSamplesInterfVolt = obj.tdsCwfFactorsVpc(iBlts) * tdsCwfSamplesTm{i};    % MULTIPLICATION
                    
                    %=====================================================
                    % CALIBRATE: TDS/BIAS interface volt --> antenna volt
                    %=====================================================
                    % APPLY TRANSFER FUNCTION
                    tempAsrSamplesVolt = bicas.utils.apply_transfer_function(...
                        dtSec(i), ...
                        tempBltsSamplesInterfVolt(:), ...
                        BiasCalibData.itf, ...
                        'enableDetrending', obj.enableDetrending);
                    asrSamplesVolt{i} = tempAsrSamplesVolt + BiasCalibData.offsetVolt;
                end

            else
                
                for i = 1:numel(tdsCwfSamplesTm)
                    asrSamplesVolt{i} = NaN * tdsCwfSamplesTm{i};
                end
            end

        end



        % ARGUMENTS
        % =========
        % tdsCwfSamplesTm : Samples.
        % iBlts           : 1..5.
        %                   NOTE: TDS does not supply any signal for BLTS 4-5. Therefore the calibration procedure is
        %                   undefined and algorithm should return only NaN. BICAS should supply NaN samples for that
        %                   case anyway though.
        % BltsSrc         : bicas.BLTS_src_dest describing where the signal comes from.
        %
        function asrSamplesVolt = calibrate_TDS_RSWF(obj, dtSec, tdsRswfSamplesTm, iBlts, BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH)
            
            % ASSERTIONS
            EJ_library.utils.assert.vector(dtSec)
            assert(iscell(tdsRswfSamplesTm))
            assert(numel(tdsRswfSamplesTm) == numel(dtSec))
            assert((1 <= iBlts) && (iBlts <= 5))
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            
            %==============================
            % Obtain calibration constants
            %==============================
            % NOTE: Low/high gain is irrelevant for TDS. Argument value arbitrary.
            BiasCalibData = obj.get_BIAS_calib_data(BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            
            asrSamplesVolt = cell(size(tdsRswfSamplesTm));   % Initialize empty output variable.
            if ismember(iBlts, [1,2,3])
                itf = @(omega) (...
                    obj.TdsRswfItfList{iBlts}.eval_linear(omega) .* ...
                    BiasCalibData.itf(omega));
                
                %====================================
                % CALIBRATE: TDS TM --> antenna volt
                %====================================
                for i = 1:numel(tdsRswfSamplesTm)
                    tempAsrSamplesVolt = bicas.utils.apply_transfer_function(...
                        dtSec(i), ...
                        tdsRswfSamplesTm{i}(:), ...
                        itf, ...
                        'enableDetrending', obj.enableDetrending);
                    
                    asrSamplesVolt{i} = tempAsrSamplesVolt + BiasCalibData.offsetVolt;
                end
            else
                for i = 1:numel(tdsRswfSamplesTm)
                    asrSamplesVolt{i} = NaN * tdsRswfSamplesTm{i};
                end
            end
            
        end



        function iCalib = get_calibration_time_L(obj, Epoch)
            iCalib = bicas.calib.get_calibration_time(Epoch, obj.Bias.epochL);
        end



        function iCalib = get_calibration_time_H(obj, Epoch)
            iCalib = bicas.calib.get_calibration_time(Epoch, obj.Bias.epochH);
        end



    end    % methods(Access=public)

    
    
    %###################################################################################################################

    
    
    methods(Access=private)
        
        
        
        % Read any single RCT file, and log it. Effectively wraps the different RCT-reading functions.
        % 
        % IMPLEMENTATION NOTE: This method exists to
        % (1) run shared code that should be run when reading any RCT (logging, algorithm for finding file),
        % (2) separate logging from the RCT-reading code, so that one can read RCTs without BICAS.
        % IMPLEMENTATION NOTE: This method is an instance method only because of needing settings for
        %   * find_RCT_by_SETTINGS_regexp, and
        %   * read_LFR_RCT.
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
                case 'LFR'      ; Rcd = bicas.RCT.read_LFR_RCT(filePath, obj.SETTINGS.get_fv('PROCESSING.RCT.LFR.EXTRAPOLATE_TF_AMOUNT_HZ'));
                case 'TDS-CWF'  ; Rcd = bicas.RCT.read_TDS_CWF_RCT(filePath);
                case 'TDS-RSWF' ; Rcd = bicas.RCT.read_TDS_RSWF_RCT(filePath);
                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal rctId="%s"', rctId);
            end
            RctCalibData = Rcd;
        end



        % Return subset of already loaded BIAS calibration data, for specified settings.
        %
        % ARGUMENTS
        % =========
        % biasHighGain : NUMERIC value: 0=Off, 1=ON, or NaN=Value not known.
        %                IMPLEMENTATION NOTE: Needs value to represent that biasHighGain is unknown.
        %                Sometimes, if biasHighGain is unknown, then it is useful to process as usual since some of the
        %                data can still be derived/calibrated, so that the caller does not need to handle the special case.
        %
        function BiasCalibData = get_BIAS_calib_data(obj, BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH)
            
            % ASSERTION
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            assert(isscalar(biasHighGain) && isnumeric(biasHighGain))
            assert(isscalar(iCalibTimeL))
            assert(isscalar(iCalibTimeH))
            
            switch(BltsSrc.category)
                case 'DC single'
                    assert(isscalar(BltsSrc.antennas))
                    BiasItfList = obj.Bias.ItfSet.DcSingle;
                    offsetVolt = obj.Bias.dcSingleOffsetsVolt(iCalibTimeH, BltsSrc.antennas);

                case 'DC diff'
                    BiasItfList = obj.Bias.ItfSet.DcDiff;
                    if     isequal(BltsSrc.antennas(:)', [1,2]);   offsetVolt = obj.Bias.DcDiffOffsets.E12Volt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [2,3]);   offsetVolt = obj.Bias.DcDiffOffsets.E23Volt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [1,3]);   offsetVolt = obj.Bias.DcDiffOffsets.E13Volt(iCalibTimeH);
                    else
                        error('BICAS:calib:Assertion:IllegalArgument', 'Illegal BltsSrc.');
                    end

                case 'AC diff'
                    if     biasHighGain == 1;   BiasItfList = obj.Bias.ItfSet.AcHighGain;   offsetVolt = 0;
                    elseif biasHighGain == 0;   BiasItfList = obj.Bias.ItfSet.AcLowGain;    offsetVolt = 0;
                    elseif isnan(biasHighGain); BiasItf     = @(omega) (omega*NaN);         offsetVolt = NaN;   % NOTE: Set TF such that data becomes NaN.
                    else
                        error('BICAS:calib:Assertion:IllegalArgument', 'Illegal argument biasHighGain=%g.', biasHighGain)
                    end

                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', ...
                        'Illegal argument BltsSrc.category=%s. Can not obtain calibration data for this type of signal.', ...
                        BltsSrc.category)
            end
            
            % CASE: Either "BiasItfList" or "BiasItf" is defined.
            if exist('BiasItfList', 'var')
                % Select (1) TF by calibration time, and (2) convert to function handle.
                BiasItf = @(omegaRps) (BiasItfList{iCalibTimeL}.eval(omegaRps));
            end
            
            BiasCalibData.itf        = BiasItf;
            BiasCalibData.offsetVolt = offsetVolt;
        end
        
        
        
        % Obtain LFR ITF, but handle the case that should never happen for actual non-NaN data (LSF F3 + BLTS 4 or 5)
        % and return a TF that only returns NaN instead. BICAS may still iterate over that combination though when
        % calibrating.
        % 
        function lfrItf = get_LFR_ITF(obj, iBlts, iLfrFreq)
            assert(ismember(iBlts,    [1:5]))
            assert(ismember(iLfrFreq, [1:4]))
            
            if (iLfrFreq == 4) && ismember(iBlts, [4,5])
                lfrItf = @(omegaRps) (omegaRps * NaN);
            else                
                tempLfrItf = obj.LfrItfTable{iLfrFreq}{iBlts};
                lfrItf = @(omegaRps) (tempLfrItf.eval_linear(omegaRps));
            end
        end



    end    % methods(Access=private)

    %###################################################################################################################

    
    
    methods(Static, Access=private)
    %methods(Static, Access=public)
        % NOTE: Public so that automatic test code can call get_calibration_time.



        % Given a sequence of Epoch values, determine for each value which calibration time index should be used. The
        % caller will have to decide which sequence of data that should be calibrated together (e.g. if calibration time
        % changes in the middle of CWF), and which Epoch values should be used to determine calibration time (e.g. first
        % Epoch value for a snapshot determines entire snapshot).
        %
        % ARGUMENTS AND RETURN VALUES
        % ===========================
        % Epoch          : Column vector with Epoch values.
        % CalibEpochList : List of monotonically increasing timestamps ("Epoch format").
        %                  In practice intended to be Bias.epochL or Bias.epochH.
        % iCalibList     : Array. iCalibList(i) = calibration time index for Epoch(i).
        %
        function [iCalib] = get_calibration_time(Epoch, CalibEpochList)
            
            % ASSERTIONS
            bicas.proc_utils.assert_Epoch(Epoch)
            bicas.proc_utils.assert_Epoch(CalibEpochList)
            % IMPLEMENTATION NOTE: Does not work if CalibEpochList is empty, since discretize behaves differently for
            % scalar second argument.
            assert(~isempty(CalibEpochList))
            
            % IMPLEMENTATION NOTE: "discretize" by itself returns NaN for Epoch values outside the outermost edges.
            % Therefore (1) must add upper edge "Inf", (2) asserts non-Nan afterwards.
            % IMPLEMENTATION NOTE: "discretize" behaves differently for scalar second argument. Adding edges at infinity hides
            % this problem. If one does not add infinities and uses a scalar edge list, then one has to treat those
            % cases manually.
            iCalib = discretize(Epoch, [CalibEpochList; Inf], 'IncludedEdge', 'left');
            assert(all(~isnan(iCalib(:))), 'Can not derive which calibration data to use for all specified timestamps.')
        end



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
        % EpochEdgeList : List of monotonically increasing timestamps ("Epoch format").
        %                 In practice intended to be Bias.epochL or Bias.epochH.
        % iEpochList    : Cell array. iEpochList{jCalib} = 1D vector of indices into Epoch for which to use jCalib as time
        %                 index into BIAS calibration data. Every vector should describe a continuous interval in the
        %                 original order. One can thus safely apply transfer functions on data selected this way.
        %                 NOTE: Intervals can be empty.
        %
%         function [iEpochList] = get_calibration_time_interval(Epoch, EpochEdgeList)
%             % PROPOSAL: Re-create as generic function that splits sorted vector into sub-arrays.
%             
%             % ASSERTIONS
%             bicas.proc_utils.assert_Epoch(Epoch)
%             bicas.proc_utils.assert_Epoch(EpochEdgeList)
%             validateattributes(Epoch,         {'numeric'}, {'increasing'})
%             validateattributes(EpochEdgeList, {'numeric'}, {'increasing'})
% 
%             % NOTE: int64(Inf) = int64 max value.
%             % NOTE: Adds Inf to edges. "discretize" assigns NaN to values outside the list of edges, which is hereby avoided.
%             jEpochCalib = discretize(Epoch, [EpochEdgeList; Inf], 'IncludedEdge', 'left');    % NaN for Epoch < epochL(1)
% 
%             % ASSERTION: All Epoch values were assigned a calibration time index.
%             assert(~any(isnan(jEpochCalib)), 'BICAS:calib:get_calibration_time_interval:Assertion:IllegalArgument', ...
%                 'Found Epoch value(s) which could not be assigned to a time interval with calibration data.')
% 
%             iEpochList = {};
%             for j = 1:numel(EpochEdgeList)
%                 iEpochList{j} = find(jEpochCalib == j);
%             end
%         end


            
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



    end    % methods(Static)

end
