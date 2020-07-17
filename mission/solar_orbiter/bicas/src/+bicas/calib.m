%
% Class for (1) loading calibration data from file, and (2) library/utility functions that calibrate data.
% An instance contains all loaded RCTs.
%
% NOTE: RCT reading functions assume that the same type of RCT (BIAS, LFR, TDS-CWF or TDS-RSWF) is identical (in all
% relevant parts) for both the RODP and ROC-SGSE pipeline.
%
%
% SHORTCOMINGS
% ============
% Does not implement parasitic capacitance yet due to lack of calibration values (at least).
%
%
% IMPLEMENTATION NOTES
% ====================
% Class is implemented with extra care, and therefore
% - has extra many assertions
% - is extra careful with identifiers with units
% - is extra careful with well defined terms/shorternings/naming conventions
% since
% (1) undiscovered calibration bugs could be considered extra bad
% (2) it is expected to be hard to detect certain bugs,
% (3) it could be hard to use automatic testing here,
% (4) to detect changing RCT formats, in particular from teams to BIAS team.
% --
% All calibration functions of measured data are assumed to accept data from all BLTS (1-5), i.e. including TDS, in
% order to reduce the number assumptions that the calling code needs to make.
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% Offset = Value (constant) that is ADDED to (not subtracted from) a measured value during the calibration process.
% CA     = Cell Array
% --
% LSF  = LFR Sampling Frequency (F0...F3)
%        NOTE: When used as an array index, 1=F0, ..., 4=F3.
% TF   = Transfer function (Z=Z(omega))
% FTF  = Forward Transfer Function = TF that describes physical input-to-output (not the reverse)
% ITF  = Inverse Transfer Function = TF that describes physical output-to-input (not the reverse)
% CTI  = CALIBRATION_TABLE_INDEX (zVar)
% CTI1 = First  value in record of zVar CALIBRATION_TABLE_INDEX.
% CTI2 = Second value in record of zVar CALIBRATION_TABLE_INDEX.
% RCTS = RCT CALIBRATION_TABLE (glob.attr)+CALIBRATION_TABLE_INDEX (zVar). S = plural
%
%
% UNITS / TYPES OF QUANTITIES
% ---------------------------
% TM      = Telemetry units (in LFR/TDS ADC), or telecommand (TC) units. Using this instead of the term "count".
% IVolt   = Interface Volt = Calibrated volt at the interface between BIAS and LFR/TDS.
% AVolt   = Antenna Volt   = Calibrated volt at the antennas, i.e. the final calibrated (measured) value, including
%           for reconstructed signals (e.g. diffs calculated from singles). May also refer to offsets and values without
%           offsets.
% AAmpere = Antenna ampere      = Calibrated ampere at the antenna.
% SAmpere = Set current ampere. Exactly proportional to bias current in TM.
% TPIV    = TM/interface volt
% IVPT    = Interface volt/TM
% AAPT    = Antenna ampere/TM
% AVPIV   = Antenna volt/interface volt
% IVPAV   = Interface volt/antenna volt
% Deg     = Degrees (angle). 1 revolution=360 degrees=2*pi radians.
% RPS     = Radians/second
% Sec     = Seconds
% 
%
% BLTS = BIAS-LFR/TDS SIGNAL
% ---------------------------
% Signals somewhere between the LFR/TDS ADCs and the non-antenna side of the BIAS demuxer
% including the BIAS transfer functions. Like BIAS_i, i=1..5, but includes various stages of calibration/non-calibration,
% including in particular
%   - TM units (inside LFR/TDS),
%   - Interface volt (at the physical boundary BIAS-LFR/TDS (BIAS_i)), and
%   - Calibrated values inside BIAS but without demuxer addition and subtraction inside BIAS (i.e. including
%     using BIAS offsets, BIAS transfer functions; volt).
% NOTE: Definition is partly created to avoid using term "BIAS_i" since it is easily confused with other things (the
% subsystem BIAS, bias currents), partly to include various stages of calibration.
%
%
% ASR = Antenna Signal Representation
% -----------------------------------
% The "physical antenna signals" which BIAS-LFR/TDS is trying to measure, or a measurement thereof. In reality, the
% terminology is:
% ASR         : Pointer to a specific physical antenna signal, e.g. V12_LF (DC diff, antenna 1-2)
% ASR samples : Samples representing a specific ASR (as opposed to BLTS)
% NOTE: There are 9 ASRs, i.e. they can refer also to signals not represented by any single BLTS, given a chosen mux
% mode (and latching relay setting).
%
%
% BIAS_i, i=1..5
% --------------
% Defined in BIAS specifications document. Equal to the physical signal at the physical boundary between BIAS and
% LFR/TDS. Unit: LFR/TDS calibrated volt. Mostly replaced by BLTS+specified unit in the code.
%
%
% REMINDER: CALIBRATION_TABLE, CALIBRATION_TABLE_INDEX
% ====================================================
% CALIBRATION_TABLE       : Global attribute
%   "Filename of the calibration table(s)."
%   "There must as many as entries than the number of calibration table files associated to the L1R file."
% CALIBRATION_TABLE_INDEX : zVariable
%   "Index of the calibration table(s) value required to generate L2 data files."
%   "Each CDF record must contain 2 elements: the first element must gives the index of the associated CALIBRATION_TABLE
%   entry (i.e., 0 for the first entry, 1 for the second, etc.). The second element must refer to the index of the value
%   to be used inside the calibration table file."
%
% Source: ROC-PRO-DAT-NTT-00006-LES_Iss01_Rev02\(ROC_Data_Products\).Draft2020-04-06.pdf
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-15
%
classdef calib < handle

% BOGIQ:
% ------
% TODO-NEED-INFO: Where does the parasitic capacitance TF fit into the calibration formulas?
% TODO-NEED-INFO: Parasitic capacitance values?
%
% PROPOSAL: Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
%
% PROPOSAL: Create general-purpose read_CDF function which handles indices correctly (1 vs many records).
% PROPOSAL: Function for permuting indices to handle dataobj's handling of 1 record-case.
%
% TODO-DECISION: How distribute the "calibration formulas/algorithms between
%   (1) calibrate_* functions, 
%   (2) functions that select relevant calibration data (get_BIAS_calib_data)
%   (2) RCT reading functions,
%   (3) classes that store TFs?
%   Ex: Invert the (rat.func., tabulated) TFs
%   Ex: Extrapolate the tabulated TFs to zero
%   Ex: Extrapolate the tabulated LFR TF to higher frequencies.
%   Ex: If one modifies the TFs before applying the inverse (hypothetical; not implemented)
%   --
%   PROPOSAL: Separate modifications/choices/code that
%       (1) only need to be done/run once during the execution: modification of calibration data,
%       (2) are done every calibration (call to calibrate_*):   exact algorithms/formulas through which data is fed
%   PROPOSAL: read_*_RCT should not modify any calibration data, just store it: Not invert TFs, not extrapolate TFs to 0 Hz.
%       CON: Must separate unmodified and modified calibration data.
%       CON: "Must" label variables that they are modified/unmodified.
%       CON-PROBLEM: No clear line for what is a modification or not.
%           NOTE: Difference between
%               (1) modification (information potentially destroyed), and
%               (2) conversion (no information destroyed; e.g. format conversion)
%           Ex: TF frequency Hz-->omega rad/s
%           Ex: TF amplitude+phase-->Z
%           Ex: Apply upper frequency cut-off to ITF, in particular analytical ITFs.
%           Ex: Extrapolate tabulated TF
%           Ex: How/where make different choices for how to calibrate?
%               (1) No calibration
%               (2) Scalar calibration (a) with/(b) without offsets
%               (3) Full calibration
%               (4) Full calibration except without parasitic capacitance.
%       TODO-DECISION: Where is it natural to modify calibration data then?!
%   PROPOSAL: General philosophy should be that calibrate_* chooses as much as possible, and thus chooses different
%             functions to call.
%
% ~DOCUMENTATION BUG?!!:
%   PROPOSAL: Abolish ASR. Define acronyms for
%       (1) Antenna signals (AC, DC, singles, diffs)
%       (2) All possible sources of signals of which (1) is a subset (a diff counts as a source).
%           PROPOSAL: Be able to use for both physical signal sources, and for where to place in dataset.
%               PRO: Can use acronym for both, and for class bicas.BLTS_src_dest.
%           PROPOSAL: PSS  = Physical Signal Source
%           PROPOSAL: PS   = Physical Signal
%           PROPOSAL: PSSD = Physical Signal Source or Destination
%           PROPOSAL: PSSR = Physical Signal Source or Representation
%           PROPOSAL: PSSR = Physical Signal Source or Dataset Representation
%
% PROPOSAL: Other name for bicas.BLTS_src_dest that does not reference BLTS.
%   PRO: Reference to BLTS is confusing.
%   PROPOSAL: Define acronym for all physical signal sources which is a superset of ASR.
%   PROPOSAL: Have different classes and acronyms for (1) physical signal sources and (2) dataset representation
%       ("BLTS src" and "BLTS dest") where (2) is in practice a subset of (1).
%
% PROPOSAL: Move TF=0 for HIGH frequencies to modify_tabulated_ITF.
%
% PROPOSAL: Initialize with RCT-data directly. Put RCT-reading outside of class.
% PROPOSAL: Calibration mode for only removing LFR LSF offsets
%   NOTE: Does not apply to TDS.
%
% PROPOSAL: Rename "semi" initialization methods ("semiconstructors")
%       read_non_BIAS_RCTs_by_regexp
%       read_non_BIAS_RCT_by_CALIBRATION_TABLE
%   to make it more clear that they are needed for initialization.
%   PROPOSAL: "init_*"
%   PROPOSAL: "complete_init_*"
%
% PROPOSAL: Function for producing combined FTFs (LFR/TDS + BIAS).
%   PRO: Useful for debugging. Can easily inspect & plot FTFs.
%   CON: Code only reads ITFs.



    %properties(Access=private)
    properties(SetAccess=private, GetAccess=public)
        
        Bias;
        BiasGain             % BIAS scalar (simplified) calibration, not in the RCTs. For debugging/testing purposes.
        HkBiasCurrent
        
        % Cell arrays. {iRct}, iRct=CALIBRATION_TABLE_INDEX(i,1), if CALIBRATION_TABLE_INDEX is used. Otherwise always
        % scalar cell array.
        LfrItfIvptTable    = {};
        tdsCwfFactorsIvpt  = {};
        TdsRswfItfIvptList = {};
        lfrLsfOffsetsTm    = [];   % EXPERIMENTAL. NOTE: Technically, the name contains a tautology (LFR+LSF).
       
        calibrationDir
        
        % Flag for whether non-BIAS data has been loaded, i.e. whether the object has been fully initialized
        % (constructor+one of two methods).
        hasLoadedNonBiasData = false;
        
        % Whether to select non-BIAS RCT using global attribute CALIBRATION_TABLE (and
        % CALIBRATION_TABLE_INDEX(iRecord,1)).
        use_CALIBRATION_TABLE_rcts        
        % Whether to use CALIBRATION_TABLE_INDEX(iRecord,2) for calibration.
        use_CALIBRATION_TABLE_INDEX2

        % Corresponds to SETTINGS key-value.
        enableDetrending
        itfHighFreqLimitFraction;
                
        % What type of calibration to use.
        allVoltageCalibDisabled
        biasOffsetsDisabled
        lfrTdsTfDisabled
        useBiasTfScalar
        
        % Needed so that it can be submitted to bicas.RCT.read_log_RCT_by_SETTINGS_regexp.
        SETTINGS
        L
    end
    
    
    
    properties(Access=private, Constant)
        % Local constant for convenience.
        NAN_TF = @(omegaRps) (omegaRps * NaN);
    end

    %###################################################################################################################

    methods(Access=public)



        % Constructor.
        %
        % IMPORTANT NOTE
        % ==============
        % The constructor only initializes BIAS calibration data. One also needs to call either
        % (1) read_non_BIAS_RCTs_by_regexp, OR
        % (2) read_non_BIAS_RCT_by_CALIBRATION_TABLE to fully initialize the object.
        %
        function obj = calib(calibrationDir, SETTINGS, L)
            % TODO-DECISION: Appropriate to use SETTINGS this way? Submit calibration data directly?
            
            % IMPLEMENTATION NOTE: Must assign obj.SETTINGS before calling methods that rely on it having been set.
            
            obj.SETTINGS         = SETTINGS;
            obj.L                = L;
            
            obj.enableDetrending               = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF_DETRENDING_ENABLED');
            obj.itfHighFreqLimitFraction       = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF_HIGH_FREQ_LIMIT_FRACTION');
            obj.calibrationDir                 = calibrationDir;
            
            obj.Bias                           = obj.read_log_RCT_by_SETTINGS_regexp('BIAS');
            obj.BiasGain.alphaIvpav            = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.ALPHA_IVPAV');
            obj.BiasGain.betaIvpav             = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.BETA_IVPAV');
            obj.BiasGain.gammaIvpav.highGain   = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.HIGH_GAIN');
            obj.BiasGain.gammaIvpav.lowGain    = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.LOW_GAIN');
            
            obj.HkBiasCurrent.offsetTm         = SETTINGS.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.OFFSET_TM');
            obj.HkBiasCurrent.gainAapt         = SETTINGS.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.GAIN_APT');
            
            obj.lfrLsfOffsetsTm                = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR.LSF_OFFSETS_TM');

            obj.allVoltageCalibDisabled        = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.DISABLE');
            obj.biasOffsetsDisabled            = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.DISABLE_OFFSETS');
            obj.lfrTdsTfDisabled               = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR_TDS.TF_DISABLED');
            settingBiasTf                      = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF');
            switch(settingBiasTf)
                case 'FULL'
                    obj.useBiasTfScalar = 0;
                case 'SCALAR'
                    obj.useBiasTfScalar = 1;
                otherwise
                    error('BICAS:calib:Assertion:ConfigurationBug', 'Illegal value for setting PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF="%s".', settingBiasTf)
            end



            %===========================================
            % Log some indicative value(s) in BIAS ITFs
            %===========================================
            DC_FREQ_HQ       = 0;
            AC_DIFF_FREQS_HZ = [0, 1000];
            nBiasEpochL = numel(obj.Bias.epochL);
            for iEpochL = 1:nBiasEpochL
                dcSingleZ = obj.Bias.ItfSet.DcSingleAvpiv{iEpochL}.eval(DC_FREQ_HQ);
                dcDiffZ   = obj.Bias.ItfSet.DcDiffAvpiv{iEpochL}.eval(DC_FREQ_HQ);
                
                % Log bias current calibration
                L.logf('debug', '(%i) BIAS current offsets (%s) [A]',         iEpochL, ...
                    strjoin(EJ_library.str.sprintf_many('% 10e', obj.Bias.Current.offsetsAAmpere(iEpochL, :)), ', '))
                L.logf('debug', '(%i) BIAS current gain    (%s) [A/TM unit]', iEpochL, ...
                    strjoin(EJ_library.str.sprintf_many('% 10e', obj.Bias.Current.gainsAapt(iEpochL, :)),      ', '))
                
                % Log transfer functions (frequency domain), selected frequencies.
                obj.log_ITF_Z(sprintf('(%i) BIAS DC single',          iEpochL), 'AVolt/IVolt', DC_FREQ_HQ,       @(omegaRps) (obj.Bias.ItfSet.DcSingleAvpiv{iEpochL}.eval(omegaRps)))
                obj.log_ITF_Z(sprintf('(%i) BIAS DC diff',            iEpochL), 'AVolt/IVolt', DC_FREQ_HQ,       @(omegaRps) (obj.Bias.ItfSet.DcDiffAvpiv{iEpochL}.eval(omegaRps)))
                
                obj.log_ITF_Z(sprintf('(%i) BIAS AC diff, low  gain', iEpochL), 'AVolt/IVolt', AC_DIFF_FREQS_HZ, @(omegaRps) (obj.Bias.ItfSet.AcLowGainAvpiv{iEpochL}.eval(omegaRps)))
                obj.log_ITF_Z(sprintf('(%i) BIAS AC diff, high gain', iEpochL), 'AVolt/IVolt', AC_DIFF_FREQS_HZ, @(omegaRps) (obj.Bias.ItfSet.AcHighGainAvpiv{iEpochL}.eval(omegaRps)))
            end
            
        end
        
        
        
        % Load non-BIAS RCTs (all types) using assumptions on filenames.
        %
        % NOTE: Can be useful for manual experimentation with calibration.
        % NOTE: Will only load one of each RCT type (no time dependence as per global attribute CALIBRATION_TABLE).
        %
        function read_non_BIAS_RCTs_by_regexp(obj, use_CALIBRATION_TABLE_INDEX2)
            % ASSERTIONS
            assert(~obj.hasLoadedNonBiasData, 'BICAS:calib:Assertion', 'Can not load non-BIAS data twice.')
            assert(isscalar(use_CALIBRATION_TABLE_INDEX2))
            
            obj.LfrItfIvptTable{1}    = obj.read_log_RCT_by_SETTINGS_regexp('LFR');
            obj.tdsCwfFactorsIvpt{1}  = obj.read_log_RCT_by_SETTINGS_regexp('TDS-CWF');
            obj.TdsRswfItfIvptList{1} = obj.read_log_RCT_by_SETTINGS_regexp('TDS-RSWF');
            
            obj.hasLoadedNonBiasData         = 1;
            obj.use_CALIBRATION_TABLE_rcts   = 0;
            obj.use_CALIBRATION_TABLE_INDEX2 = use_CALIBRATION_TABLE_INDEX2;
            
            obj.log_LFR_ITFs();
            obj.log_TDS_RSWF_ITFs();
        end


        
        % Load non-BIAS RCT(s) of ONE type (rctId) using CDF global attribute CALIBRATION_TABLE and zVar
        % CALIBRATION_TABLE_INDEX.
        %
        % IMPLEMENTATION NOTE: May load multiple RCTs (of the same type) but will only load those RCTs which are
        % actually needed, as indicated by CALIBRATION_TABLE_INDEX. This is necessary since CALIBRATION_TABLE may
        % contain unnecessary RCTs of types not recognized by BICAS (LFR's ROC-SGSE_CAL_RCT-LFR-VHF_V01.cdf
        % /2019-12-16), and which are therefore unreadable by BICAS (BICAS would crash).
        %
        %
        % ARGUMENTS
        % =========
        % ga_CALIBRATION_TABLE       : LFR/TDS RCT global attribute CALIBRATION_TABLE. 1D cell array of strings.
        % zv_CALIBRATION_TABLE_INDEX : LFR/TDS BICAS input dataset zVariable CALIBRATION_TABLE_INDEX.
        %
        %
        % RETURN VALUE
        % ============
        % rctDataList : 1D cell array, of same length as ga_CALIBRATION_TABLE. {iRct},
        %               where iRct=zv_CALIBRATION_TABLE_INDEX(i,1). Each element is the content of the corresponding RCT
        %               mentioned in ga_CALIBRATION_TABLE.
        %
        function rctDataList = read_non_BIAS_RCT_by_CALIBRATION_TABLE(obj, ...
                rctId, ga_CALIBRATION_TABLE, zv_CALIBRATION_TABLE_INDEX, use_CALIBRATION_TABLE_INDEX2)
            
            % ASSERTIONS
            assert(iscell(ga_CALIBRATION_TABLE), 'BICAS:calib:Assertion:IllegalArgument', 'ga_CALIBRATION_TABLE is not a cell array.')
            EJ_library.assert.vector(ga_CALIBRATION_TABLE)
            assert(size(zv_CALIBRATION_TABLE_INDEX, 2) == 2, 'BICAS:calib:Assertion:IllegalArgument', 'zv_CALIBRATION_TABLE_INDEX does not have two columns.')
            assert(isscalar(use_CALIBRATION_TABLE_INDEX2))

            %==================================================================================================
            % Read RCTs, but only those actually NEEDED according to zv_CALIBRATION_TABLE_INDEX
            % ---------------------------------------------------------------------------------
            % NOTE: This is important since not all RCTs in glob.attr. CALIBRATION_TABLE can be read by BICAS!
            %==================================================================================================
            iUniqueList = unique(zv_CALIBRATION_TABLE_INDEX(:,1));
            % List/cell array of data from multiple RCTs (potentially, but practically probably not).
            rctDataList = cell(numel(ga_CALIBRATION_TABLE), 1);
            for i = 1:numel(iUniqueList)
                j              = iUniqueList(i) + 1;   % NOTE: Cell array index is one greater that the stored value.
                rctDataList{j} = obj.read_log_RCT_by_filename(ga_CALIBRATION_TABLE{j}, rctId);
            end
            
            %=================
            % Store read RCTs
            %=================
            switch(rctId)
                case 'LFR'
                    obj.LfrItfIvptTable    = rctDataList;
                case 'TDS-CWF'
                    obj.tdsCwfFactorsIvpt  = rctDataList;
                case 'TDS-RSWF'
                    obj.TdsRswfItfIvptList = rctDataList;
                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal rctId="%s"', rctId)
            end
            
            obj.hasLoadedNonBiasData         = 1;
            obj.use_CALIBRATION_TABLE_rcts   = 1;
            obj.use_CALIBRATION_TABLE_INDEX2 = use_CALIBRATION_TABLE_INDEX2;
            
            obj.log_LFR_ITFs();
            obj.log_TDS_RSWF_ITFs();
        end



        function log_LFR_ITFs(obj)
            
            FREQ_HZ = 0;
            for iLfrRct = 1:numel(obj.LfrItfIvptTable)
                if ~isempty(obj.LfrItfIvptTable{iLfrRct})
                    % CASE: This index corresponds to an actually loaded RCT (some are intentionally empty).
                    for iLsf = 1:4
                        if iLsf ~= 4
                            nBltsMax = 5;
                        else
                            nBltsMax = 3;
                        end
                        
                        for iBlts = 1:nBltsMax
                            TabulatedItfIvpt = obj.LfrItfIvptTable{iLfrRct}{iLsf}{iBlts};
                            ItfIvpt          = @(omegaRps) (bicas.calib.eval_tabulated_ITF(TabulatedItfIvpt, omegaRps));
                            obj.log_ITF_Z(...
                                sprintf('LFR RCT %i, F%i, BLTS/BIAS_%i', iLfrRct, iLsf-1, iBlts), ...
                                'IVolt/TM unit', FREQ_HZ, ItfIvpt)
                        end
                    end    % for
                end    % if
            end    % for
        end
        
        
        
        function log_TDS_RSWF_ITFs(obj)
            
            FREQ_HZ = 0;
            for iTdsRswfRct = 1:numel(obj.TdsRswfItfIvptList)
                for iBlts = 1:3
                    TabulatedItfIvpt = obj.TdsRswfItfIvptList{iTdsRswfRct}{iBlts};                    
                    ItfIvpt          = @(omegaRps) (bicas.calib.eval_tabulated_ITF(TabulatedItfIvpt, omegaRps));
                    obj.log_ITF_Z(...
                        sprintf('TDS RSWF RCT %i, BLTS/BIAS_%i', iTdsRswfRct, iBlts), ...
                        'IVolt/TM unit', FREQ_HZ, ItfIvpt)
                end
            end
        end
        
        
        
        % Convert/calibrate TC bias current: TM units --> physical units.
        %
        % NOTE: This is the normal way of obtaining bias current in physical units (as opposed to HK bias current).
        %
        function biasCurrentAAmpere = calibrate_TC_bias_TM_to_bias_current(obj, biasCurrentTm, iAntenna, iCalibTimeL)
            
            %==============================
            % Obtain calibration constants
            %==============================
            offsetAAmpere = obj.Bias.Current.offsetsAAmpere(iCalibTimeL, iAntenna);
            gainAapt      = obj.Bias.Current.gainsAapt(     iCalibTimeL, iAntenna);

            % CALIBRATE
            biasCurrentAAmpere = offsetAAmpere + gainAapt .* double(biasCurrentTm);    % LINEAR FUNCTION
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
        % NOTE: Walter Puccio, IRF-U 2019-09-06: Values are measured on the order of once per second (and sent back as HK
        % even more rarely). Expect errors on the order of 5%.
        %
        % NOTE: The calibration data are NOT stored in the BIAS RCT.
        %
        % NOTE: The conversion function can be found in the BIAS specification, sections 3.4.4.{1-3} ("BIAS1" etc) under
        % "Telemetry". (Not to be confused with the corresponding telecommands.). The conversion functions are identical
        % for all three probes.
        %
        function biasCurrentAAmpere = calibrate_HK_bias_TM_to_bias_current(obj, biasCurrentTm, iAntenna)
            
            % ASSERTION: zVar HK_BIA_BIAS1/2/3's class in BIAS HK.
            % Not strictly required, but the variable has to be some integer.
            assert(isa(biasCurrentTm, 'uint16'))
            
            %===================================================================================================
            % CALIBRATE
            % ---------
            % Unsigned integer which represents ~signed integer.
            % ==> Intervals 0..0x7FFF and 0x8000...0xFFFF need to "change places".
            % ==> Need to flip bit representing sign to have one interval 0...0xFFFF with monotonic function
            %     TM-to-calibrated values.
            %===================================================================================================
            biasCurrentTm     = bitxor(biasCurrentTm, hex2dec('8000'));    % FLIP BIT
            biasCurrentAAmpere = obj.HkBiasCurrent.gainAapt(iAntenna) * ...
                (biasCurrentTm + obj.HkBiasCurrent.offsetTm(iAntenna));    % LINEAR FUNCTION
        end
        
        
        
        % Calibrate all voltages. Function will choose the more specific algorithm internally.
        function samplesCaAVolt = calibrate_voltage_all(obj, ...
                dtSec, samplesCaTm, isLfr, isTdsCwf, CalSettings, zv_CALIBRATION_TABLE_INDEX)
%                dtSec, samplesCaTm, isLfr, isTdsCwf, iBlts, BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH, iLsf, CALIBRATION_TABLE_INDEX)
            
            % ASSERTIONS
            %EJ_library.assert.struct(CalSettings, {'iBlts', 'BltsSrc', 'biasHighGain', 'iCalibTimeL', 'iCalibTimeH', 'iLsf'}, {})   % Too slow?
            assert(all(size(zv_CALIBRATION_TABLE_INDEX) == [1,2]))
            
            

            % Set cti1, cti2.
            if obj.use_CALIBRATION_TABLE_rcts
                cti1 = zv_CALIBRATION_TABLE_INDEX(1,1);
            else
                cti1 = 0;
            end
            % NOTE: NOT incrementing by one, since the variable's meaning can vary between LFR, TDS-CWF, TDS-RSWF.
            cti2 = zv_CALIBRATION_TABLE_INDEX(1,2);

            
            
            if obj.allVoltageCalibDisabled
                
                samplesCaAVolt = cell(size(samplesCaTm));
                for i = 1:numel(samplesCaTm)
                    samplesCaAVolt{i} = double(samplesCaTm{i});
                end
                
            else

                if isLfr
                    %===========
                    % CASE: LFR
                    %===========
                    samplesCaAVolt = obj.calibrate_LFR_full(dtSec, samplesCaTm, CalSettings, cti1, cti2);
                else
                    %===========
                    % CASE: TDS
                    %===========
                    if isTdsCwf
                        % CASE: TDS CWF
                        samplesCaAVolt = obj.calibrate_TDS_CWF_full(dtSec, samplesCaTm, CalSettings, cti1, cti2);
                    else
                        % CASE: TDS RSWF
                        samplesCaAVolt = obj.calibrate_TDS_RSWF_full(dtSec, samplesCaTm, CalSettings, cti1, cti2);
                    end
                end
                
            end
        end



        % ARGUMENTS
        % =========
        % samplesTm    : 1D cell array of numeric 1D arrays.
        % samplesAVolt : 1D cell array of numeric 1D arrays.
        % CalSettings  : Struct that groups together arguments
        %   .iBlts     : 1..5.
        %   .BltsSrc   : bicas.BLTS_src_dest describing where the signal comes from.
        %   ...
        %
        function samplesCaAVolt = calibrate_LFR_full(obj, dtSec, samplesCaTm, CalSettings, cti1, cti2)
            
            %EJ_library.assert.struct(CalSettings, {'iBlts', 'BltsSrc', 'biasHighGain', 'iCalibTimeL', 'iCalibTimeH', 'iLsf'}, {})   % Too slow?
            iBlts        = CalSettings.iBlts;
            BltsSrc      = CalSettings.BltsSrc;
            biasHighGain = CalSettings.biasHighGain;
            iCalibTimeL  = CalSettings.iCalibTimeL;
            iCalibTimeH  = CalSettings.iCalibTimeH;
            iLsf         = CalSettings.iLsf;
            
            % ASSERTIONS
            assert(iscell(samplesCaTm))
            EJ_library.assert.vector(samplesCaTm)
            EJ_library.assert.vector(dtSec)
            assert(numel(samplesCaTm) == numel(dtSec))
            assert((1 <= iBlts) && (iBlts <= 5))
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            assert((1 <= iLsf)  && (iLsf  <= 4), 'Illegal argument iLsf=%g.', iLsf)
            assert(cti1 >= 0, 'Illegal cti1=%g', cti1)
            % No assertion on cti2 unless used (determined later).
            
            
            
            %============================================
            % Only place to potentially make use of cti2
            %============================================
            if obj.use_CALIBRATION_TABLE_INDEX2
                % ASSERTIONS
                assert(cti2 >= 0,      'BICAS:calib:Assertion', 'Illegal cti2=%g (=zVar CALIBRATION_TABLE_INDEX(iRecord, 2))', cti2)
                assert(iLsf == cti2+1, 'BICAS:calib:Assertion', 'cti2+1=%i != iLsf=%i (before overwriting iLsf)', cti2+1, iLsf)
                
                % NOTE: Only place cti2 is used.
                iLsf = cti2 + 1;
            end



            %==============================
            % Obtain calibration constants
            %==============================
            BiasCalibData = obj.get_BIAS_calib_data(BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            if obj.lfrTdsTfDisabled
                lfrItfIvpt = @(omegaRps) (ones(omegaRps));
            else
                lfrItfIvpt = obj.get_LFR_ITF(cti1, iBlts, iLsf);
            end

            %======================================
            % Create combined ITF for LFR and BIAS
            %======================================
            itfIvpt = @(omegaRps) (...
                lfrItfIvpt(omegaRps) ...
                .* ...
                BiasCalibData.itfAvpiv(omegaRps));

            %=======================================
            % CALIBRATE: LFR TM --> TM --> AVolt
            %=======================================
            samplesCaAVolt = cell(size(samplesCaTm));
            lsfOffsetTm = obj.lfrLsfOffsetsTm(iLsf);
            for i = 1:numel(samplesCaTm)
                
                % ADD LSF OFFSET
                samplesTm = samplesCaTm{i}(:) + lsfOffsetTm;
                
                % APPLY TRANSFER FUNCTION
                tempSamplesAVolt = bicas.utils.apply_transfer_function(...
                    dtSec(i), samplesTm, itfIvpt, ...
                    'enableDetrending',        obj.enableDetrending, ...
                    'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction);

                samplesCaAVolt{i} = tempSamplesAVolt + BiasCalibData.offsetAVolt;
            end
        end



        



        %
        % ARGUMENTS
        % =========
        % See calibrate_LFR_full.
        %
        function samplesCaAVolt = calibrate_TDS_CWF_full(obj, dtSec, samplesCaTm, CalSettings, cti1, cti2)

            %EJ_library.assert.struct(CalSettings, {'iBlts', 'BltsSrc', 'biasHighGain', 'iCalibTimeL', 'iCalibTimeH'}, {'iLsf'})   % Too slow?
            iBlts        = CalSettings.iBlts;
            BltsSrc      = CalSettings.BltsSrc;
            biasHighGain = CalSettings.biasHighGain;
            iCalibTimeL  = CalSettings.iCalibTimeL;
            iCalibTimeH  = CalSettings.iCalibTimeH;
            
            % ASSERTIONS
            EJ_library.assert.vector(dtSec)
            assert(iscell(samplesCaTm))
            assert(numel(samplesCaTm) == numel(dtSec))
            assert((1 <= iBlts) && (iBlts <= 5))
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            assert(cti1 >= 0)
            
            if obj.use_CALIBRATION_TABLE_INDEX2
                % TODO? ASSERTION: cti2 = 0???
                error('BICAS:calib:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
                    'TDS-CWF calibration never uses CALIBRATION_TABLE_INDEX2.')
            end
            
            samplesCaAVolt = cell(size(samplesCaTm));   % Initialize empty output variable.
                
            if ismember(iBlts, [1,2,3])
                % CASE: BLTS 1-3 which TDS does support.
                
                %==============================
                % Obtain calibration constants
                %==============================
                % NOTE: Low/high gain is irrelevant for TDS. Argument value arbitrary.
                BiasCalibData = obj.get_BIAS_calib_data(BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
                
                for i = 1:numel(samplesCaTm)
                    
                    %===============================================
                    % CALIBRATE: TDS TM --> TDS/BIAS interface volt
                    %===============================================
                    if obj.lfrTdsTfDisabled
                        tdsFactorIvpt = 1;
                    else
                        tdsFactorIvpt = obj.tdsCwfFactorsIvpt{cti1+1}(iBlts);
                    end
                    tempSamplesIVolt = tdsFactorIvpt * samplesCaTm{i};    % MULTIPLICATION
                   
                    %=====================================================
                    % CALIBRATE: TDS/BIAS interface volt --> antenna volt
                    %=====================================================
                    % APPLY TRANSFER FUNCTION
                    tempSamplesAVolt = bicas.utils.apply_transfer_function(...
                        dtSec(i), ...
                        tempSamplesIVolt(:), ...
                        BiasCalibData.itfAvpiv, ...
                        'enableDetrending',        obj.enableDetrending, ...
                        'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction);
                    samplesCaAVolt{i} = tempSamplesAVolt + BiasCalibData.offsetAVolt;
                end

            else                
                % CASE: BLTS 4-5 which TDS does NOT support.
                
                for i = 1:numel(samplesCaTm)
                    % Always return NaN.
                    samplesCaAVolt{i} = NaN * samplesCaTm{i};
                end
            end
            
        end



        % ARGUMENTS
        % =========
        % See calibrate_LFR_full.
        %
        function samplesCaAVolt = calibrate_TDS_RSWF_full(obj, dtSec, samplesCaTm, CalSettings, cti1, cti2)
            
            %EJ_library.assert.struct(CalSettings, {'iBlts', 'BltsSrc', 'biasHighGain', 'iCalibTimeL', 'iCalibTimeH'}, {'iLsf'})   % Too slow?
            iBlts        = CalSettings.iBlts;
            BltsSrc      = CalSettings.BltsSrc;
            biasHighGain = CalSettings.biasHighGain;
            iCalibTimeL  = CalSettings.iCalibTimeL;
            iCalibTimeH  = CalSettings.iCalibTimeH;
            
            % ASSERTIONS
            EJ_library.assert.vector(dtSec)
            assert(iscell(samplesCaTm))
            assert(numel(samplesCaTm) == numel(dtSec))
            assert((1 <= iBlts) && (iBlts <= 5))
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            assert(cti1 >= 0)
            
            if obj.use_CALIBRATION_TABLE_INDEX2
                % TODO? ASSERTION: cti2 = 0???
                error('BICAS:calib:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
                    'TDS-RSWF calibration never uses CALIBRATION_TABLE_INDEX2.')
            end
            
            %==============================
            % Obtain calibration constants
            %==============================
            % NOTE: Low/high gain is irrelevant for TDS. Argument value arbitrary.
            BiasCalibData = obj.get_BIAS_calib_data(BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            
            samplesCaAVolt = cell(size(samplesCaTm));   % Initialize empty output variable.
            if ismember(iBlts, [1,2,3])
                
                %======================================
                % Create combined ITF for TDS and BIAS
                %======================================
                if obj.lfrTdsTfDisabled
                    tdsItfIvpt = @(omegaRps) (ones(omegaRps));
                else
                    tdsItfIvpt = @(omegaRps) (bicas.calib.eval_tabulated_ITF(obj.TdsRswfItfIvptList{cti1+1}{iBlts}, omegaRps));
                end

                itf = @(omegaRps) (...
                    tdsItfIvpt(omegaRps) ...
                    .* ...
                    BiasCalibData.itfAvpiv(omegaRps));
                
                %====================================
                % CALIBRATE: TDS TM --> antenna volt
                %====================================
                for i = 1:numel(samplesCaTm)
                    tempSamplesAVolt = bicas.utils.apply_transfer_function(...
                        dtSec(i), ...
                        samplesCaTm{i}(:), ...
                        itf, ...
                        'enableDetrending',        obj.enableDetrending, ...
                        'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction);
                    
                    samplesCaAVolt{i} = tempSamplesAVolt + BiasCalibData.offsetAVolt;
                end
            else
                for i = 1:numel(samplesCaTm)
                    % CASE: BLTS 4-5 which TDS does not support.
                    % Always return NaN.
                    samplesCaAVolt{i} = NaN * samplesCaTm{i};
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



        % NOTE: To be compared with read_log_RCT_by_SETTINGS_regexp.
        % NOTE: Does not need to be an instance method. Is so only to put it next to read_log_RCT_by_SETTINGS_regexp.
        function RctCalibData = read_log_RCT_by_filename(obj, filename, rctId)
            RctCalibData = obj.read_log_modify_RCT(fullfile(obj.calibrationDir, filename), rctId);
        end



        % NOTE: To be compared with read_log_RCT_by_filename.
        function RctCalibData = read_log_RCT_by_SETTINGS_regexp(obj, rctId)
            filePath     = bicas.RCT.find_RCT_by_SETTINGS_regexp(obj.calibrationDir, rctId, obj.SETTINGS, obj.L);
            RctCalibData = obj.read_log_modify_RCT(filePath, rctId);
        end


        
        % Return subset of already loaded BIAS calibration data, for specified settings.
        %
        % NOTE: May return calibration values corresponding to scalar calibration, depending on SETTINGS:
        %
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % kIvpav = Multiplication factor "k" that represents/replaces the (forward) transfer function.
            switch(BltsSrc.category)
                case 'DC single'
                    
                    BiasItfAvpiv = TF_list_2_func(obj.Bias.ItfSet.DcSingleAvpiv);    % NOTE: List of ITFs for different times.
                    kFtfIvpav    = obj.BiasGain.alphaIvpav;
                    offsetAVolt  = obj.Bias.dcSingleOffsetsAVolt(iCalibTimeH, BltsSrc.antennas);
                    
                case 'DC diff'
                    
                    BiasItfAvpiv = TF_list_2_func(obj.Bias.ItfSet.DcDiffAvpiv);
                    kFtfIvpav    = obj.BiasGain.betaIvpav;
                    if     isequal(BltsSrc.antennas(:)', [1,2]);   offsetAVolt = obj.Bias.DcDiffOffsets.E12AVolt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [2,3]);   offsetAVolt = obj.Bias.DcDiffOffsets.E23AVolt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [1,3]);   offsetAVolt = obj.Bias.DcDiffOffsets.E13AVolt(iCalibTimeH);
                    else
                        error('BICAS:calib:Assertion:IllegalArgument', 'Illegal BltsSrc.');
                    end
                    
                case 'AC diff'
                    
                    if     biasHighGain == 0
                        BiasItfAvpiv = TF_list_2_func(obj.Bias.ItfSet.AcLowGainAvpiv);
                        kFtfIvpav    = obj.BiasGain.gammaIvpav.lowGain;
                        offsetAVolt  = 0;
                    elseif biasHighGain == 1
                        BiasItfAvpiv = TF_list_2_func(obj.Bias.ItfSet.AcHighGainAvpiv);
                        kFtfIvpav    = obj.BiasGain.gammaIvpav.highGain;
                        offsetAVolt  = 0;
                    elseif isnan(biasHighGain)
                        BiasItfAvpiv = bicas.calib.NAN_TF;
                        kFtfIvpav    = NaN;
                        offsetAVolt  = NaN;
                    else
                        error('BICAS:calib:Assertion:IllegalArgument', 'Illegal argument biasHighGain=%g.', biasHighGain)
                    end

                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', ...
                        'Illegal argument BltsSrc.category=%s. Can not obtain calibration data for this type of signal.', ...
                        BltsSrc.category)
            end
            
            if obj.biasOffsetsDisabled && ~isnan(offsetAVolt)
                offsetAVolt = 0;
            end
            if obj.useBiasTfScalar
                BiasItfAvpiv = @(omegaRps) (ones(size(omegaRps)) / kFtfIvpav);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            BiasCalibData.itfAvpiv    = BiasItfAvpiv;
            BiasCalibData.offsetAVolt = offsetAVolt;



            % Convenience function: Cell array of TFs --> Time-relevant TF (as function handle)
            %
            % RFTF = Rational Function (rational_func_transform object) Transfer Function
            function Tf = TF_list_2_func(RftfList)
                Tf = @(omegaRps) (RftfList{iCalibTimeL}.eval(omegaRps));
            end
        end
        
        
        
        % Obtain LFR ITF, but handle the case that should never happen for actual non-NaN data (LSF F3 + BLTS 4 or 5)
        % and return a TF that only returns NaN instead. BICAS may still iterate over that combination though when
        % calibrating.
        % 
        function lfrItfIvpt = get_LFR_ITF(obj, cti1, iBlts, iLsf)
            assert(cti1 >= 0)
            assert(ismember(iBlts, [1:5]))
            assert(ismember(iLsf,  [1:4]))
            assert(logical(obj.hasLoadedNonBiasData))
            
            if (iLsf == 4) && ismember(iBlts, [4,5])
                lfrItfIvpt = bicas.calib.NAN_TF;
            else
                lfrItfIvpt = @(omegaRps) (bicas.calib.eval_tabulated_ITF(obj.LfrItfIvptTable{cti1+1}{iLsf}{iBlts}, omegaRps));
            end
        end



        % NOTE: Requires obj.L to be initialized.
        %
        % ARGUMENTS
        % =========
        % freqHzArray : Array of frequencies for which the ITF value should be logged.
        % TfFuncPtr   : Function pointer. Z(omegaRps).
        %
        function log_ITF_Z(obj, itfName, itfUnit, freqHzArray, itfFuncPtr)
            
            zArray = itfFuncPtr(freqHzArray);
            prefixStr = sprintf('Inverse TF, %s,', itfName);    % NOTE: Includes final comma.
            for i=1:numel(freqHzArray)
                freqHz = freqHzArray(i);
                Z      = zArray(i);
                
                inverseZValueStr = sprintf('1/%8.5f', 1/abs(Z));
                
                %======================================================================================================
                % NOTE 2020-04-30: Execution at ROC fails due to not finding function "phase" for unknown reason.
                % --------------------------------------------------------------------------------------
                % Exception.identifier = "MATLAB:UndefinedFunction"
                % Exception.message    = "Undefined function 'phase' for input arguments of type 'double'."
                % Matching MATLAB error message identifier parts (error types derived from Exception1.identifier):
                %     UntranslatableErrorMsgId : Error occurred, but code can not translate the error's MATLAB message identifier into any of BICAS's internal standard error codes.
                % MATLAB call stack:
                %     row  969, calib.m,                    calib.log_ITF_Z
                %     row  303, calib.m,                    calib.calib
                %     row   68, execute_sw_mode.m,          execute_sw_mode
                %     row  455, main.m,                     main_without_error_handling
                %     row  116, main.m,                     main
                % --------------------------------------------------------------------------------------
                % See also
                % https://se.mathworks.com/matlabcentral/answers/408657-which-toolbox-is-phase-in
                % """"phase() as a routine by itself is part of the System Identification Toolbox, in the "obsolete" category.
                % phase() is also a method of the newer iddata() class from the System Identification Toolbox.
                % But what you probably want is angle() followed by unwrap(), which is part of basic MATLAB.""""
                %
                % Therefore using "angle" instead of "phase".
                %======================================================================================================
                obj.L.logf('debug', '%-41s %4i Hz: abs(Z)=%8.5f=%12s [%s], phase(Z)=%5.1f [deg]', ...
                    prefixStr, freqHz, abs(Z), inverseZValueStr, itfUnit, rad2deg(angle(Z)))
                
                % Do not print prefix more than once.
                prefixStr = '';
            end
        end
        


        % Read any single RCT file, and log it. Effectively wraps the different RCT-reading functions.
        % 
        %
        % IMPLEMENTATION NOTES
        % ====================
        % This method exists to
        % (1) run shared code that should be run when reading any RCT (logging, modifying data),
        % (2) separate logging from the RCT-reading code, so that one can read RCTs without BICAS.
        %
        %
        % ARGUMENTS
        % =========
        % rctId : String constants representing pipeline and RCT to be read.
        %
        function RctCalibData = read_log_modify_RCT(obj, filePath, rctId)
            % PROPOSAL: Incorporate modify_LFR_data, modify_TDS_RSWF_data
            
            obj.L.logf('info', 'Reading %-4s RCT: "%s"', rctId, filePath)

            switch(rctId)
                case 'BIAS'     ; Rcd =                                  bicas.RCT.read_BIAS_RCT(    filePath);
                case 'LFR'      ; Rcd = bicas.calib.modify_LFR_data(     bicas.RCT.read_LFR_RCT(     filePath));
                case 'TDS-CWF'  ; Rcd =                                  bicas.RCT.read_TDS_CWF_RCT( filePath);
                case 'TDS-RSWF' ; Rcd = bicas.calib.modify_TDS_RSWF_data(bicas.RCT.read_TDS_RSWF_RCT(filePath));
                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', 'Illegal rctId="%s"', rctId);
            end
            RctCalibData = Rcd;
        end



    end    % methods(Access=private)

    %###################################################################################################################

    
    
    methods(Static, Access=public)


        
        % Given a sequence of Epoch values, determine for each value which calibration time index should be used. The
        % caller will have to decide which sequence of data that should be calibrated together (e.g. if calibration time
        % changes in the middle of CWF), and which Epoch values should be used to determine calibration time (e.g. first
        % Epoch value for a snapshot determines entire snapshot).
        %
        % NOTE: Method is public so that automatic test code can call get_calibration_time.
        %
        %
        % ARGUMENTS AND RETURN VALUES
        % ===========================
        % Epoch          : Column vector with Epoch values.
        % CalibEpochList : List of monotonically increasing timestamps ("Epoch format").
        %                  In practice intended to be Bias.epochL or Bias.epochH.
        % iCalib         : Array. iCalibList(i) = calibration time index for Epoch(i).
        %
        function [iCalib] = get_calibration_time(Epoch, CalibEpochList)
            
            % ASSERTIONS
            bicas.proc_utils.assert_zv_Epoch(Epoch)
            bicas.proc_utils.assert_zv_Epoch(CalibEpochList)
            % IMPLEMENTATION NOTE: Does not work if CalibEpochList is empty, since discretize behaves differently for
            % scalar second argument.
            assert(~isempty(CalibEpochList))
            
            % IMPLEMENTATION NOTE: "discretize" by itself returns NaN for Epoch values outside the outermost edges.
            % Therefore (1) must add upper edge "Inf", (2) asserts non-Nan afterwards.
            % IMPLEMENTATION NOTE: "discretize" behaves differently for scalar second argument. Adding edges at infinity hides
            % this problem. If one does not add infinities and uses a scalar edge list, then one has to treat those
            % cases manually.
            iCalib = discretize(Epoch, [CalibEpochList; Inf], 'IncludedEdge', 'left');
            assert(all(~isnan(iCalib(:))), ...
                'BICAS:calib:SWModeProcessing', ...
                'Can not derive which calibration data to use for all specified timestamps.')
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
        
        
        
        % Convert "set current" to TC/TM units.
        %
        function biasCurrentTm = calibrate_set_current_to_bias_current(currentSAmpere)
            % PROPOSAL: Make into static function.
            %   PRO: Easier to use from outside of BICAS.
            
            MAX_ABS_SET_CURRENT_SAMPERE = 60e-6;
            
            % ASSERTION
            % NOTE: max(...) ignores NaN, unless that is the only value, which then becomes the max value.
            [maxAbsSAmpere, iMax] = max(abs(currentSAmpere(:)));
            if ~(isnan(maxAbsSAmpere) || (maxAbsSAmpere <= MAX_ABS_SET_CURRENT_SAMPERE))
                
                error('BICAS:calib:Assertion:IllegalArgument', ...
                    'Argument currentSAmpere (unit: set current/ampere) contains illegally large value(s). Largest found value is %g.', ...
                    currentSAmpere(iMax))
            end
            
            biasCurrentTm = currentSAmpere * 32768/60e-6;
        end



    end    % methods(Static, Access=public)



    methods(Static, Access=private)



        function LfrItfIvptTable = modify_LFR_data(LfrItfIvptTable)
            % Modify tabulated LFR TFs.
            for iLsf = 1:numel(LfrItfIvptTable)
                for iBlts = 1:numel(LfrItfIvptTable{iLsf})
                    LfrItfIvptTable{iLsf}{iBlts} = bicas.calib.modify_tabulated_ITF(LfrItfIvptTable{iLsf}{iBlts});
                end
            end
        end
        
        
        
        function TdsRswfItfIvptList = modify_TDS_RSWF_data(TdsRswfItfIvptList)
            % Modify tabulated TDS-RSWF TFs.
            for iBlts = 1:numel(TdsRswfItfIvptList)
                TdsRswfItfIvptList{iBlts} = bicas.calib.modify_tabulated_ITF(TdsRswfItfIvptList{iBlts});
            end
        end



        % Evaluate tabulated INVERSE transfer function.
        %
        % NOTE: This function is effectively meant to specify how tabulated transfer functions should be interpreted w.r.t.
        % interpolation.
        % NOTE: Intended specifically for INVERSE transfer functions. Therefore using Z=0 for frequencies higher than
        % the table.
        %
        function Z = eval_tabulated_ITF(TabulatedTf, omegaRps)
            useTabTf = (omegaRps <= TabulatedTf.omegaRps(end));
            
            % NOTE: interp1 return NaN for values outside range.
            Z = interp1(TabulatedTf.omegaRps, TabulatedTf.Z, omegaRps, 'linear');
            Z(~useTabTf) = 0;   % Set to zero (overwrite) for values above highest tabulated frequency.
            
            % ASSERTION
            if ~all(isfinite(Z))
                % IMPLEMENTATION NOTE: Experience shows that it is useful to have an extended error message confirming
                % that the requested frequence range is outside the tabulated one, and by how much.
                errorMsg = sprintf(...
                    ['Can not evaluate tabulated transfer function for frequencies outside of the range of tabulated frequencies.\n', ...
                    'Range of frequencies for which there are tabulated Z values:\n', ...
                    '    min(TabulatedTf.omegaRps) = %g\n', ...
                    '    max(TabulatedTf.omegaRps) = %g\n', ...
                    'Range of frequencies for which evaluation (interpolation) of Z was attempted:\n', ...
                    '    min(omegaRps)     = %g\n', ...
                    '    max(omegaRps)     = %g\n'], ...
                    min(TabulatedTf.omegaRps), ...
                    max(TabulatedTf.omegaRps), ...
                    min(omegaRps), ...
                    max(omegaRps));
                
                error('BICAS:Assertion', errorMsg)
            end
        end
        
        
        
        % Modify tabulated INVERSE transfer functions, if needed.
        % Tabulated TF --> Tabulated TF
        %
        % Extrapolate to 0 Hz, if needed.
        %
        function ModifItf = modify_tabulated_ITF(Itf)
            assert(Itf.omegaRps(1) > 0)

            % NOTE: Can not just use the lowest-frequency Z value for 0 Hz since it has to be real (not complex).
            Z1       = Itf.Z(1);
            signZ0   = sign(real(Z1));
            assert(signZ0 ~= 0, ...
                'BICAS:calib:modify_tabulated_ITF:FailedToReadInterpretRCT:Assertion', ...
                'Can not extrapolate tabulated inverse transfer function (ITF) to zero Hz due to ambiguity. real(Z(1)) = 0.')
            Z0       = abs(Z1) * signZ0;   % Z value at 0 Hz.
            
            omegaRps = [0;  Itf.omegaRps(:)];
            Z        = [Z0; Itf.Z(:)       ];
            
            ModifItf = EJ_library.utils.tabulated_transform(omegaRps, Z);
        end
        
        
        
    end    % methods(Static, Access=private)

end
