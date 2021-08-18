%
% Class for
% (1) library/utility functions that calibrate data.
% (2) helper functions to find and load calibration data from file, as needed by
%     BICAS.
% An instance of this class contains
%   (1) relevant settings (loaded from SETTINGS) on how to calibrate data, and
%   (2) calibration data.
% An instance may or may not contain calibration data for __ALL__ types of
% data/RCTs depending on how it was initialized.
%
% NOTE: RCT reading functions assume that the same type of RCT (BIAS, LFR,
% TDS-CWF or TDS-RSWF) is identical (in all relevant parts) for both the RODP
% and ROC-SGSE pipeline.
%
%
% SHORTCOMINGS
% ============
% Does not implement parasitic capacitance yet due to lack of calibration values
% (at least).
%
%
% IMPLEMENTATION NOTES
% ====================
% Class is implemented with extra care, and therefore
% - has extra many assertions
% - is extra careful to include units in identifiers
% - is extra careful to use well defined terms/shortenings/naming conventions
% since
% (1) undiscovered calibration bugs could be considered extra bad,
% (2) it is expected to be hard to detect certain bugs,
% (3) it could be hard to use automatic testing here,
% (4) to detect changing RCT formats, in particular in RCTS from non-BIAS teams.
% --
% NOTE: All calibration functions of measured data are assumed to accept data
% from all BLTS (1-5), i.e. including TDS, in order to reduce the number
% assumptions that the calling code needs to make.
%
%
% DESIGN INTENT
% =============
% If needed, then this class modifies the calibration data read from RCTs, e.g.
% inverts FTFs, so that RCT-reading code (bicas.RCT) does not need to (it should
% not).
%
%
% DEFINITIONS, NAMING CONVENTIONS
% ===============================
% Offset = Value (constant) that is ADDED to (not subtracted from) a measured
%          value during the calibration process.
% CA     = Cell Array
% --
% LSF    = LFR Sampling Frequency (F0...F3)
%          NOTE: When used as a variabe (array index), 1=F0, ..., 4=F3.
% TF     = Transfer function (Z=Z(omega), i.e. in frequency domain)
% FTF    = Forward Transfer Function = TF that describes physical
%          INPUT-to-OUTPUT (not the reverse)
% ITF    = Inverse Transfer Function = TF that describes physical
%          OUTPUT-to-INPUT (not the reverse)
% CTI    = CALIBRATION_TABLE_INDEX (zVar)
% CTI1   = First  value in a record of zVar CALIBRATION_TABLE_INDEX.
% CTI2   = Second value in a record of zVar CALIBRATION_TABLE_INDEX.
% RCTS   = RCT CALIBRATION_TABLE (glob.attr)+CALIBRATION_TABLE_INDEX (zVar).
%          S = plural, RCT= RPW Calibration Table (ROC acronym).
%
%
% UNITS / TYPES OF QUANTITIES
% ---------------------------
% TM         = Telemetry units (in LFR/TDS ADC), or telecommand (TC) units. Using
%              this instead of the term "count".
% IV=ivolt   = Interface Volt = Calibrated volt at the interface between BIAS and
%              LFR/TDS.
% AV=avolt   = Antenna Volt = Calibrated volt at the antennas, i.e. the final
%              calibrated (measured) value, including for reconstructed signals
%              (e.g. diffs calculated from singles). May also refer to offsets
%              and values without offsets.
% AA=aampere = Antenna ampere = Calibrated ampere at the antenna.
% sampere    = Set current ampere. Exactly proportional to bias current in TM.
% TPIV       = TM/interface volt (=TM per interface volt)
% IVPT       = Interface volt/TM
% AAPT       = Antenna ampere/TM
% AVPIV      = Antenna volt/interface volt
% IVPAV      = Interface volt/antenna volt
% Deg        = Degrees (angle). 1 revolution=360 degrees=2*pi radians.
% RPS        = Radians/second
% Sec        = Seconds
% 
%
% BLTS = BIAS-LFR/TDS SIGNAL
% ---------------------------
% Signals somewhere between the LFR/TDS ADCs and the non-antenna side of the
% BIAS demuxer including the BIAS transfer functions. Like BIAS_i, i=1..5, but
% includes various stages of calibration/non-calibration, including in
% particular
%   - TM units (inside LFR/TDS),
%   - Interface volt (at the physical boundary BIAS-LFR/TDS (BIAS_i)), and
%   - Calibrated values inside BIAS but without demuxer addition and subtraction
%     inside BIAS (i.e. including using BIAS offsets, BIAS transfer functions;
%     volt).
% NOTE: Definition is partly created to avoid using term "BIAS_i" since it is
% easily confused with other things (the subsystem BIAS, bias currents), partly
% to include various stages of calibration.
%
%
% ASR = Antenna Signal Representation
% -----------------------------------
% The "physical antenna signals" which BIAS-LFR/TDS is trying to measure, or a
% measurement thereof. In reality, the terminology is:
% ASR         : Pointer to a specific physical antenna signal, e.g. V12_LF (DC
%               diff, antenna 1-2)
% ASR samples : Samples representing a specific ASR (as opposed to BLTS).
% NOTE: There are 9 ASRs, i.e. they can refer also to signals not represented by
% any single BLTS, given a chosen mux mode (and latching relay setting).
%
%
% BIAS_i, i=1..5
% --------------
% Defined in BIAS specifications document. Equal to the physical signal at the
% physical boundary between BIAS and LFR/TDS. Unit: LFR/TDS calibrated volt.
% Mostly replaced by BLTS+specified unit in the code.
%
%
% REMINDER: HOW CALIBRATION_TABLE & CALIBRATION_TABLE_INDEX WORK
% ==============================================================
% CALIBRATION_TABLE       : CDF L1R global attribute
%   """"Filename of the calibration table(s).""""
%   """"There must as many as entries than the number of calibration table files
%   associated to the L1R file.""""
% CALIBRATION_TABLE_INDEX : CDF L1R zVariable
%   """"Index of the calibration table(s) value required to generate L2 data
%   files.""""
%   """"Each CDF record must contain 2 elements: the first element must gives
%   the index of the associated CALIBRATION_TABLE entry (i.e., 0 for the first
%   entry, 1 for the second, etc.). The second element must refer to the index
%   of the value to be used inside the calibration table file.""""
%
% Source: ROC-PRO-DAT-NTT-00006-LES_Iss01_Rev02(ROC_Data_Products).Draft2020-04-06.pdf
%
% Summary
% -------
% CALIBRATION_TABLE{CALIBRATION_TABLE_INDEX{iRecord, 1} + 1}
%     == RCT filename
% CALIBRATION_TABLE_INDEX{iRecord, 2}
%     == Index/pointer to some calibration value(s) to use in RCT. Exact
%        interpretation depends on RCT.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-15
%
classdef cal < handle

% BOGIQ:
% ------
% TODO-NI: Where does the parasitic capacitance TF fit into the calibration formulas?
% TODO-NI: Parasitic capacitance values?
% PROPOSAL: Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
%
% ~DOCUMENTATION BUG?!!:
%   PROPOSAL: Abolish ASR. Define acronyms for
%       (1) Antenna signals (AC, DC, singles, diffs)
%       (2) All possible sources of signals of which (1) is a subset (a diff counts as a source).
%           PROPOSAL: Be able to use for both physical signal sources, and for where to place in dataset.
%               PRO: Can use acronym for both, and for class bicas.proc.L1L2.BLTS_src_dest.
%           PROPOSAL: PSS  = Physical Signal Source
%           PROPOSAL: PS   = Physical Signal
%           PROPOSAL: PSSD = Physical Signal Source or Destination
%           PROPOSAL: PSSR = Physical Signal Source or Representation
%           PROPOSAL: PSSR = Physical Signal Source or Dataset Representation
%
% PROPOSAL: Other name for bicas.proc.L1L2.BLTS_src_dest that does not reference BLTS.
%   PRO: Reference to BLTS is confusing.
%   PROPOSAL: Define acronym for all physical signal sources which is a superset of ASR.
%   PROPOSAL: Have different classes and acronyms for
%             (1) physical signal sources, and
%             (2) dataset representation
%             ("BLTS src" and "BLTS dest") where (2) is in practice a subset of (1).
%
% PROPOSAL: Assertion function for CalSettings.
%   TODO-NI: Same struct, with same fields in all cases?
%   NOTE: Function does not know which fields are actually used.
% PROPOSAL: Class for CalSettings. Mere container for fields.
% PROPOSAL: Class for init_RCT_TYPES_MAP:"Entry" structs.
%
% TODO-DEC: How distribute the calibration formulas/algorithms between
%   (1) calibrate_* functions, 
%   (2) functions that select relevant calibration data (get_BIAS_calib_data)
%   (2) RCT reading functions,
%   (3) classes that store TFs?
%   --
%   Ex: Invert the (rat.func., tabulated) TFs
%   Ex: Extrapolate the tabulated TFs to zero
%   Ex: Extrapolate the tabulated LFR TF to higher frequencies.
%   Ex: Modify AC TF at lower frequencies (change amplitude, keep phase)
%   Ex: Interpolation algorithm of tabulated TFs.
%   Ex: If one modifies the TFs before applying the inverse (hypothetical; not implemented)
%   --
%   NEED: Plot all TFs.
%   NEED: Plot all versions of any particular TF as it gets modified.
%   NEED: Plot all TFs used for a particular calibration case
%         (when calibrating using bicas.caib only, without BICAS).
%   PROPOSAL: Separate modifications/choices/code that
%       (1) only need to be done/run once during the execution: modification of calibration data,
%       (2) are done every calibration (call to calibrate_*):   exact algorithms/formulas through which data is fed
%   PROPOSAL: read_*_RCT should not modify any calibration data, just store it:
%             Not invert TFs, not extrapolate TFs to 0 Hz.
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
%       TODO-DEC: Where is it natural to modify calibration data then?!
%   PROPOSAL: General philosophy should be that calibrate_* chooses as much as
%             possible, and thus chooses different functions to call.
%
%
%
% PROPOSAL: Store all versions of TFs internally.
%   Ex: FTF, ITF, tabulated ITF with extrapolation+interpolation+modification
%   PRO: Useful for debugging. Can easily inspect & plot FTFs.
%   NOTE: BIAS & LFR RCTs contain FTFs, TDS RCT contains ITFs.
%   NOTE: Has to keep track of FTF/ITF before modifications (extrapolation
%         to 0 Hz, Z=0 for high freq.).
%   NOTE: Modification (besides inversion) happens on the combined function
%         which is not stored beforehand.
%   NOTE: The set of BIAS+LFR+TDS TFs is different from the set of TFs actually
%         used (combinations of BIAS+LFR and BIAS+TDS respectively).
%   PROPOSAL: Store LFR TFs as one 1D array of structs with fields: iLsf, iBlts, ftf, itf, ...
%       PRO: Can easily iterate over.
%       PRO: For every modification of TFs, can easily add another field for the old
%            version.
%       NOTE/CON: All structs/TFs must have the same fields if true struct array.
%
% PROPOSAL: Separate function for interpolation (no extrapolation) of tabulated TF.
%
% PROPOSAL: Do not expose internal calibration data structs and let the caller
%           specify indices. Access via methods that ask for indices.
%   PRO: Methods document the proper use of indices.
%
%
%
% PROPOSAL: Have calibrate_LFR(), calibrate_TDS_CWF(), calibrate_TDS_RSWF() return
%           function handle to a transfer function that does everything,
%           including offsets.
%   CON: Useful for debugging.
%       CON: Aren't the calibration methods available already that, if turned
%            into function handles?!!
%   PROPOSAL: Have ~special functions/methods for this, so that one does not use
%             function handles wrapped in function handles (not too many
%             anyway).
%   NOTE: CLARIFICATION: Separate TFs with offsets and calibration methods in the
%         time domain.
%   PROPOSAL: Have special function that returns transfer function for every
%             case.
%       PRO: Useful for debugging.
%           PRO: Plot
%           PRO: Save to file.
%
% PROPOSAL: Be able to return values for glob.attrs.
%       CAL_ENTITY_NAME,
%       CAL_ENTITY_AFFILIATION,
%       CAL_EQUIPMENT
%       since on per RCT (CALIBRATION_TABLE entries).



    properties(Access=private, Constant)
        
        % Local constant TF for convenience.
        NAN_TF = @(omegaRps) (omegaRps * NaN);
        
    end



    properties(SetAccess=private, GetAccess=public)
        
        %==================
        % Calibration data
        %==================
        
        % RCT calibration data
        % --------------------------------------------------
        % containers.Map: RCT Type ID --> Data
        % For BIAS, data is a struct (only one BIAS RCT is loaded).
        % For non-BIAS, data is a 1D cell array. {iRct}.
        % iRct-1 corresponds to ga. CALIBRATION_TABLE. and zv.
        % CALIBRATION_TABLE_INDEX(:,1) when those are used. May thus contain
        % empty cells for non-BIAS RCTs which should not (and can not) be
        % loaded.
        RctDataMap = containers.Map();
        
        % Non-RCT calibration data
        % ------------------------
        % BIAS scalar (simplified) calibration, not in the RCTs. For
        % debugging/testing purposes.
        BiasScalarGain
        HkBiasCurrent
        % EXPERIMENTAL. Remove functionality?!
        % NOTE: Technically, the name contains a tautology (LFR+LSF).
        lfrLsfOffsetsTm = [];
        
        
        
        %==================================================
        % Settings for what kind of calibration to perform
        %==================================================
        
        % Corresponds to SETTINGS key-value.
        dcDetrendingDegreeOf
        dcRetrendingEnabled
        acDetrendingDegreeOf
        itfHighFreqLimitFraction
        itfAcConstGainLowFreqRps
                
        % Whether to select non-BIAS RCT using global attribute
        % CALIBRATION_TABLE (and CALIBRATION_TABLE_INDEX(iRecord,1)).
        use_CALIBRATION_TABLE_rcts        
        % Whether to use CALIBRATION_TABLE_INDEX(iRecord,2) for calibration.
        use_CALIBRATION_TABLE_INDEX2

        % What type of calibration to use.
        allVoltageCalibDisabled    % Use TM values (not set to NaN).
        useBiasTfScalar
        biasOffsetsDisabled
        lfrTdsTfDisabled
        
    end
    
    
    
    %###########################################################################



    methods(Access=public)



        % Constructor.
        %
        %
        % NOTES ON INTENDED USAGE
        % =======================
        % The nominal use is that 
        % (1) the user first selects non-BIAS RCTs by initializing a
        % containers.Map using either
        %   (1a) static helper method "find_read_non_BIAS_RCTs_by_regexp", 
        %   (1b) static helper method "find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE",
        % or
        %   (1c) manually (for debugging/analysis/testing).
        % Uses that containers.Map object to call the constructor.
        %
        %
        % IMPLEMENTATION NOTE
        % ===================
        % The class (instance methods, including constructor) deliberately does
        % not itself determine which non-BIAS RCTs to read and therefore does
        % not READ the non-BIAS RCTs. This is useful since
        % ** it makes it possible to inspect & modify the RCT content before
        %    submitting it to bicas.proc.L1L2.cal
        % ** it completely separates the RCT-reading from the class
        %    (modularization)
        % ** it simplifies the constructor somewhat since it does not need to
        %    translate paths into RCT data for non-empty cells.
        %
        function obj = cal(...
                NonBiasRctDataMap, rctDir, ...
                use_CALIBRATION_TABLE_rcts, ...
                use_CALIBRATION_TABLE_INDEX2, ...
                SETTINGS, L)

            % ASSERTIONS: Arguments
            assert(isscalar(use_CALIBRATION_TABLE_INDEX2))
            assert(~ismember('BIAS', NonBiasRctDataMap.keys))
            EJ_library.assert.subset(...
                NonBiasRctDataMap.keys, ...
                bicas.proc.L1L2.cal_RCT.RCT_TYPES_MAP.keys)
            
            
            
            %====================
            % Set obj.RctDataMap
            %====================
            obj.RctDataMap = containers.Map();
            % ------------------------------------
            % Add RCT data map entry for BIAS RCT
            % ------------------------------------
            filenameRegexp = SETTINGS.get_fv(...
                bicas.proc.L1L2.cal_RCT.RCT_TYPES_MAP('BIAS').filenameRegexpSettingKey);
            filePath       = bicas.proc.L1L2.cal_RCT.find_RCT_regexp(rctDir, filenameRegexp, L);
            % IMPLEMENTATION NOTE: It has been observed that value sometimes
            % survives from previous runs, despite being an instance variable.
            % Unknown why. Therefore explicitly overwrites it.
            obj.RctDataMap('BIAS') = bicas.proc.L1L2.cal_RCT.read_RCT_modify_log(...
                'BIAS', filePath, L);
            % -------------------------------------------
            % Add RCT data map entries for non-BIAS RCTs
            % -------------------------------------------
            for rctTypeId = NonBiasRctDataMap.keys
                obj.RctDataMap(rctTypeId{1}) = NonBiasRctDataMap(rctTypeId{1});
            end
            
            
            
            %=========================================================
            % Store miscellaneous SETTINGS key values for convenience
            %=========================================================
            obj.BiasScalarGain.alphaIvpav          = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.ALPHA_IVPAV');
            obj.BiasScalarGain.betaIvpav           = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.BETA_IVPAV');
            obj.BiasScalarGain.gammaIvpav.highGain = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.HIGH_GAIN');
            obj.BiasScalarGain.gammaIvpav.lowGain  = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.LOW_GAIN');
            
            obj.HkBiasCurrent.offsetTm             = SETTINGS.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.OFFSET_TM');
            obj.HkBiasCurrent.gainAapt             = SETTINGS.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.GAIN_AAPT');
            
            obj.lfrLsfOffsetsTm                    = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR.LSF_OFFSETS_TM');

            obj.dcDetrendingDegreeOf               = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE');
            obj.dcRetrendingEnabled                = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED');
            obj.acDetrendingDegreeOf               = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE');
            obj.itfHighFreqLimitFraction           = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF_HIGH_FREQ_LIMIT_FRACTION');
            % NOTE: Converts Hz-->rad/s
            obj.itfAcConstGainLowFreqRps           = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF.AC_CONST_GAIN_LOW_FREQ_HZ') * 2*pi;
            
            obj.allVoltageCalibDisabled            = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.DISABLE');
            obj.biasOffsetsDisabled                = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED');
            obj.lfrTdsTfDisabled                   = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR_TDS.TF_DISABLED');
            
            %-------------------------
            % Set obj.useBiasTfScalar
            %-------------------------
            settingBiasTf                          = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF');
            switch(settingBiasTf)
                case 'FULL'
                    obj.useBiasTfScalar = 0;
                case 'SCALAR'
                    obj.useBiasTfScalar = 1;
                otherwise
                    error(...
                        'BICAS:Assertion:ConfigurationBug', ...
                        ['Illegal value for setting',...
                        ' PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF="%s".'], ...
                        settingBiasTf)
            end
            
            
            
            %============================
            % Store some argument values
            %============================
            obj.use_CALIBRATION_TABLE_rcts   = use_CALIBRATION_TABLE_rcts;
            obj.use_CALIBRATION_TABLE_INDEX2 = use_CALIBRATION_TABLE_INDEX2;
        end



        % Convert/calibrate TC bias current: TM units --> physical units.
        %
        % NOTE: This is the normal way of obtaining bias current in physical
        % units (as opposed to HK bias current).
        %
        function biasCurrentAAmpere = calibrate_current_TM_to_aampere(obj, ...
                biasCurrentTm, iAntenna, iCalibTimeL)
            
            %==============================
            % Obtain calibration constants
            %==============================
            offsetAAmpere = obj.RctDataMap('BIAS').Current.offsetsAAmpere(iCalibTimeL, iAntenna);
            gainAapt      = obj.RctDataMap('BIAS').Current.gainsAapt(     iCalibTimeL, iAntenna);

            % CALIBRATE
            %
            % LINEAR FUNCTION
            biasCurrentAAmpere = offsetAAmpere + gainAapt .* double(biasCurrentTm);
        end



        % Convert/calibrate diagnostic HK TM bias current values to physical
        % units. Refers to BIAS HK zVars HK_BIA_BIAS1/2/3.
        %
        %
        % NOTES
        % =====
        % IMPORTANT NOTE: The HK bias current values are measured onboard but
        % are only meant as DIAGNOSTIC values, NOT AS THE PROPER BIAS CURRENT
        % values for nominal use. Therefore the values should only be seen as
        % approximate.
        %
        % NOTE: Walter Puccio, IRF-U 2019-09-06: Values are measured on the
        % order of once per second (and sent back as HK even more rarely).
        % Expect errors on the order of 5%.
        %
        % NOTE: The calibration data are NOT stored in the BIAS RCT.
        %
        % NOTE: The conversion function can be found in the BIAS specification,
        % sections 3.4.4.{1-3} ("BIAS1" etc) under "Telemetry". (Not to be
        % confused with the corresponding telecommands.). The conversion
        % functions are identical for all three probes.
        %
        function biasCurrentAAmpere = calibrate_HK_bias_TM_to_bias_current(obj, ...
                biasCurrentTm, iAntenna)
            
            % ASSERTION: zVar HK_BIA_BIAS1/2/3's class in BIAS HK.
            % Not strictly required, but the variable has to be some integer.
            assert(isa(biasCurrentTm, 'uint16'))
            
            %=============================================================
            % CALIBRATE
            % ---------
            % Unsigned integer which represents ~signed integer.
            % ==> Intervals 0..0x7FFF and 0x8000...0xFFFF need to
            %     "change places".
            % ==> Need to flip bit representing sign to have one interval
            %     0...0xFFFF with monotonic function for TM-to-calibrated
            %     values.
            %=============================================================
            biasCurrentTm      = bitxor(biasCurrentTm, hex2dec('8000'));    % FLIP BIT
            biasCurrentAAmpere = obj.HkBiasCurrent.gainAapt(iAntenna) * ...
                (biasCurrentTm + obj.HkBiasCurrent.offsetTm(iAntenna));    % LINEAR FUNCTION
        end
        
        
        
        % Calibrate all voltages. Function will choose the more specific
        % algorithm internally.
        %
        % ARGUMENTS
        % =========
        % voltageNaN
        %       Scalar logical.
        %       Whether to set output voltages to NaN and not execute any (real)
        %       calibration.
        %       RATIONALE: This option is useful to
        %       (1) potentially speed up BICAS when it is known that
        %           data will be overwritten with fill values later.
        %       (2) avoid executing calibration algorithms when it is
        %           known that there is no calibration configuration anyway
        %           Ex: LFR zVar BW=0 ==> CALIBRATION_TABLE_INDEX(1,:) is illegal.
        %               ==> Can not calibrate.
        %           Note: This means that this function technically accepts
        %           an illegal calibration configuration when argument is set
        %           to true.
        %
        function samplesCaAVolt = calibrate_voltage_all(obj, ...
                dtSec, samplesCaTm, isLfr, isTdsCwf, CalSettings, ...
                zv_CALIBRATION_TABLE_INDEX, voltageNaN)
            
            % ASSERTIONS
            assert(isstruct(CalSettings))
%             EJ_library.assert.struct(CalSettings, {...
%                 'iBlts', 'BltsSrc', 'biasHighGain', ...
%                 'iCalibTimeL', 'iCalibTimeH', 'iLsf'}, {})   % Too slow?
            EJ_library.assert.sizes(zv_CALIBRATION_TABLE_INDEX, [1,2])
            assert(islogical(voltageNaN) && isscalar(voltageNaN))
            
            

            % Set cti1, cti2.
            if obj.use_CALIBRATION_TABLE_rcts
                cti1 = zv_CALIBRATION_TABLE_INDEX(1,1);
            else
                cti1 = 0;
            end
            % NOTE: NOT incrementing by one, since the variable's meaning can
            % vary between LFR, TDS-CWF, TDS-RSWF.
            cti2 = zv_CALIBRATION_TABLE_INDEX(1,2);

            
            
            if obj.allVoltageCalibDisabled || voltageNaN
                
                samplesCaAVolt = cell(size(samplesCaTm));
                
                for i = 1:numel(samplesCaTm)
                    if obj.allVoltageCalibDisabled
                        % CASE: Set voltages to TM values.
                        samplesCaAVolt{i} = double(samplesCaTm{i});
                    end
                    if voltageNaN
                        % CASE: Set voltages to NaN.
                        
                        % IMPLEMENTATION NOTE: Potentially overwrites TM value
                        % set in above "if" statement.
                        samplesCaAVolt{i} = nan(size(samplesCaTm{i}));
                    end
                end
                
            else

                if isLfr
                    %===========
                    % CASE: LFR
                    %===========
                    samplesCaAVolt = obj.calibrate_LFR_full(...
                        dtSec, samplesCaTm, CalSettings, cti1, cti2);
                else
                    %===========
                    % CASE: TDS
                    %===========
                    if isTdsCwf
                        % CASE: TDS CWF
                        samplesCaAVolt = obj.calibrate_TDS_CWF_full(...
                            dtSec, samplesCaTm, CalSettings, cti1, cti2);
                    else
                        % CASE: TDS RSWF
                        samplesCaAVolt = obj.calibrate_TDS_RSWF_full(...
                            dtSec, samplesCaTm, CalSettings, cti1, cti2);
                    end
                end
                
            end
        end



        % ARGUMENTS
        % =========
        % samplesTm    : 1D cell array of numeric 1D arrays.
        % samplesAVolt : 1D cell array of numeric 1D arrays.
        % CalSettings  : Struct that groups together arguments
        %   .iBlts     : Scalar integer. 1..5.
        %   .BltsSrc   : bicas.proc.L1L2.BLTS_src_dest describing where the
        %                signal comes from.
        %   ...
        %
        function samplesCaAVolt = calibrate_LFR_full(obj, ...
                dtSec, samplesCaTm, CalSettings, cti1, cti2)

            % ASSERTIONS
            EJ_library.assert.vector(samplesCaTm)
            assert(iscell(samplesCaTm))
            EJ_library.assert.vector(dtSec)
            assert(numel(samplesCaTm) == numel(dtSec))

            
            
            %=============================
            % Obtain all calibration data
            %=============================
            CalibData = obj.get_BIAS_LFR_calib_data(...
                CalSettings, cti1, cti2);

            %=======================================
            % CALIBRATE: LFR TM --> TM --> avolt
            %=======================================
            samplesCaAVolt = cell(size(samplesCaTm));
            for i = 1:numel(samplesCaTm)
                
                % ADD LSF OFFSET
                samplesTm = samplesCaTm{i}(:) + CalibData.lsfOffsetTm;
                
                % APPLY TRANSFER FUNCTION (BIAS + LFR)
                tempSamplesAVolt = bicas.tf.apply_TF_freq_modif(...
                    dtSec(i), ...
                    samplesTm, ...
                    CalibData.itfAvpt, ...
                    'detrendingDegreeOf',      CalibData.detrendingDegreeOf, ...
                    'retrendingEnabled',       CalibData.retrendingEnabled, ...
                    'tfHighFreqLimitFraction', CalibData.itfHighFreqLimitFraction);

                % ADD BIAS offset
                samplesCaAVolt{i} = tempSamplesAVolt + CalibData.BiasCalibData.offsetAVolt;
            end
        end
        


        % ARGUMENTS
        % =========
        % See calibrate_LFR_full.
        %
        function samplesCaAVolt = calibrate_TDS_CWF_full(obj, ...
                dtSec, samplesCaTm, CalSettings, cti1, cti2)

%             EJ_library.assert.struct(CalSettings, {...
%                 'iBlts', 'BltsSrc', 'biasHighGain', ...
%                 'iCalibTimeL', 'iCalibTimeH'}, {'iLsf'})   % Too slow?
            iBlts        = CalSettings.iBlts;
            BltsSrc      = CalSettings.BltsSrc;
            biasHighGain = CalSettings.biasHighGain;
            iCalibTimeL  = CalSettings.iCalibTimeL;
            iCalibTimeH  = CalSettings.iCalibTimeH;
            
            % ASSERTIONS
            EJ_library.assert.vector(dtSec)
            assert(iscell(samplesCaTm))
            assert(numel(samplesCaTm) == numel(dtSec))
            bicas.proc.L1L2.cal_utils.assert_iBlts(iBlts)
            assert(isa(BltsSrc, 'bicas.proc.L1L2.BLTS_src_dest'))
            assert(cti1 >= 0)
            
            if obj.use_CALIBRATION_TABLE_INDEX2
                % TODO? ASSERTION: cti2 = 0???
                error('BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
                    'TDS-CWF calibration never uses CALIBRATION_TABLE_INDEX2.')
            end
            
            % Initialize empty output variable.
            samplesCaAVolt = cell(size(samplesCaTm));

            if ismember(iBlts, [1,2,3])
                % CASE: BLTS 1-3 which TDS does support.
                
                %==============================
                % Obtain calibration constants
                %==============================
                % NOTE: Low/high gain is irrelevant for TDS. Argument value
                % arbitrary.
                BiasCalibData = obj.get_BIAS_calib_data(...
                    BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
                
                if obj.lfrTdsTfDisabled
                    tdsFactorIvpt = 1;
                else
                    RctList       = obj.RctDataMap('TDS-CWF');
                    tdsFactorIvpt = RctList{cti1+1}.factorsIvpt(iBlts);
                end
                
                for i = 1:numel(samplesCaTm)
                    
                    %===============================================
                    % CALIBRATE: TDS TM --> TDS/BIAS interface volt
                    %===============================================
                    % MULTIPLICATION
                    tempSamplesIVolt = tdsFactorIvpt * samplesCaTm{i};
                   
                    %=====================================================
                    % CALIBRATE: TDS/BIAS interface volt --> antenna volt
                    %=====================================================
                    % APPLY TRANSFER FUNCTION (for BIAS, but not for TDS-CWF)
                    tempSamplesAVolt = bicas.tf.apply_TF_freq_modif(...
                        dtSec(i), ...
                        tempSamplesIVolt(:), ...
                        BiasCalibData.itfAvpiv, ...
                        'detrendingDegreeOf',      obj.dcDetrendingDegreeOf, ...
                        'retrendingEnabled',       obj.dcRetrendingEnabled, ...
                        'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction);
                    
                    % ADD BIAS OFFSET
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
        function samplesCaAVolt = calibrate_TDS_RSWF_full(obj, ...
                dtSec, samplesCaTm, CalSettings, cti1, cti2)
            
%             EJ_library.assert.struct(CalSettings, {...
%                 'iBlts', 'BltsSrc', 'biasHighGain', ...
%                 'iCalibTimeL', 'iCalibTimeH'}, {'iLsf'})   % Too slow?
            iBlts        = CalSettings.iBlts;
            BltsSrc      = CalSettings.BltsSrc;
            biasHighGain = CalSettings.biasHighGain;
            iCalibTimeL  = CalSettings.iCalibTimeL;
            iCalibTimeH  = CalSettings.iCalibTimeH;
            
            % ASSERTIONS
            EJ_library.assert.vector(dtSec)
            assert(iscell(samplesCaTm))
            assert(numel(samplesCaTm) == numel(dtSec))
            bicas.proc.L1L2.cal_utils.assert_iBlts(iBlts)
            assert(isa(BltsSrc, 'bicas.proc.L1L2.BLTS_src_dest'))
            assert(cti1 >= 0)
            
            if obj.use_CALIBRATION_TABLE_INDEX2
                % TODO? ASSERTION: cti2 = 0???
                error(...
                    'BICAS:Assertion:IllegalCodeConfiguration:OperationNotImplemented', ...
                    'TDS-RSWF calibration never uses CALIBRATION_TABLE_INDEX2.')
            end
            
            %==============================
            % Obtain calibration constants
            %==============================
            % NOTE: Low/high gain is irrelevant for TDS. Argument value
            % arbitrary.
            BiasCalibData = obj.get_BIAS_calib_data(...
                BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            
            % Initialize empty output variable.
            samplesCaAVolt = cell(size(samplesCaTm));
            if ismember(iBlts, [1,2,3])
                
                %======================================
                % Create combined ITF for TDS and BIAS
                %======================================
                if obj.lfrTdsTfDisabled
                    tdsItfIvpt = @(omegaRps) (ones(omegaRps));
                else
                    RctList = obj.RctDataMap('TDS-RSWF');
                    tdsItfIvpt = RctList{cti1+1}.itfModifIvptCa{iBlts};
                end

                itfAvpt = @(omegaRps) (...
                    tdsItfIvpt(omegaRps) ...
                    .* ...
                    BiasCalibData.itfAvpiv(omegaRps));
                
                %====================================
                % CALIBRATE: TDS TM --> antenna volt
                %====================================
                % APPLY TRANSFER FUNCTION (BIAS + TDS-RSWF)
                for i = 1:numel(samplesCaTm)
                    tempSamplesAVolt = bicas.tf.apply_TF_freq_modif(...
                        dtSec(i), ...
                        samplesCaTm{i}(:), ...
                        itfAvpt, ...
                        'detrendingDegreeOf',      obj.dcDetrendingDegreeOf, ...
                        'retrendingEnabled',       obj.dcRetrendingEnabled, ...
                        'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction);
                    
                    % ADD BIAS OFFSET
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



        function iCalib = get_BIAS_calibration_time_L(obj, Epoch)
            iCalib = bicas.proc.L1L2.cal_utils.get_calibration_time(...
                Epoch, obj.RctDataMap('BIAS').epochL);
        end



        function iCalib = get_BIAS_calibration_time_H(obj, Epoch)
            iCalib = bicas.proc.L1L2.cal_utils.get_calibration_time(...
                Epoch, obj.RctDataMap('BIAS').epochH);
        end



        % Return subset of already loaded BIAS calibration data, for specified
        % settings.
        %
        % NOTE: May return calibration values corresponding to scalar
        % calibration, depending on SETTINGS:
        %
        %
        % ARGUMENTS
        % =========
        % biasHighGain
        %       NUMERIC value: 0=Off, 1=ON, or NaN=Value not known.
        %       IMPLEMENTATION NOTE: Needs value to represent that biasHighGain
        %       is unknown. Sometimes, if biasHighGain is unknown, then it is
        %       useful to process as usual since some of the data can still be
        %       derived/calibrated, so that the caller does not need to handle
        %       the special case.
        %
        function BiasCalibData = get_BIAS_calib_data(obj, ...
                BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH)
            
            % PROPOSAL: Log warning message when simultaneously biasHighGain=NaN
            % and the value is needed.
            
            % ASSERTION
            assert(isa(BltsSrc, 'bicas.proc.L1L2.BLTS_src_dest'))
            assert(isscalar(biasHighGain) && isnumeric(biasHighGain))
            assert(isscalar(iCalibTimeL))
            assert(isscalar(iCalibTimeH))
            
            BiasRct = obj.RctDataMap('BIAS');

            %###################################################################
            % kIvpav = Multiplication factor "k" that represents/replaces the
            % (forward) transfer function.
            switch(BltsSrc.category)
                case 'DC single'
                    
                    % NOTE: List of ITFs for different times.
                    biasItfAvpiv = BiasRct.ItfSet.dcSingleAvpiv{iCalibTimeL};
                    kFtfIvpav    = obj.BiasScalarGain.alphaIvpav;
                    offsetAVolt  = BiasRct.dcSingleOffsetsAVolt(...
                        iCalibTimeH, BltsSrc.antennas);
                    
                case 'DC diff'
                    
                    biasItfAvpiv = BiasRct.ItfSet.dcDiffAvpiv{iCalibTimeL};
                    kFtfIvpav    = obj.BiasScalarGain.betaIvpav;
                    if     isequal(BltsSrc.antennas(:)', [1,2]);   offsetAVolt = BiasRct.DcDiffOffsets.E12AVolt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [1,3]);   offsetAVolt = BiasRct.DcDiffOffsets.E13AVolt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [2,3]);   offsetAVolt = BiasRct.DcDiffOffsets.E23AVolt(iCalibTimeH);
                    else
                        error('BICAS:Assertion:IllegalArgument', ...
                            'Illegal BltsSrc.');
                    end
                    
                case 'AC diff'
                    
                    if     biasHighGain == 0
                        biasItfAvpiv = BiasRct.ItfSet.acLowGainAvpiv{iCalibTimeL};
                        kFtfIvpav    = obj.BiasScalarGain.gammaIvpav.lowGain;
                        offsetAVolt  = 0;
                    elseif biasHighGain == 1
                        biasItfAvpiv = BiasRct.ItfSet.acHighGainAvpiv{iCalibTimeL};
                        kFtfIvpav    = obj.BiasScalarGain.gammaIvpav.highGain;
                        offsetAVolt  = 0;
                    elseif isnan(biasHighGain)
                        % CASE: GAIN unknown when it is NEEDED for calibration.
                        biasItfAvpiv = bicas.proc.L1L2.cal.NAN_TF;
                        kFtfIvpav    = NaN;
                        offsetAVolt  = NaN;
                    else
                        error('BICAS:Assertion:IllegalArgument', ...
                            'Illegal argument biasHighGain=%g.', biasHighGain)
                    end

                otherwise
                    error('BICAS:Assertion:IllegalArgument', ...
                        ['Illegal argument BltsSrc.category=%s.', ...
                        ' Can not obtain calibration data for this type of signal.'], ...
                        BltsSrc.category)
            end
            
            if obj.biasOffsetsDisabled && ~isnan(offsetAVolt)
                % NOTE: Overwrites "offsetAVolt".
                offsetAVolt = 0;
            end
            if obj.useBiasTfScalar
                % NOTE: Overwrites "biasItfAvpiv".
                biasItfAvpiv = @(omegaRps) (ones(size(omegaRps)) / kFtfIvpav);
            end
            %###################################################################
            
            
            
            BiasCalibData.itfAvpiv    = biasItfAvpiv;
            BiasCalibData.offsetAVolt = offsetAVolt;

        end
        
        
        
        % Obtain LFR ITF, but handle the case that should never happen for
        % actual non-NaN data (LSF F3 + BLTS 4 or 5) and return a TF that only
        % returns NaN instead. BICAS may still iterate over that combination
        % though when calibrating.
        % 
        function lfrItfIvpt = get_LFR_ITF(obj, cti1, iBlts, iLsf)
            % ASSERTIONS
            assert(cti1 >= 0)
            bicas.proc.L1L2.cal_utils.assert_iBlts(iBlts)
            bicas.proc.L1L2.cal_utils.assert_iLsf(iLsf)
            
            if (iLsf == 4) && ismember(iBlts, [4,5])
                % CASE: F3 and BLTS={4,5}
                
                % NOTE: There is no tabulated LFR TF and no such combination
                % signal route, so the TF can not be returned even in principle.
                lfrItfIvpt = bicas.proc.L1L2.cal.NAN_TF;
            else
                RctDataList = obj.RctDataMap('LFR');
                
                % ASSERTION
                % IMPLEMENTATION NOTE: Anonymous function below will fail at a
                % later stage if these assertions are false. Checking for these
                % criteria here makes it easier to understand these particular
                % types of error.
                assert(numel(RctDataList) >= (cti1+1), ...
                    'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
                    ['LFR RctDataList is too small for argument cti1=%g.', ...
                    ' This could indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
                    ' value is larger than glob. attr. CALIBRATION TABLE allows.'], cti1)
                assert(~isempty(RctDataList{cti1+1}), ...
                    'BICAS:IllegalArgument:DatasetFormat:Assertion', ...
                    ['LFR RctDataList contains no RCT data corresponding', ...
                    ' to argument cti1=%g. This may indicate that', ...
                    ' a zVar CALIBRATION_TABLE_INDEX(:,1) value is wrong or', ...
                    ' that BICAS did not try to load the corresponding RCT', ...
                    ' in glob. attr. CALIBRATION_TABLE.'], cti1)
                
                lfrItfIvpt = RctDataList{cti1+1}.ItfModifIvptCaCa{iLsf}{iBlts};
            end
        end
        
        
        
        % Return calibration data for LFR+BIAS calibration.
        %
        % RATIONALE
        % =========
        % Method exists to
        % (1) simplify & clarify calibrate_LFR_full(),
        % (2) be useful for non-BICAS code to inspect the calibration data used
        %     for a particular calibration case, in particular the combined
        %     (LFR+BIAS) transfer functions.
        %     NOTE: This is also the reason why this method is public.
        %
        % IMPLEMENTATION NOTE: Return one struct instead of multiple return
        % values to make sure that the caller does not confuse the return values
        % with each other.
        function [CalData] = get_BIAS_LFR_calib_data(obj, CalSettings, cti1, cti2)
            
            % ASSERTIONS
%             EJ_library.assert.struct(CalSettings, {...
%                 'iBlts', 'BltsSrc', 'biasHighGain', ...
%                 'iCalibTimeL', 'iCalibTimeH', 'iLsf'}, {})   % Too slow?
            iBlts        = CalSettings.iBlts;
            BltsSrc      = CalSettings.BltsSrc;
            biasHighGain = CalSettings.biasHighGain;
            iCalibTimeL  = CalSettings.iCalibTimeL;
            iCalibTimeH  = CalSettings.iCalibTimeH;
            iLsf         = CalSettings.iLsf;
            
            % ASSERTIONS
            bicas.proc.L1L2.cal_utils.assert_iBlts(iBlts)
            assert(isa(BltsSrc, 'bicas.proc.L1L2.BLTS_src_dest'))
            bicas.proc.L1L2.cal_utils.assert_iLsf(iLsf)
            assert(isscalar(cti1))
            assert(cti1 >= 0, 'Illegal cti1=%g', cti1)
            % No assertion on cti2 unless used (determined later).



            %============================================
            % Only place to potentially make use of cti2
            %============================================
            if obj.use_CALIBRATION_TABLE_INDEX2
                % ASSERTIONS
                assert(isscalar(cti2), ...
                    'BICAS:IllegalArgument:Assertion', ...
                    'Argument cti2 is not scalar.')
                assert(cti2 >= 0, ...
                    'BICAS:IllegalArgument:Assertion', ...
                    ['Illegal argument cti2=%g', ...
                    ' (=zVar CALIBRATION_TABLE_INDEX(iRecord, 2))'], ...
                    cti2)
                assert(iLsf == cti2+1, ...
                    'BICAS:IllegalArgument:Assertion', ...
                    'cti2+1=%i != iLsf=%i (before overwriting iLsf)', ...
                    cti2+1, iLsf)
                
                % NOTE: Only place cti2 is used.
                iLsf = cti2 + 1;
            end
            
            

            CalData = struct();
            
            % (Inofficial) Experimental LFR calibration offsets.
            CalData.lsfOffsetTm = obj.lfrLsfOffsetsTm(CalSettings.iLsf);

            %====================================================
            % Obtain settings for bicas.tf.apply_TF_freq_modif()
            %====================================================
            if CalSettings.BltsSrc.is_AC()
                % IMPLEMENTATION NOTE: DC is (optionally) detrended via
                % bicas.tf.apply_TF_freq_modif() in the sense of a linear fit
                % being removed, TF applied, and then added back. That same
                % algorithm, or at least adding back the fit, is by its nature
                % inappropriate for non-lowpass filters, i.e. for AC. (The fit
                % can not be scaled with the 0 Hz signal amplitude)
                CalData.detrendingDegreeOf = obj.acDetrendingDegreeOf;
                CalData.retrendingEnabled  = false;   % NOTE: HARDCODED SETTING.
            else
                CalData.detrendingDegreeOf = obj.dcDetrendingDegreeOf;
                CalData.retrendingEnabled  = obj.dcRetrendingEnabled;
            end
            CalData.itfHighFreqLimitFraction = obj.itfHighFreqLimitFraction;
            
            %==============================
            % Obtain BIAS calibration data
            %==============================
            CalData.BiasCalibData = obj.get_BIAS_calib_data(...
                BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            
            %========================================
            % Obtain (official) LFR calibration data
            %========================================            
            if obj.lfrTdsTfDisabled
                CalData.lfrItfIvpt = @(omegaRps) (ones(size(omegaRps)));
            else
                CalData.lfrItfIvpt = obj.get_LFR_ITF(cti1, iBlts, iLsf);
            end
            
            %======================================
            % Create combined ITF for LFR and BIAS
            %======================================
            CalData.itfAvpt = bicas.proc.L1L2.cal_utils.create_LFR_BIAS_ITF(...
                CalData.lfrItfIvpt, ...
                CalData.BiasCalibData.itfAvpiv, ...
                BltsSrc.is_AC(), ...
                obj.itfAcConstGainLowFreqRps);
        end



    end    % methods(Access=private)



    %###########################################################################

    
    
    methods(Static, Access=public)



%         function tfZ = parasitic_capacitance_TF(tfOmega)
%             % Calculate Z(omega) values for TF representing parasitic
%             % capacitances (based on analytic function).
% 
%             % Function name? "Input capacitance"?
%             % Not read R & C from constants here? Submit as arguments?
%             capacitanceFarad =
%             impedanceOhm     =
%             
%             % Correct for a TF?!
%             % 2020-11-11, EJ: Not same as in paper note.
%             tfZ = 1 / (1 + 1j*tfOmega*capacitanceFarad*impedanceOhm);
%             
%             error('BICAS:OperationNotImplemented', 'Function not implemented Yet.')
%         end
        
        
        
        % Convert "set current" to TC/TM units.
        %
        function biasCurrentTm = calibrate_current_sampere_to_TM(currentSAmpere)
            
            % ASSERTION: Input values are within range.
            % NOTE: max(...) ignores NaN, unless that is the only value, which
            % then becomes the max value.
            [maxAbsSAmpere, iMax] = max(abs(currentSAmpere(:)));
            if ~(isnan(maxAbsSAmpere) || (maxAbsSAmpere <= EJ_library.so.constants.MAX_ABS_SAMPERE))
                
                error('BICAS:Assertion:IllegalArgument', ...
                    ['Argument currentSAmpere (unit: set current/ampere)', ...
                    ' contains illegally large value(s).', ...
                    ' Largest found value is %g.'], ...
                    currentSAmpere(iMax))
            end
            
            biasCurrentTm = currentSAmpere * EJ_library.so.constants.TM_PER_SAMPERE;
        end
        
        
        
    end    % methods(Static, Access=public)
    
    

end
