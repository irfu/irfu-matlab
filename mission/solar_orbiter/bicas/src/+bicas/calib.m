%
% Class for
% (1) library/utility functions that calibrate data.
% (2) helper functions to find and load calibration data from file, as needed by
%     BICAS.
% An instance of this class may or may not contain calibration data for all
% types of data/RCTs depending on how it was initialized.
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
% - is extra careful with identifiers with units
% - is extra careful with well defined terms/shorternings/naming conventions
% since
% (1) undiscovered calibration bugs could be considered extra bad
% (2) it is expected to be hard to detect certain bugs,
% (3) it could be hard to use automatic testing here,
% (4) to detect changing RCT formats, in particular in RCTS from non-BIAS teams.
% --
% All calibration functions of measured data are assumed to accept data from all
% BLTS (1-5), i.e. including TDS, in order to reduce the number assumptions that
% the calling code needs to make.
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
%          NOTE: When used as an array index, 1=F0, ..., 4=F3.
% TF     = Transfer function (Z=Z(omega))
% FTF    = Forward Transfer Function = TF that describes physical input-to-output (not the reverse)
% ITF    = Inverse Transfer Function = TF that describes physical output-to-input (not the reverse)
% CTI    = CALIBRATION_TABLE_INDEX (zVar)
% CTI1   = First  value in record of zVar CALIBRATION_TABLE_INDEX.
% CTI2   = Second value in record of zVar CALIBRATION_TABLE_INDEX.
% RCTS   = RCT CALIBRATION_TABLE (glob.attr)+CALIBRATION_TABLE_INDEX (zVar).
%          S = plural
%
%
% UNITS / TYPES OF QUANTITIES
% ---------------------------
% TM       = Telemetry units (in LFR/TDS ADC), or telecommand (TC) units. Using
%            this instead of the term "count".
% IV=ivolt = Interface Volt = Calibrated volt at the interface between BIAS and
%            LFR/TDS.
% AV=avolt = Antenna   Volt = Calibrated volt at the antennas, i.e. the final
%            calibrated (measured) value, including for reconstructed signals
%            (e.g. diffs calculated from singles). May also refer to offsets and
%            values without offsets.
% aampere  = Antenna ampere = Calibrated ampere at the antenna.
% sampere  = Set current ampere. Exactly proportional to bias current in TM.
% TPIV     = TM/interface volt
% IVPT     = Interface volt/TM
% AAPT     = Antenna ampere/TM
% AVPIV    = Antenna volt/interface volt
% IVPAV    = Interface volt/antenna volt
% Deg      = Degrees (angle). 1 revolution=360 degrees=2*pi radians.
% RPS      = Radians/second
% Sec      = Seconds
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
% CALIBRATION_TABLE       : Global attribute
%   """"Filename of the calibration table(s).""""
%   """"There must as many as entries than the number of calibration table files
%   associated to the L1R file.""""
% CALIBRATION_TABLE_INDEX : zVariable
%   """"Index of the calibration table(s) value required to generate L2 data
%   files.""""
%   """"Each CDF record must contain 2 elements: the first element must gives
%   the index of the associated CALIBRATION_TABLE entry (i.e., 0 for the first
%   entry, 1 for the second, etc.). The second element must refer to the index
%   of the value to be used inside the calibration table file.""""
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
% PROPOSAL: Store both FTFs and ITFs, despite that FTFs are not used for calibration directly.
%   NOTE: BIAS & LFR RCTs contain FTFs, TDS RCT contains ITFs.
%   NOTE: Has to keep track of FTF/ITF before modifications (extrapolation to 0 Hz, Z=0 for high freq.).
%   PRO: Useful for debugging. Can easily inspect & plot FTFs.
%
% PROPOSAL: Assertion function for CalSettings.
%   TODO-NI: Same struct, with same fields in all cases?
%   NOTE: Function does not know which fields are actually used.



    properties(Access=private, Constant)
        
        % containers.Map: RCT Type ID --> Info about RCT type.
        % Its keys defines the set of RCT Type ID strings.
        RCT_TYPES_MAP       = bicas.calib.init_RCT_Types_Map();
        
        % Local constant for convenience.
        NAN_TF              = @(omegaRps) (omegaRps * NaN);
        
        READING_RCT_PATH_LL = 'info';    % LL = Log Level
        RCT_DATA_LL         = 'debug';
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
        BiasScalarGain             % BIAS scalar (simplified) calibration, not in the RCTs. For debugging/testing purposes.
        HkBiasCurrent
        lfrLsfOffsetsTm    = [];   % EXPERIMENTAL. NOTE: Technically, the name contains a tautology (LFR+LSF).
        
        
        
        %==================================================
        % Settings for what kind of calibration to perform
        %==================================================
        
        % Corresponds to SETTINGS key-value.
        enableDetrending
        itfHighFreqLimitFraction;
                
        % Whether to select non-BIAS RCT using global attribute CALIBRATION_TABLE (and
        % CALIBRATION_TABLE_INDEX(iRecord,1)).
        use_CALIBRATION_TABLE_rcts        
        % Whether to use CALIBRATION_TABLE_INDEX(iRecord,2) for calibration.
        use_CALIBRATION_TABLE_INDEX2

        % What type of calibration to use.
        allVoltageCalibDisabled    % Use TM values (not set to NaN).
        useBiasTfScalar
        biasOffsetsDisabled
        lfrTdsTfDisabled
        
    end
    
    
    
    %###################################################################################################################

    methods(Access=public)



        % Constructor.
        %
        %
        % IMPLEMENTATION NOTE
        % ===================
        % The class (instance methods, including constructor) deliberately does
        % not itself determine which non-BIAS RCTs to read. This is intended so
        % that the user can select non-BIAS RCTs by initializing a
        % containers.Map using either
        % (1) static helper method "find_read_non_BIAS_RCTs_by_regexp", 
        % (2) static helper method "find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE",
        %     or
        % (3) manually (for debugging/analysis/testing).
        % --
        % The class (instance methods, including constructor) deliberately does
        % not READ the non-BIAS RCTs. This is useful since
        % ** it makes it possible to inspect & modify the RCT content before
        %    submitting it to bicas.calib
        % ** it completely separates the RCT-reading from the class
        %    (modularization)
        % ** it simplifies the constructor somewhat since it does not need to
        %    translate paths into RCT data for non-empty cells.
        %
        function obj = calib(NonBiasRctDataMap, rctDir, use_CALIBRATION_TABLE_rcts, use_CALIBRATION_TABLE_INDEX2, SETTINGS, L)
            % PROPOSAL: Abolish second init methods:
            %       read_non_BIAS_RCTs_by_regexp
            %       read_non_BIAS_RCT_by_CALIBRATION_TABLE
            %   PROPOSAL: Initialize completely using constructor.
            %       NOTE: Can not initialize in execute_sw_mode (since does not know LFR/TDS since does not know s/w mode).
            %       NOTE: Will have to submit rctDir to prodFunc or to swmode_defs (analogous to SETTINGS, L).
            %       PRO: Can abolish obj.hasLoadedNonBiasData .
            %   PROPOSAL: Initialize using only ONE second init method.
            %   PROPOSAL: (1) Two static methods (two alternate ways) to initialize 1D cell array of non-BIAS RCTs to
            %             load, corresponding to those RCTs in ga. CALIBRATION_TABLE that should be loaded, leaving
            %             ~empty cells for those that should not (and can not) be loaded.
            %             (2) Cell array is passed to calib init method (primary or secondary), combined with flag for
            %             whether calib should read zv CALIBRATION_TABLE_INDEX.
            %       PROBLEM: zv BW only applies to LFR, not TDS.
            %       CON: Requires knowledge of zv BW and CALIBRATION_TABLE_INDEX which is then used in both static calib
            %            methods (above) and the calibration code.
            %           CON: Relatively harmless since reading too few or too many RCTs will lead to error.
            %   PROPOSAL: Dynamic RCT loading. Load RCT first time it is requested.
            %       NOTE: Must still be aware of how to load: using CALIBRATION_TABLE_INDEX+CALIBRATION_TABLE, or using
            %             pattern matching/regexp..
            %       NOTE: Must connect calib fields with rctTypeId and which RCT to load.
            %           PROPOSAL: containers.Map: rctTypeId --> (data, bicas.RCT reading method, SETTINGS key for filename regexp.).
            %       PRO: Knowledge of BW and CALIBRATION_TABLE_INDEX+CALIBRATION_TABLE in one place.
            %       CON: Code will not immediately know whether all non-BIAS RCTs can be loaded. May crash after a delay.

            % ASSERTIONS
            assert(isscalar(use_CALIBRATION_TABLE_INDEX2))
            
            filenameRegexp = SETTINGS.get_fv(bicas.calib.RCT_TYPES_MAP('BIAS').filenameRegexpSettingKey);
            filePath       = bicas.RCT.find_RCT_regexp(rctDir, filenameRegexp, L);
            % IMPLEMENTATION NOTE: It has been observed that value sometimes survives from previous runs, despite being
            % an instance variable. Unknown why. Therefore explicitly overwrites it.
            obj.RctDataMap = containers.Map();
            obj.RctDataMap('BIAS')                 = bicas.calib.read_log_modify_RCT('BIAS', filePath, L);
            
            EJ_library.assert.subset(NonBiasRctDataMap.keys, bicas.calib.RCT_TYPES_MAP.keys)
            for rctTypeId = NonBiasRctDataMap.keys
                obj.RctDataMap(rctTypeId{1}) = NonBiasRctDataMap(rctTypeId{1});
            end
            
            
            obj.BiasScalarGain.alphaIvpav          = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.ALPHA_IVPAV');
            obj.BiasScalarGain.betaIvpav           = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.BETA_IVPAV');
            obj.BiasScalarGain.gammaIvpav.highGain = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.HIGH_GAIN');
            obj.BiasScalarGain.gammaIvpav.lowGain  = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.LOW_GAIN');
            
            obj.HkBiasCurrent.offsetTm             = SETTINGS.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.OFFSET_TM');
            obj.HkBiasCurrent.gainAapt             = SETTINGS.get_fv('PROCESSING.CALIBRATION.CURRENT.HK.GAIN_AAPT');
            
            obj.lfrLsfOffsetsTm                    = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR.LSF_OFFSETS_TM');

            
            
            obj.enableDetrending                   = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF_DETRENDING_ENABLED');
            obj.itfHighFreqLimitFraction           = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF_HIGH_FREQ_LIMIT_FRACTION');
            
            obj.allVoltageCalibDisabled            = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.DISABLE');
            obj.biasOffsetsDisabled                = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.DISABLE_OFFSETS');
            obj.lfrTdsTfDisabled                   = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.LFR_TDS.TF_DISABLED');
            
            % Set obj.useBiasTfScalar.
            settingBiasTf                          = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF');
            switch(settingBiasTf)
                case 'FULL'
                    obj.useBiasTfScalar = 0;
                case 'SCALAR'
                    obj.useBiasTfScalar = 1;
                otherwise
                    error(...
                        'BICAS:calib:Assertion:ConfigurationBug', ...
                        'Illegal value for setting PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF="%s".', ...
                        settingBiasTf)
            end
            
            obj.use_CALIBRATION_TABLE_rcts   = use_CALIBRATION_TABLE_rcts;
            obj.use_CALIBRATION_TABLE_INDEX2 = use_CALIBRATION_TABLE_INDEX2;
        end



        % Convert/calibrate TC bias current: TM units --> physical units.
        %
        % NOTE: This is the normal way of obtaining bias current in physical units (as opposed to HK bias current).
        %
        function biasCurrentAAmpere = calibrate_current_TM_to_aampere(obj, biasCurrentTm, iAntenna, iCalibTimeL)
            
            %==============================
            % Obtain calibration constants
            %==============================
            offsetAAmpere = obj.RctDataMap('BIAS').Current.offsetsAAmpere(iCalibTimeL, iAntenna);
            gainAapt      = obj.RctDataMap('BIAS').Current.gainsAapt(     iCalibTimeL, iAntenna);

            % CALIBRATE
            biasCurrentAAmpere = offsetAAmpere + gainAapt .* double(biasCurrentTm);    % LINEAR FUNCTION
        end



        % Convert/calibrate diagnostic HK TM bias current values to physical units.
        % Refers to BIAS HK zVars HK_BIA_BIAS1/2/3.
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
        %
        % ARGUMENTS
        % =========
        % voltageNaN : Scalar logical. Whether to set voltages to NaN and not
        %              executing any calibration.
        %              RATIONALE: This option is useful to
        %              (1) potentially speed up BICAS when it is known that
        %              data will be overwritten with fill values later.
        %              (2) avoid executing calibration algorithms when it is
        %              known that there is no calibration configuration anyway
        %              Ex: LFR zVar BW=0 => CALIBRATION_TABLE_INDEX(1,:) illegal.
        %              ==> Can not calibrate.
        %              Note: This means that this function technically accepts
        %              an illegal calibration configuration when argument is set
        %              to true.
        %
        function samplesCaAVolt = calibrate_voltage_all(obj, ...
                dtSec, samplesCaTm, isLfr, isTdsCwf, CalSettings, ...
                zv_CALIBRATION_TABLE_INDEX, voltageNaN)
            
            % ASSERTIONS
            %EJ_library.assert.struct(CalSettings, {'iBlts', 'BltsSrc', 'biasHighGain', 'iCalibTimeL', 'iCalibTimeH', 'iLsf'}, {})   % Too slow?
            EJ_library.assert.sizes(zv_CALIBRATION_TABLE_INDEX, [1,2])
            assert(islogical(voltageNaN) && isscalar(voltageNaN))
            
            

            % Set cti1, cti2.
            if obj.use_CALIBRATION_TABLE_rcts
                cti1 = zv_CALIBRATION_TABLE_INDEX(1,1);
            else
                cti1 = 0;
            end
            % NOTE: NOT incrementing by one, since the variable's meaning can vary between LFR, TDS-CWF, TDS-RSWF.
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
                        % (Potentially overwrite TM value.)
                        samplesCaAVolt{i} = nan(size(samplesCaTm{i}));
                    end
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
        %   .BltsSrc   : bicas.BLTS_src_dest describing where the signal comes
        %                from.
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
            bicas.calib.assert_iBlts(iBlts)
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            bicas.calib.assert_iLsf(iLsf)
            assert(isscalar(cti1))
            assert(cti1 >= 0, 'Illegal cti1=%g', cti1)
            % No assertion on cti2 unless used (determined later).



            %============================================
            % Only place to potentially make use of cti2
            %============================================
            if obj.use_CALIBRATION_TABLE_INDEX2
                % ASSERTIONS
                assert(isscalar(cti2), 'BICAS:calib:IllegalArgument:Assertion', 'Argument cti2 is not scalar.')
                assert(cti2 >= 0,      'BICAS:calib:IllegalArgument:Assertion', 'Illegal argument cti2=%g (=zVar CALIBRATION_TABLE_INDEX(iRecord, 2))', cti2)
                assert(iLsf == cti2+1, 'BICAS:calib:IllegalArgument:Assertion', 'cti2+1=%i != iLsf=%i (before overwriting iLsf)', cti2+1, iLsf)
                
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
            % CALIBRATE: LFR TM --> TM --> avolt
            %=======================================
            samplesCaAVolt = cell(size(samplesCaTm));
            lsfOffsetTm = obj.lfrLsfOffsetsTm(iLsf);
            for i = 1:numel(samplesCaTm)
                
                % ADD LSF OFFSET
                samplesTm = samplesCaTm{i}(:) + lsfOffsetTm;
                
                % APPLY TRANSFER FUNCTION
                tempSamplesAVolt = bicas.utils.apply_TF_freq(...
                    dtSec(i), samplesTm, itfIvpt, ...
                    'enableDetrending',        obj.enableDetrending, ...
                    'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction);

                samplesCaAVolt{i} = tempSamplesAVolt + BiasCalibData.offsetAVolt;
            end
        end



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
            bicas.calib.assert_iBlts(iBlts)
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
                        RctList       = obj.RctDataMap('TDS-CWF');
                        tdsFactorIvpt = RctList{cti1+1}.factorsIvpt(iBlts);
                    end
                    tempSamplesIVolt = tdsFactorIvpt * samplesCaTm{i};    % MULTIPLICATION
                   
                    %=====================================================
                    % CALIBRATE: TDS/BIAS interface volt --> antenna volt
                    %=====================================================
                    % APPLY TRANSFER FUNCTION
                    tempSamplesAVolt = bicas.utils.apply_TF_freq(...
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
            bicas.calib.assert_iBlts(iBlts)
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
                    RctList = obj.RctDataMap('TDS-RSWF');
                    tdsItfIvpt = @(omegaRps) (bicas.calib.eval_tabulated_ITF(...
                        RctList{cti1+1}.ItfIvptList{iBlts}, omegaRps));
                end

                itf = @(omegaRps) (...
                    tdsItfIvpt(omegaRps) ...
                    .* ...
                    BiasCalibData.itfAvpiv(omegaRps));
                
                %====================================
                % CALIBRATE: TDS TM --> antenna volt
                %====================================
                for i = 1:numel(samplesCaTm)
                    tempSamplesAVolt = bicas.utils.apply_TF_freq(...
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
            iCalib = bicas.calib.get_calibration_time(Epoch, obj.RctDataMap('BIAS').epochL);
        end



        function iCalib = get_calibration_time_H(obj, Epoch)
            iCalib = bicas.calib.get_calibration_time(Epoch, obj.RctDataMap('BIAS').epochH);
        end



    end    % methods(Access=public)

    
    
    %###################################################################################################################



    methods(Access=private)



        % Return subset of already loaded BIAS calibration data, for specified
        % settings.
        %
        % NOTE: May return calibration values corresponding to scalar
        % calibration, depending on SETTINGS:
        %
        %
        % ARGUMENTS
        % =========
        % biasHighGain : NUMERIC value: 0=Off, 1=ON, or NaN=Value not known.
        %                IMPLEMENTATION NOTE: Needs value to represent that
        %                biasHighGain is unknown. Sometimes, if biasHighGain is
        %                unknown, then it is useful to process as usual since
        %                some of the data can still be derived/calibrated, so
        %                that the caller does not need to handle the special
        %                case.
        %
        function BiasCalibData = get_BIAS_calib_data(obj, BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH)
            
            % ASSERTION
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            assert(isscalar(biasHighGain) && isnumeric(biasHighGain))
            assert(isscalar(iCalibTimeL))
            assert(isscalar(iCalibTimeH))
            
            BiasRct = obj.RctDataMap('BIAS');

            %###################################################################
            % kIvpav = Multiplication factor "k" that represents/replaces the
            % (forward) transfer function.
            switch(BltsSrc.category)
                case 'DC single'
                    
                    BiasItfAvpiv = TF_list_2_func(BiasRct.ItfSet.DcSingleAvpiv);    % NOTE: List of ITFs for different times.
                    kFtfIvpav    = obj.BiasScalarGain.alphaIvpav;
                    offsetAVolt  = BiasRct.dcSingleOffsetsAVolt(iCalibTimeH, BltsSrc.antennas);
                    
                case 'DC diff'
                    
                    BiasItfAvpiv = TF_list_2_func(BiasRct.ItfSet.DcDiffAvpiv);
                    kFtfIvpav    = obj.BiasScalarGain.betaIvpav;
                    if     isequal(BltsSrc.antennas(:)', [1,2]);   offsetAVolt = BiasRct.DcDiffOffsets.E12AVolt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [1,3]);   offsetAVolt = BiasRct.DcDiffOffsets.E13AVolt(iCalibTimeH);
                    elseif isequal(BltsSrc.antennas(:)', [2,3]);   offsetAVolt = BiasRct.DcDiffOffsets.E23AVolt(iCalibTimeH);
                    else
                        error('BICAS:calib:Assertion:IllegalArgument', ...
                            'Illegal BltsSrc.');
                    end
                    
                case 'AC diff'
                    
                    if     biasHighGain == 0
                        BiasItfAvpiv = TF_list_2_func(BiasRct.ItfSet.AcLowGainAvpiv);
                        kFtfIvpav    = obj.BiasScalarGain.gammaIvpav.lowGain;
                        offsetAVolt  = 0;
                    elseif biasHighGain == 1
                        BiasItfAvpiv = TF_list_2_func(BiasRct.ItfSet.AcHighGainAvpiv);
                        kFtfIvpav    = obj.BiasScalarGain.gammaIvpav.highGain;
                        offsetAVolt  = 0;
                    elseif isnan(biasHighGain)
                        BiasItfAvpiv = bicas.calib.NAN_TF;
                        kFtfIvpav    = NaN;
                        offsetAVolt  = NaN;
                    else
                        error('BICAS:calib:Assertion:IllegalArgument', ...
                            'Illegal argument biasHighGain=%g.', biasHighGain)
                    end

                otherwise
                    error('BICAS:calib:Assertion:IllegalArgument', ...
                        ['Illegal argument BltsSrc.category=%s.', ...
                        ' Can not obtain calibration data for this type of signal.'], ...
                        BltsSrc.category)
            end
            
            if obj.biasOffsetsDisabled && ~isnan(offsetAVolt)
                % NOTE: Overwrites "offsetAVolt".
                offsetAVolt = 0;
            end
            if obj.useBiasTfScalar
                % NOTE: Overwrites "BiasItfAvpiv".
                BiasItfAvpiv = @(omegaRps) (ones(size(omegaRps)) / kFtfIvpav);
            end
            %###################################################################
            
            
            
            BiasCalibData.itfAvpiv    = BiasItfAvpiv;
            BiasCalibData.offsetAVolt = offsetAVolt;



            % Convenience function: Cell array of TFs --> Time-relevant TF (as function handle)
            %
            % RFTF = Rational Function (rational_func_transform object) Transfer Function
            function Tf = TF_list_2_func(RftfList)
                Tf = @(omegaRps) (RftfList{iCalibTimeL}.eval(omegaRps));
            end
        end
        
        
        
        % Obtain LFR ITF, but handle the case that should never happen for
        % actual non-NaN data (LSF F3 + BLTS 4 or 5) and return a TF that only
        % returns NaN instead. BICAS may still iterate over that combination
        % though when calibrating.
        % 
        function lfrItfIvpt = get_LFR_ITF(obj, cti1, iBlts, iLsf)
            % ASSERTIONS
            assert(cti1 >= 0)
            bicas.calib.assert_iBlts(iBlts)
            bicas.calib.assert_iLsf(iLsf)
            
            if (iLsf == 4) && ismember(iBlts, [4,5])
                lfrItfIvpt = bicas.calib.NAN_TF;
            else
                RctDataList = obj.RctDataMap('LFR');
                
                % ASSERTION
                % IMPLEMENTATION NOTE: Anonymous function below will fail at a
                % later stage if these assertions are false. Checking for these
                % criteria here makes it easier to understand these particular
                % types of error.
                assert(numel(RctDataList) >= (cti1+1), ...
                    'BICAS:calib:IllegalArgument:DatasetFormat:Assertion', ...
                    ['LFR RctDataList is too small for argument cti1=%g.', ...
                    ' This could indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
                    ' value is larger than g.attr. CALIBRATION TABLE allows.'], cti1)
                assert(~isempty(RctDataList{cti1+1}), ...
                    'BICAS:calib:IllegalArgument:DatasetFormat:Assertion', ...
                    ['LFR RctDataList contains no RCT data corresponding to argument cti1=%g.', ...
                    ' This may indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
                    ' value is wrong or that BICAS did not try to load the corresponding', ...
                    ' RCT in g.attr. CALIBRATION_TABLE.'], cti1)
                
                lfrItfIvpt = @(omegaRps) (bicas.calib.eval_tabulated_ITF(...
                    RctDataList{cti1+1}.ItfIvptTable{iLsf}{iBlts}, omegaRps));
            end
        end



    end    % methods(Access=private)

    %###################################################################################################################

    
    
    methods(Static, Access=public)



        % Load all non-BIAS RCTs (all types) using assumptions on filenames.
        %
        % NOTE: Can be useful for manual experimentation with calibration.
        % NOTE: Necessary when processing L1-->L2 (inofficially) since L1 does
        %       not have CALIBRATION_TABLE+CALIBRATION_TABLE_INDEX.
        % NOTE: Will only load one of each RCT type (no potential time
        %       dependence as per global attribute CALIBRATION_TABLE) and
        %       requires user to not use CALIBRATION_TABLE_INDEX.
        %
        % IMPLEMENTATION NOTE: BICAS only needs one non-BIAS RCT type at a time.
        % However, it is useful to be able to initialize bicas.calib so that it
        % can simultanteously calibrate all kinds of data for debugging
        % purposes. Therefore loads ALL non-BIAS RCT types.
        %
        % RETURN VALUE
        % ============
        % RctDataMap : containers.Map. One key per non-BIAS RCT type ID. Value =
        %              1x1 cell array with RCT data.
        %              IMPLEMENTATION NOTE: Returns containers.Map to provide
        %              the same interface to bicas.calib constructor as
        %              bicas.calib.find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE.
        % 
        function RctDataMap = find_read_non_BIAS_RCTs_by_regexp(rctDir, SETTINGS, L)
            
            RctDataMap = containers.Map();
            
            for rctTypeId = {'LFR', 'TDS-CWF', 'TDS-RSWF'}
                
                settingKey     = bicas.calib.RCT_TYPES_MAP(rctTypeId{1}).filenameRegexpSettingKey;
                filenameRegexp = SETTINGS.get_fv(settingKey);
                filePath       = bicas.RCT.find_RCT_regexp(rctDir, filenameRegexp, L);
                RctDataList    = {bicas.calib.read_log_modify_RCT(rctTypeId{1}, filePath, L)};
                
                % NOTE: Placing all non-BIAS RCT data inside 1x1 cell arrays so
                % that they are stored analogously with when using ga.
                % CALIBRATION_TABLE.
                RctDataMap(rctTypeId{1}) = RctDataList;
            end
        end



        % Load non-BIAS RCT(s) of ONE type (rctTypeId) using CDF global
        % attribute CALIBRATION_TABLE and zVars CALIBRATION_TABLE_INDEX and BW.
        %
        % IMPLEMENTATION NOTE: May load multiple RCTs (of the same RCT type) but
        % will only load those RCTs which are actually needed, as indicated by
        % CALIBRATION_TABLE_INDEX and BW. This is necessary since
        % CALIBRATION_TABLE may reference unnecessary RCTs of types not
        % recognized by BICAS (LFR's ROC-SGSE_CAL_RCT-LFR-VHF_V01.cdf
        % /2019-12-16), and which are therefore unreadable by BICAS (BICAS would
        % crash).
        %
        %
        % ARGUMENTS
        % =========
        % ga_CALIBRATION_TABLE       : LFR/TDS RCT global attribute
        %                              CALIBRATION_TABLE. 1D cell array of
        %                              strings.
        % zv_CALIBRATION_TABLE_INDEX : LFR/TDS BICAS input dataset zVariable
        %                              CALIBRATION_TABLE_INDEX.
        % zv_BW                      : Either
        %                               (1) [] (as for TDS data), or
        %                               (2) LFR input dataset zVariable BW.
        %
        %
        % RETURN VALUE
        % ============
        % RctDataMap : containers.Map with
        %               keys   = non-BIAS RCT type ID.
        %               values = 1D cell. Non-empty indices {iRct}
        %              come from zv_CALIBRATION_TABLE_INDEX(i,1). Each element
        %              is the content of the corresponding RCT mentioned in
        %              ga_CALIBRATION_TABLE.
        %              IMPLEMENTATION NOTE: Returns containers.Map to provide
        %              the same interface to bicas.calib constructor as
        %              bicas.calib.find_read_non_BIAS_RCTs_by_regexp.
        %
        function RctDataMap = find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE(...
                rctDir, rctTypeId, ga_CALIBRATION_TABLE, zv_CALIBRATION_TABLE_INDEX, zv_BW, L)
            
            % ASSERTION
            assert(iscell(ga_CALIBRATION_TABLE))
            
            if isempty(zv_BW)
                % ASSERTION
                nCt = EJ_library.assert.sizes(...
                    ga_CALIBRATION_TABLE,       [-1, 1], ...
                    zv_CALIBRATION_TABLE_INDEX, [-2, 2]);
                
                iCtArray = unique(zv_CALIBRATION_TABLE_INDEX(:, 1));     % CT = CALIBRATION_TABLE
                
            else
                % ASSERTIONS
                nCt = EJ_library.assert.sizes(...
                    ga_CALIBRATION_TABLE,       [-1, 1], ...
                    zv_CALIBRATION_TABLE_INDEX, [-2, 2], ...
                    zv_BW,                      [-2, 1]);
                assert(all(ismember(zv_BW, [0,1])))
                
                iCtArray = unique(zv_CALIBRATION_TABLE_INDEX(logical(zv_BW), 1));
            end
            
            
            
            % Cell array of paths to RCTs of the same RCT type.
            RctDataList = cell(nCt, 1);
            
            % NOTE: Iterate over those entries in CALIBRATION_TABLE that should be considered, NOT all indices. May
            % therefore legitimately leave some cells in cell array empty.
            for i = 1:numel(iCtArray)
                j              = iCtArray(i) + 1;   % NOTE: Cell array index is one greater than the stored value.
                filePath       = fullfile(rctDir, ga_CALIBRATION_TABLE{j});
                RctDataList{j} = bicas.calib.read_log_modify_RCT(rctTypeId, filePath, L);
            end
            
            RctDataMap = containers.Map();
            RctDataMap(rctTypeId) = RctDataList;
        end


        
        % Given a sequence of Epoch values, determine for each value which
        % calibration time index should be used. The caller will have to decide
        % which sequence of data that should be calibrated together (e.g. if
        % calibration time changes in the middle of CWF), and which Epoch values
        % should be used to determine calibration time (e.g. first Epoch value
        % for a snapshot determines entire snapshot).
        %
        % NOTE: Method is public so that automatic test code can call
        % get_calibration_time.
        %
        %
        % ARGUMENTS AND RETURN VALUES
        % ===========================
        % Epoch          : Column vector with Epoch values.
        % CalibEpochList : List of monotonically increasing timestamps ("Epoch
        %                  format"). In practice intended to be Bias.epochL or
        %                  Bias.epochH.
        % iCalib         : Array. iCalibList(i) = calibration time index for
        %                  Epoch(i).
        %
        function [iCalib] = get_calibration_time(Epoch, CalibEpochList)
            
            % ASSERTIONS
            bicas.proc_utils.assert_zv_Epoch(Epoch)
            bicas.proc_utils.assert_zv_Epoch(CalibEpochList)
            % IMPLEMENTATION NOTE: Does not work if CalibEpochList is empty,
            % since discretize behaves differently for scalar second argument.
            assert(~isempty(CalibEpochList))
            
            % IMPLEMENTATION NOTE: "discretize" by itself returns NaN for Epoch
            % values outside the outermost edges. Therefore (1) must add upper
            % edge "Inf", (2) asserts non-Nan afterwards.
            % IMPLEMENTATION NOTE: "discretize" behaves differently for scalar
            % second argument. Adding edges at infinity hides this problem. If
            % one does not add infinities and uses a scalar edge list, then one
            % has to treat those cases manually.
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
        function biasCurrentTm = calibrate_current_sampere_to_TM(currentSAmpere)
            
            % ASSERTION
            % NOTE: max(...) ignores NaN, unless that is the only value, which
            % then becomes the max value.
            [maxAbsSAmpere, iMax] = max(abs(currentSAmpere(:)));
            if ~(isnan(maxAbsSAmpere) || (maxAbsSAmpere <= EJ_library.so.constants.MAX_ABS_SAMPERE))
                
                error('BICAS:calib:Assertion:IllegalArgument', ...
                    ['Argument currentSAmpere (unit: set current/ampere) contains illegally large value(s).', ...
                    ' Largest found value is %g.'], ...
                    currentSAmpere(iMax))
            end
            
            biasCurrentTm = currentSAmpere * EJ_library.so.constants.TM_PER_SAMPERE;
        end



    end    % methods(Static, Access=public)



    methods(Static, Access=private)
        
        
        
        % Code to initialize hard-coded static constant.
        function RctTable = init_RCT_Types_Map()
            RctTable = containers.Map();
            
            NOTHING_FUNC = @(S) (S);
            
            RctTable('BIAS')     = RCT_entry(@bicas.RCT.read_BIAS_RCT,     @bicas.calib.modify_BIAS_RCT_data,     @bicas.calib.log_BIAS_RCT,      'PROCESSING.RCT_REGEXP.BIAS');
            RctTable('LFR')      = RCT_entry(@bicas.RCT.read_LFR_RCT,      @bicas.calib.modify_LFR_RCT_data,      @bicas.calib.log_LFR_RCTs,      'PROCESSING.RCT_REGEXP.LFR');
            RctTable('TDS-CWF')  = RCT_entry(@bicas.RCT.read_TDS_CWF_RCT,  NOTHING_FUNC,                          @bicas.calib.log_TDS_CWF_RCTs,  'PROCESSING.RCT_REGEXP.TDS-LFM-CWF');
            RctTable('TDS-RSWF') = RCT_entry(@bicas.RCT.read_TDS_RSWF_RCT, @bicas.calib.modify_TDS_RSWF_RCT_data, @bicas.calib.log_TDS_RSWF_RCTs, 'PROCESSING.RCT_REGEXP.TDS-LFM-RSWF');
            
            %###################################################################
            function Entry = RCT_entry(readRctFunc, modifyRctFunc, logRctFunc, filenameRegexpSettingKey)
                Entry = struct(...
                    'readRctFunc',              readRctFunc, ...      % Pointer to function that reads one RCT.
                    'modifyRctFunc',            modifyRctFunc, ...    % Pointer to function that modifies data for one RCT.
                    'logRctFunc',               logRctFunc, ...       % Pointer to function that logs data for one RCT.
                    'filenameRegexpSettingKey', filenameRegexpSettingKey);
            end
            %###################################################################
        end
        
        
        
        function RctData = modify_BIAS_RCT_data(RctData)
            FtfSet = RctData.FtfSet;
            
            % ASSERTIONS
            nTime = EJ_library.assert.sizes(...
                FtfSet.DcSingleAvpiv,   [-1, 1], ...
                FtfSet.DcDiffAvpiv,     [-1, 1], ...
                FtfSet.AcLowGainAvpiv,  [-1, 1], ...
                FtfSet.AcHighGainAvpiv, [-1, 1]);
            
            % NOTE: Derive ITFs.
            ItfSet = [];
            for iTf = 1:nTime
                % INVERT: FTF --> ITF
                ItfSet.DcSingleAvpiv{  iTf} = FtfSet.DcSingleAvpiv{  iTf}.inverse();
                ItfSet.DcDiffAvpiv{    iTf} = FtfSet.DcDiffAvpiv{    iTf}.inverse();
                ItfSet.AcLowGainAvpiv{ iTf} = FtfSet.AcLowGainAvpiv{ iTf}.inverse();
                ItfSet.AcHighGainAvpiv{iTf} = FtfSet.AcHighGainAvpiv{iTf}.inverse();
            end
            
            RctData.ItfSet = ItfSet;
        end
        
        
            
        function RctData = modify_LFR_RCT_data(RctData)
            
            FtfTpivTable = RctData.FtfTpivTable;
            
            % Read LFR FTFs, derive ITFs and modify them.
            ItfIvptTable = {};
            for iLsf = 1:numel(FtfTpivTable)
                
                ItfIvptTable{end+1} = {};
                for iBlts = 1:numel(FtfTpivTable{iLsf})

                    % INVERT: FTF --> ITF
                    ItfIvpt = FtfTpivTable{iLsf}{iBlts}.inverse();
                    % MODIFY ITF
                    ItfIvptTable{iLsf}{iBlts} = bicas.calib.modify_tabulated_ITF(ItfIvpt);
                    
                end                
            end
            
            RctData.ItfIvptTable = ItfIvptTable;            
        end
        
        
        
        function RctData = modify_TDS_RSWF_RCT_data(RctData)
            
            % Modify tabulated TDS-RSWF TFs.
            for iBlts = 1:numel(RctData.ItfIvptList)
                RctData.ItfIvptList{iBlts} = bicas.calib.modify_tabulated_ITF(RctData.ItfIvptList{iBlts});
            end
            
        end



        % Evaluate tabulated INVERSE transfer function.
        %
        % NOTE: This function is effectively meant to specify how tabulated
        % transfer functions should be interpreted w.r.t. interpolation.
        % NOTE: Intended specifically for INVERSE transfer functions. Therefore
        % using Z=0 for frequencies higher than the table.
        %
        function Z = eval_tabulated_ITF(TabulatedItf, omegaRps)
            useTabTf = (omegaRps <= TabulatedItf.omegaRps(end));
            
            % NOTE: interp1 returns NaN for values outside range.
            Z = interp1(TabulatedItf.omegaRps, TabulatedItf.Z, omegaRps, 'linear');
            Z(~useTabTf) = 0;   % Set to zero (overwrite) for values above highest tabulated frequency.
            
            % ASSERTION
            if ~all(isfinite(Z))
                % IMPLEMENTATION NOTE: Experience shows that it is useful to
                % have an extended error message confirming that the requested
                % frequence range is outside the tabulated one, and by how much.
                errorMsg = sprintf(...
                    ['Can not evaluate tabulated transfer function for', ...
                    ' frequencies outside of the range of tabulated frequencies.\n', ...
                    'Range of frequencies for which there are tabulated Z values:\n', ...
                    '    min(TabulatedTf.omegaRps) = %g\n', ...
                    '    max(TabulatedTf.omegaRps) = %g\n', ...
                    'Range of frequencies for which evaluation (interpolation) of Z was attempted:\n', ...
                    '    min(omegaRps)     = %g\n', ...
                    '    max(omegaRps)     = %g\n'], ...
                    min(TabulatedItf.omegaRps), ...
                    max(TabulatedItf.omegaRps), ...
                    min(omegaRps), ...
                    max(omegaRps));
                
                error('BICAS:Assertion', errorMsg)
            end
        end
        
        
        
        % Modify TABULATED INVERSE transfer functions, if needed.
        % Tabulated TF --> Tabulated TF
        %
        % Extrapolate to 0 Hz, if needed.
        %
        % NOTE: Does NOT remove high frequencies.
        %
        function ModifItf = modify_tabulated_ITF(Itf)
            assert(Itf.omegaRps(1) > 0)
            assert(isa(Itf, 'EJ_library.utils.tabulated_transform'))

            % NOTE: Can not just use the lowest-frequency Z value for 0 Hz since
            % it has to be real (not complex).
            Z1       = Itf.Z(1);
            signZ0   = sign(real(Z1));
            assert(signZ0 ~= 0, ...
                'BICAS:calib:modify_tabulated_ITF:FailedToReadInterpretRCT:Assertion', ...
                ['Can not extrapolate tabulated inverse transfer function', ...
                ' (ITF) to zero Hz due to ambiguity. real(Z(1)) = 0.'])
            Z0       = abs(Z1) * signZ0;   % Z value at 0 Hz.
            
            omegaRps = [0;  Itf.omegaRps(:)];
            Z        = [Z0; Itf.Z(:)       ];
            
            ModifItf = EJ_library.utils.tabulated_transform(omegaRps, Z);
        end
        


        % Log some indicative value(s) for a BIAS RCT.
        % NOTE: Does not log file path. Caller is assumed to do that.
        function log_BIAS_RCT(RctData, L)
            
            DC_FREQ_HZ       = 0;
            AC_DIFF_FREQS_HZ = [0, 1000];
            LL               = bicas.calib.RCT_DATA_LL;
            
            for iEpochL = 1:numel(RctData.epochL)
                
                L.logf(LL, 'Below values are used for data beginning %s:', EJ_library.cdf.tt2000_to_UTC_str(RctData.epochL(iEpochL)))
                
                % Log bias current calibration
                L.logf(LL, '    BIAS current offsets: %s [aampere]',         bicas.calib.vector_string('% 10e', RctData.Current.offsetsAAmpere(iEpochL, :)))
                L.logf(LL, '    BIAS current gain   : %s [aampere/TM unit]', bicas.calib.vector_string('% 10e', RctData.Current.gainsAapt(     iEpochL, :)))
                
                % Log transfer functions (frequency domain), selected frequencies.
                log_TF('    BIAS ITF DC single',          DC_FREQ_HZ,       RctData.ItfSet.DcSingleAvpiv)
                log_TF('    BIAS ITF DC diff',            DC_FREQ_HZ,       RctData.ItfSet.DcDiffAvpiv)
                
                log_TF('    BIAS ITF AC diff, low  gain', AC_DIFF_FREQS_HZ, RctData.ItfSet.AcLowGainAvpiv)
                log_TF('    BIAS ITF AC diff, high gain', AC_DIFF_FREQS_HZ, RctData.ItfSet.AcHighGainAvpiv)
            end
            
            % NOTE: Must work for multiple CDF records.
            dcDiffOffsetsAVolt = [...
                RctData.DcDiffOffsets.E12AVolt, ...
                RctData.DcDiffOffsets.E13AVolt, ...
                RctData.DcDiffOffsets.E23AVolt];
            EJ_library.assert.sizes(dcDiffOffsetsAVolt, [NaN, 3]);            
            for iEpochH = 1:numel(RctData.epochH)
                L.logf(LL, 'Below values are used for data beginning %s:', ...
                    EJ_library.cdf.tt2000_to_UTC_str(RctData.epochH(iEpochH)))
                
                L.logf(LL, '    BIAS DC single voltage offsets ( V1, V2, V3): %s [avolt]', ...
                    bicas.calib.vector_string('%g', RctData.dcSingleOffsetsAVolt(iEpochH, :)))
                L.logf(LL, '    BIAS DC diff   voltage offsets (E12,E13,E23): %s [avolt]', ...
                    bicas.calib.vector_string('%g', dcDiffOffsetsAVolt(iEpochH)))

            end
                
            %###################################################################
            % Nested utility function.
            % NOTE: Impicitly function of iEpochL, L, LL.
            function log_TF(name, freqArray, ItfList)
                bicas.calib.log_TF(...
                    LL, name, 'avolt/ivolt', freqArray, ...
                    @(omegaRps) (ItfList{iEpochL}.eval(omegaRps)), L);
            end
            %###################################################################
        end


        
        % Analogous to log_BIAS_RCT.
        function log_LFR_RCTs(RctData, L)
            FREQ_HZ = 0;
            
            % CASE: This index corresponds to an actually loaded RCT (some are
            % intentionally empty).
            for iLsf = 1:4
                if iLsf ~= 4   nBltsMax = 5;
                else           nBltsMax = 3;
                end
                
                for iBlts = 1:nBltsMax
                    TabulatedItfIvpt = RctData.ItfIvptTable{iLsf}{iBlts};
                    ItfIvpt          = @(omegaRps) (bicas.calib.eval_tabulated_ITF(TabulatedItfIvpt, omegaRps));
                    bicas.calib.log_TF(...
                        bicas.calib.RCT_DATA_LL, ...
                        sprintf('LFR ITF, F%i, BLTS/BIAS_%i', iLsf-1, iBlts), ...
                        'ivolt/TM unit', FREQ_HZ, ItfIvpt, L)
                end
            end    % for
            
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_CWF_RCTs(RctData, L)
            
            L.logf(bicas.calib.RCT_DATA_LL, 'TDS CWF calibration factors: %s [ivolt/TM]', ...
                bicas.calib.vector_string('%g', RctData.factorsIvpt));
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_RSWF_RCTs(RctData, L)
            
            FREQ_HZ = 0;
            
            for iBlts = 1:3
                TabulatedItfIvpt = RctData.ItfIvptList{iBlts};
                ItfIvpt          = @(omegaRps) (bicas.calib.eval_tabulated_ITF(TabulatedItfIvpt, omegaRps));
                bicas.calib.log_TF(...
                    bicas.calib.RCT_DATA_LL, ...
                    sprintf('TDS RSWF ITF, BLTS/BIAS_%i', iBlts), ...
                    'ivolt/TM unit', FREQ_HZ, ItfIvpt, L)
            end
        end



        % Read any single RCT file, and log it. Effectively wraps the different
        % RCT-reading functions.
        % 
        %
        % IMPLEMENTATION NOTES
        % ====================
        % This method exists to
        % (1) run shared code that should be run when reading any RCT (logging,
        %     modifying data),
        % (2) separate the logging from the RCT-reading code, so that external
        %     code can read RCTs without BICAS.
        %
        %
        % ARGUMENTS
        % =========
        % rctTypeId : String constants representing pipeline and RCT to be read.
        %
        function RctData = read_log_modify_RCT(rctTypeId, filePath, L)
            
            L.logf(bicas.calib.READING_RCT_PATH_LL, 'Reading RCT (rctTypeId=%s): "%s"', rctTypeId, filePath)
            
            readRctFunc   = bicas.calib.RCT_TYPES_MAP(rctTypeId).readRctFunc;
            modifyRctFunc = bicas.calib.RCT_TYPES_MAP(rctTypeId).modifyRctFunc;
            logRctFunc    = bicas.calib.RCT_TYPES_MAP(rctTypeId).logRctFunc;
            
            RctData = readRctFunc(filePath);
            RctData = modifyRctFunc(RctData);
            logRctFunc(RctData, L);
        end
        
        
        
        % ARGUMENTS
        % =========
        % freqHzArray : Array of frequencies for which the ITF value should be logged.
        % TfFuncPtr   : Function pointer. Z(omegaRps).
        %
        function log_TF(logLevel, tfName, tfUnit, freqHzArray, tfFuncPtr, L)
            
            zArray = tfFuncPtr(freqHzArray);
            for i=1:numel(freqHzArray)
                freqHz = freqHzArray(i);
                Z      = zArray(i);
                
                inverseZValueStr = sprintf('1/%10.5f', 1/abs(Z));
                
                %======================================================================================================
                % NOTE 2020-04-30: Execution at ROC fails due to not finding
                % function "phase" for unknown reason.
                % --------------------------------------------------------------------------------------
                % Exception.identifier = "MATLAB:UndefinedFunction"
                % Exception.message    = "Undefined function 'phase' for input arguments of type 'double'."
                % Matching MATLAB error message identifier parts (error types derived from Exception1.identifier):
                %     UntranslatableErrorMsgId : Error occurred, but code can not translate the error's MATLAB message
                %     identifier into any of BICAS's internal standard error codes.
                % MATLAB call stack:
                %     row  969, calib.m,                    calib.log_ITF_Z
                %     row  303, calib.m,                    calib.calib
                %     row   68, execute_sw_mode.m,          execute_sw_mode
                %     row  455, main.m,                     main_without_error_handling
                %     row  116, main.m,                     main
                % --------------------------------------------------------------------------------------
                % See also
                % https://se.mathworks.com/matlabcentral/answers/408657-which-toolbox-is-phase-in
                % """"phase() as a routine by itself is part of the System
                % Identification Toolbox, in the "obsolete" category. phase() is
                % also a method of the newer iddata() class from the System
                % Identification Toolbox. But what you probably want is angle()
                % followed by unwrap(), which is part of basic MATLAB.""""
                %
                % Therefore using function "angle" instead of "phase".
                %======================================================================================================
                assert(numel(tfName) <= 31)    % Check that string is not too long for neat printouts.
                L.logf(logLevel, '%-31s %4i Hz: abs(Z)=%8.5f=%12s [%s], phase(Z)=% 6.1f [deg]', ...
                    tfName, freqHz, abs(Z), inverseZValueStr, tfUnit, rad2deg(angle(Z)))
            end    % for
        end
        
        
        
        % Utility function for creating string representing 1D vector.
        % Ex: '(3.1416, 2.7183, 1.6180)'
        function s = vector_string(pattern, v)
            assert(~isempty(v))
            EJ_library.assert.vector(v)
            s = sprintf('(%s)', strjoin(EJ_library.str.sprintf_many(pattern, v), ', '));
        end
        
        
        
        function assert_iBlts(iBlts)
            assert(ismember(iBlts, [1:5]), 'BICAS:calib:IllegalArgument:Assertion', 'Illegal value iBlts=%g', iBlts)
        end
        
        
        
        function assert_iLsf(iLsf)
            assert(ismember(iLsf,  [1:4]), 'BICAS:calib:IllegalArgument:Assertion', 'Illegal value iLsf=%g.', iLsf)
        end
        
        
        
    end    % methods(Static, Access=private)

end
