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
% TPIV     = TM/interface volt (=TM per interface volt)
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
% PROPOSAL: Add TF for (arbitrary) capacitance. (Needed for ~debugging/testing.)
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
% PROPOSAL: Assertion function for CalSettings.
%   TODO-NI: Same struct, with same fields in all cases?
%   NOTE: Function does not know which fields are actually used.
%
% TODO-DECISION: How distribute the calibration formulas/algorithms between
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
%   NEED: Plot all TFs used for a particular calibration case (when calibrating using bicas.caib only, without BICAS).
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
% PROPOSAL: Store all versions of TFs internally.
%   Ex: FTF, ITF, tabulated ITF with extrapolation+interpolation+modification
%   NOTE: Modification (besides inversion) happens on the combined function
%   which is not stored beforehand.
%   NOTE: The set of BIAS+LFR+TDS TFs is different from the set of TFs actually
%           used (combinations of BIAS+LFR and BIAS+TDS respectively).
%
% PROPOSAL: Separate function for interpolation (no extrapolation) of tabulated TF.
%
% PROPOSAL: Store both FTFs and ITFs, despite that FTFs are not used for calibration directly.
%   NOTE: BIAS & LFR RCTs contain FTFs, TDS RCT contains ITFs.
%   NOTE: Has to keep track of FTF/ITF before modifications (extrapolation to 0 Hz, Z=0 for high freq.).
%   PRO: Useful for debugging. Can easily inspect & plot FTFs.
%
% PROPOSAL: Do not expose internal calibration data structs and let the caller specify indices. Access via
%   methods that ask for indices.
%   PRO: Methods document the proper use of indices.
%
% PROPOSAL: Move out static functions, at least public static functions.
%   PROPOSAL: Class calib_utils?
%
% PROPOSAL: Store LFR TFs as one 1D array of structs with fields: iLsf, iBlts, ftf, itf, ...
%   PRO: Can easily iterate over.
%   PRO: For every modification of TFs, can easily add another field for the old
%       version.
%       NOTE/CON: All structs/TFs must have the same fields if true struct array.
%
% PROPOSAL: Have calibrate_LFR, calibrate_TDS_CWF, calibrate_TDS_RSWF return
%           function handle to a transfer function that does everything,
%           including offsets.
%   CON: Useful for debugging.
%       CON: Aren't the calibration methods available already that, if turned
%            into function handles?!!
%   PROPOSAL: Have ~special functions/methods for this, so that one does not use
%             function handles wrapped in function handles (not too many
%             anyway).
%   NOTE: CLARIFICATION: Separate TFs with offsets and calibration metods in the
%         time domain.
%   PROPOSAL: Have special function that returns transfer function for every
%             case.
%       PRO: Useful for debugging.
%           PRO: Plot
%           PRO: Save to file.



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
        % ------------------------
        % BIAS scalar (simplified) calibration, not in the RCTs. For
        % debugging/testing purposes.
        BiasScalarGain             
        HkBiasCurrent
        % EXPERIMENTAL. NOTE: Technically, the name contains a tautology
        % (LFR+LSF).
        lfrLsfOffsetsTm    = [];
        
        
        
        %==================================================
        % Settings for what kind of calibration to perform
        %==================================================
        
        % Corresponds to SETTINGS key-value.
        dcDetrendingDegreeOf
        dcRetrendingEnabled
        acDetrendingDegreeOf
        itfHighFreqLimitFraction;
                
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

            % ASSERTIONS
            assert(isscalar(use_CALIBRATION_TABLE_INDEX2))
            assert(~ismember('BIAS', NonBiasRctDataMap.keys))
            EJ_library.assert.subset(NonBiasRctDataMap.keys, bicas.calib.RCT_TYPES_MAP.keys)
            
            %=======================
            % Assign obj.RctDataMap
            %=======================
            obj.RctDataMap = containers.Map();
            % ----------------------------------
            % Create RCT data map entry for BIAS
            % ----------------------------------
            filenameRegexp = SETTINGS.get_fv(bicas.calib.RCT_TYPES_MAP('BIAS').filenameRegexpSettingKey);
            filePath       = bicas.RCT.find_RCT_regexp(rctDir, filenameRegexp, L);
            % IMPLEMENTATION NOTE: It has been observed that value sometimes
            % survives from previous runs, despite being an instance variable.
            % Unknown why. Therefore explicitly overwrites it.
            obj.RctDataMap('BIAS') = bicas.calib.read_RCT_modify_log('BIAS', filePath, L);
            % -------------------------------------
             % Add RCT data map entry for non-BIAS
             % -------------------------------------
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



            obj.dcDetrendingDegreeOf               = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE');
            obj.dcRetrendingEnabled                = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED');
            obj.acDetrendingDegreeOf               = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE');
            obj.itfHighFreqLimitFraction           = SETTINGS.get_fv('PROCESSING.CALIBRATION.TF_HIGH_FREQ_LIMIT_FRACTION');
            
            obj.allVoltageCalibDisabled            = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.DISABLE');
            obj.biasOffsetsDisabled                = SETTINGS.get_fv('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED');
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
        % NOTE: This is the normal way of obtaining bias current in physical
        % units (as opposed to HK bias current).
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
        function biasCurrentAAmpere = calibrate_HK_bias_TM_to_bias_current(obj, biasCurrentTm, iAntenna)
            
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
            biasCurrentTm     = bitxor(biasCurrentTm, hex2dec('8000'));    % FLIP BIT
            biasCurrentAAmpere = obj.HkBiasCurrent.gainAapt(iAntenna) * ...
                (biasCurrentTm + obj.HkBiasCurrent.offsetTm(iAntenna));    % LINEAR FUNCTION
        end
        
        
        
        % Calibrate all voltages. Function will choose the more specific
        % algorithm internally.
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
            assert(isstruct(CalSettings))
            %EJ_library.assert.struct(CalSettings, {'iBlts', 'BltsSrc', 'biasHighGain', 'iCalibTimeL', 'iCalibTimeH', 'iLsf'}, {})   % Too slow?
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
            bicas.calib_utils.assert_iBlts(iBlts)
            assert(isa(BltsSrc, 'bicas.BLTS_src_dest'))
            bicas.calib_utils.assert_iLsf(iLsf)
            assert(isscalar(cti1))
            assert(cti1 >= 0, 'Illegal cti1=%g', cti1)
            % No assertion on cti2 unless used (determined later).



            %============================================
            % Only place to potentially make use of cti2
            %============================================
            if obj.use_CALIBRATION_TABLE_INDEX2
                % ASSERTIONS
                assert(isscalar(cti2), ...
                    'BICAS:calib:IllegalArgument:Assertion', ...
                    'Argument cti2 is not scalar.')
                assert(cti2 >= 0, ...
                    'BICAS:calib:IllegalArgument:Assertion', ...
                    'Illegal argument cti2=%g (=zVar CALIBRATION_TABLE_INDEX(iRecord, 2))', cti2)
                assert(iLsf == cti2+1, ...
                    'BICAS:calib:IllegalArgument:Assertion', ...
                    'cti2+1=%i != iLsf=%i (before overwriting iLsf)', cti2+1, iLsf)
                
                % NOTE: Only place cti2 is used.
                iLsf = cti2 + 1;
            end



            %==============================
            % Obtain calibration constants
            %==============================
            BiasCalibData = obj.get_BIAS_calib_data(BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            biasItfAvpiv = BiasCalibData.itfAvpiv;
            if obj.lfrTdsTfDisabled
                lfrItfIvpt = @(omegaRps) (ones(size(omegaRps)));
            else
                lfrItfIvpt = obj.get_LFR_ITF(cti1, iBlts, iLsf);
            end
            
            % TEST: Change sign of BIAS TF.
%             if 0
%                 biasItfAvpiv = @(omegaRps) (-biasItfAvpiv(omegaRps));
%             end

            %======================================
            % Create combined ITF for LFR and BIAS
            %======================================
            itfIvpt = @(omegaRps) (bicas.calib_utils.multiply_TFs(...
                omegaRps, lfrItfIvpt, biasItfAvpiv));



            % TEST
%             if 0
%                 fHz = 100;
%                 zLimit = itfIvpt(fHz * 2*pi);
%                 %zLimit = 1;
%                 itfIvpt = @(omegaRps) (bicas.calib_utils.TF_LF_constant_abs_Z(itfIvpt, omegaRps, fHz, zLimit));
%             end
            
            % TEST DEBUG: Invert ITF (again) to FTF.
            %itfIvpt = @(omegaRps) (1./itfIvpt(omegaRps));
            
            

            if BltsSrc.is_AC()
                % IMPLEMENTATION NOTE: DC is (optionally) detrended via
                % bicas.utils.apply_TF_freq_modif in the sense of a linear fit
                % being removed, TF applied, and then added back. That same
                % algorithm is inappropriate for non-lowpass filters.
                detrendingDegreeOf = obj.acDetrendingDegreeOf;
                retrendingEnabled  = 0;
            else
                detrendingDegreeOf = obj.dcDetrendingDegreeOf;
                retrendingEnabled  = obj.dcRetrendingEnabled;
            end
            
            %=======================================
            % CALIBRATE: LFR TM --> TM --> avolt
            %=======================================
            samplesCaAVolt = cell(size(samplesCaTm));
            lsfOffsetTm = obj.lfrLsfOffsetsTm(iLsf);
            for i = 1:numel(samplesCaTm)
                
                % ADD LSF OFFSET
                samplesTm = samplesCaTm{i}(:) + lsfOffsetTm;
                
                % APPLY TRANSFER FUNCTION (BIAS + LFR)
                tempSamplesAVolt = bicas.utils.apply_TF_freq_modif(...
                    dtSec(i), samplesTm, itfIvpt, ...
                    'detrendingDegreeOf',      detrendingDegreeOf, ...
                    'retrendingEnabled',       retrendingEnabled, ...
                    'tfHighFreqLimitFraction', obj.itfHighFreqLimitFraction);

                % ADD BIAS offset
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
            bicas.calib_utils.assert_iBlts(iBlts)
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
                    tempSamplesAVolt = bicas.utils.apply_TF_freq_modif(...
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
            bicas.calib_utils.assert_iBlts(iBlts)
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
            % NOTE: Low/high gain is irrelevant for TDS. Argument value
            % arbitrary.
            BiasCalibData = obj.get_BIAS_calib_data(...
                BltsSrc, biasHighGain, iCalibTimeL, iCalibTimeH);
            
            samplesCaAVolt = cell(size(samplesCaTm));   % Initialize empty output variable.
            if ismember(iBlts, [1,2,3])
                
                %======================================
                % Create combined ITF for TDS and BIAS
                %======================================
                if obj.lfrTdsTfDisabled
                    tdsItfIvpt = @(omegaRps) (ones(omegaRps));
                else
                    RctList = obj.RctDataMap('TDS-RSWF');
                    %tdsItfIvpt = @(omegaRps) (bicas.calib_utils.eval_tabulated_TF(...
                    %    RctList{cti1+1}.ItfModifIvptList{iBlts}, omegaRps));
                    tdsItfIvpt = RctList{cti1+1}.itfModifIvptCa{iBlts};
                end

                itf = @(omegaRps) (...
                    tdsItfIvpt(omegaRps) ...
                    .* ...
                    BiasCalibData.itfAvpiv(omegaRps));
                
                %====================================
                % CALIBRATE: TDS TM --> antenna volt
                %====================================
                % APPLY TRANSFER FUNCTION (BIAS + TDS-RSWF)
                for i = 1:numel(samplesCaTm)
                    tempSamplesAVolt = bicas.utils.apply_TF_freq_modif(...
                        dtSec(i), ...
                        samplesCaTm{i}(:), ...
                        itf, ...
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



        function iCalib = get_calibration_time_L(obj, Epoch)
            iCalib = bicas.calib.get_calibration_time(...
                Epoch, obj.RctDataMap('BIAS').epochL);
        end



        function iCalib = get_calibration_time_H(obj, Epoch)
            iCalib = bicas.calib.get_calibration_time(...
                Epoch, obj.RctDataMap('BIAS').epochH);
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



            % Convenience function:
            % Cell array of TFs --> Time-relevant TF (as function handle)
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
            bicas.calib_utils.assert_iBlts(iBlts)
            bicas.calib_utils.assert_iLsf(iLsf)
            
            if (iLsf == 4) && ismember(iBlts, [4,5])
                % CASE: F3 and BLTS={4,5}
                
                % NOTE: There is no tabulated LFR TF and no such combination
                % signal route, so the TF can not be returned even in principle.
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
                    ' value is larger than glob. attr. CALIBRATION TABLE allows.'], cti1)
                assert(~isempty(RctDataList{cti1+1}), ...
                    'BICAS:calib:IllegalArgument:DatasetFormat:Assertion', ...
                    ['LFR RctDataList contains no RCT data corresponding to argument cti1=%g.', ...
                    ' This may indicate that a zVar CALIBRATION_TABLE_INDEX(:,1)', ...
                    ' value is wrong or that BICAS did not try to load the corresponding', ...
                    ' RCT in glob. attr. CALIBRATION_TABLE.'], cti1)
                
%                 lfrItfIvpt = @(omegaRps) (bicas.calib_utils.eval_tabulated_TF(...
%                     RctDataList{cti1+1}.ItfModifIvptTable{iLsf}{iBlts}, omegaRps));
                lfrItfIvpt = RctDataList{cti1+1}.ItfModifIvptCaCa{iLsf}{iBlts};
            end
        end



    end    % methods(Access=private)

    %###################################################################################################################

    
    
    methods(Static, Access=public)



        % Load all non-BIAS RCTs (all types) using assumptions on filenames.
        %
        % NOTES
        % =====
        % NOTE: Can be useful for manual experimentation with calibration.
        % NOTE: Necessary when processing L1-->L2 (inofficially) since L1 does
        %       not have CALIBRATION_TABLE+CALIBRATION_TABLE_INDEX.
        % NOTE: Will only load ONE of each RCT type (no potential RCT time
        %       dependence as per global attribute CALIBRATION_TABLE) and
        %       requires user to not use CALIBRATION_TABLE_INDEX.
        %
        % IMPLEMENTATION NOTE: BICAS only needs one non-BIAS RCT type at a time.
        % However, it is useful to be able to initialize bicas.calib so that it
        % can simultanteously calibrate all kinds of data for debugging
        % purposes. Therefore loads ALL non-BIAS RCT types.
        %
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
                RctDataList    = {bicas.calib.read_RCT_modify_log(rctTypeId{1}, filePath, L)};
                
                % NOTE: Placing all non-BIAS RCT data inside 1x1 cell arrays so
                % that they are stored analogously with when using ga.
                % CALIBRATION_TABLE.
                RctDataMap(rctTypeId{1}) = RctDataList;
            end
        end



        % Load non-BIAS RCT(s) of ONE type (rctTypeId) using CDF global
        % attribute CALIBRATION_TABLE and zVars CALIBRATION_TABLE_INDEX and BW.
        %
        % IMPLEMENTATION NOTE
        % ===================
        % May load MULTIPLE RCTs (of the same RCT type) but will only load those
        % RCTs which are actually needed, as indicated by
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
            
            % NOTE: Iterate over those entries in CALIBRATION_TABLE that should
            % be considered, NOT all indices. May therefore legitimately leave
            % some cells in cell array empty.
            for i = 1:numel(iCtArray)
                % NOTE: Cell array index is one greater than the stored value.
                j              = iCtArray(i) + 1;
                filePath       = fullfile(rctDir, ga_CALIBRATION_TABLE{j});
                RctDataList{j} = bicas.calib.read_RCT_modify_log(rctTypeId, filePath, L);
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
        %
        function RctTypesMap = init_RCT_Types_Map()
            RctTypesMap = containers.Map();
            
            % NOTE: Not TF.
            MODIFY_TDS_CWF_DATA_FUNC = @(S) (S);
            
            RctTypesMap('BIAS')     = entry(@bicas.RCT.read_BIAS_RCT,     @bicas.calib.modify_BIAS_RCT_data,     @bicas.calib.log_BIAS_RCT,      'PROCESSING.RCT_REGEXP.BIAS');
            RctTypesMap('LFR')      = entry(@bicas.RCT.read_LFR_RCT,      @bicas.calib.modify_LFR_RCT_data,      @bicas.calib.log_LFR_RCTs,      'PROCESSING.RCT_REGEXP.LFR');
            RctTypesMap('TDS-CWF')  = entry(@bicas.RCT.read_TDS_CWF_RCT,  MODIFY_TDS_CWF_DATA_FUNC,              @bicas.calib.log_TDS_CWF_RCTs,  'PROCESSING.RCT_REGEXP.TDS-LFM-CWF');
            RctTypesMap('TDS-RSWF') = entry(@bicas.RCT.read_TDS_RSWF_RCT, @bicas.calib.modify_TDS_RSWF_RCT_data, @bicas.calib.log_TDS_RSWF_RCTs, 'PROCESSING.RCT_REGEXP.TDS-LFM-RSWF');
            
            %###################################################################
            function Entry = entry(readRctFunc, modifyRctFunc, logRctFunc, filenameRegexpSettingKey)
                Entry = struct(...
                    'readRctFunc',              readRctFunc, ...      % Pointer to function that reads one RCT.
                    'modifyRctFunc',            modifyRctFunc, ...    % Pointer to function that modifies data for one RCT.
                    'logRctFunc',               logRctFunc, ...       % Pointer to function that logs data for one RCT.
                    'filenameRegexpSettingKey', filenameRegexpSettingKey);
            end
            %###################################################################
        end
        
        
        
        function RctData = modify_BIAS_RCT_data(RctData)
            FtfRctSet = RctData.FtfSet;
            
            % Change name of field (sic!).
            % (There are many fields which are just kept untouched by this
            % function.)
            RctData = rmfield(RctData, 'FtfSet');
            RctData.FtfRctSet = FtfRctSet;
            
            % ASSERTIONS
            nTime = EJ_library.assert.sizes(...
                FtfRctSet.DcSingleAvpiv,   [-1, 1], ...
                FtfRctSet.DcDiffAvpiv,     [-1, 1], ...
                FtfRctSet.AcLowGainAvpiv,  [-1, 1], ...
                FtfRctSet.AcHighGainAvpiv, [-1, 1]);
            
            % NOTE: Derive ITFs.
            ItfSet = [];
            for iTf = 1:nTime
                % INVERT: FTF --> ITF
                ItfSet.DcSingleAvpiv{  iTf} = FtfRctSet.DcSingleAvpiv{  iTf}.inverse();
                ItfSet.DcDiffAvpiv{    iTf} = FtfRctSet.DcDiffAvpiv{    iTf}.inverse();
                ItfSet.AcLowGainAvpiv{ iTf} = FtfRctSet.AcLowGainAvpiv{ iTf}.inverse();
                ItfSet.AcHighGainAvpiv{iTf} = FtfRctSet.AcHighGainAvpiv{iTf}.inverse();
            end
            
            RctData.ItfSet = ItfSet;
        end
        
        
            
        function RctData2 = modify_LFR_RCT_data(RctData1)
            
            FtfRctTpivCaCa = RctData1.FtfTpivTable;
            
            % Read LFR FTFs, derive ITFs and modify them.
            itfModifIvptCaCa = {};
            for iLsf = 1:numel(FtfRctTpivCaCa)
                
                itfModifIvptCaCa{end+1} = {};
                for iBlts = 1:numel(FtfRctTpivCaCa{iLsf})

                    % INVERT: tabulated FTF --> tabulated ITF
                    ItfIvpt = FtfRctTpivCaCa{iLsf}{iBlts}.inverse();
                    
                    % MODIFY tabulated ITF
                    ItfModifIvpt = bicas.calib_utils.extrapolate_tabulated_TF_to_zero_Hz(ItfIvpt);
                    
                    % MODIFY tabulated ITF --> Function TF
                    %
                    % NOTE: Can not blindly forbid extrapolation (beyond the
                    % extrapolation to 0 Hz already done above) by setting value
                    % outside table=NaN (which deliberately triggers error
                    % elsewhere). LFR's tabulated TFs do in theory cover
                    % frequencies up to the Nyquist frequency, but in practice,
                    % the sampling frequency varies slightly. This means that
                    % when the Nyquist frequency varies slightly and sometimes
                    % it exceeds the tabulated frequencies.
                    % Ex: solo_L1R_rpw-lfr-surv-cwf-e-cdag_20201102_V01.cd
                    %
                    % 2020-11-06: LFR tables (RCT):
                    % F0=24576 Hz: f={  12--12288} [Hz]
                    % F1= 4096 Hz: f={0.01-- 2048} [Hz]
                    % F2=  256 Hz: f={0.01--  128} [Hz]
                    % F3=   16 Hz: f={0.01--    8} [Hz]
                    VALUE_OUTSIDE_TABLE = 0;
                    %VALUE_OUTSIDE_TABLE = NaN;   % Does not work. See comments above.
                    itfModifIvpt = @(omegaRps) (bicas.calib_utils.eval_tabulated_TF(...
                        ItfModifIvpt, omegaRps, VALUE_OUTSIDE_TABLE));
                    clear ItfModifIvpt
                    
                    itfModifIvptCaCa{iLsf}{iBlts} = itfModifIvpt;
                end
            end
            
            RctData2 = [];
            % NOTE: RctData.FtfRctTpivCaCa is still kept (for debugging).
            RctData2.FtfRctTpivCaCa   = FtfRctTpivCaCa;    % Just copied.
            RctData2.ItfModifIvptCaCa = itfModifIvptCaCa;
        end
        
        
        
        function RctData2 = modify_TDS_RSWF_RCT_data(RctData1)
            RctData2 = [];
            
            % Modify tabulated TDS-RSWF TFs.
            for iBlts = 1:numel(RctData1.ItfIvptList)
                % NOTE: Overwriting.
                
                ItfRctIvpt = RctData1.ItfIvptList{iBlts};
                
                % Store tabulated ITF EXACTLY AS THEY ARE in the RCT (before
                % modification).
                % NOTE: Struct field does not need to be
                % pre-initialized/pre-allocated.
                RctData2.ItfRctIvptCa{iBlts} = ItfRctIvpt;
                
                % MODIFY __tabulated__ ITF
                % (Does NOT wrap function handle in function handle.)
                ItfModifIvpt = bicas.calib_utils.extrapolate_tabulated_TF_to_zero_Hz(ItfRctIvpt);
                
                % MODIFY tabulated ITF --> function ITF
                %
                % NOTE: Use zero outside of tabulated frequencies (beyond
                % already made extrapolation). TDS-RSWF data requires
                % extrapolation.
                itfModifIvpt = @(omegaRps) (bicas.calib_utils.eval_tabulated_TF(ItfModifIvpt, omegaRps, 0));
                
                
                RctData2.itfModifIvptCa{iBlts} = itfModifIvpt;
                    
            end
            
        end
        


        % Log some indicative value(s) for a BIAS RCT.
        %
        % NOTE: Does not log file path. Caller is assumed to do that.
        function log_BIAS_RCT(RctData, L)
            
            DC_FREQ_HZ       = [0];   % Single & diffs.
            AC_DIFF_FREQS_HZ = [0, 1000];
            LL               = bicas.calib.RCT_DATA_LL;
            
            %=====================
            % Iterate over EpochL
            %=====================
            for iEpochL = 1:numel(RctData.epochL)
                
                L.logf(LL, 'Below values are used for data beginning %s:', ...
                    EJ_library.cdf.TT2000_to_UTC_str(RctData.epochL(iEpochL)))
                
                % Log bias current calibration
                L.logf(LL, '    BIAS current offsets: %s [aampere]',         bicas.calib_utils.vector_string('% 10e', RctData.Current.offsetsAAmpere(iEpochL, :)))
                L.logf(LL, '    BIAS current gain   : %s [aampere/TM unit]', bicas.calib_utils.vector_string('% 10e', RctData.Current.gainsAapt(     iEpochL, :)))
                
                % Log transfer functions (frequency domain), selected frequencies.
                L.logf(LL, ...
                    '    Note: Not logging the exact RCT BIAS TFs (FTFs; RctData.FtfRctSet) since the inversion is trivial.')
                log_TF('    BIAS ITF DC single',          DC_FREQ_HZ,       RctData.ItfSet.DcSingleAvpiv)
                log_TF('    BIAS ITF DC diff',            DC_FREQ_HZ,       RctData.ItfSet.DcDiffAvpiv)
                log_TF('    BIAS ITF AC diff, low  gain', AC_DIFF_FREQS_HZ, RctData.ItfSet.AcLowGainAvpiv)
                log_TF('    BIAS ITF AC diff, high gain', AC_DIFF_FREQS_HZ, RctData.ItfSet.AcHighGainAvpiv)
            end
            
            %=====================
            % Iterate over EpochH
            %=====================
            % NOTE: Must work for multiple CDF records.
            dcDiffOffsetsAVolt = [...
                RctData.DcDiffOffsets.E12AVolt, ...
                RctData.DcDiffOffsets.E13AVolt, ...
                RctData.DcDiffOffsets.E23AVolt];
            EJ_library.assert.sizes(dcDiffOffsetsAVolt, [NaN, 3]);            
            for iEpochH = 1:numel(RctData.epochH)
                L.logf(LL, 'Below values are used for data beginning %s:', ...
                    EJ_library.cdf.TT2000_to_UTC_str(RctData.epochH(iEpochH)))
                
                L.logf(LL, '    BIAS DC single voltage offsets ( V1, V2, V3): %s [avolt]', ...
                    bicas.calib_utils.vector_string('%g', RctData.dcSingleOffsetsAVolt(iEpochH, :)))
                L.logf(LL, '    BIAS DC diff   voltage offsets (E12,E13,E23): %s [avolt]', ...
                    bicas.calib_utils.vector_string('%g', dcDiffOffsetsAVolt(iEpochH)))

            end
                
            %###################################################################
            % Nested utility function.
            % NOTE: Impicitly function of iEpochL, L, LL.
            function log_TF(name, freqArray, ItfList)
                bicas.calib_utils.log_TF_function_handle(...
                    LL, name, 'avolt/ivolt', freqArray, ...
                    @(omegaRps) (ItfList{iEpochL}.eval(omegaRps)), L);
            end
            %###################################################################
        end


        
        % Analogous to log_BIAS_RCT.
        function log_LFR_RCTs(RctData, L)
            % NOTE: Frequencies may go outside of tabulated data.
            %FREQ_HZ = [0, 1, 5];
            FREQ_HZ = [0, 1, 100];
            
            % CASE: This index corresponds to an actually loaded RCT (some are
            % intentionally empty).
            for iLsf = 1:4
                if iLsf ~= 4   nBltsMax = 5;
                else           nBltsMax = 3;
                end
                
                for iBlts = 1:nBltsMax
                    
                    itfNamePrefix = sprintf('LFR, F%i, BLTS/BIAS_%i', iLsf-1, iBlts);
                    
                    itfName = sprintf('%s FTF (as in RCT)', itfNamePrefix);
                    bicas.calib_utils.log_TF_tabulated(...
                        bicas.calib.RCT_DATA_LL, ...
                        itfName, ...
                        RctData.FtfRctTpivCaCa{iLsf}{iBlts}, ...
                        L);
                    
                    %TabulatedItfIvpt = RctData.ItfModifIvptTable{iLsf}{iBlts};
                    %ItfIvpt          = @(omegaRps) (bicas.calib_utils.eval_tabulated_TF(TabulatedItfIvpt, omegaRps));
                    itfIvpt          = RctData.ItfModifIvptCaCa{iLsf}{iBlts};
                    itfName = sprintf('%s ITF (modif., interp.)', itfNamePrefix);
                    
                    bicas.calib_utils.log_TF_function_handle(...
                        bicas.calib.RCT_DATA_LL, ...
                        itfName, ...
                        'ivolt/TM unit', FREQ_HZ, itfIvpt, L)
                    
                end
            end    % for
            
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_CWF_RCTs(RctData, L)
            
            L.logf(bicas.calib.RCT_DATA_LL, ...
                'TDS CWF calibration factors: %s [ivolt/TM]', ...
                bicas.calib_utils.vector_string('%g', RctData.factorsIvpt));
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_RSWF_RCTs(RctData, L)
            % TODO: Log tabulated TFs.
            
            FREQ_HZ = 0;
            
            for iBlts = 1:3
                itfNamePrefix = sprintf('TDS RSWF, BLTS/BIAS_%i, ITF', iBlts);
                
                bicas.calib_utils.log_TF_tabulated(...
                    bicas.calib.RCT_DATA_LL, ...
                    sprintf('%s (as in RCT)', itfNamePrefix), ...
                    RctData.ItfRctIvptCa{iBlts}, ...
                    L);
                
%                 bicas.calib_utils.log_TF_tabulated(...
%                     bicas.calib.RCT_DATA_LL, ...
%                     sprintf('%s (modified)', itfNamePrefix), ...
%                     RctData.ItfModifIvptList{iBlts}, ...
%                     L);
                
                bicas.calib_utils.log_TF_function_handle(...
                    bicas.calib.RCT_DATA_LL, ...
                    sprintf('%s (modif., interp.)', itfNamePrefix), ...
                    'ivolt/TM unit', FREQ_HZ, RctData.itfModifIvptCa{iBlts}, L)
            end
        end



        % Read any single RCT file, modify the content as required (in practice
        % extrapolate TFs) and log it. Effectively wraps the different
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
        function RctData = read_RCT_modify_log(rctTypeId, filePath, L)
            
            L.logf(bicas.calib.READING_RCT_PATH_LL, ...
                'Reading RCT (rctTypeId=%s): "%s"', rctTypeId, filePath)
            
            readRctFunc   = bicas.calib.RCT_TYPES_MAP(rctTypeId).readRctFunc;
            modifyRctFunc = bicas.calib.RCT_TYPES_MAP(rctTypeId).modifyRctFunc;
            logRctFunc    = bicas.calib.RCT_TYPES_MAP(rctTypeId).logRctFunc;
            
            RctData = readRctFunc(filePath);
            RctData = modifyRctFunc(RctData);
            logRctFunc(RctData, L);
        end
        
        

    end    % methods(Static, Access=private)
    
    

end
