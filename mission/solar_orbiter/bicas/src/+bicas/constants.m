%
% Hard-coded constants which are needed
% -- for error handling
% -- early, before regular settings are initialized,
% and which  thus need to be initialized independent of settings and in a way
% which is unlikely to trigger errors.
%
% NOTE: This file contains the authoritative definitions of the meaning of error
% codes that should (maybe) be used in documentation.
%
%
% VARIABLE NAMING CONVENTIONS
% ===========================
% TBW   = To Bash Wrapper.
% EMIDP = (MATLAB) Error Message Identifier Part. One of the colon-separated
%         parts of the MException .identifier string field (error message ID).
%         NOTE: "Component" has a special meaning in the context of error
%         message IDs. Therefore uses the term "part" instead.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-07-09, as a replacement for the FUNCTION
% error_safe_constant created 2016-06-02.
%
classdef constants   % < handle
    % PROPOSAL: Error category for bad input datasets (both science and HK).
    %   PRO: Has similar for RCTs.
    %
    % PROPOSAL: Move N_MIN_SAMPLES_PER_DWNS_BIN to settings?
    %
    % PROPOSAL: Log all constants.
    
    
    
    properties(Constant)

        % Permissible string values returned by "version('-release')" when using
        % the correct MATLAB version.
        %
        % NOTE: BICAS originally required MATLAB R2016a but no longer does.
        % NOTE: ROC only needs MATLAB R2019b. Source:
        % https://gitlab.obspm.fr/ROC/RCS/BICAS/issues/2#note_10804
        PERMITTED_MATLAB_VERSIONS         = {'2019b'};

        % Path to default config file relative to BICAS's directory root. Note
        % that this is also implicitly the constant for the default config file
        % filename.
        DEFAULT_CONFIG_FILE_RELATIVE_PATH = fullfile('config', 'bicas.conf');
        
        % MATLAB stdout prefix to signal to bash wrapper that the log message
        % should be passed on to STDOUT (without the prefix).
        STDOUT_PREFIX_TBW                 = 'STDOUT: ';
        
        % MATLAB stdout prefix to signal to bash wrapper that the log message
        % should be passed on to LOG FILE (without the prefix).
        LOG_FILE_PREFIX_TBW               = 'LOG FILE: ';

        % Information to "interpret" and "translate" captured exceptions
        % --------------------------------------------------------------
        % containers.Map with
        %   key   = Any one of the colon-separated parts of a MATLAB error
        %           message identifier string (see "error" function).
        %   value = Struct with fields representing a type of error:
        %       .errorCode
        %           The error code/number to be returned from BICAS' main
        %           function.
        %           IMPORTANT NOTE: A MATLAB error message identifier may match
        %           multiple "error types" (keys). The error-handling code
        %           (try-catch) should decide whether every message identifier
        %           should be used to identify only one error type if there are
        %           multiple ones to choose from.
        %       .description
        %           English human-readable text describing the error. Implicitly
        %           defines what kinds of errors this error code should cover.
        %
        %
        EMIDP_2_INFO = bicas.constants.init_EMIDP_2_INFO;
        
        
        
        SWD_METADATA = bicas.constants.init_swd_metadata();
        
        
        
        % Regular expression which the CLI name of a s/w mode must satisfy.
        %
        % The RCS ICD 00037, iss1rev2, draft 2019-07-11, section 5.3 seems to
        % imply this regex for S/W mode CLI parameters: ^[A-Za-z][\\w-]+$
        % NOTE: Only one backslash in MATLAB regex as opposed to in the RCS ICD.
        %
        % NOTE: Must not begin with "--" since it could be confused with CLI
        % options, but the RCS ICD constraints already ensure this.
        %
        % NOTE: help regexp: "\w    A word character [a-z_A-Z0-9]"
        %
        SW_MODE_CLI_OPTION_REGEX = '[A-Za-z][\w-]+';
        
        
        
        % The RCS ICD 00037 iss1rev2 draft 2019-07-11, section 3.1.2.3 only
        % permits these characters (and only lowercase!).
        % This regexp only describes the "option body", i.e. not the preceding
        % "--".
        SIP_CLI_OPTION_BODY_REGEX = '[a-z0-9_]+';
        
        
        
        % Field values = Legal RCS NSO IDs used in the RCS NSO XML file.
        % Field names can be used as constants for those strings inside BICAS.
        %
        % IMPLEMENTATION NOTE: Specified as struct so that the struct can
        % simultaneously be used to
        % ** compile a list of legal NSO IDs in NSO table file
        % ** reference specific constants (fields) throughout BICAS without
        %    hardcoding the actual NSO IDs in multiple places.
        %
        % IMPLEMENTATION NOTE: One does not want to use the RCS NSO ID string
        % constants directly inside the code, in case of typos.
        % 
        NSOID = struct(...
            'TEST_QF0',                'TEST_QUALITY_FLAG_0', ...
            'TEST_UFV',                'TEST_UFV', ...            
            'TEST_THRUSTER_FIRING',    'TEST_THRUSTER_FIRING', ...
            'PARTIAL_SATURATION',      'PARTIAL_SATURATION', ...
            'FULL_SATURATION',         'FULL_SATURATION', ...
            'THRUSTER_FIRING',         'THRUSTER_FIRING');
        
        % Define the bits in L2_QUALITY_BITMASK (L2QBM).
        % Intended for bit operations.
        L2QBM_PARTIAL_SATURATION = uint16(1);
        L2QBM_FULL_SATURATION    = uint16(2);
        
        
        
        % Minimum number of non-downsampled records per downsampled record.
        % NOTE: Does not consider whether samples are fill values.
        %
        % IMPLEMENTATION NOTE: Not perfect solution since includes fill values,
        % but easy to implement as at least a temporary algorithm.
        N_MIN_SAMPLES_PER_DWNS_BIN = 3;
        


        GA_MODS = bicas.constants.init_GA_MODS();

    end    % properties(Constant)
    
    

    methods(Static, Access=private)
        
        
        
        % Various S/W descriptor (SWD) release data for the entire software (not
        % specific outputs)
        % ----------------------------------------------------------------------
        %
        function MAP = init_swd_metadata()
            MAP = containers.Map();
            
            IRF_LONG_NAME = 'Swedish Institute of Space Physics (IRF)';
            %
            MAP('SWD.identification.project')     = 'ROC';
            MAP('SWD.identification.name')        = 'BIAS Calibration Software (BICAS)';
            MAP('SWD.identification.identifier')  = 'BICAS';
            MAP('SWD.identification.description') = ...
                ['Calibration software meant to', ...
                ' (1) calibrate electric field L2 data from electric L1R LFR and TDS (LFM) data, and', ...
                ' (2) calibrate bias currents from L1 data.'];
            
            % 2020-11-24: Latest document version is 01/04.
            MAP('SWD.identification.icd_version') = '1.4';
            
            % ROC-GEN-SYS-NTT-00019-LES, "ROC Engineering Guidelines for External
            % Users":
            % """"""""
            % 2.2.3 RCS versioning
            % The RCS version must be a unique number sequence identifier “X.Y.Z”,
            % where “X” is an integer indicating the release (major changes, not
            % necessarily retro-compatible), “Y” is an integer indicating the issue
            % (minor changes, necessarily retro-compatible) and “Z” is an integer
            % indicating a revision (e.g., bug correction).
            % """"""""
            %
            %  ROC-PRO-PIP-ICD-00037-LES, "RPW Calibration Software Interface
            %  Control Document", 01/04:
            % """"""""
            % "version" : Current version of the S/W. The RCS version shall be a unique number
            % sequence identifier “X.Y.Z”, where “X” is an integer indicating the release (major
            % changes, not necessarily retro-compatible), “Y” is an integer indicating the issue (minor
            % changes, necessarily retro-compatible) and “Z” is an integer indicating a revision (e.g.,
            % bug correction). The first stable release of software (S/W) must have its major number
            % “X” equals to 1, its minor number “Y” equals to 0 and its revision number “Z” equals
            % to 0 (i.e., “1.0.0”). S/W preliminary versions (e.g., alpha, beta, etc.) must have their
            % version number “X” equals to 0 and must not have a character as a prefix/suffix
            % (“0.Y.Zb” for the 0.Y.Z beta version for instance). In all cases, any change in the S/W
            % must lead to update the version number.
            % """"""""
            MAP('SWD.release.version')            = '5.0.0';
            MAP('SWD.release.date')               = '2021-02-02T15:15:00Z';
            MAP('SWD.release.author')             = 'Erik P G Johansson, BIAS team, IRF';
            MAP('SWD.release.contact')            = 'erjo@irfu.se';
            MAP('SWD.release.institute')          = IRF_LONG_NAME;   % Full name or abbreviation?
            % 'Various updates and refactoring; close to complete support for
            % LFR & TDS datasets (but untested); Removed ROC-SGSE_* dataset
            % support.' 'Almost-complete support for LFR & TDS datasets
            % (voltages) with transfer functions (partially tested).'
            % /Earlier version
            %
            % SWD.release.modification = ...
            % ['Modified default settings: inverted transfer function', ...
            % ' cutoff at 0.8*omega_Nyquist, duplicate bias current gives error'];
            % /v3.1.1
%             MAP('SWD.release.modification')       = ...
%                 ['Non-Standard Operations (NSO) table for setting QUALITY_FLAG, L2_QUALITY_BITMASK (new)', ...
%                 '; Set glob.attr. Datetime, OBS_ID, SOOP_TYPE, TIME_MIN, TIME_MAX', ...
%                 '; Modified default setting: PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY=ERROR'];   % v4.0.0
%             MAP('SWD.release.modification')       = ...
%                 ['Non-Standard Operations (NSO) table for thruster firings up until 2020-12-05', ...
%                 '; Zero order AC detrending, no AC re-trending', ...
%                 '; Combined transfer functions (LFR+BIAS) modified during', ...
%                 ' execution to have constant gain for low freqs. for AC data', ...
%                 '; Bugfix: glob.attr. Datetime, OBS_ID, SOOP_TYPE', ...
%                 '; Bugfix: not crashing when reading CURRENT datasets with one bias setting'];   % v4.1.0
            MAP('SWD.release.modification')       = ...
                ['Non-Standard Operations (NSO) table for thruster firings & saturation up until 2021-01-26', ...
                '; Better output CDF standards compliance', ...
                '; Using new master CDFs']; % v5.0.0
            MAP('SWD.release.source')             = 'https://github.com/irfu/irfu-matlab/commits/SOdevel';
            % Appropriate branch? "master" instead?
            %
            % Relative path to BICAS executable. See RCS ICD.
            MAP('SWD.environment.executable')     = 'roc/bicas';
            
            

            %======================
            % ASSERTIONS: SETTINGS
            %======================
            EJ_library.assert.castring_regexp(MAP('SWD.release.version'), '[0-9]+\.[0-9]+\.[0-9]+')
            EJ_library.assert.castring_regexp(MAP('SWD.release.date'),    '20[1-3][0-9]-[01][0-9]-[0-3][0-9]T[0-2][0-9]:[0-5][0-9]:[0-6][0-9]Z')
            % Validate S/W release version
            % ----------------------------
            % RCS ICD 00037, iss1rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
            % NOTE: It is hard to thoroughly follow the description, but the end result should be under
            % release-->version-->pattern (not to be confused with release_dataset-->version--pattern).
            EJ_library.assert.castring_regexp(MAP('SWD.release.version'), '(\d+\.)?(\d+\.)?(\d+)')

        end
        
        

        function MAP = init_EMIDP_2_INFO()
            % NOTE: The RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section
            % 3.4.3 specifies
            %   error code 0 : No error
            %   error code 1 : Every kind of error (!)
            
            MAP = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            MAP('NoError')                      = init_struct(0, 'No error');
            MAP('BadMatlabVersion')             = init_struct(1, 'Using the wrong MATLAB version.');
            MAP('UntranslatableErrorMsgId')     = init_struct(1, 'Error occurred, but code can not translate the error''s MATLAB message identifier into any of BICAS''s internal standard error codes.');
            MAP('MatlabCodeErrorHandlingError') = init_struct(1, 'The MATLAB code''s own error handling failed.');
            MAP('CLISyntax')                    = init_struct(1, 'Can not interpret command-line interface (CLI) arguments syntax.');
            MAP('PathNotFound')                 = init_struct(1, 'A specified directory or file does not exist.');
            MAP('OperationNotImplemented')      = init_struct(1, 'Execution has reached a portion of the code that has not been implemented yet.');
            MAP('Assertion')                    = init_struct(1, 'Detected an internal state that should never be possible in a bug-free code that receives correct inputs.');
            MAP('IllegalArgument')              = init_struct(1, 'An argument to an internal function had an illegal value.');
            MAP('SWModeProcessing')             = init_struct(1, 'Error in s/w mode processing (processing datasets).');
            MAP('DatasetFormat')                = init_struct(1, 'Error when interpreting (official CDF) input datasets, including master CDF files.');
            MAP('IllegalCodeConfiguration')     = init_struct(1, 'Bad hard-coded configuration (or possibly configurable setting but should not be), e.g. constants, S/W descriptor. This should ideally indicate a pure code bug, i.e. it is not triggered by certain user-controlled input.');
            MAP('CannotInterpretConfigFile')    = init_struct(1, 'Can not interpret the content of the configuration file. This implies a problem with the syntax.');
            MAP('ConfigurationBug')             = init_struct(1, 'Trying to configure BICAS in an illegal way.');
            MAP('FailedToReadInterpretRCT')     = init_struct(1, 'Can not interpret the content of the calibration file (RCT) file, e.g. because the RCT contains invalid calibration values.');
            MAP('FailedToReadInterpretNsOps')   = init_struct(1, 'Can not read or interpret the content of the non-standard operations file.');
            MAP('CannotFindRegexMatchingRCT')   = init_struct(1, 'Can not find any matching calibration file to read. No file matches regular expression.');
            
            % IMPLEMENTATION NOTE: Using a nested function merely to keep the
            % function call short.
            function ErrorTypeInfo = init_struct(errorCode, errorDescription)
                ErrorTypeInfo = struct(...
                    'errorCode',   errorCode, ...
                    'description', errorDescription);
            end
            
        end
        
        
        
        % Initialize data structure containing the CDF global attribute (GA)
        % MODS.
        %
        % RCS ICD 01/05 draft 2021-05-04, Table 4, on MODS:
        % =================================================
        % """"It shall be at least one entry for each release of
        % the RCS software that brings significant change in
        % the data content. Entry format shall be
        % “YYYY-MM-DD :: change #1 short description | change
        % #2 short description”,
        % where YYYY, MM and DD are respectively the year,
        % month and day of the new RCS release. Then followed
        % by the change descriptions, which shall be separated
        % by the pipe character (“|”)""""
        %
        function MAP = init_GA_MODS()
            % PROPOSAL: Exclude VHT since not produced by BICAS.
            % PROPOSAL: Use function to create MODS entries (strings).
            %   CON: Makes strings less readable.
            %   CON: Makes overlap between similar entries harder.
            %   PRO: Enforces correct format.
            %   PROPOSAL: %create_entry = @(dateStr, varargin) ([dateStr, ' :: ', join(varargin, ' | ')]);
            %
            % PROPOSAL: Utility function to add entry for selected DATASET_IDs.
            %       add_entry(entryStr, datasetIdsCa).
            %       Verify that specified DATASET_IDs are valid keys.
            %   PRO: Less use of temporary string variables. ==> Less risk of
            %        mistakes.
            %   CON: DATASET_IDs are long for being repeated hardcoded constants.
            %   PROPOSAL: Define lists of datasets if adding to recurring sets of
            %             DATASET_IDs.
            
            MAP = containers.Map('KeyType', 'char', 'ValueType', 'Any');
            
            %===================================================================
            % Empty lists of MODS entries for different (groups of) DATASET_IDs
            %===================================================================
            L2_LFR_SURV_SBM12_CWF = {};
            L2_LFR_SWF           = {};
            L2_TDS_CWF           = {};
            L2_TDS_RSWF          = {};
            %
            L3_DENSITY     = {};
            L3_DENSITY_10S = {};
            L3_EFIELD      = {};
            L3_EFIELD_10S  = {};
            L3_SCPOT       = {};
            L3_SCPOT_10S   = {};
            %
            L3_VHT = {};
            
            
            
            %##############################################################
            % ACTUAL MODS ENTRIES, ADDED FOR ONLY THE RELEVANT DATASET_IDs
            %##############################################################
            
            %===================================================================
            % L2
            % --
            % NOTE: L2 dates should be taken from ROC's BICAS git repo commits.
            %===================================================================
            % BICAS v1.0.0 : No MODS needed.
            
            % BICAS v2.0.1
            s = '2020-05-18 :: Bias currents bugfixed to be correct unit.';
            L2_LFR_SURV_SBM12_CWF{end+1} = s;
            L2_LFR_SWF{end+1}  = s;
            L2_TDS_CWF{end+1}  = s;
            L2_TDS_RSWF{end+1} = s;
            
            % BICAS v3.0.0
            sTds = ...
                ['2020-07-07 :: Bias currents changed to nA (not ampere).', ...
                ' | Ignoring frequencies above high-frequency cutoff at 0.7 times Nyquist frequency.'];
            sLfr = [sTds, ' | Hereafter copying LFR L1 zVar BW.'];            
            L2_LFR_SURV_SBM12_CWF{end+1} = sLfr;
            L2_LFR_SWF{end+1}  = sLfr;
            L2_TDS_CWF{end+1}  = sTds;
            L2_TDS_RSWF{end+1} = sTds;
            
            % BICAS v3.1.0
            s = ['2020-09-01 :: Bugfix to handle LFR L1 zVar BW=0.', ...
                ' | Crude sweep removal.', ...
                ' | Preliminary setting of QUALITY_FLAG (max 2).'];
            L2_LFR_SURV_SBM12_CWF{end+1} = s;
            L2_LFR_SWF{end+1}  = s;
            L2_TDS_CWF{end+1}  = s;
            L2_TDS_RSWF{end+1} = s;
          
            % BICAS v3.1.1
            s = ['2020-09-15 :: ', ...
                'Ignoring frequencies above high-frequency cutoff at 0.8', ...
                ' (instead of 0.7) times Nyquist frequency.'];
            L2_LFR_SURV_SBM12_CWF{end+1} = s;
            L2_LFR_SWF{end+1}  = s;
            L2_TDS_CWF{end+1}  = s;
            L2_TDS_RSWF{end+1} = s;
            
            % BICAS v4.0.0
            s = ['2020-10-07 :: ', ...
                'Uses table to set zVars QUALITY_FLAG and L2_QUALITY_BITMASK.'];
            L2_LFR_SURV_SBM12_CWF{end+1} = s;
            L2_LFR_SWF{end+1}  = s;
            L2_TDS_CWF{end+1}  = s;
            L2_TDS_RSWF{end+1} = s;

            % BICAS v4.1.0
            sTds = ['2020-12-07 :: ', ...
                'Set QUALITY_FLAG and L2_QUALITY_BITMASK based on tabulated thruster firings.'];
            sLfr = [sTds, ...
                ' | Bugfixed AC detrending that only removes mean and does not add linear component (mostly SWF).', ...
                ' | Inverting AC using artificial constant gain for low frequencies to not amplify noise.'];
            L2_LFR_SURV_SBM12_CWF{end+1} = sLfr;
            L2_LFR_SWF{end+1}  = sLfr;
            L2_TDS_CWF{end+1}  = sTds;
            L2_TDS_RSWF{end+1} = sTds;
            
            % BICAS v5.0.0 (already delivered): No new L2 MODS entries (if excluding NSOPS update).
            
            
            
            %===================================================================
            % L3 DENSITY+EFIELD+SCPOT (not VHT)
            % ---------------------------------
            % NOTE: L3 dates are effectively determined by dataset deliveries to
            % ROC were generated.
            %===================================================================
            % NOTE: No MODS entries for initial delivery.
            % NOTE: Including VHT, since VHT uses the same BICAS functions for
            % writing datasets (including
            % bicas.derive_output_dataset_GlobalAttributes).
            
            % Delivery 1: ~2021-01-29
            % NOTE: No entries, but the date is needed for determining MODS
            % between delivery 1 and 2.
            
            % Delivery 2: ~2021-02-16
            % NOTE: Master CDFs updated according to feedback. ==> No MODS.
            % psp2ne.m updated ==> DENSITY
            s = '2021-02-16 :: Updated algorithm for density.';
            L3_DENSITY{end+1}     = s;
            L3_DENSITY_10S{end+1} = s;
            
            % Delivery 3: ~2021-04-09
            % vdccal.m updated ==> EFIELD updated.
            s = '2021-04-09 :: Updated antenna scaling of E_z.';
            L3_EFIELD{end+1}      = s;
            L3_EFIELD_10S{end+1}  = s;
            
            
            
            %=================================================================
            % L3 VHT
            % ------
            % NOTE: L3 dates are effectively determined by deliveries to ROC.
            %=================================================================
            
            % Delivery 1: 2021-04-27 (Generation_time)
            % NOTE: No entries, but the date is needed for determining MODS
            % between delivery 1 and 2.
            
            % 2020-05-05: There has not been any second delivery, and therefore
            % no MODS.
            
            % L3_VHT{end+1} = 
            
            
            
            %=====================================================
            % Assign MODS entries to final data structure
            % -------------------------------------------
            % NOTE: Same list for all three SURV-CWF DATASET_IDs.
            %=====================================================
            MAP('SOLO_L2_RPW-LFR-SBM1-CWF-E') = L2_LFR_SURV_SBM12_CWF;
            MAP('SOLO_L2_RPW-LFR-SBM2-CWF-E') = L2_LFR_SURV_SBM12_CWF;
            MAP('SOLO_L2_RPW-LFR-SURV-CWF-E') = L2_LFR_SURV_SBM12_CWF;
            MAP('SOLO_L2_RPW-LFR-SURV-SWF-E') = L2_LFR_SWF;
            MAP('SOLO_L2_RPW-TDS-LFM-CWF-E')  = L2_TDS_CWF;
            MAP('SOLO_L2_RPW-TDS-LFM-RSWF-E') = L2_TDS_RSWF;
            
            % NOTE: INOFFICIAL DATASET. Sensible not strictly required.
            % NOTE: Formal parent dataset(s) might be changed due to
            % reorganizing s/w mode, which could change the technically correct
            % value.
            MAP('SOLO_L2_RPW-LFR-SURV-CWF-E-1-SECOND') = L2_LFR_SURV_SBM12_CWF;
            
            MAP('SOLO_L3_RPW-BIA-DENSITY')            = L3_DENSITY;
            MAP('SOLO_L3_RPW-BIA-DENSITY-10-SECONDS') = L3_DENSITY_10S;
            MAP('SOLO_L3_RPW-BIA-EFIELD')             = L3_EFIELD;
            MAP('SOLO_L3_RPW-BIA-EFIELD-10-SECONDS')  = L3_EFIELD_10S;
            MAP('SOLO_L3_RPW-BIA-SCPOT')              = L3_SCPOT;
            MAP('SOLO_L3_RPW-BIA-SCPOT-10-SECONDS')   = L3_SCPOT_10S;
            MAP('SOLO_L3_RPW-BIA_VHT')                = L3_VHT;
            
            
            
            % ASSERTIONS
            for keyCa = MAP.keys
                MODS = MAP(keyCa{1});
                for i = 1:numel(MODS)
                    s = MODS{i};
                    
                    EJ_library.assert.castring_regexp(s, ...
                        '20[1-9][0-9]-[0-1][0-9]-[0-3][0-9] :: [-=_|.() a-zA-Z0-9]+')
                    
                    % No more than one whitespace per occurrence.
                    assert(~contains(s, '  '))
                    
                    % All pipes surrounded by whitespace.
                    iPipes1 = strfind(s, '. | ') + 2;
                    iPipes2 = strfind(s,  '|');
                    assert(isequal(iPipes1, iPipes2))
                end
                
                % ASSERT: All strings are unique.
                EJ_library.assert.castring_set(MODS)
                
            end
            
        end
        
        
        
    end    % methods(Static)

end
