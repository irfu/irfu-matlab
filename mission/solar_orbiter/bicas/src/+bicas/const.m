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
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-09, as a replacement for the FUNCTION
% error_safe_constant created 2016-06-02.
%
classdef const
    % PROPOSAL: Error category for bad input datasets (both science and HK).
    %   PRO: Has similar for RCTs.
    %
    % PROPOSAL: Move N_MIN_OSR_SAMPLES_PER_BIN to settings?
    %
    % PROPOSAL: Log all constants.
    %   CON: Not straightforward/easy to log all constants since they use
    %        "non-primitive" data structures.
    %       Ex: SWD_METADATA = bicas.const.init_SWD_metadata();
    %           Map-->Strings
    %       Ex: EMIDP_2_INFO = bicas.const.init_EMIDP_2_INFO;
    %           Map-->Struct
    %       Ex: GA_MODS_DB = bicas.const.init_GA_MODS_DB();
    %           Custom objects.
    %       PROPOSAL: Only log those which are easy.
    % PROPOSAL: Move bicas.const.N_BLTS to solo.hwzv.const.



    properties(Constant)

        % Permissible string values returned by "version('-release')" when using
        % the correct MATLAB version.
        %
        % NOTE: BICAS originally required MATLAB R2016a but no longer does.
        % NOTE: ROC only needs MATLAB R2019b. Source:
        %       https://gitlab.obspm.fr/ROC/RCS/BICAS/issues/2#note_10804
        % NOTE: Added MATLAB 2020a in order to make it possible to run test code
        %       on GitHub CI (earliest supported MATLAB version).
        %       https://github.com/matlab-actions/setup-matlab/
        PERMITTED_MATLAB_VERSIONS         = {'2019b', '2020a'};

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
        EMIDP_2_INFO = bicas.const.init_EMIDP_2_INFO();



        SWD_METADATA = bicas.const.init_SWD_metadata();



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
        SWM_CLI_OPTION_REGEX = '[A-Za-z][\w-]+';



        % The RCS ICD 00037 iss1rev2 draft 2019-07-11, section 3.1.2.3 only
        % permits these characters (and only lowercase!).
        % This regexp only describes the "option body", i.e. not the preceding
        % "--".
        SIP_CLI_OPTION_BODY_REGEX = '[a-z0-9_]+';



        % Field values = Legal RCS NSOIDs used in the RCS NSO XML file.
        % Field names can be used as constants for those strings inside BICAS.
        %
        % IMPLEMENTATION NOTE: Specified as struct so that the struct can
        % simultaneously be used to
        % (1) compile a complete list of legal NSOIDs in the NSO table file
        % (2) reference specific constants (fields) throughout BICAS without
        %     hardcoding the actual NSOIDs in multiple places.
        %
        % IMPLEMENTATION NOTE: One does not want to use the RCS NSOID string
        % constants directly inside the code, in case of typos.
        %
        NSOID = struct(...
            'PARTIAL_SATURATION', 'PARTIAL_SATURATION', ...
            'FULL_SATURATION',    'FULL_SATURATION', ...
            'THRUSTER_FIRING',    'THRUSTER_FIRING');

        % Define the bits (bitmasks) in L2_QUALITY_BITMASK and
        % L3_QUALITY_BITMASK. Intended for bit operations.
        %
        % NOTE: The definitions here must be consistent with the definitions in
        % the corresponding CDF skeletons.
        %
        % NOTE: L3_QUALITY_BITMASK bits might (maybe) assign different meanings
        % to the same bit in different L3 data sets in the future (DENSITY.
        % EFIELD, SCPOT).
        L2QBM_PARTIAL_SATURATION = uint16(1);
        L2QBM_FULL_SATURATION    = uint16(2);
        L3QBM_BAD_DENSITY        = uint16(1);



        % Highest allowed QUALITY_FLAG values for various conditions.
        QUALITY_FLAG_CAP_PARTIAL_SATURATION = uint8(1);
        QUALITY_FLAG_CAP_FULL_SATURATION    = uint8(0);
        QUALITY_FLAG_CAP_THRUSTER_FIRING    = uint8(1);


            
        % Minimum number of non-FV OSR records per bin/DSR record.
        %
        % NOTE: Currently used for both L2 downsampled and L3 downsampled.
        % /2021-05-24
        N_MIN_OSR_SAMPLES_PER_BIN = 3;
        
        
        
        % Absolute min & max for ZV QUALITY_FLAG, according to the definition in
        % external metadata standards.
        QUALITY_FLAG_MIN = uint8(0);
        QUALITY_FLAG_MAX = uint8(3);
        
        % Number of BLTSs.
        N_BLTS = 5;



        GA_MODS_DB = bicas.const.init_GA_MODS_DB();

    end    % properties(Constant)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        function MAP = init_EMIDP_2_INFO()
            % NOTE: The RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section
            % 3.4.3 specifies
            %   error code 0 : No error
            %   error code 1 : Every kind of error (!)

            MAP = containers.Map('KeyType', 'char', 'ValueType', 'any');

            MAP('NoError')                      = init_struct(0, ...
                'No error');
            MAP('BadMatlabVersion')             = init_struct(1, ...
                'Using the wrong MATLAB version.');
            MAP('UntranslatableErrorMsgId')     = init_struct(1, ...
                ['Error occurred, but code can not translate the error''s', ...
                ' MATLAB message identifier into any of BICAS''s internal', ...
                ' standard error codes.']);
            MAP('MatlabCodeErrorHandlingError') = init_struct(1, ...
                'The MATLAB code''s own error handling failed.');
            MAP('CLISyntax')                    = init_struct(1, ...
                'Can not interpret command-line interface (CLI) arguments syntax.');
            MAP('PathNotFound')                 = init_struct(1, ...
                'A specified directory or file does not exist.');
            MAP('PathNotAvailable')             = init_struct(1, ...
                ['A specified file can not be created since the path', ...
                ' matches a pre-existing file/directory.']);
            MAP('CanNotOpenFile')               = init_struct(1, ...
                'Can not open a file for reading/writing.');
            MAP('OperationNotImplemented')      = init_struct(1, ...
                ['Execution has reached a portion of the code that has', ...
                ' not been implemented yet.']);
            MAP('Assertion')                    = init_struct(1, ...
                ['Detected an internal state that should never be possible', ...
                ' in a bug-free code that receives correct inputs.']);
            MAP('IllegalArgument')              = init_struct(1, ...
                'An argument to an internal function had an illegal value.');
            MAP('SWMProcessing')             = init_struct(1, ...
                'Error in s/w mode processing (processing datasets).');
            MAP('DatasetFormat')                = init_struct(1, ...
                ['Error when interpreting (official CDF) input datasets,', ...
                ' including master CDF files.']);
            MAP('IllegalCodeConfiguration')     = init_struct(1, ...
                ['Bad hard-coded configuration (or possibly configurable', ...
                ' setting but should not be), e.g. constants, S/W', ...
                ' descriptor. This should ideally indicate a pure code', ...
                ' bug, i.e. it is not triggered by certain user-controlled input.']);
            MAP('CannotInterpretConfigFile')    = init_struct(1, ...
                ['Can not interpret the content of the configuration file.', ...
                ' This implies a problem with the syntax.']);
            MAP('ConfigurationBug')             = init_struct(1, ...
                'Trying to configure BICAS in an illegal way.');
            MAP('FailedToReadInterpretRCT')     = init_struct(1, ...
                ['Can not interpret the content of the calibration file', ...
                ' (RCT) file, e.g. because the RCT contains invalid', ...
                ' calibration values.']);
            MAP('FailedToReadInterpretNsOps')   = init_struct(1, ...
                ['Can not read or interpret the content of the non-standard', ...
                ' operations file.']);
            MAP('CannotFindRegexMatchingRCT')   = init_struct(1, ...
                ['Can not find any matching calibration file to read. No', ...
                ' file matches regular expression.']);

            % IMPLEMENTATION NOTE: Using a nested function merely to keep the
            % function call short.
            function ErrorTypeInfo = init_struct(errorCode, errorDescription)
                ErrorTypeInfo = struct(...
                    'errorCode',   errorCode, ...
                    'description', errorDescription);
            end

        end



        % Various S/W descriptor (SWD) release data for the entire software (not
        % specific outputs)
        % ----------------------------------------------------------------------
        %
        function MAP = init_SWD_metadata()
            MAP = containers.Map();

            IRF_LONG_NAME = 'Swedish Institute of Space Physics (IRF)';
            %
            MAP('SWD.identification.project')     = 'ROC';
            MAP('SWD.identification.name')        = 'BIAS Calibration Software (BICAS)';
            MAP('SWD.identification.identifier')  = 'BICAS';
            MAP('SWD.identification.description') = ...
                ['Calibration software meant to be run at LESIA/ROC to', ...
                ' (1) calibrate electric field L2 data from', ...
                    ' electric L1R LFR and TDS-LFM data, and', ...
                ' (2) calibrate bias currents from L1R data.', ...
                ' Also has undocumented support for processing', ...
                ' L1 (not L1R) to L2, and BIAS L2 to L3 data', ...
                ' (both disabled by default).'];

            % 2021-09-21: Latest RCS ICD version is 01/05.
            MAP('SWD.identification.icd_version') = '1.4';

            % ROC-GEN-SYS-NTT-00019-LES, "ROC Engineering Guidelines for External
            % Users":
            % """"""""
            % 2.2.3 RCS versioning
            % The RCS version must be a unique number sequence identifier
            % “X.Y.Z”, where “X” is an integer indicating the release (major
            % changes, not necessarily retro-compatible), “Y” is an integer
            % indicating the issue (minor changes, necessarily retro-compatible)
            % and “Z” is an integer indicating a revision (e.g., bug
            % correction).
            % """"""""
            %
            %  ROC-PRO-PIP-ICD-00037-LES, "RPW Calibration Software Interface
            %  Control Document", 01/04:
            % """"""""
            % "version" : Current version of the S/W. The RCS version shall be a
            % unique number sequence identifier “X.Y.Z”, where “X” is an integer
            % indicating the release (major changes, not necessarily
            % retro-compatible), “Y” is an integer indicating the issue (minor
            % changes, necessarily retro-compatible) and “Z” is an integer
            % indicating a revision (e.g., bug correction). The first stable
            % release of software (S/W) must have its major number “X” equals to
            % 1, its minor number “Y” equals to 0 and its revision number “Z”
            % equals to 0 (i.e., “1.0.0”). S/W preliminary versions (e.g.,
            % alpha, beta, etc.) must have their version number “X” equals to 0
            % and must not have a character as a prefix/suffix (“0.Y.Zb” for the
            % 0.Y.Z beta version for instance). In all cases, any change in the
            % S/W must lead to update the version number.
            % """"""""
            MAP('SWD.release.version')            = '7.0.0';
            MAP('SWD.release.date')               = '2023-02-17T11:10:00Z';
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
%             MAP('SWD.release.modification')       = ...
%                 ['Non-Standard Operations (NSO) table for thruster firings', ...
%                 ' & saturation up until 2021-01-26', ...
%                 '; Better output CDF standards compliance', ...
%                 '; Using new master CDFs']; % v5.0.0
%             MAP('SWD.release.modification')       = ...
%                 ['Non-Standard Operations (NSO) table for thruster firings', ...
%                 ' & saturation up until 2021-09-21', ...
%                 '; Setting zVar attributes', ...
%                 ' SCALEMIN & SCALEMAX using zVar min & max values', ...
%                 '; Using new L2 master CDFs', ...
%                 '; Less excessive log messages for processing L1/L1R-->L2 LFR SWF', ...
%                 '; Permits HK and science input datasets to not overlap at', ...
%                 ' all in order to salvage some LFR DC data.']; % v6.0.0
%             MAP('SWD.release.modification')       = [...
%                 'Non-Standard Operations (NSO) table for thruster firings', ...
%                 ' & saturation up until 2022-09-03', ...
%                 '; Bugfix: Use LFR''s R0/R1/R2 for splitting time into', ...
%                 ' time intervals', ...
%             ]; % v6.0.1
%             MAP('SWD.release.modification')       = [...
%                 'Non-Standard Operations (NSO) table for thruster firings', ...
%                 ' updated for up until 2022-12-17', ...
%                 '; Bugfix: Corrected algorithm/formula for calculating Ez_SRF (vdccal.m).', ...
%                 '; Recalculated and added E-field calibration data.', ...
%                 ' Now covers 2020-02-28 to 2022-12-03.'
%             ]; % v6.0.2
            MAP('SWD.release.modification')       = [...
                'Non-Standard Operations (NSO) table for thruster firings', ...
                ' updated for up until 2023-02-05', ...
                '; Uses new master CDFs.'
            ]; % v7.0.0
            MAP('SWD.release.source')             = 'https://github.com/irfu/irfu-matlab/commits/SOdevel';
            % Appropriate branch? "master" instead?
            %
            % Relative path to BICAS executable. See RCS ICD.
            MAP('SWD.environment.executable')     = 'roc/bicas';



            %======================
            % ASSERTIONS: SETTINGS
            %======================
            irf.assert.castring_regexp(MAP('SWD.release.version'), ...
                '[0-9]+\.[0-9]+\.[0-9]+')
            irf.assert.castring_regexp(MAP('SWD.release.date'), ...
                '20[1-3][0-9]-[01][0-9]-[0-3][0-9]T[0-2][0-9]:[0-5][0-9]:[0-6][0-9]Z')
            % Validate S/W release version
            % ----------------------------
            % RCS ICD 00037, iss1rev2, Section 5.3 S/W descriptor file
            % validation scheme implies this regex.
            % NOTE: It is hard to thoroughly follow the description, but the end
            % result should be under release-->version-->pattern (not to be
            % confused with release_dataset-->version--pattern).
            irf.assert.castring_regexp(MAP('SWD.release.version'), '(\d+\.)?(\d+\.)?(\d+)')

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
        function Gmdb = init_GA_MODS_DB()
            % PROPOSAL: Exclude VHT since not produced by BICAS proper.
            %   CON: Production uses BICAS infrastrucutre for writing datasets.
            %
            % PROPOSAL: DSI lists as public constants.
            %   PROPOSAL: solo.hwzv.const.
            %   PRO: Could be used by functions for classifying DSIs.
            %       Ex: bicas.classify_BICAS_L1_L1R_to_L2_DSI().
            %       CON: Not if want to be really general, e.g. accounting for
            %            ROC-SGSE/SOLO distinctions.
            % PROPOSAL: Setting L2 and L3 in separate (sub)functions.

            %====================================================
            % Lists of commonly used GROUPS of DSIs
            % --------------------------------------------------
            % NOTE: Groups are allowed to overlap.
            % NOTE: Only include OFFICIAL datasets.
            %====================================================
            L2_LFR_DSIs = {...
                'SOLO_L2_RPW-LFR-SBM1-CWF-E', ...
                'SOLO_L2_RPW-LFR-SBM2-CWF-E', ...
                'SOLO_L2_RPW-LFR-SURV-CWF-E', ...
                'SOLO_L2_RPW-LFR-SURV-SWF-E'};
            L2_TDS_DSIs = {...
                'SOLO_L2_RPW-TDS-LFM-CWF-E', ...
                'SOLO_L2_RPW-TDS-LFM-RSWF-E'};
            % NOTE: Only include OFFICIAL L2 datasets.
            L2_LFR_TDS_DSIs = [L2_LFR_DSIs, L2_TDS_DSIs];

            L3_DENSITY_DSIs = {...
                'SOLO_L3_RPW-BIA-DENSITY', ...
                'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS'};
            L3_EFIELD_SCPOT_DSIs = {...
                'SOLO_L3_RPW-BIA-EFIELD', ...
                'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS', ...
                'SOLO_L3_RPW-BIA-SCPOT', ...
                'SOLO_L3_RPW-BIA-SCPOT-10-SECONDS'};
            % DES = Density+EField+Scpot (not VHT).
            L3_DES_DSIs = [ ...
                L3_DENSITY_DSIs, ...
                L3_EFIELD_SCPOT_DSIs];

            % All L3.
            L3_DSIs = [L3_DES_DSIs, {'SOLO_L3_RPW-BIA-VHT'}];

            %======================================================
            % Initialize empty data structure for all MODS entries)
            %======================================================
            % NOTE: Includes UNOFFICIAL DATASETS. Not strictly required.
            % NOTE: Formal parent dataset(s) might be changed due to
            %       reorganizing SWM, which could change the technically
            %       correct value.
            DSI_CA = [...
                L2_LFR_TDS_DSIs, L3_DSIs ...
                %{'SOLO_L2_RPW-LFR-SURV-CWF-E-1-SECOND'}...
            ]';
            Gmdb = bicas.gamods.Database(DSI_CA);

            
            
            %##############################################################
            % ACTUAL MODS ENTRIES, ADDED FOR ONLY THE RELEVANT DSIs
            %##############################################################

            %===================================================================
            % L2: At most one entry per BICAS version
            % ---------------------------------------
            % NOTE: L2 dates should be taken from ROC's BICAS git repo commits
            %       since those represent deliveries to ROC.
            %===================================================================
            % L3 DENSITY+EFIELD+SCPOT (not VHT): At most one entry per delivery
            % -----------------------------------------------------------------
            % NOTE: L3 dates are effectively determined by when dataset
            %       deliveries to ROC were generated.
            % --
            % NOTE:
            % (1) BICAS is used for generating L3 at IRF (not ROC), and
            % (2) BICAS version numbers are only updated when delivering to ROC.
            % Therefore,
            % (1) the MODS BICAS version numbers do not exactly specify the
            %     BICAS version, only the previous official version, and
            % (2) the dates may conflict with the combinations of BICAS version
            %     and date for other MODS entries.
            %===================================================================
            % L3 VHT
            % ------
            % NOTE: L3 dates are effectively determined by deliveries to ROC.
            % NOTE: Including VHT, since VHT uses the same BICAS functions for
            %       writing datasets (including
            %       bicas.derive_output_dataset_GAs).
            %===================================================================

            % BICAS v1.0.0 : No MODS needed since there are no changes compared
            %                to an earlier version.



            Gmdb.add_version_entry(L2_LFR_TDS_DSIs, ...
                bicas.gamods.VersionEntry('2020-05-18', '2.0.1', ...
                    {'Bias currents bugfixed to be correct unit.'}))



            GmveTds = bicas.gamods.VersionEntry('2020-07-07', '3.0.0', ...
                {'Bias currents changed to nA (not ampere).', ...
                 'Ignoring frequencies above high-frequency cutoff at 0.7 times Nyquist frequency.'});
            GmveLfr = GmveTds.add_comments({'Hereafter copying LFR L1 zVar BW.'});
            Gmdb.add_version_entry(L2_LFR_DSIs, GmveLfr)
            Gmdb.add_version_entry(L2_TDS_DSIs, GmveTds)
            clear GmveLfr GmveTds
            


            GmveTds = bicas.gamods.VersionEntry('2020-09-01', '3.1.0', {...
                    'Crude sweep removal based on mux mode.', ...
                    'Preliminary setting of QUALITY_FLAG (max 2).'});
            GmveLfr = GmveTds.add_comments({'Bugfix to handle LFR L1 zVar BW=0.'});
            Gmdb.add_version_entry(L2_LFR_DSIs, GmveLfr)
            Gmdb.add_version_entry(L2_TDS_DSIs, GmveTds)
            clear GmveLfr GmveTds



            Gmdb.add_version_entry(L2_LFR_TDS_DSIs, ...
                bicas.gamods.VersionEntry('2020-09-15', '3.1.1', {...
                    ['Ignoring frequencies above high-frequency cutoff at 0.8', ...
                    ' (instead of 0.7) multiplied by Nyquist frequency.']}))



            Gmdb.add_version_entry(L2_LFR_TDS_DSIs, ...
                bicas.gamods.VersionEntry('2020-10-07', '4.0.0', ...
                {'Uses table to set zVars QUALITY_FLAG and L2_QUALITY_BITMASK.'}))



            GmveTds = bicas.gamods.VersionEntry('2020-12-07', '4.1.0', {...
                ['Set QUALITY_FLAG and L2_QUALITY_BITMASK based on', ...
                ' tabulated thruster firings.']...
            });
            GmveLfr = GmveTds.add_comments({...
                ['Bugfixed AC detrending that only removes mean and does', ...
                ' not add linear component (mostly SWF).'], ...
                ['Inverting AC using artificial constant gain for low', ...
                ' frequencies to not amplify noise.']...
            });
            Gmdb.add_version_entry(L2_LFR_DSIs, GmveLfr)
            Gmdb.add_version_entry(L2_TDS_DSIs, GmveTds)
            clear GmveLfr GmveTds



            % L3 delivery 1: ~2021-01-29
            % NOTE: No entries, but the date is needed for determining MODS
            % between delivery 1 and 2.



            % No new L2 MODS entries (if excluding NSOPS update).
            Gmdb.add_version_entry(L2_TDS_DSIs, bicas.gamods.VersionEntry(...
                '2021-02-02', '5.0.0', ...
                {['Cap QUALITY_FLAG<=1 for tabulated thruster firings up', ...
                ' until 2021-01-26.']}))



            % L3 delivery 2: ~2021-02-16
            % NOTE: Master CDFs updated according to feedback. ==> No MODS.
            % psp2ne.m updated ==> DENSITY
            Gmdb.add_version_entry(...
                {'SOLO_L3_RPW-BIA-DENSITY', ...
                 'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS'}, ...
                bicas.gamods.VersionEntry('2021-02-16', '5.0.0', ...
                    {'Updated algorithm for density.'}))



            % L3 delivery 3: ~2021-04-09
            % vdccal.m updated ==> EFIELD updated.
            Gmdb.add_version_entry(...
                {'SOLO_L3_RPW-BIA-EFIELD', ...
                 'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS'}, ...
                bicas.gamods.VersionEntry('2021-04-09', '5.0.0', ...
                    {'Updated antenna scaling of E_z.'}))



            % VHT delivery 1: 2021-04-27 (Generation_time)
            % NOTE: No entries, but the date is needed for determining MODS for
            % delivery 2 (i.e. determine modifications between delivery 1 and 2).
            % 2023-02-07: There has not been any second delivery, and therefore
            % no MODS.



            % NOTE: Not included since it does not affect any already existent
            % datasets: "Salvage LFR DC data when HK does not overlap with
            % science anywhere in dataset."
            Gmdb.add_version_entry(L2_LFR_TDS_DSIs, ...
                bicas.gamods.VersionEntry('2021-09-21', '6.0.0', ...
                    {'Set zVar attributes SCALEMIN & SCALEMAX using data min & max.', ...
                    ['Cap QUALITY_FLAG<=1 for tabulated thruster firings up', ...
                    ' until 2021-09-11.']}))



            GmveTds = bicas.gamods.VersionEntry('2022-09-15', '6.0.1', ...
                {['Cap QUALITY_FLAG<=1 for tabulated thruster firings up', ...
                ' until 2022-09-03.']});
            GmveLfr = GmveTds.add_comments(...
                {'Bugfix: Use LFR''s R0/R1/R2 for splitting into time intervals.'});
            Gmdb.add_version_entry(L2_LFR_DSIs, GmveLfr)
            Gmdb.add_version_entry(L2_TDS_DSIs, GmveTds)
            clear GmveLfr GmveTds



            % BICAS v6.0.2
            Gmdb.add_version_entry(L2_LFR_TDS_DSIs, ...
                bicas.gamods.VersionEntry('2022-12-17', '6.0.2', ...
                    {['Cap QUALITY_FLAG<=1 for tabulated thruster firings up', ...
                    ' until 2022-12-17.']}))



            % L3 delivery 4: ~2022-12-20
            Gmdb.add_version_entry(...
                {'SOLO_L3_RPW-BIA-EFIELD', ...
                 'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS'}, ...
                 bicas.gamods.VersionEntry('2022-12-20', '6.0.2', ...
                    {'Bugfix: Updated formula for E_z.', ...
                     'New E field calibration data.'}))



            % BICAS v7.0.0
            Gmdb.add_version_entry(DSI_CA, ...
                bicas.gamods.VersionEntry('2023-02-17', '7.0.0', ...
                {['Updated all CDF skeletons to correct values for', ...
                ' GAs APPLICABLE and Data_type and correct usage of', ...
                ' zVar attributes DELTA_PLUS_VAR and DELTA_MINUS_VAR.']}))



            % FUTURE YET-TO-BE-DETERMINED BICAS VERSION
%             Gmdb.add_version_entry(L3_DSIs, ...
%                 bicas.gamods.VersionEntry('202x-xx-xx', 'x.y.z', ...
%                 {['Updated L3 CDF skeletons to remove LFR from', ...
%                 ' GAs Dataset_ID, Descriptor, and SKELETON_PARENT.']}))
%             Gmdb.add_version_entry(L2_LFR_TDS_DSIs, ...
%                 bicas.gamods.VersionEntry('202x-xx-xx', 'x.y.z', ...
%                 {['Support demultiplexer latching relay setting changing over time.']}))
%             Gmdb.add_version_entry(L3_DENSITY_DSIs, ...
%                 bicas.gamods.VersionEntry('202x-xx-xx', 'x.y.z', ...
%                 {['Add zVar L3_QUALITY_BITS with density quality bit.']}))
        end    % init_GA_MODS_DB



    end    % methods(Static, Access=private)



end
