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
  %       Ex: bicas.const.swdmd.SWD_METADATA;
  %           Map-->Strings
  %       Ex: bicas.const.EMIDP_2_INFO;
  %           Map-->Struct
  %       Ex: bicas.const.gamods.GA_MODS_DB;
  %           Custom objects.
  %       PROPOSAL: Only log those which are easy.
  % PROPOSAL: Move bicas.const.N_BLTS to solo.hwzv.const.
  %
  % PROPOSAL: Derive lists of datasets using
  %           bicas.classify_BICAS_L1_L1R_to_L2_DSI() or reverse.
  %
  % PROPOSAL: Lists for both official+unofficial DSIs.
  %   NOTE: See bicas.const.gamods.init_GA_MODS_DB().
  %   TODO-DEC: Notation?



  %###########
  %###########
  % CONSTANTS
  %###########
  %###########
  properties(Constant)
    % NOTE: MATLAB versions are specified as the value returned by
    % version('-release'). Note that is omits "R".

    % MATLAB version which shall be used when OFFICIALLY running BICAS, in
    % particular at LIRA which only suports one particular MATLAB version.
    % NOTE: Not to be confused with bicas.const.PERMITTED_MATLAB_VERSIONS_CA.
    OFFICIAL_MATLAB_VERSION      = '2024a';

    % Permissible string values when running BICAS in general, both officially
    % and unofficially (value returned by version('-release')).
    %
    % NOTE: This value is more permissable than
    %       bicas.const.OFFICIAL_MATLAB_VERSION since it is sometimes useful to
    %       run BICAS for other than the official version.
    % NOTE: BICAS originally required MATLAB R2016a.
    %       This was later changed to only require MATLAB R2019b.
    %       Source: https://gitlab.obspm.fr/ROC/RCS/BICAS/issues/2#note_10804
    %       Official MATLAB version again was later changed to MATLAB R2024a.
    % NOTE: Added MATLAB 2023b since it is currently (2025-07-31) the latest
    %       MATLAB version running on brain, spis, anna (IRFU servers). This
    %       should be abolished eventually when beforementioned IRFU servers
    %       support MATLAB 2024a.
    % NOTE: Additionally added MATLAB R2025b since (1) this MATLAB version
    %       supports AI (in the GUI) which is useful for development, and (2) it
    %       is installed by primary BICAS author Erik P G Johansson
    %       (2025-11-24).
    PERMITTED_MATLAB_VERSIONS_CA = {'2023b', bicas.const.OFFICIAL_MATLAB_VERSION, '2025b'};

    % Path to "source directory" (the directory under which all source code is
    % located) relative to BICAS's directory root.
    SOURCE_DIR_RPATH             = "src";

    % Path to "config directory" (the directory where the default config file is
    % located, if any) relative to BICAS's directory root.
    DEFAULT_CONFIG_DIR_RPATH     = "config";

    DEFAULT_CONFIG_FILENAME      = 'bicas.conf';

    % MATLAB stdout prefix to signal to bash wrapper that the log message
    % should be passed on to STDOUT (without the prefix).
    STDOUT_PREFIX_TBW            = 'STDOUT: ';

    % MATLAB stdout prefix to signal to bash wrapper that the log message
    % should be passed on to LOG FILE (without the prefix).
    LOG_FILE_PREFIX_TBW          = "LOG FILE: ";

    SWD_FILENAME                 = "descriptor.json";

    BRVF_FILENAME                = "bias_rct_validity.json";

    DEFAULT_NSO_TABLE_RPATH      = fullfile("data", "solo_ns_ops.xml");

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
    EMIDP_2_INFO = bicas.const.init_EMIDP_2_INFO();



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



    % Number of BLTSs.
    N_BLTS = 5;



    % Minimum number of non-FV OSR records per bin/DSR record.
    %
    % NOTE: Currently used for both L2 downsampled and L3 downsampled.
    % /2021-05-24
    N_MIN_OSR_SAMPLES_PER_BIN = 3;

    % Threshold used for determining what counts as a data gap for CWF data.
    MAX_SAMPLE_GAP_RATIO = 2.0



    %===============
    % Lists of DSIs
    %===============
    % NOTE: Only official datasets.
    L2_LFR_DSI_CA = {...
      'SOLO_L2_RPW-LFR-SBM1-CWF-E'; ...
      'SOLO_L2_RPW-LFR-SBM2-CWF-E'; ...
      'SOLO_L2_RPW-LFR-SURV-CWF-E'; ...
      'SOLO_L2_RPW-LFR-SURV-SWF-E'};
    L2_TDS_DSI_CA = {...
      'SOLO_L2_RPW-TDS-LFM-CWF-E'; ...
      'SOLO_L2_RPW-TDS-LFM-RSWF-E'};
    % Only L1/L1R->L2 datasets
    L2_LFR_TDS_DSI_CA = [...
      bicas.const.L2_LFR_DSI_CA; ...
      bicas.const.L2_TDS_DSI_CA];
    L2_CWF_DSI_CA = {
      'SOLO_L2_RPW-LFR-SBM1-CWF-E'; ...
      'SOLO_L2_RPW-LFR-SBM2-CWF-E'; ...
      'SOLO_L2_RPW-LFR-SURV-CWF-E'; ...
      'SOLO_L2_RPW-TDS-LFM-CWF-E'}
    L2_SWF_DSI_CA = {
      'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
      'SOLO_L2_RPW-TDS-LFM-RSWF-E'}
    %
    L3_DENSITY_DSI_CA = {...
      'SOLO_L3_RPW-BIA-DENSITY'; ...
      'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS'};
    L3_EFIELD_DSI_CA = {...
      'SOLO_L3_RPW-BIA-EFIELD'; ...
      'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS'};
    L3_SCPOT_DSI_CA = {...
      'SOLO_L3_RPW-BIA-SCPOT'; ...
      'SOLO_L3_RPW-BIA-SCPOT-10-SECONDS'};
    % All L3, excluding VHT (OSR+DSR).
    L3_DENSITY_EFIELD_SCPOT_DSI_CA = [ ...
      bicas.const.L3_DENSITY_DSI_CA; ...
      bicas.const.L3_EFIELD_DSI_CA; ...
      bicas.const.L3_SCPOT_DSI_CA];

    RCT_DSI = 'SOLO_CAL_RPW-BIAS';



    NAN_TF = @(omegaRps) (omegaRps * NaN);



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

      MAP('NoError')                           = init_struct(0, ...
        'No error');
      MAP('BadMatlabVersion')                  = init_struct(1, ...
        'Using the wrong MATLAB version.');
      MAP('UntranslatableErrorMsgId')          = init_struct(1, ...
        ['Error occurred, but code can not translate the error''s', ...
        ' MATLAB message identifier into any of BICAS''s internal', ...
        ' standard error codes.']);
      MAP('MatlabCodeErrorHandlingError')      = init_struct(1, ...
        'The MATLAB code''s own error handling failed.');
      MAP('CLISyntax')                         = init_struct(1, ...
        'Can not interpret command-line interface (CLI) arguments syntax.');

      MAP('OperationNotImplemented')           = init_struct(1, ...
        ['Execution has reached a portion of the code that has', ...
        ' not been implemented yet.']);
      MAP('Assertion')                         = init_struct(1, ...
        ['Detected an internal state that should never be possible', ...
        ' in a bug-free code that receives correct inputs.']);
      MAP('IllegalArgument')                   = init_struct(1, ...
        'An argument to an internal function had an illegal value.');

      % Generic file/path error
      MAP('PathNotFound')                      = init_struct(1, ...
        'A specified directory or file does not exist.');
      MAP('PathNotAvailable')                  = init_struct(1, ...
        ['A specified file can not be created since the path', ...
        ' matches a pre-existing file/directory.']);
      MAP('CanNotOpenFile')                    = init_struct(1, ...
        'Can not open a file for reading/writing.');

      MAP('SWMProcessing')                     = init_struct(1, ...
        'Error in s/w mode processing (processing datasets).');

      MAP('DatasetFormat')                     = init_struct(1, ...
        ['Error when interpreting (official CDF) input datasets,', ...
        ' including master CDF files.']);
      MAP('IllegalOutputCdfFormatVersion')     = init_struct(1, ...
        ['The generated CDF has a CDF format version that is not', ...
        ' permitted by settings']);

      % Configuration/settings file error
      MAP('IllegalCodeConfiguration')          = init_struct(1, ...
        ['Bad hard-coded configuration (or possibly configurable', ...
        ' setting but should not be), e.g. constants, S/W', ...
        ' descriptor. This should ideally indicate a pure code', ...
        ' bug, i.e. it is not triggered by certain user-controlled input.']);
      MAP('IllegalOverridingSettingValueType') = init_struct(1, ...
        ['Overring setting value with data type that is not', ...
        ' allowed for that setting']);
      MAP('CannotInterpretConfigFile')         = init_struct(1, ...
        ['Can not interpret the content of the configuration file.', ...
        ' This implies a problem with the syntax.']);
      MAP('ConfigurationBug')                  = init_struct(1, ...
        'Trying to configure BICAS in an illegal way.');

      % RCT error
      MAP('FailedToReadInterpretRCT')          = init_struct(1, ...
        ['Can not interpret the content of the calibration file', ...
        ' (RCT) file, e.g. because the RCT contains invalid', ...
        ' calibration values.']);
      MAP('CannotFindRegexMatchingRCT')        = init_struct(1, ...
        ['Can not find any matching calibration file to read. No', ...
        ' file matches regular expression.']);

      % NSO table error
      MAP('FailedToReadInterpretNsoTable')     = init_struct(1, ...
        ['Can not read or interpret the content of the non-standard', ...
        ' operations (NSO) table.']);

      % IMPLEMENTATION NOTE: Using a nested function merely to keep the
      % function call short.
      function ErrorTypeInfo = init_struct(errorCode, errorDescription)
        % PROPOSAL: Replace with class.

        ErrorTypeInfo = struct(...
          'errorCode',   errorCode, ...
          'description', errorDescription);
      end

    end



  end    % methods(Static, Access=private)



end
