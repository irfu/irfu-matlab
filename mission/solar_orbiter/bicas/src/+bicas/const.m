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
  % PROPOSAL: Split up in multiple files.
  %   NOTE: There is already bicas.proc.L1L2.const.
  %   --
  %   PRO: Too large file. ~545 rows  /2025-07-30
  %     PRO: init_EMIDP_2_INFO() is ~80 rows.
  %       PRO: Is expected to be quite constant.
  %   --
  %   PROBLEM: Functions used for setting constants may need to use constants
  %            themselves. ==> Need to avoid cyclic dependence.
  %            ==> Affects splitting.
  %     Ex: init_QRC_constants() does not truly (?) use bicas.const.* values
  %         but could conceivably do in the future.
  %     PROPOSAL: See them as separate "const" files, used only for setting
  %               the corresponding constants.
  %       PROPOSAL: Modules bicas.const.*.
  %         Ex: bicas.const.qrc:     init_QRC_constants()
  %
  %
  %
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

    % MATLAB version which shall be used when officially running BICAS (value
    % returned by version('-release')).
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
    % NOTE: Added MATLAB 2023b since it is currently (2024-05-28) the latest
    %       MATLAB version running on brain, spis, anna (IRFU servers). This
    %       should be abolished eventually when beforementioned IRFU servers
    %       support MATLAB 2024a.
    PERMITTED_MATLAB_VERSIONS_CA = {'2023b', bicas.const.OFFICIAL_MATLAB_VERSION};

    % Path to "config directory" (the directory where the default config file is
    % located, if any) relative to BICAS's directory root.
    DEFAULT_CONFIG_DIR_RPATH     = 'config';

    DEFAULT_CONFIG_FILENAME      = 'bicas.conf';

    % MATLAB stdout prefix to signal to bash wrapper that the log message
    % should be passed on to STDOUT (without the prefix).
    STDOUT_PREFIX_TBW            = 'STDOUT: ';

    % MATLAB stdout prefix to signal to bash wrapper that the log message
    % should be passed on to LOG FILE (without the prefix).
    LOG_FILE_PREFIX_TBW          = 'LOG FILE: ';

    SWD_FILENAME                 = 'descriptor.json';

    BRVF_FILENAME                = 'bias_rct_validity.json';

    DEFAULT_NSO_TABLE_RPATH      = fullfile('data', 'solo_ns_ops.xml');

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



    % Constants for common values
    % ---------------------------
    % NOTE: Useful for e.g. testing.
    LxQBM_NONE       = uint16(0);
    % Absolute min & max for ZV QUALITY_FLAG, according to the definition in
    % external metadata standards.
    QUALITY_FLAG_MIN = uint8(0);
    QUALITY_FLAG_MAX = uint8(4);

    % QUALITY_FLAG value for both
    % (1) global saturation's "full saturation" QRC, and
    % (2) channel saturation QRCs.
    QUALITY_FLAG_SATURATION = uint8(0);

    % QRC-related constants initialized with code.
    Q = bicas.const.init_QRC_constants();



    % Minimum number of non-FV OSR records per bin/DSR record.
    %
    % NOTE: Currently used for both L2 downsampled and L3 downsampled.
    % /2021-05-24
    N_MIN_OSR_SAMPLES_PER_BIN = 3;



    % Number of BLTSs.
    N_BLTS = 5;



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
    L3_EFIELD_SCPOT_DSI_CA = {...
      'SOLO_L3_RPW-BIA-EFIELD'; ...
      'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS'; ...
      'SOLO_L3_RPW-BIA-SCPOT'; ...
      'SOLO_L3_RPW-BIA-SCPOT-10-SECONDS'};
    % DES = Density+EField+Scpot (OSR+DSR; not VHT).
    L3_DENSITY_EFIELD_SCPOT_DSI_CA = [ ...
      bicas.const.L3_DENSITY_DSI_CA; ...
      bicas.const.L3_EFIELD_SCPOT_DSI_CA];

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



    % Function for the initializing QRCSM with all QRCSs for channel saturation
    % in L2 CDFs and L3 CDFs.
    %
    function [L2Qrcsm, L3Qrcsm] = init_L2_L3_CHANNEL_SATURATION_QRCSM()
      L2Qrcsm = bicas.proc.QrcSettingsMap();
      L3Qrcsm = bicas.proc.QrcSettingsMap();

      %====================
      % CHANNEL SATURATION
      %====================
      % NOTE: DC/AC diffs use the same QRCIDs and settings.
      % NOTE: Channel saturation QRCIDs technically represent saturation on the
      %       corresponding ZVs, not the specified ASR SSIDs. This is only
      %       important if adding automatic saturation detection for non-ASRs
      %       (which is highly unlikely).
      function add_L2_ch_sat_QRCS(qrcid, l2qbmBit)
        L2Qrcs = bicas.proc.QrcSettingL2(...
          QUALITY_FLAG      =bicas.const.QUALITY_FLAG_SATURATION, ...
          L2_QUALITY_BITMASK=l2qbmBit);
        L2Qrcsm.add(qrcid, L2Qrcs);
      end

      add_L2_ch_sat_QRCS("SATURATION_ZV_V1",   1)
      add_L2_ch_sat_QRCS("SATURATION_ZV_V2",   2)
      add_L2_ch_sat_QRCS("SATURATION_ZV_V3",   4)
      add_L2_ch_sat_QRCS("SATURATION_ZV_V12",  8)
      add_L2_ch_sat_QRCS("SATURATION_ZV_V13", 16)
      add_L2_ch_sat_QRCS("SATURATION_ZV_V23", 32)

      L3Qrcsm.add("SATURATION_ZV_V1",  bicas.proc.QrcSettingL3(vdcFvIndexAr=1))
      L3Qrcsm.add("SATURATION_ZV_V2",  bicas.proc.QrcSettingL3(vdcFvIndexAr=2))
      L3Qrcsm.add("SATURATION_ZV_V3",  bicas.proc.QrcSettingL3(vdcFvIndexAr=3))
      L3Qrcsm.add("SATURATION_ZV_V12", bicas.proc.QrcSettingL3(edcFvIndexAr=1))
      L3Qrcsm.add("SATURATION_ZV_V13", bicas.proc.QrcSettingL3(edcFvIndexAr=2))
      L3Qrcsm.add("SATURATION_ZV_V23", bicas.proc.QrcSettingL3(edcFvIndexAr=3))

      assert(isequal(L2Qrcsm.qrcidAr, L3Qrcsm.qrcidAr))
    end



    % Function for initializing QRCSMs containing all QRCSs for producing
    % (1) all (official) L2 CDFs, and
    % (3) all L3 density CDFs.
    %
    % NOTE: The quality bit definitions here must be consistent with the
    % definitions in the corresponding CDF skeletons. Definitions in general
    % must be consistent with documentation.
    %
    function Q = init_L2_L3_L3Density_QRCSM()

      [Q.L2_CHANNEL_SATURATION_QRCSM, ...
        Q.L3_CHANNEL_SATURATION_QRCSM] = ...
        bicas.const.init_L2_L3_CHANNEL_SATURATION_QRCSM();

      L2Qrcsm        = bicas.proc.QrcSettingsMap();
      L3Qrcsm        = bicas.proc.QrcSettingsMap();
      L3DensityQrcsm = bicas.proc.QrcSettingsMap();

      L2Qrcsm.add_QRCSM(Q.L2_CHANNEL_SATURATION_QRCSM);
      L3Qrcsm.add_QRCSM(Q.L3_CHANNEL_SATURATION_QRCSM);

      %=================
      % Local constants
      %=================
      % No circular dependence?!! bicas.proc.L1L2.const <-> bicas.const
      S = bicas.proc.L1L2.const.C.SSID_DICT;
      % Global saturation quality variable scheme:
      % NOTE: L2QBM_BIT_PARTIAL_SATURATION is used for two different QRCIDs.
      L2QBM_BIT_PARTIAL_SATURATION = uint16(1);
      L2QBM_BIT_FULL_SATURATION    = uint16(2);



      %====================
      % PARTIAL_SATURATION
      %====================
      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG      =uint8(1), ...
        L2_QUALITY_BITMASK=L2QBM_BIT_PARTIAL_SATURATION);
      L2Qrcsm.add("PARTIAL_SATURATION", Qrcs)

      %=================
      % FULL_SATURATION
      %=================
      % NOTE: Also set PARTIAL saturation bit when FULL
      % saturation. /YK 2020-10-02.
      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG      =bicas.const.QUALITY_FLAG_SATURATION, ...
        L2_QUALITY_BITMASK=L2QBM_BIT_FULL_SATURATION + L2QBM_BIT_PARTIAL_SATURATION);
      L2Qrcsm.add("FULL_SATURATION", Qrcs);

      %=================
      % THRUSTER_FIRING
      %=================
      % NOTE: There will be an L1 QUALITY_BITMASK bit for thruster firings in
      % the future according to
      % https://confluence-lesia.obspm.fr/display/ROC/RPW+Data+Quality+Verification
      % Therefore(?) not setting any bit in L2_QUALITY_BITMASK.
      % (YK 2020-11-03 did not ask for any to be set.)
      Qrcs = bicas.proc.QrcSettingL2(QUALITY_FLAG=uint8(1));
      L2Qrcsm.add("THRUSTER_FIRING", Qrcs);

      %=============
      % BIAS_HW_OFF
      %=============
      % In practice, when LFR ZV "BW" says that BIAS h/w is off.
      % NOTE: Delete all values, both voltage & current.
      Qrcs = bicas.proc.QrcSettingL2(...
        voltageFvSsidAr=S.values, ...
        currentFvIantAr=[1:3]');
      L2Qrcsm.add("BIAS_HW_OFF", Qrcs);

      %=======
      % SWEEP
      %=======
      Qrcs = bicas.proc.QrcSettingL2(...
        voltageFvSsidAr=S.values, ...   % Blank all SSIDs.
        currentFvIantAr=[1:3]');
      L2Qrcsm.add("SWEEP", Qrcs);

      %==============
      % ANTx_FAILING
      %==============
      QUALITY_FLAG_ANTx_FAILING = bicas.const.QUALITY_FLAG_MAX;   % TEMPORARY. Good value not decided.
      L2QBM_ANTx_FAILING        = bicas.const.LxQBM_NONE;         % TEMPORARY. Good value not decided.
      %QUALITY_FLAG_ANTx_FAILING =bicas.const.QUALITY_FLAG_MIN;    % TEST
      %L2QBM_ANTx_FAILING        = uint16(65535);   % TEST

      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG      =QUALITY_FLAG_ANTx_FAILING, ...
        L2_QUALITY_BITMASK=L2QBM_ANTx_FAILING, ...
        voltageFvSsidAr =S(["DC_V1" "DC_V12" "DC_V13" "AC_V12" "AC_V13"]'));
      L2Qrcsm.add("ANT1_FAILING", Qrcs);

      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG      =QUALITY_FLAG_ANTx_FAILING, ...
        L2_QUALITY_BITMASK=L2QBM_ANTx_FAILING, ...
        voltageFvSsidAr =S(["DC_V2" "DC_V12" "DC_V23" "AC_V12" "AC_V23"]'));
      L2Qrcsm.add("ANT2_FAILING", Qrcs);

      Qrcs = bicas.proc.QrcSettingL2(...
        QUALITY_FLAG     =QUALITY_FLAG_ANTx_FAILING, ...
        L2_QUALITY_BITMASK=L2QBM_ANTx_FAILING, ...
        voltageFvSsidAr =S(["DC_V3" "DC_V13" "DC_V23" "AC_V13" "AC_V23"]'));
      L2Qrcsm.add("ANT3_FAILING", Qrcs);



      %=============
      % BAD_DENSITY
      %=============
      % NOTE: L3_QUALITY_BITMASK currently only exists for L3 DENSITY datasets
      % (i.e. not for other L3 datasets). If L3_QUALITY_BITMASK is used for other
      % L3 datasets in the future, then the bits may have different meanings for
      % those datasets.
      Qrcs = bicas.proc.QrcSettingL3Density(...
        QUALITY_FLAG      =uint8(1), ...
        L3_QUALITY_BITMASK=uint16(1));
      L3DensityQrcsm.add("BAD_DENSITY", Qrcs);



      %===============
      % REMOVE_EFIELD
      %===============
      % Remove EFIELD output from solo.vdccal() (i.e. not by removing input to
      % solo.vdccal()).
      Qrcs = bicas.proc.QrcSettingL3(efieldFvIndexAr=[1 2 3]');
      L3Qrcsm.add("REMOVE_EFIELD", Qrcs);



      Q.L2_QRCSM         = L2Qrcsm;
      Q.L3_QRCSM         = L3Qrcsm;
      Q.L3_DENSITY_QRCSM = L3DensityQrcsm;
    end



    function Q = init_QRC_constants()
      Q = bicas.const.init_L2_L3_L3Density_QRCSM();

      Q.CHANNEL_SATURATION_QRCID_AR = string(Q.L2_CHANNEL_SATURATION_QRCSM.qrcidAr);
      Q.SATURATION_QRCID_AR = [...
        Q.CHANNEL_SATURATION_QRCID_AR; ...
        "PARTIAL_SATURATION"; ...
        "FULL_SATURATION"];

      % All legal QRCIDs, or all kinds of processing. This defines the set of
      % legal QRCIDs, including ones that can be used in the NSO table file.
      % IMPLEMENTATION NOTE: This is required for asserting QRCIDs in the NSO
      % table file.
      Q.ALL_QRCID_AR = [...
        Q.L2_QRCSM.qrcidAr; ...
        Q.L3_QRCSM.qrcidAr; ...
        Q.L3_DENSITY_QRCSM.qrcidAr];
    end



  end    % methods(Static, Access=private)



end
