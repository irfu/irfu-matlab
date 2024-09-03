%
% Create a BSO by creating a bicas.Settings object and
% (1) define the set of permitted/existing settings keys, and
% (2) set all settings keys to their initial default values.
% Note: Does NOT make the object write-only.
%
% NOTE: Slightly deceiving name, since it defines which keys are permitted.
%
%
% NAMING CONVENTIONS
% ==================
% CUR : CURRENT (type of data, dataset)
% ENABLE(D)/DISABLE(D) always at the end of a setting key.
%
%
% NOTES ON BSO KEY NAMING CONVENTION
% ==================================
% Some constants
%   (1) correspond exactly to fields in the (JSON) S/W descriptor, and
%   (2) are unlikely to be used for anything else.
% These are labeled with a prefix "SWD." and the remainder is in lowercase
% (because that is what they are in the S/W descriptor).
%
%
% RETURN VALUE
% ============
% Bso
%       bicas.Settings object with all BICAS settings pre-defined with default
%       values. Can not define more settings, but can modify settings.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2018-01-24
%
function Bso = create_default_BSO()
% PROPOSAL: PROCESSING.CALIBRATION.CURRENT.HK.DISABLE      : Whether to calibrate HK current or use HK TM.
%                                                            Not which data to use (HK or TC).
%           PROCESSING.CALIBRATION.CURRENT.SOURCE = TC, HK : Which data to use.
%
% PROPOSAL: Setting keys should use cased version of ZVs and glob.attrs..
%   Ex: Epoch, (GA) Test_id, (GA) Dataset_ID.
%
% PROBLEM: Setting values "ERROR", "WARNING" are identical to the ICD-specified log row prefixes.
%   ==> Problems with grepping log files.
%   NOTE: Want setting value convention to be consistent with other settings values.
%       Ex: CORRECT, SORT, FULL, SCALAR, ROUND
%   PROPOSAL: Keep as is and grep log files using surrounding characters.
%       Ex: " ERROR "
%   PROPOSAL: Use lower case "error", "warning"
%   PROPOSAL: Use shortenings: "ERR", "WARN", "E", "W"
%
% PROPOSAL: Separate function for validating/asserting settings.
%   NOTE: Must be done AFTER all settings have been set.
%   PROPOSAL: Do every time settings are set, i.e. for default values,
%       config file values, CLI argument values.
%
% PROPOSAL: Abolish settings/functionality:
%   OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY
%   OUTPUT_CDF.EMPTY_NONNUMERIC_ZV_POLICY
%   OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE                     -- ALREADY ABOLISHED
%   OUTPUT_CDF.WRITE_FILE_DISABLED                          -- ALREADY ABOLISHED
%   PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY
%   PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY
%   --
%   PRO: Functionality appears to be obsolete.
%   PRO: Default ERROR has been used for a long time without raising exception.
%
%
%
% =========================
% BOGIQ: SETTING KEY NAMING
% =========================
% PROPOSAL: INPUT_CDF.* : Settings that apply to ALL input datasets.
% PROPOSAL: Only INPUT_CDF.ALL.* apply to all input datasets.
%   PROPOSAL: OUTPUT_CDF.ALL
%
% PROPOSAL: Use naming convention for settings keys for testing ONLY:
%
% PROPOSAL: Relevant setting keys should always be on the form ENABLE, never DISABLE.
%   PRO: More consistent.
%   CON: Less clear what is a deviation from the default.
% PROPOSAL: Relevant setting keys should always be on format "USE" (not "ENABLE", "DISABLE").
%   CON: Sounds bad to put "USE" at the end of settings key.
% PROPOSAL: Always use either DISABLED/ENABLED or DISABLE/ENABLE.
%   TODO-DEC: Which?
%
% PROPOSAL: Need (settings name) terminology for temporary
%           "bugfixes"/corrections due to bugs outside of BICAS,
%   Ex: Bugs datasets (bugs in other RCS, ROC's pipeline).
%   Ex: ~Corrupted data (different from bugs?)
%   Ex: PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY
%   NOTE: Need something is compatible with different course of action, not just permit or error.
%       Ex: PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY = ERROR, ROUND, PERMIT
%   PROPOSAL: Always include name of zVar.
%   PROPOSAL: "correction"
%   PROPOSAL: "mitigation"
%   PROPOSAL: "bugfix"
%       CON: Sounds like a bug in BICAS which it is not.
%   PROPOSAL: workaround
%   PROPOSAL: behaviour
%   PROPOSAL: action
%   PROPOSAL: ~anomaly
%
% PROPOSAL: Policy on usage of dataset levels in settings keys.
%   Ex: For now, L1R refers to algorithms to use WHEN processing L1R as input
%       Ex: PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY
%       NOTE: Does not have to refer to the L1R datasets as such.
%           Ex: L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS
%           Ex: L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2
%   Ex: For now, L2 refers to algorithms to use when processing L2 as output.
%       Ex: PROCESSING.L2.REMOVE_DATA.MUX_MODES
%   --
%   NEED: Specify whether refers to input or output data (not necessarily datasets).
%       Ex: Distinguish processing L1/L1R-->L2, L2-->L3.
%       PROPOSAL: Some kind of prefix before data level.
%           PROPOSAL: IN_L1R, OUT_L2 etc.
%           PROPOSAL: INPUT_CDF.L1R, OUTPUT_CDF.L2 etc.
%   NEED: Distinguish processing (science) data L1-->L2, L1R-->L2.
%   NEED: Distinguish input L1 voltage and L1 current.
%   NEED: Distinguish
%       (1) datasets of specific level, and
%       (2) algorithms that run when using given input level.
%   PROPOSAL: When speaking of datasets, use
%       (1) INPUT_CDF.<level>  : How to interpret, read datasets
%       (2) OUTPUT_CDF.<level> : How to output, write datasets.
%       PROBLEM: How distinguish from processing?
%
% TODO-DEC: How specify units for seconds? "_S" or "_SEC"?
%   NOTE: Is currently inconsistent.


S = bicas.Settings();

% The MATLAB command (e.g. path) to use to launch MATLAB for BICAS.
% NOTE: Only the value in the BICAS config file is actually used. The normal
% priority order for how BSO values are being obtained does apply here but
% does not matter since the value is only used by the bash wrapper script
% for launching MATLAB.
S.define_setting('MATLAB_COMMAND', '');



% Prefix used to identify the subset of stdout that should actually be
% passed on as stdout by the bash launcher script.
%S.define_setting('STDOUT_PREFIX',               'STDOUT: ');
% NOTE: Analogous LOG_PREFIX is hard-coded for safety.

% When logging contents of matrix/vector, maximum number of unique values
% printed before switching to shorter representation (min-max range)
S.define_setting('LOGGING.MAX_NUMERIC_UNIQUES_PRINTED', 5);
% When logging contents of TT2000 vector (in practice zVar Epoch), maximum
% number of unique TT2000 values printed before switching to shorter
% representation (min-max range).
S.define_setting('LOGGING.MAX_TT2000_UNIQUES_PRINTED', 2);

% Enable unofficial (to ROC) support for S/W modes
% ------------------------------------------------
% Enable s/w modes for processing LFR & TDS datasets L1-->L2 in addition to
% the official support for L1R. LFR_TDS refers to LFR/TDS input datasets, as
% opposed to L1 current datasets.
S.define_setting('SWM.L1-L2_ENABLED',          0);
% Enable s/w mode for processing L2 LFR-CWF-E to L2 LFR-CWF-E-1-SECONDS.
S.define_setting('SWM.L2-L2_CWF-DSR_ENABLED',  0);
% Enable s/w mode for processing L2-->L3 datasets.
S.define_setting('SWM.L2-L3_ENABLED',          0);



%####################
% ENV_VAR_OVERRIDE.*
%####################
% Variables which, if non-empty, are used to override the corresponding
% environment variables.
% ROC_RCS_CAL_PATH : Defined in RCS ICD. Path to dir. with calibration files.
S.define_setting('ENV_VAR_OVERRIDE.ROC_RCS_CAL_PATH',    '');
% ROC_RCS_MASTER_PATH : Defined in RCS ICD. Path to dir. with master CDF
%                       files.
S.define_setting('ENV_VAR_OVERRIDE.ROC_RCS_MASTER_PATH', '');



%###########################################
% INPUT_CDF.*
% ------------
% Settings that apply to ALL input datasets
%###########################################

% The epoch for zVar ACQUISITION_TIME.
% The time in UTC at which ACQUISITION_TIME is [0,0].
% Format: Year-month-day
%         -hour-minute-second
%         -millisecond-microsecond(0-999)-nanoseconds(0-999)
%
% PROPOSAL: Store the value returned by
%           spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) instead?
S.define_setting('INPUT_CDF.ACQUISITION_TIME_EPOCH_UTC', [2000,01,01, 12,00,00, 000,000,000]);

% NOTE: Requires INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY = non-error.
S.define_setting('INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED', 1)
% NOTE: See INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED
S.define_setting('INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY',  'WARNING')    % WARNING, ERROR

S.define_setting('INPUT_CDF.USING_GA_NAME_VARIANT_POLICY',  'WARNING')    % WARNING, ERROR

% Require input CDF Global Attribute "DSI" to match the expected
% value.
S.define_setting('INPUT_CDF.GA_DSI_MISMATCH_POLICY',        'WARNING')    % ERROR, WARNING
S.define_setting('INPUT_CDF.GA_PARENTS_MISMATCH_POLICY',    'WARNING')    % ERROR, WARNING

% NOTE: This modification applies BEFORE
% PROCESSING.HK.USE_ZV_ACQUISITION_TIME and therefore always applies to zVar
% Epoch.
% NOTE: Only check for increasing, not monotonically.
% NOTE: There is a known, mitigatable anomaly in SOLO_L1_RPW_BIA-CURRENT
% which duplicated settings (same timestamp, same bias setting on same
% antenna) which would be triggered by an assertion on an assert on
% monotonically increasing timestamps.
S.define_setting('INPUT_CDF.NON-INCREMENTING_ZV_EPOCH_POLICY', 'ERROR')      % ERROR, WARNING, SORT

S.define_setting('INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY', 'ERROR')    % ERROR, REMOVE_DUPLICATES



%############################################
% OUTPUT_CDF.*
% ------------
% Settings that apply to ALL output datasets
%############################################
% What BICAS should do when there is a pre-existing file on a output dataset
% file path.
% NOTE: Not known if the RCS ICD says anything about what should be the
% default, or what ROC thinks it should be.
S.define_setting('OUTPUT_CDF.PREEXISTING_OUTPUT_FILE_POLICY', 'WARNING');    % ERROR, WARNING.

% "There must be one entry for each entry in the CALIBRATION_TABLE
% attribute"  /RCS ICD 1.6
% NOTE: None of these disabled settings are actually used in the code.
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_ENTITY_NAME.BIAS',        'BIAS team')
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_ENTITY_NAME.LFR',         'LFR team')
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_ENTITY_NAME.TDS',         'TDS team')
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_ENTITY_AFFILIATION.BIAS', IRF_LONG_NAME)
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_ENTITY_AFFILIATION.LFR',  'Laboratoire de Physique des Plasmas (LPP)')      % Should be checked.
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_ENTITY_AFFILIATION.TDS',  'Institute of Atmospheric Physics AS CR (IAP)')   % Should be checked.
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_EQUIPMENT.BIAS', 'BIAS')   % Abolish?
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_EQUIPMENT.LFR',  'LFR')    % Abolish?
% S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.CAL_EQUIPMENT.TDS',  'TDS')    % Abolish?



% What to do with zVariables which are still empty after copying data into
% the master CDF. This indicates that something is wrong, either in the
% master CDF or in the processing.
S.define_setting('OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY',    'ERROR');   % ERROR, WARNING, USE_FILLVAL
% Ex: Non-numeric ACQUISITION_TIME_UNITS in (master?)
%     SOLO_L2_RPW-LFR-SBM1-CWF-E_V05.cdf is empty
% Ex: VDC_LABEL etc can be empty due to ROC bug updating skeletons.
S.define_setting('OUTPUT_CDF.EMPTY_NONNUMERIC_ZV_POLICY', 'ERROR');   % ERROR, WARNING



% NOTE: ACQUSITION_TIME_UNITS being empty in the master CDF requires value
% 0/false.
S.define_setting('OUTPUT_CDF.write_dataobj.strictEmptyZvClass',                1)

% Whether the size per record of an empty (0 records) output DF zVar has to
% be in agreement with the master CDF's size per record.
% NOTE: ACQUSITION_TIME_UNITS being empty in the master CDF requires value
% 0/false.
S.define_setting('OUTPUT_CDF.write_dataobj.strictEmptyNumericZvSizePerRecord', 1)

% Whether the size per record of an output CDF zVar has to be in agreement
% with the master CDF's size per record. Disabling is useful if the master
% CDF has not been updated in this regard only.
%
% 2020-10-22: Disabled due to L2 ZVs IBIAS1/2/3 mismatch between
% master CDFs and BICAS code. Thus only affects SWF/RSWF, not CWF.
%   Master CDFs: Multiple samples/snapshot.
%   BICAS code:  1 sample/snapshot.
% 2021-02-02: Skeletons fixed in L2 skeletons V12. Can now enable.
S.define_setting('OUTPUT_CDF.write_dataobj.strictNumericZvSizePerRecord',      1)   % 0/false, 1/true.

% Permitted CDF versions as a reg.expr.. Entire CDF version string must match
% the reg.expr..
% CDF format version 3.9 is required by ROC (Solo?). /2024-07-24
%
% NOTE: If irfu-matlab and BICAS (as delivered to ROC) use different CDF format
% versions, then BICAS should be called in irfu-matlab in such a way that this
% setting is overriden. (The default value should be what is required by
% ROC/SolO, as with all other settings.)
S.define_setting('OUTPUT_CDF.FORMAT_VERSION_REGEXP', '3\.9\.[0-9]+')



%##############
% PROCESSING.*
%##############
% Whether to use ACQUISITION_TIME instead of Epoch for HK.
% NOTE: This change happens AFTER
% INPUT_CDF.NON-INCREMENTING_ZV_EPOCH_POLICY.
% NOTE: Setting created so that HK can use ACQUISITION_TIME for
% interpolating its data to SCI time. Not trivial (but doable) to generalize
% to SCI data (voltages) since the naming implies using this for all data
% use, not just the HK-SCI interpolation. Such generalization should ideally
% be made when reading the dataset, but then code which treats datasets as
% generic, has to dentify which dataset is SCI. Should not be worth the
% effort.
S.define_setting('PROCESSING.HK.USE_ZV_ACQUISITION_TIME',    0)

% How to react to HK not overlapping with SCI.
% NOTE: Switch is shared for LFR & TDS, but WARNING only(?) makes sense for
% LFR, since some data can be salvaged in the event of non-overlap for LFR
% (using LFR BDM), but not for TDS.
S.define_setting('PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY',       'WARNING')    % WARNING, ERROR
% NOTE: "WARNING": Will lead to using nearest interpolation.
S.define_setting('PROCESSING.HK.TIME_NOT_SUPERSET_OF_SCI_POLICY',  'WARNING')    % WARNING, ERROR
S.define_setting('PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY', 'WARNING')    % WARNING, ERROR

% Quick ~BUGFIX for bad values in zv SAMPLING_RATE in L1R TDS-LFM-RSWF
% datasets. Abolish?
S.define_setting('PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY', 'ERROR')   % WARNING, ERROR, CORRECT

% Mitigation: How to handle that LFR ZVs QUALITY_FLAG QUALITY_BITMASK are
% empty.
% Needed for test files
%   ROC-SGSE_L1R_RPW-LFR-SBM1-CWF-E_4129f0b_CNE_V02.cdf
%   ROC-SGSE_L1R_RPW-LFR-SBM2-CWF-E_6b05822_CNE_V02.cdf
%   These are not used any more. /2020-10-07
S.define_setting('PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY', 'ERROR')   % ERROR, USE_FILL_VALUE

% ~BUGFIX for bug in L1/L1R TDS-LFM RSWF datasets.
% TDS has bugfixed. /2019-12-19
% PROPOSAL: Rename.
S.define_setting('PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY', 'ERROR')   % ERROR, WARNING, ROUND

%============================================================================
% Where to obtain the BDM (mux mode)
% ----------------------------------
%
% BIAS HK data
% ------------
% Contains BDM using its own Epoch (typically ~30 s time resolution?), which
% means that ~interpolation to SCI data is necessary, which means that the
% effective BDM value can briefly be wrong. NOTE: BIAS HK may potentially
% NOT cover the same time range as SCI data at all, and then the BDM can be
% really wrong (e.g. when using different versions).
%
% LFR SCI data (L1/L1R)
% ---------------------
% Contains a ZV for BDM using the same Epoch as the data.
% NOTE: Th LFR L1/L1R BDM may be available when the BIAS HK BDM is not.
% NOTE: The relevant TDS datasets do not contain BDM.
%============================================================================
S.define_setting('PROCESSING.LFR.MUX_MODE_SOURCE', 'LFR_SCI')    % BIAS_HK, LFR_SCI, BIAS_HK_LFR_SCI



%============================================================================
% Settings for when to remove data by setting it to fill value
% ------------------------------------------------------------
% "L2" refers to output datasets. Both voltage and current data. In
% practice, this functionality is there as a temporary solution for removing
% sweeps.
%============================================================================
%S.define_setting('PROCESSING.L2.REMOVE_DATA.MUX_MODES', [1,2,3,4,5,6,7])
S.define_setting('PROCESSING.L2.REMOVE_DATA.MUX_MODES', zeros(0, 1))

% Unit: S = Seconds
% Lower number since using LFR BDM (mux mode; unless configured not to),
% which has same cadence as science data.
% See PROCESSING.LFR.MUX_MODE_SOURCE.
S.define_setting('PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S',  0)

% Higher number since using BIAS HK for TDS, which means that the BDM is
% known with a lower time resolution.
S.define_setting('PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S', 30)

% Lowest zVar QUALITY_FLAG value that may be used for deriving L3 DENSITY,
% EFIELD, and SCPOT data. Affects both OSR and DSR.
S.define_setting('PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN',    2)

% Lowest zVar QUALITY_FLAG value that may be used for deriving downsampled
% SOLO_L2_RPW-LFR-SURV-CWF-E-1-SECOND data.
% NOTE: This does not affect the corresponding OSR dataset and is therefore
% not entirely analogous to PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN.
S.define_setting('PROCESSING.L2-CWF-DSR.ZV_QUALITY_FLAG_MIN',  2)

% Maximum value for zVar QUALITY_FLAG in output datasets.
% YK        2020-08-31: Use cap "2=Survey data, possibly not publication-quality"
% YK, Slack 2024-01-19: Use cap "3=Good for publication, subject to PI approval"
%
% NOTE: Used for both L2 and L3 datasets.
S.define_setting('PROCESSING.ZV_QUALITY_FLAG_MAX', 3)

% Path to NSO table file. Relative to BICAS root.
S.define_setting('PROCESSING.NSO_TABLE.FILE.RELATIVE_PATH', fullfile('data', 'solo_ns_ops.xml'))
% Path to NSO table file for debugging purposes.
% If non-empty, then it overrides PROCESSING.NSO_TABLE.FILE.RELATIVE_PATH.
% Can be set to absolute path. Intended for testing.
S.define_setting('PROCESSING.NSO_TABLE.FILE.OVERRIDE_PATH', '')

%-------------------------------------------------------------------------------
% Configuration parameters for the automatic detection of saturation performed
% for L1R-->L2 processing.
% ------------------------------------------------------------------------------
% NOTE: The measured VDC3 potential is stuck at 0.8 V due to apparent h/w
% failure on 2023-11-13T23:35. This is relevant since the returned value is very
% constant and falsely look like a regular saturated value, but at a much lower
% value. Saturation detection is not supposed to detect this value, but the
% existence of this constant VDC3, and its consequences for other reconstructed
% values in datasets such as correlated diffs, may mistakenly be interpreted as
% saturation if one is not careful.
%-------------------------------------------------------------------------------
% TODO-DEC: Too high AC diff thresholds?
S.define_setting('PROCESSING.SATURATION.CWF_SLIDING_WINDOW_LENGTH_SEC',            60.0);
S.define_setting('PROCESSING.SATURATION.TSF_FRACTION_THRESHOLD',                    0.5);
S.define_setting('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.SINGLE',         40.0);
S.define_setting('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.DC.DIFF',            2.0);
S.define_setting('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.LOW_GAIN',   0.3);
S.define_setting('PROCESSING.SATURATION.HIGHER_THRESHOLD_AVOLT.AC.DIFF.HIGH_GAIN',  0.3/20);

%--------------------------------------------------------------------------
% Settings for autodetecting sweeps (so that they can be excluded from L2)
%--------------------------------------------------------------------------
% NOTE: This is intended as a temporary functionality. The long-term solution
% should be to use L1/L1R QUALITY_BITMASK. It has been created so that sweeps
% can still be removed while BIAS is commanded to use BDM=4 ("mux=4") for bulk
% data.
%--------------------------------------------------------------------------
% PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC:
% Before this time: SBDA is used for detecting sweeps (BDM=4).
% After this time:  SCDA is used for detecting sweeps
%                   (moving window + varying HK currents).
% Format: Year-month-day
%         -hour-minute-second
%         -millisecond-microsecond(0-999)-nanoseconds(0-999)
%
% NOTE: BDM changed from 0 to 4 for bulk science data on/at:
% (1) Xavier Bonnin e-mail 2023-12-22: "BIAS is set to MUX_4 on-board after
%     Dec. 25."
% (2) According to BIAS HK: 2023-12-25T23:29:10 (+/-30 s).
% (3) According to SOLO_L1R_RPW-LFR-SURV-CWF-E: between about
%     2023-12-25T23:28:21 and 2023-12-25T23:28:44.
% However, a test with multiple BDMs (mux modes) ran on 2023-12-16 so it is
% still worth NOT setting PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC to after
% that.
%
% NOTE: Might be that SCDA window length=3 pts is a bit too short (or possibly
% diff minimum=500 TM is too large), but that it is saved by window margin 120
% s. In combination, the SCDA settings seem to be good enough though.
% /2024-02-01
%
% NOTE: Empirically, sweeps are surrounded by small data gaps, 1-4 min long(?).
%-------------------------------------------------------------------------------
S.define_setting('PROCESSING.L2.DETECT_SWEEPS.SBDA.END_UTC', [2023, 12, 16, 0, 0, 0, 0, 0, 0])
% SCDA window length. Unit: Data points/HK CDF records.
S.define_setting('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_LENGTH_PTS', 3)
% SCDA threshold for HK bias current difference between min and max within a
% window. If the value exceeds this value, then the interval is labelled as
% sweeping. Unit: TM units
% NOTE: Empirically, fluctuations around constant bias current (on a single
% channel) is on the order of 30-40 TM units (2024-01-01).
S.define_setting('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MINMAX_DIFF_MINIMUM_TM', 500)
% Amount of margin to add around regions labelled as sweeps by the SCDA. The
% sweeps autodetection works on BIAS HK which has a lower time resolution, and
% may therefore be incorrect at the beginning and end of a labelled region.
S.define_setting('PROCESSING.L2.DETECT_SWEEPS.SCDA.WINDOW_MARGIN_SEC', 120)



%============================================================================
% PROCESSING.RCT_REGEXP.*
% Regular expressions for RCT filenames
% -------------------------------------
%
% NOTES
% -----
% ** (a) When reading *L1* (voltage) datasets, then regular expression are
%        needed for identifying all necessary RCTs.
%    (b) When reading *L1R* (voltage) datasets, then regular expressions are
%        only needed for identifying the BIAS RCT.
% ** BIAS & TDS have previously not followed the correct filenaming
%    convention but does now (2020-11-20).
% ** LFR do not seem to follow the filenaming convenction (2020-11-20)
%    NOTE: LFR RCTs use 2+6+6 digits in the timestamps (they add
%    seconds=2_digits).
% ** Only the last filename in a sorted list of matching filenames will
%    actually be used (BICAS algorithm).
%
%
% OFFICIAL DOCUMENTATION ON RCT FILENAMING CONVENTION
% ---------------------------------------------------
% The RCT filenaming has changed over time and de facto filenaming is confusing.
% One may therefore be forced to set these reqular expression to match filenames
% which do not follow the official filenaming conventions.
% See bicas.tools.rct.create_RCT() for RCT filenaming conventions.
%
%
% EXAMPLES OF DE FACTO RCT FILENAMES (2019 SEPT + LATER)
% ------------------------------------------------------
% BIAS:
%       ROC-SGSE_CAL_RCT-BIAS_V201803211625.cdf   (old impl convention)
%       ROC-SGSE_CAL_RPW_BIAS_V201908231028.cdf   (new impl convention,
%                                                  closer to documentation)
%           SOLO_CAL_RCT-BIAS_V201901141146.cdf   (old impl convention)
%           SOLO_CAL_RPW_BIAS_V202004062127.cdf   (almost correct)
%           solo_CAL_rpw-bias_20220210-20990101_V01.cdf
% LFR:
%       ROC-SGSE_CAL_RCT-LFR-BIAS_V20180724165443.cdf
%           SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf
%           solo_CAL_rpw-lfr-bias_20200210-20990101_V01.cdf
% TDS:
%           SOLO_CAL_RCT-TDS-LFM-CWF-E_V20190128.cdf
%           SOLO_CAL_RCT-TDS-LFM-RSWF-E_V20190128.cdf
%           SOLO_CAL_RPW-TDS-LFM-CWF-E_V20200512000000.cdf
%           SOLO_CAL_RPW-TDS-LFM-RSWF-E_V20200512000000.cdf
%           (Two types of calibration files, but only RODP versions)
%
%============================================================================
CDF_SUFFIX_REGEXP = '\.(cdf|CDF)';
%S.define_setting('PROCESSING.RCT_REGEXP.BIAS',         ['SOLO_CAL_RPW-BIAS_V20[0-9]{10}',                    CDF_SUFFIX_REGEXP]);   % Old illegal filenaming convention
S.define_setting('PROCESSING.RCT_REGEXP.BIAS',         ['solo_CAL_rpw-bias_[0-9]{8}-[0-9]{8}_V[0-9][0-9]+', CDF_SUFFIX_REGEXP]);

% NOTE: 2024-09-02: TDS still uses old/illegal RCT filenaming convention.
% S.define_setting('PROCESSING.RCT_REGEXP.LFR',
% ['SOLO_CAL_RCT-LFR-BIAS_V20[0-9]{12}',                   CDF_SUFFIX_REGEXP]);   % Old illegal filenaming convention.
S.define_setting('PROCESSING.RCT_REGEXP.LFR',          ['solo_CAL_rpw-lfr-bias_20200210-20990101_V[0-9][0-9]+', CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.TDS-LFM-CWF',  ['SOLO_CAL_RPW-TDS-LFM-CWF-E_V20[0-9]{12}',              CDF_SUFFIX_REGEXP]);   % Old illegal filenaming convention.
S.define_setting('PROCESSING.RCT_REGEXP.TDS-LFM-RSWF', ['SOLO_CAL_RPW-TDS-LFM-RSWF-E_V20[0-9]{12}',             CDF_SUFFIX_REGEXP]);   % Old illegal filenaming convention.



% CALIBRATION_TABLE_INDEX2 = Second value in zVar CALIBRATION_TABLE_INDEX
% (in every record), that contains an index to calibration data inside a
% given RCT.
% "L1R" refers to when using L1R datasets as input, as opposed to L1.
S.define_setting('PROCESSING.L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS',      1)
S.define_setting('PROCESSING.L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2',    1)
S.define_setting('PROCESSING.L1R.TDS.CWF.USE_GA_CALIBRATION_TABLE_RCTS',  1)
% CALIBRATION_TABLE_INDEX is not set for TDS. Therefore no such setting for TDS.
S.define_setting('PROCESSING.L1R.TDS.RSWF.USE_GA_CALIBRATION_TABLE_RCTS', 1)



%============================================================================
% Calibration constants for the "scalar" calibration mode
% -------------------------------------------------------
% Unit: IVPAV = Interface volt per antenna volt.
%
% Calibration constants that are used instead of the corresponding BIAS
% transfer functions.
% NOTE: These values do not influence the nominal, "full" calibration. They
%       are entirely separate.
% NOTE: The sign should preferably be consistent with the BIAS transfer
%       functions, i.e. positive values as of 2020-04-27.
% NOTE: There are no equivalent (alternative) scalar values to replace the
%       LFR & TDS transfer functions.
%============================================================================
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.ALPHA_IVPAV',           1/17);
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.BETA_IVPAV',               1);
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.HIGH_GAIN',  100);
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.GAIN.GAMMA_IVPAV.LOW_GAIN',     5);

%============================================================================
% Constants for calibrating bias currents from the HK bias currents
% -----------------------------------------------------------------
% Values taken from BIAS specifications, 01/16, Section 3.4.4.1-3. Not to be
% confused with registers which set the bias command.
%
% NOTE: This is a non-standard way of deriving the bias currents.
%============================================================================
% NOTE: OFFSET_TM value is added to the TM value (not the ampere value).
S.define_setting('PROCESSING.CALIBRATION.CURRENT.HK.OFFSET_TM', -hex2dec('56C0') * [1,1,1])
S.define_setting('PROCESSING.CALIBRATION.CURRENT.HK.GAIN_AAPT', -0.008198754     * [1,1,1])



%===============================================================
% Disable/simplify different parts of the calibration algorithm
%===============================================================
% Disable all voltage calibration. Output dataset data contain TM units.
% BIAS demultiplexer addition/subtraction of BLTS necessary to derive
% antenna signals is still done though.
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.DISABLE',               0);
% Whether to disable BIAS offsets.
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.OFFSETS_DISABLED', 0);
% Whether to use transfer functions or scalar multiplication for calibration
% of signals between antennas and BIAS-LFR/TDS interface. It does not affect
% the LFR/TDS transfer functions.
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF',              'FULL');    % SCALAR, FULL

S.define_setting('PROCESSING.CALIBRATION.TF.METHOD',             'FFT')   % FFT, kernel
%S.define_setting('PROCESSING.CALIBRATION.TF.METHOD',             'kernel')   % FFT, kernel
%S.define_setting('PROCESSING.CALIBRATION.TF.KERNEL.EDGE_POLICY', 'zeros')   % zeros, cyclic, mirror
S.define_setting('PROCESSING.CALIBRATION.TF.KERNEL.EDGE_POLICY', 'mirror')   % zeros, cyclic, mirror
S.define_setting('PROCESSING.CALIBRATION.TF.KERNEL.HANN_WINDOW_ENABLED', false)   % false, true



%===========================================================================
% De-/re-trending
% ---------------
% Whether/how to use:
% (1a) De-trending (remove fit) before applying TF.
% (1b) Which degree of fit to use.
% (2)  Re-trending (add back fit previously removed) after applying TF.
% --
% NOTE: There is no re-trending for AC since detrending is ill-defined for
% non-lowpass filters, which is due to that the linear fit needs to be
% scaled by the TF at 0 Hz, when er-trending which does not make sense for
% non-lowpass filters.
%
% *_FIT_DEGREE >= 0 : Enable de-trending using a polynomial fit of this
%                     degree. 0=constant, 1=linear, and so on.
% *_FIT_DEGREE  < 0 : Disable de-trending (requires disabled re-trending).
%
% YK 2020-11-02: Detrend AC data, but do not add linear fit back.
%===========================================================================
S.define_setting('PROCESSING.CALIBRATION.TF.DC_DE-TRENDING_FIT_DEGREE', 1)
S.define_setting('PROCESSING.CALIBRATION.TF.DC_RE-TRENDING_ENABLED',    true)
S.define_setting('PROCESSING.CALIBRATION.TF.AC_DE-TRENDING_FIT_DEGREE', 0)



% --------------------------------------------
% Frequency above which the ITF is set to zero
% --------------------------------------------
% Expressed as a fraction of the Nyquist frequency (half the sampling
% frequency; 1 sample/s = 1 Hz).
% inf = No limit.
% YK 2020-09-15: Set inverted transfer function to zero for
% omega>0.8*omega_Nyquist (not 0.7).
S.define_setting('PROCESSING.CALIBRATION.TF.HIGH_FREQ_LIMIT_FRACTION',  0.8)

% When using AC, the combined LFR+BIAS ITFs are modified to have constant
% gain between 0 Hz and this frequency. The gain used is taken from this
% frequency. Phase remains unchanged.
% NOTE: Setting it to 0 Hz or lower effectively disables the functionality.
% NOTE: This does not set gain at 0 Hz to zero, but AC de-trending
% influences that.
% NOTE: "BIAS specifications", Section 2.3.2.4 specifies (AC) "a high pass
% filter at 7 Hz".
S.define_setting('PROCESSING.CALIBRATION.TF.AC_CONST_GAIN_LOW_FREQ_HZ', 7)

% Whether to split a time series into shorter time series, separated by fill
% values, before de-trending and applying the (modified) TF. This avoids
% applying TF to fill values which avoids destroying non-fill value data.
S.define_setting('PROCESSING.CALIBRATION.TF.FV_SPLITTING.ENABLED',     true)
% Minimum number of samples in a time series (after splitting).
% NOTE: Limit does not apply if there was no splitting (for "backward
% compatibility").
S.define_setting('PROCESSING.CALIBRATION.TF.FV_SPLITTING.MIN_SAMPLES', 128)

% Whether to disable LFR/TDS transfer functions (but still potentially use
% the BIAS transfer functions). This effectively means that TM voltage
% corresponds to interface volt.
% NOTE: This useful for separately using bicas.proc.L1L2.cal.Cal for analyzing
% BIAS standalone calibration tables (BSACT).
S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.LFR_TDS.TF_DISABLED',  0);



S.disable_define();

Bso = S;

end
