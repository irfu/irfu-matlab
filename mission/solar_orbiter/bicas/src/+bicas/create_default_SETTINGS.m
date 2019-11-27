%
% Create settings object and
% (1) define the set of permitted/existing settings keys, and
% (2) set all settings keys to their initial default values.
%
% NOTE: This function does not declare any global SETTINGS object and does not rely on such one being already defined.
% NOTE: Slightly deceiving name, since it defines which keys are permitted.
%
%
% NOTES ON SETTINGS KEY NAMING CONVENTION
% =======================================
% Some constants (1) correspond exactly to fields in the (JSON) S/W descriptor, and (2) are unlikely to be used for
% anything else. These are labeled with a prefix "SWD." and the remainder is in lowercase (because that is what they are
% in the S/W descriptor).
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-24
%
function SETTINGS = create_default_SETTINGS()
% PROPOSAL: Move STDOUT_PREFIX to ~constants (not overridable), error_safe_constants?
% PROPOSAL: Setting for latching relay? Setting for only mode 0?
% PROPOSAL: Better name for *.CALIBRATION.SCALAR.*. Neither indicates simple calibration, nor BIAS calibration.
%   PROPOSAL: BIAS_SIMPLE

S = bicas.settings();

% The MATLAB command (e.g. path) to use to launch MATLAB for BICAS.
% NOTE: Only the value in the BICAS config file is actually used. The normal priority order for how SETTINGS values are
% being obatined does apply here.
S.define_setting('MATLAB_COMMAND', '');



% Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher script.
S.define_setting('STDOUT_PREFIX',                  'STDOUT: ');
% NOTE: Analogous LOG_PREFIX is hard-coded for safety.

% Parameters influencing how JSON objects are printed with function JSON_object_str.
S.define_setting('JSON_OBJECT_STR.INDENT_SIZE',     4);



S.define_setting('INPUT_CDF_ASSERTIONS.STRICT_DATASET_ID', 0);    % Require input CDF Global Attribute "DATASET_ID" to match the expected value.
S.define_setting('INPUT_CDF_ASSERTIONS.MATCHING_TEST_ID',  0);    % Require Test_id to be identical for all input CDF datasets.

S.define_setting('OUTPUT_CDF.SET_TEST_ID',                 1);    % Set CDF GlobalAttribute "Test_id". ROC DFMD says that it should really be set by ROC.
S.define_setting('OUTPUT_CDF.DATA_VERSION',                '01'); % Set CDF GlobalAttribute "Data_version". ROC DFMD says it should be updated in a way which can not be automatized?!!! Set here for now.

% What to do with zVariables which are still empty after copying data into the master CDF.
% This indicates that something is wrong, either in the master CDF or in the processing.
%
S.define_setting('OUTPUT_CDF.EMPTY_NUMERIC_ZVARIABLES_SET_TO_FILL', 0);
% Ex: Non-numeric ACQUISITION_TIME_UNITS in SOLO_L2_RPW-LFR-SBM1-CWF-E_V05.cdf is empty
S.define_setting('OUTPUT_CDF.EMPTY_NONNUMERIC_ZVARIABLES_IGNORE',   1);



% Value that shows up in output dataset GlobalAttributes.Calibration_version.
% Value that is used to set the output dataset GlobalAttribute "Calibration_version". String value. TEMPORARY SOLUTION.
S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.Calibration_version', '0.2; Only (the wrong) proportionality constants i.e. no voltage offset tables, no transfer functions; No bias currents');

% Variables, if non-empty, are used to override the corresponding environment variables.
S.define_setting('ENV_VAR_OVERRIDE.ROC_PIP_NAME',        '');   % ROC_PIP_NAME        defined in RCS ICD. Which pipeline to run, "RGTS" or "RODP".
S.define_setting('ENV_VAR_OVERRIDE.ROC_RCS_CAL_PATH',    '');   % ROC_RCS_CAL_PATH    defined in RCS ICD. Path to dir. with calibration files.
S.define_setting('ENV_VAR_OVERRIDE.ROC_RCS_MASTER_PATH', '');   % ROC_RCS_MASTER_PATH defined in RCS ICD. Path to dir. with master CDF files.

% Whether to enable (make visible & accessible to the user) certain s/w modes.
%S.define_setting('SW_MODES.ENABLE_INPUT_L2R',   0);    % Enable OLD s/w modes which accept L2R input datasets.
S.define_setting('SW_MODES.ENABLE_TDS',         1);    % Enable     s/w modes which accept TDS input datasets. NOTE: Not implemented.

S.define_setting('LOGGING.MAX_UNIQUES_PRINTED', 5);    % When logging contents of matrix/vector, maximum number of unique values printed before switching to shorter representation (min-max range)



% The epoch for ACQUISITION_TIME.
% The time in UTC at which ACQUISITION_TIME is [0,0].
% Year-month-day-hour-minute-second-millisecond-mikrosecond(0-999)-nanoseconds(0-999)
% PROPOSAL: Store the value returned by spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) instead?
S.define_setting('PROCESSING.ACQUISITION_TIME_EPOCH_UTC',                   [2000,01,01, 12,00,00, 000,000,000]);

S.define_setting('PROCESSING.USE_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION', 0);



%===========================================================================================================
% Various S/W descriptor release data for the entire software (not specific outputs)
% ----------------------------------------------------------------------------------
% EXCEPTION TO VARIABLE NAMING CONVENTION: Field names are used for constructing the JSON object struct and
% can therefore NOT follow variable naming conventions without modifying other code.
%===========================================================================================================
S.define_setting('SWD.identification.project',     'ROC');
S.define_setting('SWD.identification.name',        'BIAS Calibration Software (BICAS)');
S.define_setting('SWD.identification.identifier',  'BICAS');
S.define_setting('SWD.identification.description', ...
    '(Incomplete) calibration software meant to (1) calibrate electric field L2 data from electric L1R LFR and TDS (LFM) data, and (2) calibrate bias currents.');
S.define_setting('SWD.identification.icd_version', '1.2');   % Technically wrong. In reality iss1rev2, draft 2019-07-11.

S.define_setting('SWD.release.version',            '0.2.0');
S.define_setting('SWD.release.date',               '2019-09-20');
S.define_setting('SWD.release.author',             'Erik P G Johansson, BIAS team');
S.define_setting('SWD.release.contact',            'erik.johansson@irfu.se');
%S.define_setting('SWD.release.institute',          'IRF-U');
S.define_setting('SWD.release.institute',          'Swedish Institute of Space Physics (IRF)');   % Full name or abbreviation?
S.define_setting('SWD.release.modification',       'Various updates and refactoring; Incomplete but ICD compliant support for TDS data sets & RODP pipeline, L2R-->L1R, updated ICD compliance (pre-release).');
S.define_setting('SWD.release.source',             'https://github.com/irfu/irfu-matlab/commits/SOdevel');    % Appropriate branch?
%S.define_setting('SWD.release.source',             'https://github.com/irfu/irfu-matlab/commits/master');
%
S.define_setting('SWD.environment.executable',     'roc/bicas');   % Relative path to BICAS executable. See RCS ICD.


%=======================================================================================================================
% Regular expressions for the filenames of RCTs
% ---------------------------------------------
% RCT filenaming convention is described in:	ROC-PRO-DAT-NTT-00006-LES, 1/1 draft, Sect 4.3.2-3.
%
% IMPLEMENTATION NOTE: RCT filenaming is implemented as settings since filenaming seems likely to change.
% (1) LFR & TDS do not seem to follow the filenaming convenction
% (2) BIAS has (previously at least) not followed the filenaming convention.
% (3) it is uncertain how it can be applied to BIAS RCTs (which receiver should the BIAS RCT specify when BIAS uses the
% same RCT for both LFR & TDS data?).
%
% NOTE: LFR RCTs use 2+6+6 digits instead of 2+6 in the timestamps (they add seconds=2 digits). NOTE: TDS RCTs use 2+6
% digits instead of 10 in the timestamps (the have no time of day, only date)
%
% Examples of RCT filenames (2019 Sept)
% -------------------------------------
% BIAS:
%       ROC-SGSE_CAL_RCT-BIAS_V201803211625.cdf   (old implemented convention)
%       ROC-SGSE_CAL_RPW_BIAS_V201908231028.cdf   (new implemented convention, closer to documentation)
%           SOLO_CAL_RCT-BIAS_V201901141146.cdf   (old implemented convention)
% LFR:
%       ROC-SGSE_CAL_RCT-LFR-BIAS_V20180724165443.cdf
%           SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf
% TDS:
%           SOLO_CAL_RCT-TDS-LFM-CWF-E_V20190128.cdf
%           SOLO_CAL_RCT-TDS-LFM-RSWF-E_V20190128.cdf
%           (Two types of calibration files, but only RODP versions)
%
% NOTE: Only the last filename in a sorted list of matching filenames will actually be used.
%=======================================================================================================================
CDF_SUFFIX_REGEXP = '\.(cdf|CDF)';
%S.define_setting('PROCESSING.RCT_REGEXP.RGTS.BIAS',         ['ROC-SGSE_CAL_RCT-BIAS_V20[0-9]{10}',          CDF_SUFFIX_REGEXP]);   % Old BIAS convention
S.define_setting('PROCESSING.RCT_REGEXP.RGTS.BIAS',         ['ROC-SGSE_CAL_RPW_BIAS_V20[0-9]{10}',          CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.RODP.BIAS',         [    'SOLO_CAL_RPW_BIAS_V20[0-9]{10}',          CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.RGTS.LFR',          ['ROC-SGSE_CAL_RCT-LFR-BIAS_V20[0-9]{12}',      CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.RODP.LFR',          [    'SOLO_CAL_RCT-LFR-BIAS_V20[0-9]{12}',      CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.RGTS.TDS-LFM-CWF',  ['ROC-SGSE_CAL_RCT-TDS-LFM-CWF-E_V20[0-9]{6}',  CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.RODP.TDS-LFM-CWF',  [    'SOLO_CAL_RCT-TDS-LFM-CWF-E_V20[0-9]{6}',  CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.RGTS.TDS-LFM-RSWF', ['ROC-SGSE_CAL_RCT-TDS-LFM-RSWF-E_V20[0-9]{6}', CDF_SUFFIX_REGEXP]);
S.define_setting('PROCESSING.RCT_REGEXP.RODP.TDS-LFM-RSWF', [    'SOLO_CAL_RCT-TDS-LFM-RSWF-E_V20[0-9]{6}', CDF_SUFFIX_REGEXP]);

S.define_setting('PROCESSING.RCT.LFR.EXTRAPOLATE_TF_AMOUNT_HZ', 1);



%====================================================================================================================
% Define constants relating to interpreting LFR datasets
% ------------------------------------------------------
% F0, F1, F2, F3: Frequencies with which samples are taken. The variables names (F[0-3]) follow LFR's naming scheme.
% Unit: Hz.
% Only used for sequences of sample on the same CDF record.
%====================================================================================================================
S.define_setting('PROCESSING.LFR.F0_HZ', 24576);  % = 6 * 4096
S.define_setting('PROCESSING.LFR.F1_HZ',  4096);
S.define_setting('PROCESSING.LFR.F2_HZ',   256);
S.define_setting('PROCESSING.LFR.F3_HZ',    16);

%========================================================
% Constants for how the "simple demuxer" calibrates data
%========================================================
%S.define_setting('PROCESSING.CALIBRATION.SCALAR.ALPHA',           1/17);
%S.define_setting('PROCESSING.CALIBRATION.SCALAR.BETA',               1);
%S.define_setting('PROCESSING.CALIBRATION.SCALAR.GAMMA_HIGH_GAIN',  100);
%S.define_setting('PROCESSING.CALIBRATION.SCALAR.GAMMA_LOW_GAIN',     5);

%=================================================================================
% Constants for using HK bias currents for deriving/calibrating the bias currents
% NOTE: This is a non-standard way of deriving the bias currents.
% NOTE 2019-09-12: THIS HAS NOT BEEN IMPLEMENTED IN THE CODE YET EXCEPT FOR CALIBRATION FUNCTION.
%=================================================================================
%S.define_setting('PROCESSING.USE_UNCALIBRATED_BIAS_CURRENTS_FROM_HK', 0);
% NOTE: OFFSET_TM value is added to the TM value (not the ampere value).
S.define_setting('PROCESSING.CALIBRATION.HK_BIAS_CURRENT.OFFSET_TM', -hex2dec('56C0') * [1,1,1])
S.define_setting('PROCESSING.CALIBRATION.HK_BIAS_CURRENT.GAIN_APC',  -0.008198754     * [1,1,1])

S.define_setting('PROCESSING.CALIBRATION.DISABLE_CALIBRATION', 0)
S.define_setting('PROCESSING.CALIBRATION.DETRENDING_ENABLED',  0)

S.disable_define();

SETTINGS = S;

end
