%
% Create settings object with
% (1) all permitted key defined, and
% (2) with all keys set to default values.
%
% NOTE: This function does not declare any global SETTINGS object and does not rely on such one being already defined.
% NOTE: Slightly deceiving name, since it defines which keys are permitted.
%
%
% NOTES ON SETTINGS KEY NAMING CONVENTION
% =======================================
% Some constants (1) correspond exactly to fields in the (JSON) S/W descriptor, and (2) are unlikely to be used for
% anything else. These are labeled with a prefix "SWD_" and the last part is in lowercase (because that is what they are
% in the S/W descriptor). Other variables which do not have the prefix may also be used for the S/W descriptor too but
% they are probably more unambiguous in their meaning.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2018-01-24
%
function SETTINGS = create_default_SETTINGS()
% PROPOSAL: Rename SIMPLE_DEMUXER. SIMPLE_CALIBRATION?
% PROPOSAL: Move STDOUT_PREFIX, LOG_PREFIX to ~constants (not overridable). error_safe_constants?


%-------------------------------------------------------------------------------------
% Values common to multiple settings
% ----------------------------------
% Only used INDIRECTLY and only INTERNALLY to set the values of the "real" constants.
%
% TODO-DECISION: Better to define separate constants, rather than have some settings copy other settings?
%-------------------------------------------------------------------------------------
D = [];
D.AUTHOR_NAME  = 'Erik P G Johansson';
D.AUTHOR_EMAIL = 'erik.johansson@irfu.se';
D.INSTITUTE    = 'IRF-U';
%D.SWD_OUTPUT_RELEASE_VERSION = '01';  % For the S/W descriptor output CDFs' release version. Unknown what a sensible value is.



S = bicas.settings();

S.define_setting('MATLAB_COMMAND', '');

S.define_setting('AUTHOR_NAME',  D.AUTHOR_NAME);
S.define_setting('AUTHOR_EMAIL', D.AUTHOR_EMAIL);
S.define_setting('INSTITUTE',    D.INSTITUTE);

% Value that shows up in EOut dataset GlobalAttributes.Calibration_version.
% String value. TEMPORARY SOLUTION.
S.define_setting('CALIBRATION_VERSION', '0.1; Only proportionality constants i.e. no voltage offset tables, no transfer functions; No bias currents');



%===========================================================================================================
% Various S/W descriptor release data for the entire software (not specific outputs)
% ----------------------------------------------------------------------------------
% EXCEPTION TO VARIABLE NAMING CONVENTION: Field names are used for constructing the JSON object struct and
% can therefore NOT follow variable naming conventions without modifying other code.
%===========================================================================================================
S.define_setting('SWD_IDENTIFICATION.project',     'ROC-SGSE');
S.define_setting('SWD_IDENTIFICATION.name',        'BICAS');
S.define_setting('SWD_IDENTIFICATION.identifier',  'ROC-SGSE-BICAS');
S.define_setting('SWD_IDENTIFICATION.description', 'BIAS Calibration Software (BICAS) which derives the BIAS L2S input signals (plus some) from the BIAS L2R output signals.');
%
S.define_setting('SWD_RELEASE.version',            '0.1.0');
S.define_setting('SWD_RELEASE.date',               '2018-01-22');
S.define_setting('SWD_RELEASE.author',             D.AUTHOR_NAME);
S.define_setting('SWD_RELEASE.contact',            D.AUTHOR_EMAIL);
S.define_setting('SWD_RELEASE.institute',          D.INSTITUTE);
S.define_setting('SWD_RELEASE.modification',       'No modification (initial release)');
%
S.define_setting('SWD_ENVIRONMENT.executable',     'roc/bicas');     % Relative path to BICAS executable. See RCS ICD.



% Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher script.
S.define_setting('STDOUT_PREFIX',                  'STDOUT: ');
% NOTE: Analogous LOG_PREFIX is hard-coded for safety.

% Parameters influencing how JSON objects are printed with function JSON_object_str.
S.define_setting('JSON_OBJECT_STR.INDENT_SIZE',     4);
S.define_setting('JSON_OBJECT_STR.VALUE_POSITION', 15);

% The epoch for ACQUISITION_TIME.
% The time in UTC at which ACQUISITION_TIME is [0,0].
% Year-month-day-hour-minute-second-millisecond-mikrosecond(0-999)-nanoseconds(0-999)
% PROPOSAL: Store the value returned by spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) instead?
S.define_setting('ACQUISITION_TIME_EPOCH_UTC',                   [2000,01,01, 12,00,00, 000,000,000]);

S.define_setting('INPUT_CDF_ASSERTIONS.STRICT_DATASET_ID',       0);    % Require input CDF Global Attribute "DATASET_ID"       to match the expected value.
S.define_setting('INPUT_CDF_ASSERTIONS.STRICT_SKELETON_VERSION', 1);    % Require input CDF Global Attribute "Skeleton_version" to match the expected value.
S.define_setting('INPUT_CDF_ASSERTIONS.MATCHING_TEST_ID',        0);    % Require Test_id to be identical for all input CDF datasets.
S.define_setting('OUTPUT_CDF.SET_TEST_ID',                       1);    % Set CDF GlobalAttribute "Test_id". ROC DFMD says that it should really be set by ROC.
S.define_setting('OUTPUT_CDF.DATA_VERSION',                      '01'); % Set CDF GlobalAttribute "Data_version". ROC DFMD says it should be updated in a way which can not be automatized?!!! Set here for now.

S.define_setting('PROCESSING.USE_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION', 1);

S.define_setting('PROCESSING.ROC_PIP_NAME_OVERRIDE',        '');   % If set, override environment variable ROC_PIP_NAME        defined in RCS ICD. Which pipeline to run, "RGTS" or "RODP".
S.define_setting('PROCESSING.ROC_RCS_CAL_PATH_OVERRIDE',    '');   % If set, override environment variable ROC_RCS_CAL_PATH    defined in RCS ICD. Path to dir. with calibration files.
S.define_setting('PROCESSING.ROC_RCS_MASTER_PATH_OVERRIDE', '');   % If set, override environment variable ROC_RCS_MASTER_PATH defined in RCS ICD. Path to dir. with master CDF files.

S.define_setting('SW_MODES.ENABLE_INPUT_L2R', 1);
S.define_setting('SW_MODES.ENABLE_TDS',       0);

% zVariables which are still empty after copying data into the master CDF assigned a correctly sized array
% with fill values. This should only be necessary for S/W modes with incomplete processing.
S.define_setting('OUTPUT_CDF.EMPTY_ZVARIABLES_SET_TO_FILL', 0);

S.define_setting('LOGGING.MAX_UNIQUES_PRINTED',             5);    % When logging contents of matrix/vector, maximum number of unique values printed before switching to shorter representation (min-max range)



%=====================================================================
% Define constants relating to interpreting LFR datasets
% ------------------------------------------------------
% F0, F1, F2, F3: Frequencies with which samples are taken. Unit: Hz. Names are LFR's naming.
%=====================================================================
S.define_setting('LFR.F0', 24576);  % = 6 * 4096
S.define_setting('LFR.F1',  4096);
S.define_setting('LFR.F2',   256);
S.define_setting('LFR.F3',    16);

%========================================================
% Constants for how the "simple demuxer" calibrates data
%========================================================
S.define_setting('SIMPLE_DEMUXER.ALPHA',           1/17);
S.define_setting('SIMPLE_DEMUXER.BETA',               1);
S.define_setting('SIMPLE_DEMUXER.GAMMA_HIGH_GAIN',  100);
S.define_setting('SIMPLE_DEMUXER.GAMMA_LOW_GAIN',     5);   % NOTE/POSSIBLE BUG: Uncertain which value is high-gain, and low-gain.

S.disable_define();

SETTINGS = S;

end
