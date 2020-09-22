%
% Create settings object and
% (1) define the set of permitted/existing settings keys, and
% (2) set all settings keys to their initial default values.
%
% NOTE: This function does not declare any global SETTINGS object and does not rely on such one being already defined.
% NOTE: Slightly deceiving name, since it defines which keys are permitted.
%
%
% NAMING CONVENTIONS
% ==================
% ZV  : zVariable
% GA  : Global Attribute (in CDF file)
% CUR : CURRENT (type of data, dataset)
%
%
% NOTES ON SETTINGS KEY NAMING CONVENTION
% =======================================
% Some constants (1) correspond exactly to fields in the (JSON) S/W descriptor, and (2) are unlikely to be used for
% anything else. These are labeled with a prefix "SWD." and the remainder is in lowercase (because that is what they are
% in the S/W descriptor).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2018-01-24
%
function SETTINGS = create_default_SETTINGS()
    % PROPOSAL: Setting for latching relay? Setting for only mode 0?
    %
    % PROPOSAL: Need (settings name) terminology for temporary "bugfixes"/corrections due to bugs outside of BICAS,
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
    % PROPOSAL: PROCESSING.CALIBRATION.CURRENT.HK.DISABLE      : Whether to calibrate HK current or use HK TM. Not which data to use (HK or TC).
    %           PROCESSING.CALIBRATION.CURRENT.SOURCE = TC, HK : Which data to use.
    %
    % PROPOSAL: Abolish INPUT_CDF.HK.MOVE_TIME_TO_SCI.
    % PROPOSAL: Naming convention for settings keys for testing ONLY:
    %
    % PROPOSAL: Setting keys should used cased version of zVars and glob.attrs..
    %   Ex: Epoch, (GA) Test_id, (GA) Dataset_ID.
    %
    % PROPOSAL: Setting keys should always be on the form ENABLE, never DISABLE.
    %   PRO: More consistent.
    %   CON: Less clear what is a deviation from the default.
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
    % PROPOSAL: INPUT_CDF.* : Settings that apply to ALL input datasets.
    % PROPOSAL: Only INPUT_CDF.ALL.* apply to all input datasets.
    %
    % PROPOSAL: Policy on usage of dataset levels in settings keys.
    %   Ex: For now, L1R refers to algorithms to use WHEN processing L1R as input
    %       Ex: PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY
    %       NOTE: Does not have to refer to the L1R datasets as such.
    %           Ex: L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS
    %           Ex: L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2
    %   Ex: For now, L2 refers to algorithms to use when processing L2 as output.
    %       Ex: PROCESSING.L2.REMOVE_DATA.MUX_MODES
    %   NEED: Specify whether refers to input or output data (not necessarily datasets).
    %       Ex: Distinguish processing L1R-->L2, L2-->L3.
    %   NEED: Distinguish processing (science) data L1-->L2, L1R-->L2
    %   NEED: Distinguish
    %       (1) datasets of specific level, and
    %       (2) algorithms that run when using given input level.
    %   PROPOSAL: Some kind of prefix before data level.
    %       PROPOSAL: IN_L1R, OUT_L2 etc.
    %   PROPOSAL: when speaking of datasets, use
    %       (1) INPUT_CDF.<level>  : How to interpret, read datasets
    %       (2) OUTPUT_CDF.<level> : How to output, write datasets.
    %       PROBLEM: How distinguish from processing?
    %
    % PROPOSAL: Setting name change SW_MODES.L1_LFR_TDS_ENABLED--> SW_MODES.L1-L2_LFR_TDS_ENABLED
    %   NOTE: Name change likely influences BICAS testing code and pipeline. Should therefore only be implemented at the
    %         right time.
    %
    % PROPOSAL: Separate function for validating/asserting settings.
    %   NOTE: Must be done AFTER all settings have been set.
    %   PROPOSAL: Do every time settings are set, i.e. for default values,
    %       config file values, CLI argument values.



    S = bicas.settings();

    % The MATLAB command (e.g. path) to use to launch MATLAB for BICAS.
    % NOTE: Only the value in the BICAS config file is actually used. The normal priority order for how SETTINGS values are
    % being obtained does apply here but does not matter since the value is only used by the bash wrapper script for
    % launching MATLAB.
    S.define_setting('MATLAB_COMMAND', '');
    
    
    
    % Prefix used to identify the subset of stdout that should actually be passed on as stdout by the bash launcher
    % script.
    %S.define_setting('STDOUT_PREFIX',               'STDOUT: ');
    % NOTE: Analogous LOG_PREFIX is hard-coded for safety.
    
    % Parameters influencing how JSON objects are printed with function JSON_object_str.
    S.define_setting('JSON_OBJECT_STR.INDENT_SIZE', 4);
    
    % When logging contents of matrix/vector, maximum number of unique values printed before switching to shorter
    % representation (min-max range)
    S.define_setting('LOGGING.MAX_NUMERIC_UNIQUES_PRINTED', 5);
    % When logging contents of TT2000 vector (in practise zVar Epoch), maximum number of unique TT2000 values printed
    % before switching to shorter representation (min-max range).
    S.define_setting('LOGGING.MAX_TT2000_UNIQUES_PRINTED', 2);
    
    % Enable inofficial (to ROC) support for S/W modes
    % ------------------------------------------------
    % Enable s/w modes for processing LFR & TDS datasets L1-->L2 in addition to the official support
    % for L1R. LFR_TDS refers to LFR/TDS input datasets, as opposed to L1 current datasets.
    S.define_setting('SW_MODES.L1_LFR_TDS_ENABLED', 0);
    % Enable s/w modes for processing L2-->L3 datasets.
    S.define_setting('SW_MODES.L2-L3_ENABLED',      0);
    
    
    
    %###########################################################################################################
    % SWD.*
    % Various S/W descriptor (SWD) release data for the entire software (not specific outputs)
    % ----------------------------------------------------------------------------------------
    % EXCEPTION TO VARIABLE NAMING CONVENTION: Field names are used for constructing the JSON object struct and
    % can therefore NOT follow variable naming conventions without modifying other code.
    %
    % ROC-GEN-SYS-NTT-00019-LES, "ROC Engineering Guidelines for External Users":
    % """"""""
    % 2.2.3 RCS versioning
    % The RCS version must be a unique number sequence identifier “X.Y.Z”, where “X” is an
    % integer indicating the release (major changes, not necessarily retro-compatible), “Y” is an
    % integer indicating the issue (minor changes, necessarily retro-compatible) and “Z” is an
    % integer indicating a revision (e.g., bug correction).
    % """"""""
    %###########################################################################################################
    IRF_LONG_NAME = 'Swedish Institute of Space Physics (IRF)';
    %
    S.define_setting('SWD.identification.project',     'ROC');
    S.define_setting('SWD.identification.name',        'BIAS Calibration Software (BICAS)');
    S.define_setting('SWD.identification.identifier',  'BICAS');
    S.define_setting('SWD.identification.description', ...
        ['Calibration software meant to', ...
        ' (1) calibrate electric field L2 data from electric L1R LFR and TDS (LFM) data, and', ...
        ' (2) calibrate bias currents from L1 data.']);
    S.define_setting('SWD.identification.icd_version', '1.2');   % Technically wrong. In reality iss1rev2, draft 2019-07-11.
    S.define_setting('SWD.release.version',            '3.1.1');
    S.define_setting('SWD.release.date',               '2020-09-15T11:22:00Z');
    S.define_setting('SWD.release.author',             'Erik P G Johansson, BIAS team, IRF');
    S.define_setting('SWD.release.contact',            'erjo@irfu.se');
    S.define_setting('SWD.release.institute',          IRF_LONG_NAME);   % Full name or abbreviation?
    % 'Various updates and refactoring; close to complete support for LFR & TDS datasets (but untested); Removed ROC-SGSE_* dataset support.'
    % 'Almost-complete support for LFR & TDS datasets (voltages) with transfer functions (partially tested).'
    S.define_setting('SWD.release.modification',       'Modified default settings: inverted transfer function cutoff at 0.8*omega_Nyquist, duplicate bias current gives error');   % 3.1.1
    S.define_setting('SWD.release.source',             'https://github.com/irfu/irfu-matlab/commits/SOdevel');    % Appropriate branch? "master" instead?
    %
    S.define_setting('SWD.environment.executable',     'roc/bicas');   % Relative path to BICAS executable. See RCS ICD.
    % NOTE: See also OUTPUT_CDF.GLOBAL_ATTRIBUTES.Calibration_version.
    
    
    
    %####################
    % ENV_VAR_OVERRIDE.*
    %####################
    % Variables which, if non-empty, are used to override the corresponding
    % environment variables.
    S.define_setting('ENV_VAR_OVERRIDE.ROC_RCS_CAL_PATH',    '');   % ROC_RCS_CAL_PATH    defined in RCS ICD. Path to dir. with calibration files.
    S.define_setting('ENV_VAR_OVERRIDE.ROC_RCS_MASTER_PATH', '');   % ROC_RCS_MASTER_PATH defined in RCS ICD. Path to dir. with master CDF files.
    
    
    
    %######################################################
    % INPUT_CDF.*
    %######################################################
    
    % The epoch for ACQUISITION_TIME.
    % The time in UTC at which ACQUISITION_TIME is [0,0].
    % Year-month-day-hour-minute-second-millisecond-mikrosecond(0-999)-nanoseconds(0-999)
    % PROPOSAL: Store the value returned by spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) instead?
    S.define_setting('INPUT_CDF.ACQUISITION_TIME_EPOCH_UTC',                       [2000,01,01, 12,00,00, 000,000,000]);
    
    % NOTE: Requires INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY = non-error.
    S.define_setting('INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED', 1)
    % NOTE: See INPUT_CDF.LFR.BOTH_SYNCHRO_FLAG_AND_TIME_SYNCHRO_FLAG_WORKAROUND_ENABLED
    S.define_setting('INPUT_CDF.USING_ZV_NAME_VARIANT_POLICY',     'WARNING')    % WARNING, ERROR
    
    S.define_setting('INPUT_CDF.USING_GA_NAME_VARIANT_POLICY',     'WARNING')    % WARNING, ERROR
    
    % Require input CDF Global Attribute "DATASET_ID" to match the expected value.
    S.define_setting('INPUT_CDF.GA_DATASET_ID_MISMATCH_POLICY',    'WARNING')    % ERROR, WARNING
    S.define_setting('INPUT_CDF.GA_PROVIDER_MISMATCH_POLICY',      'WARNING')    % ERROR, WARNING
    
    % NOTE: This modification applies BEFORE
    % PROCESSING.HK.USE_ZV_ACQUISITION_TIME and therefore always applies to zVar
    % Epoch.
    % NOTE: Only check for increasing, not monotonically.
    % NOTE: There is a known, mitigatable anomaly in SOLO_L1_RPW_BIA-CURRENT
    % which duplicated settings (same timestamp, same bias setting on same
    % antenna) which would be triggered by an assertion on an assert on
    % monotonically increasing timestamps.
    S.define_setting('INPUT_CDF.NON-INCREMENTING_ZV_EPOCH_POLICY',                   'ERROR')      % ERROR, WARNING, SORT
    
    S.define_setting('INPUT_CDF.CUR.DUPLICATE_BIAS_CURRENT_SETTINGS_POLICY', 'ERROR')    % ERROR, REMOVE_DUPLICATES
    
    % Whether to replace PAD VALUES values with NaN internally.
    % IMPORTANT NOTE: Refers to CDF PAD VALUES, NOT CDF FILL VALUES!
    % NOTE: SOLO_L1_RPW-BIA-CURRENT_V06.skt uses pad value=zero (BUG). Therefore
    % useful.
    S.define_setting('INPUT_CDF.REPLACE_PAD_VALUE_DISABLED',       1)            % 0/false, 1/true.
    
    % List of zVar names for which alternate fill value should be used when the
    % zVars are loaded and interpreted.
    %S.define_setting('INPUT_CDF.OVERRIDE_FILL_VALUE.ZV_NAMES',     {'IBIAS_1', 'IBIAS_2', 'IBIAS_3'})
    S.define_setting('INPUT_CDF.OVERRIDE_FILL_VALUE.ZV_NAMES',     cell(0,1))
    % Alternate fill value to use.
    S.define_setting('INPUT_CDF.OVERRIDE_FILL_VALUE.FILL_VALUE',   single(-1e31))

    % For testing, when HK and SCI time are completely different and do not
    % overlap (though HK time still has to cover a larger interval than SCI).
    % Adds/subtracts HK time so that the first HK timestamp equals the first SCI
    % timestamp.
    S.define_setting('INPUT_CDF.HK.MOVE_TIME_TO_SCI',          0)
    
    
    
    %############################################
    % OUTPUT_CDF.*
    % ------------
    % Settings that apply to ALL output datasets
    %############################################
    % Flag to disable writing output files AFTER PROCESSING. Useful for
    % debugging.
    S.define_setting('OUTPUT_CDF.WRITE_FILE_DISABLED',            0)
    % What BICAS should do when there is a pre-existing file on a output dataset
    % file path.
    % NOTE: Not known if the RCS ICD says anything about what should be the
    % default, or what ROC thinks it should be.
    S.define_setting('OUTPUT_CDF.PREEXISTING_OUTPUT_FILE_POLICY', 'WARNING');    % ERROR, WARNING.
    % Disable processing, but generate empty output files. Useful for debugging
    % code that calls BICAS many times (batch processing) and when dataset
    % content is unimportant since it speeds up BICAS.
    S.define_setting('OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE',       0)
    
    % Value that shows up in output dataset
    % GlobalAttributes.Calibration_version. Value that is used to set the output
    % dataset GlobalAttribute "Calibration_version". String value.
    S.define_setting('OUTPUT_CDF.GLOBAL_ATTRIBUTES.Calibration_version', ...
        '1.1; Voltages: Using combined BIAS and LFR/TDS transfer functions (freq. dependent), BIAS offsets. Calibrates currents.');
    
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
    % Ex: Non-numeric ACQUISITION_TIME_UNITS in (master?) SOLO_L2_RPW-LFR-SBM1-CWF-E_V05.cdf is empty
    % Ex: VDC_LABEL etc can be empty due to ROC bug updating skeletons.
    S.define_setting('OUTPUT_CDF.EMPTY_NONNUMERIC_ZV_POLICY', 'ERROR');   % ERROR, WARNING



    % NOTE: ACQUSITION_TIME_UNITS being empty in the master CDF requires value
    % 0/false. Abolish?
    S.define_setting('OUTPUT_CDF.write_CDF_dataobj.strictEmptyZvClass',                 0)
    % Whether the size per record of an empty (0 records) output DF zVar has to
    % be in agreement with the master CDF's size per record.
    % NOTE: ACQUSITION_TIME_UNITS being empty in the master CDF requires value
    % 0/false. Abolish?
    S.define_setting('OUTPUT_CDF.write_CDF_dataobj.strictEmptyNumericZvSizePerRecord',  0)
    % Whether the size per record of an output CDF zVar has to be in agreement
    % with the master CDF's size per record. This is useful if the master CDF
    % has not been updated in this regard only.
    S.define_setting('OUTPUT_CDF.write_CDF_dataobj.strictNumericZvSizePerRecord',       0)
    
    
    
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
    
    S.define_setting('PROCESSING.HK.SCI_TIME_NONOVERLAP_POLICY',       'ERROR')      % WARNING, ERROR
    % NOTE: "WARNING": Will lead to using nearest interpolation.
    S.define_setting('PROCESSING.HK.TIME_NOT_SUPERSET_OF_SCI_POLICY',  'WARNING')    % WARNING, ERROR
    S.define_setting('PROCESSING.CUR.TIME_NOT_SUPERSET_OF_SCI_POLICY', 'WARNING')    % WARNING, ERROR
    
    % Quick ~BUGFIX for bad values in zv SAMPLING_RATE in L1R TDS-LFM-RSWF
    % datasets. Abolish?
    S.define_setting('PROCESSING.L1R.TDS.RSWF_ZV_SAMPLING_RATE_255_POLICY', 'ERROR')   % WARNING, ERROR, CORRECT
    
    % Mitigation: How to handle that LFR zVars QUALITY_FLAG QUALITY_BITMASK are
    % empty.
    S.define_setting('PROCESSING.L1R.LFR.ZV_QUALITY_FLAG_BITMASK_EMPTY_POLICY', 'USE_FILL_VALUE')   % ERROR, USE_FILL_VALUE
    
    % ~BUGFIX for bug in L1/L1R TDS-LFM RSWF datasets.
    % TDS has bugfixed. /2019-12-19
    % PROPOSAL: Rename.
    S.define_setting('PROCESSING.TDS.RSWF.ILLEGAL_ZV_SAMPS_PER_CH_POLICY', 'ERROR')   % ERROR, WARNING, ROUND
    
    %============================================================================
    % Settings for when to remove data by setting it to fill value
    % ------------------------------------------------------------
    % "L2" refers to output datasets. Both voltage and current data. In
    % practise, this functionality is there as a temporary solution for removing
    % sweeps.
    %============================================================================
    S.define_setting('PROCESSING.L2.REMOVE_DATA.MUX_MODES', [1,2,3,4,5,6,7])
    % Unit: S = Seconds
    % Lower number since using LFR mux mode (unless configured not to), which
    % has same cadence as science data.
    % See PROCESSING.LFR.MUX_MODE_SOURCE.
    S.define_setting('PROCESSING.L2.LFR.REMOVE_DATA.MUX_MODE.MARGIN_S',  0)    
    % Higher number since using BIAS HK for TDS, which means that the mux mode
    % is known with a lower time resolution.
    S.define_setting('PROCESSING.L2.TDS.REMOVE_DATA.MUX_MODE.MARGIN_S', 30)

    %============================================================================
    % Where to obtain the mux mode
    % ----------------------------
    % BIAS HK data
    % ------------
    % Contains mux mode using its own Epoch (typically ~30 s time resolution?),
    % which means that ~interpolation to SCI data is necessary, which means that
    % the effective mux mode value can briefly be wrong. BIAS HK may also
    % potentially not cover the same time range as SCI data at all, and then the
    % mux mode can be really wrong (e.g. when using different versions).
    %
    % LFR SCI data (L1/L1R)
    % ---------------------
    % Contains a zVar for mux mode using the same Epoch as the data.
    %
    % NOTE: The relevant TDS datasets do not contain mux mode.
    %============================================================================
    S.define_setting('PROCESSING.LFR.MUX_MODE_SOURCE', 'LFR_SCI')    % BIAS_HK, LFR_SCI

    %====================================================================
    % TEMPORARY SOLUTION
    % Maximum value for zVar QUALITY_FLAG in output datasets.
    % YK 2020-08-31: Use 2=Survey data, possibly not publication-quality
    %====================================================================
    S.define_setting('PROCESSING.ZV_QUALITY_FLAG_MAX', 2)
    
    % Path to RCS NSO file. Relative to BICAS root.
    S.define_setting('PROCESSING.RCS_NSO_FILE.RELATIVE_PATH', fullfile('data', 'solo_ns_ops.xml'))
    
    
    
    %===================================================================================================================
    % PROCESSING.RCT_REGEXP.*
    % Regular expressions for the filenames of RCTs
    % ---------------------------------------------
    %
    % OFFICIAL DOCUMENTATION ON RCT FILENAME CONVENTION
    % =================================================
    % RCT filenaming convention is described in ROC-PRO-DAT-NTT-00006-LES. This document refers to the RODP.
    %
    % Version 1/1:
    % """"""""
    %   4.3.2 RCT data versioning convention
    %
    %   The version of the RCT CDF data file must be the local date and time of creation of the file,
    %   in the format: “YYYYMMDDHHNN”, where “YYYY”, “MM”, “DD”, “HH” and “NN” are
    %    respectively the 4-digits year, 2-digits month, 2-digits day, 2-digits hours, 2-digits minutes of
    %    the file creation.
    %    In the RCT filename, the version number must appear with the “V” prefix (e.g.,
    %    “V202210122359”.
    %
    %
    %   4.3.3 RCT file naming convention
    %
    %   The RCT shall comply the following file naming convention:
    %   SOLO_CAL_RPW-[receiver]_[free-field]_[Version].cdf
    %   Where [receiver] is the name of the receiver in uppercase characters (i.e., “TDS” or
    %   “LFR”) of the corresponding RPW L1R dataset, [free-field] is a field that can be used to
    %   specify the content of the file (e.g., “BIAS-F0”) and [Version] is the version of the
    %   calibration table file (see previous section).
    %   Note that this RCT naming convention is not fully compliant with the SOC definition [AD1]. /.../
    % """"""""
    % Version 1/2, draft:
    % Section 2.2.6.3-4: Slightly different filenaming convention:
    % """"""""
    %   2.2.6.3 File naming
    %   The CAL file must comply the following naming convention:
    %   SOLO_CAL_[Descriptor]_[free-field]_V[CALIBRATION_VERSION].cdf
    %
    %   Where [Descriptor], [free-field] and [CALIBRATION_VERSION] correspond
    %   respectively to the short value in the “Descriptor”, “Free_field” and
    %   “CALIBRATION_VERSION” global attributes (see section 2.2.6.6).
    %   N.B. The CAL file naming convention is not fully compliant with the SOC definition [AD1]. /.../
    % """"""""
    %
    % Xavier Bonnin, 2020-07-02:
    %   Wrong:   SOLO_CAL_RPW_BIAS_V202004062127.cdf   # NOTE: No minus!
    %   Correct: SOLO_CAL_RPW-BIAS_V202004062127.cdf   # NOTE: Minus!
    %
    %
    % RATIONALE
    % =========
    % RCT filenaming is implemented as settings since filenaming seems likely to
    % change.
    % (1) LFR & TDS do not seem to follow the filenaming convenction
    % (2) BIAS has (previously at least) not followed the filenaming convention.
    % (3) it is uncertain how it (doc version 1/1) can be applied to BIAS RCTs
    % (which receiver should the BIAS RCT specify when BIAS uses the same RCT
    % for both LFR & TDS data?).
    %
    % NOTE: LFR RCTs use 2+6+6 digits in the timestamps (they add seconds=2 digits).
    % NOTE: TDS RCTs use 2+6+0 digits in the timestamps (the have no time of day, only date)
    %
    % Examples of de facto RCT filenames (2019 Sept + later)
    % ------------------------------------------------------
    % BIAS:
    %       ROC-SGSE_CAL_RCT-BIAS_V201803211625.cdf   (old implemented convention)
    %       ROC-SGSE_CAL_RPW_BIAS_V201908231028.cdf   (new implemented convention, closer to documentation)
    %           SOLO_CAL_RCT-BIAS_V201901141146.cdf   (old implemented convention)
    %           SOLO_CAL_RPW_BIAS_V202004062127.cdf   (almost correct)
    % LFR:
    %       ROC-SGSE_CAL_RCT-LFR-BIAS_V20180724165443.cdf
    %           SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf
    % TDS:
    %           SOLO_CAL_RCT-TDS-LFM-CWF-E_V20190128.cdf
    %           SOLO_CAL_RCT-TDS-LFM-RSWF-E_V20190128.cdf
    %           (Two types of calibration files, but only RODP versions)
    %
    % NOTE: Only the last filename in a sorted list of matching filenames will
    % actually be used.
    %===================================================================================================================
    CDF_SUFFIX_REGEXP = '\.(cdf|CDF)';
    S.define_setting('PROCESSING.RCT_REGEXP.BIAS',         ['SOLO_CAL_RPW-BIAS_V20[0-9]{10}',          CDF_SUFFIX_REGEXP]);
    S.define_setting('PROCESSING.RCT_REGEXP.LFR',          ['SOLO_CAL_RCT-LFR-BIAS_V20[0-9]{12}',      CDF_SUFFIX_REGEXP]);
    S.define_setting('PROCESSING.RCT_REGEXP.TDS-LFM-CWF',  ['SOLO_CAL_RCT-TDS-LFM-CWF-E_V20[0-9]{6}',  CDF_SUFFIX_REGEXP]);
    S.define_setting('PROCESSING.RCT_REGEXP.TDS-LFM-RSWF', ['SOLO_CAL_RCT-TDS-LFM-RSWF-E_V20[0-9]{6}', CDF_SUFFIX_REGEXP]);
    

    
    % CALIBRATION_TABLE_INDEX2 = Second value in zVar CALIBRATION_TABLE_INDEX
    % (in every record), that contains an index to calibration data inside a
    % given RCT.
    % "L1R" refers to when using L1R datasets as input, as opposed to L1.
    S.define_setting('PROCESSING.L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS',      1)
    S.define_setting('PROCESSING.L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2',    1)
    S.define_setting('PROCESSING.L1R.TDS.CWF.USE_GA_CALIBRATION_TABLE_RCTS',  1)
    % CALIBRATION_TABLE_INDEX is not set for TDS. Therefore no such setting for TDS.
    S.define_setting('PROCESSING.L1R.TDS.RSWF.USE_GA_CALIBRATION_TABLE_RCTS', 1)
    
    
    
    %======================================================================
    % EXPERIMENTAL: LFR sampling frequency-dependent offsets.
    %
    % Values obtained from manually fitting F0,F1,F2 (not F3) snapshots in
    % ROC-SGSE_L1R_RPW-LFR-SURV-SWF-E_59e82ff_CNE_V02.cdf.
    % NOTE: Values are relative as the absolute level is not known.
    % NOTE: Might be that LFR offsets also depend on BLTS.
    % NOTE: Has not set any value for F3.
    %======================================================================
    %S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.LFR.LSF_OFFSETS_TM', [-638, -610, 0, 0])
    S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.LFR.LSF_OFFSETS_TM', [0, 0, 0, 0])
    
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
    S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.DISABLE',              0);
    % Whether to disable BIAS offsets.
    S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.DISABLE_OFFSETS', 0);
    % Whether to use transfer functions or scalar multiplication for calibration
    % of signals between antennas and BIAS-LFR/TDS interface. It does not affect
    % the LFR/TDS transfer functions.
    S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.BIAS.TF',              'FULL');    % SCALAR, FULL
    % Whether to use de-trending before applying transfer functions.
    S.define_setting('PROCESSING.CALIBRATION.TF_DETRENDING_ENABLED',        1)
    % Frequency above which the ITF is set to zero.
    % Expressed as a fraction of the Nyquist frequency (half the sampling
    % frequency; 1 sample/s = 1 Hz).
    % inf = No limit.
    % YK 2020-09-15: Set inverted transfer function to zero for
    % omega>0.8*omega_Nyquist (not 0.7).
    S.define_setting('PROCESSING.CALIBRATION.TF_HIGH_FREQ_LIMIT_FRACTION',  0.8)
    
    % Whether to disable LFR/TDS transfer functions (but still potentially use
    % the BIAS transfer functions). This effectively means that TM voltage
    % corresponds to interface volt.
    % NOTE: This useful for separately using bicas.calib for analyzing BIAS
    % standalone calibration tables (BSACT).
    S.define_setting('PROCESSING.CALIBRATION.VOLTAGE.LFR_TDS.TF_DISABLED',  0);
    
    
    
    S.disable_define();
    
    SETTINGS = S;
    
end
