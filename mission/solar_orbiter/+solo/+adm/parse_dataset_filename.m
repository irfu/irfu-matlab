%
% Given a filename, interpret it as a SO dataset and parse it.
%
% Primarily meant to:
% ** Distinguish between datasets and non-datasets.
% ** Identify DATASET_ID using filenames, so that the files can be further
%    classfied using that, and without trying to read CDF global attributes.
% ** Centralize the parsing of filenames, the understanding of filenaming
%    conventions (as far as relevant).
%
%
% NOTES
% =====
% NOTE: Meant to cover DATASET_IDs widely.
%   ** Old DATASET_IDs beginning with ROC-SGSE, in order to cover old filenames.
%   ** Datasets unrelated to BIAS processing, including other instruments.
%   ** RCTs / CAL level.
%   ** Lower case filenames.
%   ** Recognize -cdag
% NOTE: Lowercase DATASET_ID in filenames are recognized, and returned as
%       uppercase DATASET_ID.
% NOTE: Does not work on RCTs (technically has no DATASET_ID).
% NOTE: Should be possible to use together with
%       solo.adm.create_dataset_filename().
% NOTE/BUG: Can not handle e.g.
%       solo_L1_swa-eas2-NM3D_20201027T000007-20201027T030817_V01.cdf
%       since it has mixed case in the descriptor (not level).
% NOTE: RCTs counts as datasets in this context since the conform to the same
%       filenaming convention and are described in SOL-SGS-TN-0009, "Metadata
%       Definition for Solar Orbiter Science Data", 01/06, though using its own
%       level "CAL".
%
%
% RECOGNIZED FILENAMING CONVENTION
% ================================
% 1    = Basename suffix
% 2    = Additions to filename that
% LCUC = LowerCase or UpperCase
%
%
% RECOGNIZED FILENAMING CONVENTION: EXAMPLES
% ==========================================
%
% Datasets from ground testing:
% -----------------------------
%   solo_L1_rpw-tds-lfm-rswf_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
%   ROC-SGSE_HK_RPW-BIA_19850de_CNE_V02.cdf
%
% Datasets from flight
% --------------------
%   solo_HK_rpw-bia_20200301_V01.cdf                   # NOTE: No -cdag.
%   solo_L2_rpw-lfr-surv-cwf-e-cdag_20200213_V01.cdf   # NOTE: -cdag.
%   solo_L1_rpw-bia-sweep-cdag_20200307T053018-20200307T053330_V01.cdf
%       NOTE: Time interval, not day.
%   solo_L1_rpw-bia-current-cdag_20200401T000000-20200421T000000_V01.cdf
%       NOTE: Should eventually be phased out for currents (not sweeps).
%             /EJ+XB-mail 2020-05-27
%   solo_L1_rpw-bia-current-cdag_20200301-20200331_V01.cdf
%       NOTE: Future replacement for currents (not sweeps).
%             /EJ+XB-mail 2020-05-27
%
%
% UN-RECOGNIZED FILENAMING CONVENTIONS: EXAMPLES
% ==============================================
% Master CDFs:
%   SOLO_HK_RPW-BIA_V01.cdf
%   SOLO_L2_RPW-LFR-SURV-CWF-E_V04.cdf
% Draft Master CDFs:
%   SOLO_L1_RPW-BIA-CURRENT_V01.Draft.cdf
%       (Probably not official; probably ROC ad hoc)
% Summary plots (which BTW are NOT datasets!):
%   solo_L3_rpw-lfr-surv-cwf-e_20200423_V01.png
%   solo_L3_rpw-lfr-surv-swf-e_20200423_V01.png
%
%
% IMPLEMENTATION NOTES
% ====================
% Designed to be easy to
% ** understand which and how filenaming conventions are recognized
% ** edit filenaming conventions
% ** add filenaming conventions
% ** modify the behaviour upon partial recognition, i.e. when to signal
%   ** filename is not a dataset
%   ** filename seems like dataset, but with corrupted filename (not
%      implemented)
% ** make it possible to completely reconstruct the filename from the return
%    values.
%    This may change though due to e.g.
%       -- permitting both file suffix cases,
%       -- permitting own extensions to filenaming conventions.
%
%
% ARGUMENTS
% =========
% filename
%       Filename. To be identified as a dataset filename it must have the
%       form <DATASET_ID; upper-/lowercase>_<arbitrary>.cdf.
%
%
% RETURN VALUES
% =============
% Result
%       If not a recognizable dataset filename, then: []
%       If     a recognizable dataset filename, then:
%           Struct with a varying set of fields, depending on the filenaming
%           convention the filename adheres to.
%       Fields always present:
%       .datasetId            : DATASET_ID. (Always uppercase.)
%       .isCdag               : Logical. Whether or not the file is a CDAG
%                               (DATASET_ID in filename is appended with
%                               "-CDAG"/"-cdag").
%       .dsicdagCase          : String constant describing the case of
%                               DATASET_ID+CDAG: 'upper', 'lower'.
%       .versionStr           : String (not number). Excludes "V".
%       .unoffExtension       : []       : There is no unofficial basename
%                                          extension.
%                             : 1x1 cell : {1} = String. The arbitrary string
%                                          part of unofficial basename extension.
%       .fnDatasetIdCdag      : DATASET_ID + CDAG as present in filename (FN),
%                               including case. In practice meant to be
%                               interpreted as dataset glob.attr.
%                               Logical_source, which should include -CDAG when
%                               present (for now).
%       .dateVec1
%       .dateVec2
%       .timeIntervalFormat   : String constant specifying the time interval
%                               format.
%       Fields sometimes present
%           .timeIntervalStr   : Equals "Datetime" in dataset specifications.
%           + varying fields corresponding to content in filename.
%       NOTE: dateVec* may be either 1x3 or 1x6.
% --
%
%
% VARIABLE NAMING CONVENTION
% ==========================
% DSICDAG = String containing the combination of DATASET_ID (DSI) and
%           (optionally) the suffix "-CDAG".
%
%
% ~BUG
% ====
% filename == "solo_L3_rpw-bia-density_20240101_V01.cdf  solo_L3_rpw-bia-density_20240201_V01.cdf"
% is falsely recognized as a legal dataset filename.
% This is due to the unofficial extension being interpreted as non-empty.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-12-17.
%
function R = parse_dataset_filename(filename)
%
% PROPOSAL: MATLAB package for filenames: solo.adm.dsfn (DataSet FileNames)
% PROPOSAL: Class with static methods for filenames.
%   CON: Too large.
%   CON: Can not define class for argument and return value in it.
%
% PROPOSAL: Return version NUMBER, not string.
%   CON: Harder to adapt to changing versioning scheme.
%       Ex: V2.3.4, V2_3_4
%       CON: Unlikely to happen
%   PRO: Often want to do this anyway.
%       Ex: Find latest version.
%
% PROPOSAL: Change name dateVec --> timeVec
%   PROPOSAL: MATLAB seems to use "date vector", and "time vector" for vector of timestamps.
%
% PROPOSAL: Reorganize to be more similar to python erikpgohansson.so.parse_dataset_filename?
%   NOTE: Partly done.
%   NOTE: python function
%       * does not cover LES, CNES and time-less filename formats.
%       * always returns same "struct" format, including time format.
%       * can not reconstruct filename.
%   PROPOSAL: Just ONE topmost regex_str_parts. Then regex_str_parts on parts.
%       CON: Impossible(?) to reduce number of top-level case due to special
%            cases for LES and CNES test files.
%
% PROPOSAL: Use time interval string to denote time in both
%           parse_dataset_filename() and create_dataset_filename(). No time vectors.
%           Have instead pre-existing separate functions
%               parse_time_interval_str(), and
%               create_time_interval_str()
%           be globally available.
%   PRO: No redundant information (that can contradict) in R.
%       PRO: Does not need to ignore but permit time interval string in
%           create_dataset_filename().
%   PRO: Naturally puts the varying time vectors in separate struct.
%       PRO: R does not change in format as much. Fewer cases.
%   PRO: Can chose to parse time interval string in different ways depending on needs:
%       * Reversible way (varying format)
%       * Non-reversible way that describes nominal time coverage
%         (start-stop; one constant format)
% PROPOSAL: Convert parse_time_interval_str() and create_time_interval_str() (in
%           solo.adm.create_dataset_filename()) to public functions which are
%           inverses of each other (are already?) and can be tested separately.
% PROPOSAL: Separate return struct for time vectors. Time interval string in
%           return value.
%
% PROPOSAL: Abolish dsicdagCase. Should be regarded as part of the
%           respective filenaming conventions.
% PROPOSAL: Abolish unoffExtension.
%   PRO: Causes some filenames to unexpectedly be intepreted as compliant.
%     Ex: "solo_L3_rpw-bia-density_20240101_V01.cdf  solo_L3_rpw-bia-density_20240201_V01.cdf"
%
% PROPOSAL: Separate (globally available) parse/create functions for every
%           separate filenaming scheme.
%   CON: Harder to reuse similarities.
%       Ex: DATASET_ID, version
%       Ex: Unofficial filenaming scheme
%   PRO: Cleaner functions.
%       PRO: Consistent return format.
%           CON: Not true for time vectors which reflect format of time
%                interval string.
%       PRO: Can easily have different requirements for different filenaming conventions
%           Ex: case (always the same for a given convention)
%               ==> Can eliminate dsicdagCase.
%           Ex: isCdag (only for official; not LES and CNES filenaming).
%   PROPOSAL: One class "fn" with shared code
%       Ex: utility functions (used by top-level functions)
%       Ex: regexp constants
%   PROPOSAL: One top-level function that handles all filenaming conventions
%             together. Should assign(?) string for identified convention.
%
% PROPOSAL: Return value for basename without IRFU-internal filenaming
%           extension.
%
% PROPOSAL: Refactor to return class.
%   PROBLEM: solo.adm.parse_dataset_filename_many() adds field "path" to return
%            value and passes it on in its own return value.
%   PROBLEM: Class should work with solo.adm.create_dataset_filename().
%     PROBLEM: How handle fields which
%          (1) are returned from parsing, but are simultaneously
%          (2) redundant when creating filenames.
%       Ex: fnDatasetIdCdag
%     PROBLEM: How handle fields which are only optional in certain combinations.
%       Either
%           {'dateVec1', 'dateVec2', 'lesTestStr'}
%           {'cneTestStr'}
%           {'dateVec'}
%           {'dateVec1', 'dateVec2'}
%
%     PRO: Having two timestamps can not make distinction between a single date
%          and an explicit time interval in filename.
%       Ex: "20240101" vs "20240101T000000-20240102T000000"
%     PROPOSAL: Variable for exact sub-naming convention used for time interval.
%       Ex: DAY, DAY_TO_DAY, SECOND_TO_SECOND
%   PRO: More rigorous.
%   PROBLEM: Must always have the same set of fields. ==> Backwards incompatibility.
%   PROPOSAL: Different functions for different naming conventions. One function
%             for official naming convention (+RPW's "-cdag") which (1) does not
%             return all values, and (2) can be used in most places anyway.
%     CON: Does not help much with backwards incompatibility.
%     CON: Too much overlap between implementations.
%       CON-PROPOSAL: Use parsing function for shared implemention for shared
%                     parts of naming conventions.
%
% PROPOSAL: Replace date vectors with datetime objects (UTC).
%   CON: Currently using the length of date vectors to specify the filename
%        format (time interval format).
%        In particular, with datetime only, create_dataset_filename() would
%        not know which time interval format to use!
%        YYYYMMDD, YYYYMMDD-YYYYMMDD, or YYYYMMDThhmmss-YYYYMMDThhmmss.
%       PROPOSAL: Separate argument for time interval format.
%           PROPOSAL: String constant.
%   NOTE: Would need assertion on hour=minute=second=0 for YYYYMMDD
%         format.
%
% PROPOSAL: Use field names identical to the terms used in specifications (RCS
%           ICD, SOL-SGS-TN-0009).
%   Ex: Datetime, descriptor
%   CON: Terms do not follow variable naming conventions.

NO_MATCH_RETURN_VALUE = [];
UNUSED_DATE_VECTOR    = [0, 0, 0, 0, 0, 0];



% NOTE: Parse from the END.
[~, trueBasename, n] = irf.str.read_token(filename, -1, '\.cdf');
if n == -1
  R = NO_MATCH_RETURN_VALUE;
  return
end



%==========================
% Parse DATA_SET_ID + CDAG
%==========================
% IMPLEMENTATION NOTE: Could use separate "-CDAG" regexp but then have to
% use negative lookahead to exclude "-CDAG" from the preceding DATASET_ID
% (due to maximal munch). Could then use irf.str.read_token()
% again for "-CDAG", but would then have to have separate cases of lowercase
% and uppercase anyway. Might also be slower because of negative lookahead.
%
% NOTE: Lowercase DATASET_ID+CDAG always have uppercase dataset level.
[fnDatasetIdCdag, str, n] = irf.str.read_token(trueBasename, 1, ...
  '(SOLO|ROC-SGSE)_(HK|L1|L1R|L2|L3|CAL)_[A-Z0-2-]*', ...
  '(solo|roc-sgse)_(HK|L1|L1R|L2|L3|CAL)_[a-z0-2-]*');
switch(n)
  case 1
    R.dsicdagCase = 'upper';
    R.isCdag      = strcmp(fnDatasetIdCdag(end-4:end), '-CDAG');
  case 2
    R.dsicdagCase = 'lower';
    R.isCdag      = strcmp(fnDatasetIdCdag(end-4:end), '-cdag');
  otherwise
    R = NO_MATCH_RETURN_VALUE;
    return
end
if R.isCdag
  R.datasetId = upper(fnDatasetIdCdag(1:end-5));
else
  R.datasetId = upper(fnDatasetIdCdag);
end
R.fnDatasetIdCdag = fnDatasetIdCdag;



% Recurring regular expressions.
VERSION_RE     = 'V[0-9][0-9]+';       % NOTE: Permits more than two digits.
LES_TESTSTR_RE = 'les-[0-9a-f]{7,7}';
CNE_TESTSTR_RE = '[0-9a-f]{7,7}_CNE';

% Unofficial extension, arbitrary string part amended to the end of the
% legal basename
% --------------------------------------------------------------------------
% IMPORTANT NOTE: Important to require initial separator string between
% official basename and arbitrary part of extension. Can otherwise not
% rigorously distinguish (1) some ground-test datasets from (2) some
% non-CDAG in-flight datasets.
% IMPORTANT NOTE:
%   Official filenames do not use this separator string.
%   ==> Must permit UNOFF_EXTENSION_RE to match empty string.
%   ==> Can make irf.str.regexp_str_parts mistakenly match
%       UNOFF_EXTENSION_RE to an EMPTY STRING (when there is an unofficial
%       extension).
%   ==> Not "perfect match".
%   ==> Function classifies filename as not dataset filename.
% Therefore, must express UNOFF_EXTENSION_RE such that it will try to match
% a non-empty string FIRST, and an empty string SECOND. THE ORDER MATTERS.
UNOFF_EXTENSION_RE = '(\..*|)';



%===========================================================================
% Different types of supported filenaming conventions
% ---------------------------------------------------
% solo_L2_rpw-lfr-surv-cwf-e-cdag_20200213_V01.cdf
% solo_L1_rpw-bia-current-cdag_20200401T000000-20200421T000000_V01.cdf
% solo_L1_rpw-bia-current-cdag_20200301-20200331_V01.cdf
% SOLO_L2_RPW-LFR-SURV-CWF-E_V04.cdf
% ROC-SGSE_HK_RPW-BIA_19850de_CNE_V02.cdf
% solo_L1_rpw-tds-lfm-rswf_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
%
% NOTE: Checks for different naming conventions in assumed order of
% decreasing likelyhood of being used.
% ==> Potential speedup.
%===========================================================================

TIME_INTERVAL_STR_RE = '[0-9T-]{8,31}';

%=================================================================
% Standard filenaming convention (incl. CDAG; tested for earlier)
%=================================================================
[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
  {'_', TIME_INTERVAL_STR_RE, '_', VERSION_RE, UNOFF_EXTENSION_RE}, ...
  'permit non-match');
if perfectMatch
  R                 = parse_time_interval_str(R, subStrCa{2});
  R.timeIntervalStr = subStrCa{2};
  R.versionStr      = ver_2_versionStr(subStrCa{4});
  R.unoffExtension  = unoff_extension_RE_to_str(subStrCa{5});
  return
end

%=============================
% "LES" filenaming convention
%=============================
[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
  {'_', TIME_INTERVAL_STR_RE, '_', VERSION_RE, '_', ...
  LES_TESTSTR_RE, UNOFF_EXTENSION_RE}, 'permit non-match');
if perfectMatch
  R                 = parse_time_interval_str(R, subStrCa{2});
  R.timeIntervalStr = subStrCa{2};
  R.versionStr      = ver_2_versionStr(         subStrCa{4});
  R.lesTestStr      =                           subStrCa{6};
  R.unoffExtension  = unoff_extension_RE_to_str(subStrCa{7});
  return
end

%===============================
% "CNES" filenameing convention
%===============================
[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
  {'_', CNE_TESTSTR_RE, '_', VERSION_RE, UNOFF_EXTENSION_RE}, ...
  'permit non-match');
if perfectMatch
  R.dateVec1           = UNUSED_DATE_VECTOR;
  R.dateVec2           = UNUSED_DATE_VECTOR;
  R.timeIntervalFormat = 'NO_TIME_INTERVAL';
  %
  R.cneTestStr         =                           subStrCa{2};
  R.versionStr         = ver_2_versionStr(         subStrCa{4});
  R.unoffExtension     = unoff_extension_RE_to_str(subStrCa{5});
  return
end

% CASE: Could not match filename to anything.

R = NO_MATCH_RETURN_VALUE;
end    % parse_dataset_filename()



function versionStr = ver_2_versionStr(s)
versionStr = s(2:end);   % NOTE: No conversion to number.
end



% "unoff_extension_RE" = The part of string that matches regular expression
% (RE).
function unoffExtension = unoff_extension_RE_to_str(s)
if isempty(s)
  unoffExtension = [];
else
  unoffExtension = {s(2:end)};
end
end



% Returns different fields and formats depending on format of time interval
% string.
%
% NOTE: Adds on to existing struct to avoid having to use
% irf.ds.add_struct_to_struct() (slow?).
%
function R = parse_time_interval_str(R, s)
% yyyymmdd (8 digits).
DATE_RE     = '20[0-9][0-9][01][0-9][0-3][0-9]';

% NOTE: DATETIME_RE not same as glob.attr. Datetime, but component of.
% yyyymmddThhmmss (8+1+6=15 digits/T)
DATETIME_RE = '20[0-9]{6,6}T[0-9]{6,6}';



[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(s, ...
  {DATE_RE}, 'permit non-match');
if perfectMatch
  R.dateVec1 = date_str_2_dateVec(subStrCa{1});
  R.dateVec2 = datevec(datetime(R.dateVec1) + caldays(1));
  R.timeIntervalFormat = 'DAY';
  return
end

[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(s, ...
  {DATE_RE, '-', DATE_RE}, 'permit non-match');
if perfectMatch
  R.dateVec1 = date_str_2_dateVec(subStrCa{1});
  R.dateVec2 = date_str_2_dateVec(subStrCa{3});
  R.dateVec2 = datevec(datetime(R.dateVec2) + caldays(1));
  R.timeIntervalFormat = 'DAY_TO_DAY';
  return
end

[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(s, ...
  {DATETIME_RE, '-', DATETIME_RE}, 'permit non-match');
if perfectMatch
  R.dateVec1 = date_time_str_2_dateVec(subStrCa{1});
  R.dateVec2 = date_time_str_2_dateVec(subStrCa{3});
  R.timeIntervalFormat = 'SECOND_TO_SECOND';
  return
end

error('Can not interpret time interval string "%s".', s)
end



% Utility function
function dateVec = date_time_str_2_dateVec(s)
dateVec = str2double({s(1:4), s(5:6), s(7:8), s(10:11), s(12:13), s(14:15)});

% NOTE: Is not a check on filename, but on implementation. read_token()
% should guarantee that strings can be parsed as numbers.
%assert(~any(isnan(dateVec)))
end



% Utility function
function dateVec = date_str_2_dateVec(s)
dateVec = str2double({s(1:4), s(5:6), s(7:8)});
dateVec(4:6) = [0, 0, 0];

% NOTE: Is not a check on filename, but on implementation. read_token should
% guarantee that strings can be parsed as numbers.
%assert(~any(isnan(dateVec)))
end

