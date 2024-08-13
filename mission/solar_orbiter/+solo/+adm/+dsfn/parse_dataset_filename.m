%
% Given a filename, parse it as a SolO dataset filename if possible.
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
%       solo.adm.dsfn.create_dataset_filename().
% NOTE/BUG: Can not handle e.g.
%       solo_L1_swa-eas2-NM3D_20201027T000007-20201027T030817_V01.cdf
%       since it has mixed case in the descriptor (not level).
% NOTE: RCTs counts as datasets in this context since the conform to the same
%       filenaming convention and are described in SOL-SGS-TN-0009, "Metadata
%       Definition for Solar Orbiter Science Data", 01/06, though using its own
%       level "CAL".
%
%
% IMPLEMENTATION NOTES
% ====================
% Designed to be easy to:
% * understand which and how filenaming conventions are recognized
% * edit filenaming conventions
% * add filenaming conventions
% * modify the behaviour upon partial recognition, i.e. when to signal
%   * filename is not a dataset
%   * filename seems like dataset, but with corrupted filename (not
%      implemented)
% * make it possible to completely reconstruct the filename from the return
%   values.
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
% R
%       If not a recognizable dataset filename, then: []
%       If     a recognizable dataset filename, then:
%       Struct with a varying set of fields, depending on the filenaming
%       convention the filename adheres to.
%       Fields always present:
%       .datasetId
%           DATASET_ID. (Always uppercase.)
%       .isCdag
%           Logical. Whether or not the file is a CDAG
%           (DATASET_ID in filename is appended with "-CDAG"/"-cdag").
%       .versionNbr
%           Version number.
%       .dateVec1
%       .dateVec2
%       .timeIntervalFormat
%           String constant specifying the time interval format.
%       Fields sometimes present
%           Varying fields corresponding to content in filename.
% S
%       If not a recognizable dataset filename, then: []
%       If     a recognizable dataset filename, then:
%       Struct with information which is technically redundant when combined
%       with return value "R", but which can be useful to have for certain
%       applications.
%       .timeIntervalStr
%       .filenameDsiCdag
%             The combination of DSI and optionally "-cdag" as found in
%             the filename, including case.
%             In practice meant to be interpreted as dataset glob.attr.
%             "Logical_source", which should include -CDAG when present (for
%             now at least).
%       .versionStr
%           Version number as a string the way it is represented in
%           the filename. Excludes "V".
%
%
% VARIABLE NAMING CONVENTION
% ==========================
% DSICDAG = String containing the combination of DATASET_ID (DSI) and
%           (optionally) the suffix "-CDAG".
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-12-17.
%
function [R, S] = parse_dataset_filename(filename)
%
% PROPOSAL: Change name dateVec --> timeVec
%   PROPOSAL: MATLAB seems to use "date vector", and "time vector" for vector of timestamps.
%
% PROPOSAL: Replace date vectors with datetime objects (UTC).
%   CON: CNES and LES filenaming *might* not have used UTC.
%        ==> Better to use date vectors.
%     CON-PROPOSAL: Use datetime without UTCLeapSeconds.
%
% PROPOSAL: Reorganize to be more similar to python erikpgohansson.solo.utils.parse_dataset_filename()?
%   NOTE: Partly done.
%   NOTE: Python function (not MATLAB).
%       * Does not cover LES, CNES and time-less filename formats.
%       * Always returns same "struct" format, including time format.
%       * Can not reconstruct filename.
%   PROPOSAL: Just ONE topmost regex_str_parts. Then regex_str_parts on parts.
%       CON: Impossible(?) to reduce number of top-level case due to special
%            cases for LES and CNES test files.
%
% PROPOSAL: Use time interval string to denote time in both
%           parse_dataset_filename() and create_dataset_filename(). No timestamps.
%           Caller should instead use pre-existing separate functions
%               solo.adm.dsfn.parse_time_interval_str()
%               solo.adm.dsfn.create_time_interval_str()
%   PRO: No redundant (time) information (that can contradict) in R.
%       PRO: Does not need to ignore but permit time interval string in
%           create_dataset_filename().
%   CON: Turns the filename parsing into a two stage process. ==> Bad for
%        quickly identifying filenames as datasets/non-datasets.
% PROPOSAL: Separate return struct for time vectors. Time interval string in
%           return value.
%
% PROPOSAL: Separate (globally available) parse/create functions for every
%           separate filenaming scheme.
%   CON: Harder to reuse similarities.
%       Ex: DATASET_ID, version, time interval string
%       CON-PROPOSAL: Use shared RE constants, functions.
%   CON: Requires one return value class (future) for every filenaming convention.
%     CON: Requires one pair of parse/create fuctions per filenaming convention already.
%   CON: Can not build composite function for simultaneously handling all
%        filenaming conventions (as the current create/parse functions do).
%     CON: Filenaming conventions are too dissimilar for doing that anyway: Have
%          different variables.
%       Ex: "CNE" filenaming convention has hash, has no time interval string.
%       Ex: "LES" filenaming convention has hash.
%       Ex: Inflight filenaming has isCdag (has no hash).
%   CON: solo.adm.dsfn.parse_dataset_filename_many() would only apply to one
%        filenaming convention.
%     CON-PROPOSAL: Can call functions for multiple filenaming conventions.
%       PRO: Functions would be made for making this easy.
%     NOTE: solo.adm.dsfn.parse_dataset_filename_many() is only used by
%           solo.adm.paths_to_DSMD_array().
%     CON-PROPOSAL: Merge solo.adm.dsfn.parse_dataset_filename_many() into
%                   solo.adm.paths_to_DSMD_array().
%       CON: No need for a replacement of return value + ".path" field. DSMD is
%            the "replacement".
%   PRO: Cleaner functions.
%       PRO: Consistent return format.
%           CON: Not true for time vectors which reflect format of time
%                interval string.
%       PRO: Can easily have different requirements for different filenaming
%            conventions.
%           Ex: isCdag (only for official; not LES and CNES filenaming).
%   PROPOSAL: One MATLAB package per filenaming convention.
%   PROPOSAL: One class "fn" (filenames) with shared code
%       Ex: utility functions (used by top-level functions)
%       Ex: regexp constants
%   PROPOSAL: One top-level function that handles all filenaming conventions
%             together. Should assign(?) string for identified convention.
%
% PROPOSAL: Refactor to return class.
%   PROBLEM: solo.adm.dsfn.parse_dataset_filename_many() adds field ".path" to
%            return value and passes it on in its own return value.
%   PROBLEM: Class should work with solo.adm.dsfn.create_dataset_filename().
%     PROBLEM: How handle fields which
%          (1) are returned from parsing, but are simultaneously
%          (2) redundant when creating filenames.
%       Ex: filenameDsiCdag
%
%   PRO: More rigorous.
%   PROBLEM: Must always have the same set of fields. ==> Backwards incompatibility.
%   PROPOSAL: Different functions for different naming conventions. One function
%             for official naming convention (+RPW's "-cdag") which (1) does not
%             return all values, and (2) can be used in most places anyway.
%     CON: Does not help much with backwards incompatibility.
%     CON: Too much overlap between implementations.
%       CON-PROPOSAL: Use parsing function for shared implemention for shared
%                     parts of naming conventions.
%   PROPOSAL: Return value "S" converted into methods which derive
%             the corresponding values.
%     ~CON: If the class contains non-redundant information, then it has to
%           re-derive timeIntervalStr, filenameDsiCdag, (future versionStr).
%
% PROPOSAL: Return time interval string separately from struct. -- IMPLEMENTED
%   PRO: Used by bicas.ga.derive_output_dataset_GAs().
% PROPOSAL: Return file basename separately from struct.
%   PRO: Used by bicas.ga.derive_output_dataset_GAs().
% PROPOSAL: Return filenameDsiCdag separately from struct. -- IMPLEMENTED
%   PRO: Used by bicas.ga.derive_output_dataset_GAs(), solo.db_list_files___UTEST.
%
% PROPOSAL: Use field names identical to the terms used in specifications (RCS
%           ICD, SOL-SGS-TN-0009).
%   Ex: Datetime, descriptor
%   CON: Terms do not follow own variable naming conventions.
%
% PROPOSAL: datasetId --> DSID
%
% PROPOSAL: Forbid longer version strings beginning with zero.
% PROPOSAL: Forbid uppercase "-CDAG". Should never happen.
% PROPOSAL: Only permit "-cdag"" for inflight datasets.

NO_MATCH_RETURN_VALUE = [];
UNUSED_DATE_VECTOR    = [0, 0, 0, 0, 0, 0];



% NOTE: Parse from the END.
[~, trueBasename, n] = irf.str.read_token(filename, -1, '\.cdf');
if n == -1
  R = NO_MATCH_RETURN_VALUE;
  S = NO_MATCH_RETURN_VALUE;
  return
end



%=========================
% Parse DATASET_ID + CDAG
%=========================
% IMPLEMENTATION NOTE: Could use separate "-CDAG" regexp but then have to
% use negative lookahead to exclude "-CDAG" from the preceding DATASET_ID
% (due to maximal munch). Could then use irf.str.read_token()
% again for "-CDAG", but would then have to have separate cases of lowercase
% and uppercase anyway. Might also be slower because of negative lookahead.
%
% NOTE: Lowercase DATASET_ID+CDAG always have uppercase dataset level.
[filenameDsicdag, str, n] = irf.str.read_token(trueBasename, 1, ...
  '(SOLO|ROC-SGSE)_(HK|L1|L1R|L2|L3|CAL)_[A-Z0-2-]*', ...
  '(solo|roc-sgse)_(HK|L1|L1R|L2|L3|CAL)_[a-z0-2-]*');
switch(n)
  case 1
    dsicdagUppercase = true;
    R.isCdag         = strcmp(filenameDsicdag(end-4:end), '-CDAG');
  case 2
    dsicdagUppercase = false;
    R.isCdag         = strcmp(filenameDsicdag(end-4:end), '-cdag');
  otherwise
    R = NO_MATCH_RETURN_VALUE;
    S = NO_MATCH_RETURN_VALUE;
    return
end
if R.isCdag
  R.datasetId = upper(filenameDsicdag(1:end-5));
else
  R.datasetId = upper(filenameDsicdag);
end
S.filenameDsiCdag = filenameDsicdag;



% Recurring regular expressions.
VERSION_RE     = 'V[0-9][0-9]+';       % NOTE: Permits more than two digits.
LES_TESTSTR_RE = 'les-[0-9a-f]{7,7}';
CNE_TESTSTR_RE = '[0-9a-f]{7,7}_CNE';



%===========================================================================
% NOTE: Checks for different naming conventions in assumed order of
% decreasing likelyhood of being used.
% ==> Potential speedup.
%===========================================================================

TIME_INTERVAL_STR_RE = '[0-9T-]{8,31}';

%===========================================================================
% Standard in-flight filenaming convention (incl. CDAG; tested for earlier)
%===========================================================================
[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
  {'_', TIME_INTERVAL_STR_RE, '_', VERSION_RE}, ...
  'permit non-match');
if perfectMatch & ~dsicdagUppercase
  S.timeIntervalStr = subStrCa{2};

  [R.dateVec1, R.dateVec2, R.timeIntervalFormat] = ...
    solo.adm.dsfn.parse_time_interval_str(S.timeIntervalStr);
  S.versionStr      = version_RE_match_to_versionStr(subStrCa{4});
  R.versionNbr      = str2double(S.versionStr);
  return
end

%=============================
% "LES" filenaming convention
%=============================
[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
  {'_', TIME_INTERVAL_STR_RE, '_', VERSION_RE, '_', LES_TESTSTR_RE}, ...
  'permit non-match');
if perfectMatch & ~dsicdagUppercase
  assert(~R.isCdag)

  S.timeIntervalStr = subStrCa{2};

  [R.dateVec1, R.dateVec2, R.timeIntervalFormat] = ...
    solo.adm.dsfn.parse_time_interval_str(S.timeIntervalStr);
  S.versionStr      = version_RE_match_to_versionStr(subStrCa{4});
  R.versionNbr      = str2double(S.versionStr);
  R.lesTestStr      = subStrCa{6};
  return
end

%=============================
% "CNE" filenaming convention
%=============================
[subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
  {'_', CNE_TESTSTR_RE, '_', VERSION_RE}, ...
  'permit non-match');
if perfectMatch & dsicdagUppercase
  assert(~R.isCdag)

  S.timeIntervalStr    = NO_MATCH_RETURN_VALUE;

  R.dateVec1           = UNUSED_DATE_VECTOR;
  R.dateVec2           = UNUSED_DATE_VECTOR;
  R.timeIntervalFormat = 'NO_TIME_INTERVAL';
  %
  R.cneTestStr         = subStrCa{2};
  S.versionStr         = version_RE_match_to_versionStr(subStrCa{4});
  R.versionNbr         = str2double(S.versionStr);
  return
end

%============================================
% CASE: Could not match filename to anything
%============================================
R = NO_MATCH_RETURN_VALUE;
S = NO_MATCH_RETURN_VALUE;
end    % parse_dataset_filename()



function versionStr = version_RE_match_to_versionStr(s)
versionStr = s(2:end);   % NOTE: No conversion to number.
end
