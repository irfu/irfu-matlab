%
% Create dataset filename according to filenaming convention.
%
%
% NOTES
% =====
% * Function is the inverse of solo.adm.dsfn.parse_dataset_filename(). See that
%   function for different dataset filename formats.
% * Function uses the exact set of fields to determine which naming convention
%   to follow, not just the field values.
% * Will not (yet) work for dataset filenames which uppercase outside of archiving level
%   Ex: solo_L1_swa-eas2-NM3D_20201027T000007-20201027T030817_V01.cdf
%
%
% RATIONALE
% =========
% Useful for creating filenames to be output from BICAS during IRF-internal
% testing, so that they can be analogous with the input datasets. The
% functionality is in a way overkill, but the function is easy to write and
% useful to have together with solo.adm.dsfn.parse_dataset_filename,
% since having both, and having them be inverses of each other make both of them
% easy to automatically test, and make it easy to test that they are compatible
% with each other.
%
%
% ARGUMENTS
% =========
% R
%       Struct with varying set of fieldnames, depending on which filenaming
%       convention to follow. Identical to the one in
%       solo.adm.dsfn.parse_dataset_filename(), except timeIntervalStr
%       which is optional.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-24.
%
function filename = create_dataset_filename(R)
% NOTE: See comments for solo.adm.dsfn.parse_dataset_filename().
%
% PROPOSAL: Ability to add IRF-internal basename suffix.
%   Ex: solo_L2_rpw-tds-lfm-rswf-e_20200410_V01.suffix1.cdf
%   PROPOSAL: Ability to add sequence of such suffixes.
%       Ex: solo_L2_rpw-tds-lfm-rswf-e_20200410_V01.suffix1.suffix2.cdf
%   PROPOSAL: <basename><extension>.cdf
%       PRO: Covers all reasonable extensions.
%       CON: Not rigorously parsable if one official basename is an initial substring of another official
%       basename.
%           NOTE: Version string always last except for old ground-testing datasets, except
%           (1) some non-CDAG inflight datasets which filename format is an initial substring of
%           (2) some ground-test dataset filenaming format.
%               Ex: solo_L1_rpw-tds-lfm-rswf_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
%               Ex: solo_L1_rpw-bia-current_20200401T000000-20200421T000000_V01.cdf
%               NOTE: Unlikely to matter in practice, since using different DATASET_IDs.
%   --
%   PROPOSAL: Use special character(s) to always separate extension from official basename.
%       PROPOSAL: Use character/string never used by ROC.
%           PROPOSAL: <basename>.<extension>.cdf
%               CON: Period could cause trouble on some platforms.
%                   CON: Seems OK on Mac, Linux.
%           PROPOSAL: <basename>___<extension>.cdf
%               CON: Waste of space.
%       PROPOSAL: <basename>_<extension>.cdf
%   --
%   TODO-DEC: How store unofficial extension?
%       NOTE: Want to distinguish absence of extension from empty string represented in extension.
%       PROPOSAL: [] = No extension
%                 1x1 cell with string = Extension.
%       PROPOSAL: 0x1 cell = No extension
%                 1x1 cell with string = Extension.
%           NOTE: Why that empty size?

% ASSERTION
% NOTE: Useful, so that later code does not need to be as rigorous when
% determining case.
irf.assert.struct(R, ...
  { ...
    'isCdag', 'datasetId', 'versionStr', 'unoffExtension' ...
    'dateVec1', 'dateVec2', 'timeIntervalFormat' ...
  }, ...
  { ...
    'fnDatasetIdCdag', 'timeIntervalStr', ...
    'cneTestStr', 'lesTestStr' ...
  })



% Convert R.unoffExtension to actual string to add to basename
% ------------------------------------------------------------
% IMPLEMENTATION NOTE: Thorougly check R.unoffExtension format since it
% often gets wrong.
if isempty(R.unoffExtension) && ~iscell(R.unoffExtension)
  unoffExtension = '';
elseif isscalar(R.unoffExtension) && iscell(R.unoffExtension) && ischar(R.unoffExtension{1})
  unoffExtension = ['.', R.unoffExtension{1}];
else
  error('Illegal R.unoffExtension. Must be (1) empty non-cell or (2) 1x1 cell array of string.')
end


assert(isrow(R.dateVec1) & (numel(R.dateVec1) == 6))
assert(isrow(R.dateVec2) & (numel(R.dateVec2) == 6))
irf.assert.castring(R.timeIntervalFormat)
assert(islogical(R.isCdag) & isscalar(R.isCdag))
dateVec1           = R.dateVec1;
dateVec2           = R.dateVec2;
timeIntervalFormat = R.timeIntervalFormat;
dsi                = R.datasetId;
isCdag             = R.isCdag;



% Remove all mandatory field names
% --------------------------------
% So that they do not need to be checked for when identifying filenaming
% convention from set of fields.
assert(ischar(R.versionStr));
versionStr = R.versionStr;
R = rmfield(R, {...
  'datasetId', 'versionStr', 'unoffExtension', 'isCdag', ...
  'dateVec1', 'dateVec2', 'timeIntervalFormat' ...
 });

% Remove some optional field names
% --------------------------------
% "timeIntervalStr" is not used by this function. Only permitted for
% compatibility with solo.adm.dsfn.parse_dataset_filename().
if isfield(R, 'timeIntervalStr')
  R = rmfield(R, 'timeIntervalStr');
end
if isfield(R, 'fnDatasetIdCdag')
  R = rmfield(R, 'fnDatasetIdCdag');
end



fnCa = fieldnames(R);

% NOTE: Ignores R.timeIntervalStr (if present) but derives it instead.
timeIntervalStr = solo.adm.dsfn.create_time_interval_str(dateVec1, dateVec2, timeIntervalFormat);



if sets_equal(fnCa, {'lesTestStr'}) && ismember(timeIntervalFormat, {'DAY_TO_DAY', 'SECOND_TO_SECOND'})
  %=============================
  % "LES" filenaming convention
  %=============================
  dsicdagStr = get_cased_DSICDAG(dsi, false, isCdag);
  filename = sprintf('%s_%s_V%02s_%s%s.cdf', ...
    dsicdagStr, timeIntervalStr, versionStr, R.lesTestStr, unoffExtension);

elseif sets_equal(fnCa, {'cneTestStr'}) && strcmp(timeIntervalFormat, {'NO_TIME_INTERVAL'})
  %===============================
  % "CNES" filenameing convention
  %===============================
  % Ex: ROC-SGSE_HK_RPW-BIA_19850de_CNE_V02.cdf
  dsicdagStr = get_cased_DSICDAG(dsi, true, isCdag);
  filename = sprintf('%s_%s_V%02s%s.cdf', ...
    dsicdagStr, R.cneTestStr, versionStr, unoffExtension);

elseif sets_equal(fnCa, {}) && ismember(timeIntervalFormat, {'DAY', 'DAY_TO_DAY', 'SECOND_TO_SECOND'})
  %===========================================
  % In-space filename convention (incl. CDAG)
  %===========================================
  dsicdagStr = get_cased_DSICDAG(dsi, false, isCdag);
  filename = sprintf('%s_%s_V%02s%s.cdf', ...
    dsicdagStr, timeIntervalStr, versionStr, unoffExtension);

else

  error('Illegal argument R.')
end

end    % create_dataset_filename()



% dsicdagStr
%       DSI+optional CDAG, but with the case to be used in the filename.
function [dsicdagStr] = get_cased_DSICDAG(dsi, dsicdagUppercase, isCdag)

if isCdag
  cdagStr = '-CDAG';
else
  cdagStr = '';
end

if dsicdagUppercase
    cdagStr = upper(cdagStr);
    dsiStr  = upper(dsi);
else
    cdagStr = lower(cdagStr);
    dsiStr  = lowercase_DSI(dsi);
end

dsicdagStr = [dsiStr, cdagStr];

end



function lowercaseDsi = lowercase_DSI(dsi)

[sourceName, level, descriptor] = solo.adm.disassemble_DATASET_ID(dsi);

% NOTE: The "level" is always uppercase.
lowercaseDsi = sprintf('%s_%s_%s', lower(sourceName), level, lower(descriptor));

end



function equal = sets_equal(set1, set2)
equal = isempty(setxor(set1, set2));
end
