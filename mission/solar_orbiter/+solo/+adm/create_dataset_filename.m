%
% Create dataset filename according to filenaming convention.
%
%
% NOTES
% =====
% * Function is the inverse of solo.adm.parse_dataset_filename(). See that
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
% useful to have together with solo.adm.parse_dataset_filename,
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
%       solo.adm.parse_dataset_filename(), except timeIntervalStr
%       which is optional.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-24.
%
function filename = create_dataset_filename(R)
% NOTE: See comments for solo.adm.parse_dataset_filename().
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
  {'isCdag', 'datasetId', 'dsicdagCase', 'versionStr', 'unoffExtension'}, ...
  {'fnDatasetIdCdag', 'timeIntervalStr', ...
  'dateVec', 'dateVec1', 'dateVec2', ...
  'cneTestStr', 'lesTestStr'})



if R.isCdag
  cdagStr = '-CDAG';
else
  cdagStr = '';
end

switch(R.dsicdagCase)
  case 'upper'
    cdagStr   = upper(cdagStr);
    datasetId = upper(R.datasetId());
  case 'lower'
    cdagStr   = lower(cdagStr);
    datasetId = lowercase_DATASET_ID(R.datasetId());
  otherwise
    error('Illegal R.dsicdagCase.')
end


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



% Remove all mandatory field names
% --------------------------------
% So that they do not need to be checked for when identifying filenaming
% convention from set of fields.
assert(ischar(R.versionStr));
versionStr = R.versionStr;
R = rmfield(R, {'datasetId', 'versionStr', 'unoffExtension', 'isCdag', 'dsicdagCase'});

% Remove some optional field names
% --------------------------------
% "timeIntervalStr" is not used by this function. Only permitted for
% compatibility with solo.adm.parse_dataset_filename().
if isfield(R, 'timeIntervalStr')
  R = rmfield(R, 'timeIntervalStr');
end
if isfield(R, 'fnDatasetIdCdag')
  R = rmfield(R, 'fnDatasetIdCdag');
end



fnCa = fieldnames(R);

% NOTE: Ignores R.timeIntervalStr (if present) but derives it instead.
timeIntervalStr = create_time_interval_str(R);



if sets_equal(fnCa, {'dateVec1', 'dateVec2', 'lesTestStr'})
  filename = sprintf('%s%s_%s_V%02s_%s%s.cdf', ...
    datasetId, cdagStr, timeIntervalStr, versionStr, R.lesTestStr, unoffExtension);

elseif sets_equal(fnCa, {'cneTestStr'})
  % ROC-SGSE_HK_RPW-BIA_19850de_CNE_V02.cdf
  filename = sprintf('%s%s_%s_V%02s%s.cdf', ...
    datasetId, cdagStr, R.cneTestStr, versionStr, unoffExtension);

elseif sets_equal(fnCa, {'dateVec'}) || sets_equal(fnCa, {'dateVec1', 'dateVec2'})
  filename = sprintf('%s%s_%s_V%02s%s.cdf', ...
    datasetId, cdagStr, timeIntervalStr, versionStr, unoffExtension);

  %     elseif sets_equal(fnCa, {})
  %
  %         % SOLO_L2_RPW-LFR-SURV-CWF-E_V04.cdf
  %         % NOTE: No timeIntervalStr.
  %         filename = sprintf('%s%s_V%02s%s.cdf', ...
  %         datasetId, cdagStr, versionStr, unoffExtension);
else

  error('Illegal argument R.')
end

end



function s = create_time_interval_str(R)
  function b = is_FV(fnName, nComp)
    b = isfield(R, fnName) && (numel(R.(fnName)) == nComp);
  end

dv  = isfield(R, 'dateVec');
dv1 = isfield(R, 'dateVec1');
dv2 = isfield(R, 'dateVec2');

if     is_FV('dateVec', 3) && ~dv1 && ~dv2
  s = sprintf(...
    '%04i%02i%02i', ...
    R.dateVec);

elseif ~dv && is_FV('dateVec1', 3) && is_FV('dateVec2', 3)
  s = sprintf(...
    '%04i%02i%02i-%04i%02i%02i', ...
    R.dateVec1(:), R.dateVec2(:));

elseif ~dv && is_FV('dateVec1', 6) && is_FV('dateVec2', 6)
  s = sprintf(...
    '%04i%02i%02iT%02i%02i%02i-%04i%02i%02iT%02i%02i%02i', ...
    R.dateVec1(:), R.dateVec2(:));

elseif ~dv && ~dv1 && ~dv2
  s = [];

else
  error('Can not interpret combination of date & time field values.')
end
end



function lcDatasetId = lowercase_DATASET_ID(datasetId)

[sourceName, level, descriptor] = solo.adm.disassemble_DATASET_ID(datasetId);

lcDatasetId = sprintf('%s_%s_%s', lower(sourceName), level, lower(descriptor));
end



function equal = sets_equal(set1, set2)
equal = isempty(setxor(set1, set2));
end
