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
% * Will not (yet) work for dataset filenames which uppercase outside of
%   the archiving level.
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
%       solo.adm.dsfn.parse_dataset_filename().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-24.
%
function filename = create_dataset_filename(R)
% NOTE: See comments for solo.adm.dsfn.parse_dataset_filename().

% ASSERTION
% NOTE: Useful, so that later code does not need to be as rigorous when
% determining case.
irf.assert.struct(R, ...
  { ...
    'isCdag', 'datasetId', 'versionNbr' ...
    'dateVec1', 'dateVec2', 'timeIntervalFormat' ...
  }, ...
  { ...
    'cneTestStr', 'lesTestStr' ...
  })



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
assert(isnumeric(R.versionNbr));
versionNbr = R.versionNbr;
R = rmfield(R, {...
  'datasetId', 'versionNbr', 'isCdag', ...
  'dateVec1', 'dateVec2', 'timeIntervalFormat' ...
 });



fnCa = fieldnames(R);

timeIntervalStr = solo.adm.dsfn.create_time_interval_str(dateVec1, dateVec2, timeIntervalFormat);



if sets_equal(fnCa, {'lesTestStr'}) && ismember(timeIntervalFormat, {'DAY_TO_DAY', 'SECOND_TO_SECOND'})
  %=============================
  % "LES" filenaming convention
  %=============================
  assert(~isCdag)
  dsicdagStr = get_cased_DSICDAG(dsi, false, isCdag);
  filename = sprintf('%s_%s_V%02i_%s.cdf', ...
    dsicdagStr, timeIntervalStr, versionNbr, R.lesTestStr);

elseif sets_equal(fnCa, {'cneTestStr'}) && strcmp(timeIntervalFormat, {'NO_TIME_INTERVAL'})
  %===============================
  % "CNES" filenameing convention
  %===============================
  % Ex: ROC-SGSE_HK_RPW-BIA_19850de_CNE_V02.cdf
  assert(~isCdag)
  dsicdagStr = get_cased_DSICDAG(dsi, true, isCdag);
  filename = sprintf('%s_%s_V%02i.cdf', ...
    dsicdagStr, R.cneTestStr, versionNbr);

elseif sets_equal(fnCa, {}) && ismember(timeIntervalFormat, {'DAY', 'DAY_TO_DAY', 'SECOND_TO_SECOND'})
  %===========================================
  % In-space filename convention (incl. CDAG)
  %===========================================
  dsicdagStr = get_cased_DSICDAG(dsi, false, isCdag);
  filename = sprintf('%s_%s_V%02i.cdf', ...
    dsicdagStr, timeIntervalStr, versionNbr);

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
% NOTE: There does not appear to be a generic function for this.
equal = isempty(setxor(set1, set2));
end
