%
% Function to create OFFICIAL basename (no suffix) for summary plot according to
% the official filenaming convention.
%
%
% ARGUMENTS
% =========
% versionNbr
%       Should be dataset version number (not plot file version).
%
%
% Xavier Bonnin e-mail 2020-04-10
% ===============================
% """"""""
% We should stay as much as possible compliant with the definition in the "Solar
% Orbiter data definition document" (see SOL-SGS-TN-0009-MetadataStandard-2.4
% enclosed).
%
% In this case SP file should be labelled as level 3 data products (L3) and
% should have the form:
%
% 	solo_L3_<descriptor>_<datetime>_V<version>_<Free_field>.<extension>
%
% Where:
% 	<descriptor> should be the same than the parent L1/L2 CDF used to generate
% 	the SP file (e.g. "rpw-tnr-surv", "rpw-tds-surv-rswf-e", etc.),
% 	<datetime> is the time of file (can be a daily or time range file, see
% 	SOL-SGS-TN-0009-MetadataStandard-2.4 for convention),
% 	 <version> is the version of the file (same convention then for CDF,
% 	 starting at 01 and iterating at each new version),
% 	<Free_field> this field is not used for other SolO RPW data products for
% 	now, but it can be used here to specify that it is a summary plots
% 	(otherwise it can be left empty and the file format allows to distinguish
% 	between CDF and image files)
% 	<extension> the format of the plot file should be as much as possible
% 	homogenous over Solar Orbiter data. Should be "PNG" or "JPEG2000" as
% 	proposed by SOC.
% """""""" /Xavier Bonnin e-mail 2020-04-10
% """"
% The filename (including the version "VXX") is the version of the input L2 file
% used to generate the plots.
% """" /Xavier Bonnin e-mail 2020-09-01
% --
% NOTE: Filenaming convention description only allows dates as timestamps
% (yyyymmdd)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-08.
%
function basename = create_SP_basename(srcDatasetId, dateVec3, versionNbr)
% Ex: solo_L3_rpw-tnr-surv_20200315_V01.png
%
% PROPOSAL: Create corresponding parse_summary_plot_filename().
%   PRO: May be useful for algorithms for selecting datasets when batch
%        processing.
%   PROBLEM: Future SPs may produce multiple SPs per dataset. ==> Unclear
%            return format, complexity.
% PROPOSAL: Use (future) assertion function on DATASET_ID.

% ASSERTIONS
irf.assert.castring(srcDatasetId)
assert(isnumeric(dateVec3))
assert(numel(dateVec3) == 3)
assert(isscalar(versionNbr))
assert(isnumeric(versionNbr))
assert(versionNbr >= 1)

% NOTE: Output should have uppercase "L3", despite most of basename being
% lower case.
modifDatasetId = regexprep(lower(srcDatasetId), '^solo_(hk|l2)_', 'solo_L3_');

% ASSERTION: DATASET_ID was modified.
assert(~strcmp(srcDatasetId, modifDatasetId), ...
  'Can not handle datasetId="%s".', srcDatasetId)

% NOTE: Uppercase "V" (version) despite most of basename being lower case.
basename = sprintf('%s_%04g%02g%02g_V%02g', ...
  modifDatasetId, dateVec3(:), versionNbr);
end
