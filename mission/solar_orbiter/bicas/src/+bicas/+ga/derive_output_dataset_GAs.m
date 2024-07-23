%
% Function for dynamically deriving GAs for a specific output
% dataset given the GAs of multiple input datasets.
%
% NOTE: Some of the GA values determined here are
%   (1) unique for this particular output dataset,
%   (2) common for all output datasets for the current s/w mode,
%   (3) common for alla output datasets.
%
%
% SOOP_TYPE, Datetime, OBS_ID
% ===========================
% XB on RCS telecon 2020-09-17: SOOP_TYPE, Datetime, OBS_ID should be taken from
% L1 (not HK, unless implicit that it should).
% --
% Global attributes Datetime, OBS_ID, SOOP_TYPE appear to be present in BICAS
% input L1R datasets, CURRENT datasets, and BIAS HK datasets. Not true for old
% SBM1 datasets (at least).
% Exception: OBS_ID is not in BIAS HK. /2020-09-17
%--
% GAs "Provider", "Parent_version" have been abolished according to RCS ICD
% 01/07, and SOL-SGS-TN-0009, 02/06.
%
%
% ARGUMENTS
% =========
% InputDatasetsMap
%       NOTE: This function does not really need all of InputDatasetsMap as
%       input (contains ZVs) but the function uses that input argument since
%       it is easily accessible where this function is called.
% OutputDataset
%       Class bicas.OutputDataset
% outputFilename
%       Output dataset filename. Could potentially be used for deriving
%       Glob.attrs. Datetime (time interval string), Data_version,
%       (DSI).
%
%
% RETURN VALUE
% ============
% OutGaSubset
%       Struct where each field name corresponds to a CDF global atttribute.
%       NOTE: Deviates from the usual variable naming conventions.
%       GlobalAttributesSubset field names have
%               the exact names of CDF GAs.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-03-31, based on two functions broken out of
% execute_SWM().
%
function OutGaSubset = derive_output_dataset_GAs(...
  InputDatasetsMap, OutputDataset, outputFilename, outputDsi, ...
  Bso, L)

% PROPOSAL: Automatic test code.
%
% PROPOSAL: Create class for storing GAs.
%   PRO: Can detect accidental overwriting/reuse of keys.
%
% PROPOSAL: Move from setting constants in skeletons to setting them here.
%   TODO-NI: Allowed by ROC?!
%   NOTE: Only applies to GAs. Skeletons also set ZVAs and available ZVs.
%   Ex: Acknowledgements, Instrument, Instrument_type, APPLICABLE.
%   PRO: Easier to set overlapping constants in code.
%   CON: Can not inspect differences/similarities by diffing skeletons.
%   CON: Rare to change them.
%   CON/PROBLEM: Risk of setting values in skeletons without realizing they are
%                overridden by BICAS. ==> Confusion ==> Wasted time.
%     PROPOSAL: Set to special human-readable value in skeleton.
%       Ex: "Value overwritten by RCS."

% ASSERTIONS
irf.assert.struct(OutputDataset.Ga, ...
  {'OBS_ID', 'SOOP_TYPE'}, {'Misc_calibration_versions'})

[~, level, ~] = solo.adm.disassemble_DATASET_ID(outputDsi);



OutGaSubset = OutputDataset.Ga;

OutGaSubset.Parents       = get_GA_Parents(InputDatasetsMap);

% IMPLEMENTATION NOTE: SPICE_KERNELS should be set also for L3, but this has not
% yet been implemented in skeletons.
if strcmp(level, 'L2')
  OutGaSubset.SPICE_KERNELS = get_GA_SPICE_KERNELS(InputDatasetsMap);
end

OutGaSubset.Software_name    = bicas.const.SWD_METADATA('SWD.identification.identifier');
OutGaSubset.Software_version = bicas.const.SWD_METADATA('SWD.release.version');

% BUG? Assigns local time, not UTC!!! ROC DFMD does not mention time zone.
OutGaSubset.Generation_date  = char(datetime("now","Format","uuuu-MM-dd'T'HH:mm:ss"));

% NOTE: Parsing OUTPUT dataset filename to set some GAs.
[logicalFileId, logicalSource, dataVersionStr, timeIntervalStr] = ...
  parse_dataset_filename(outputFilename);

% Ex: Logical_file_id="solo_L1_rpw-tds-surv-hist2d_20220301_V01"
OutGaSubset.Logical_file_id  = logicalFileId;

% Logical_source:
% NOTE: Overwrites skeleton value. Can otherwise not handle -cdag.
% NOTE: Could in principle be set by assuming
%       lowercase(ga_DATASET_ID) = ga_Logical_source if not for -cdag and.
% Ex: Logical_source="solo_L1_rpw-tds-surv-hist2d"
OutGaSubset.Logical_source   = logicalSource;   % Override skeleton.
OutGaSubset.Data_version     = dataVersionStr;
OutGaSubset.Datetime         = timeIntervalStr;
% OutGaSubset.Dataset_ID       = outputDsi; % Override skeleton. Wise?

% IMPLEMENTATION NOTE: Unclear if it is wise to overwrite GAs (1) Dataset_ID and
% (2) Logical_source. In principle, the skeletons should contain the correct
% values. In principle, the ideal solution is to assert that GAs Dataset_ID
% and Logical_source in the master CDF are the expected ones, but it is
% unclear if this fits well with how master CDFs are currently loaded and if
% one did, should one not do so for other values as well?
% NOTE: Could almost overwrite "Descriptor" too (can be derived from
% Logical_source or DATASET_ID), but it also includes human-readable text.



%===============================================================================
% "Metadata Definition for Solar Orbiter Science Data", SOL-SGS-TN-0009, 2/4:
%   "TIME_MIN   The date and time of the beginning of the first acquisition
%               for the data contained in the file"
%   "TIME_MAX   The date and time of the end of the last acquisition for the
%               data contained in the file"
%
% NOTE: The implementation does not consider the integration time of each
%       sample.
% NOTE: ZvsSubset.Epoch is already asserted to be increasing.
%
% NOTE: The exact format has historically been unclear and confused. Previously
%       interpreted as being Julian date and now as "ISO" (interpreted as
%       YYYY-MM-DDThh:mm:ss.xyz...")
%       ------------------------------------------------------------------------
%       https://gitlab.obspm.fr/ROC/DataPool/-/issues/16#note_19377
%       check_cdf_istp.solo_L3_rpw-bia.txt:
%       """"
%     	Global attribute TIME_MAX is of type CDF_DOUBLE.
% 	        Datatypes other than CDF_CHAR may be problematic.
%     	Global attribute TIME_MIN is of type CDF_DOUBLE.
% 	        Datatypes other than CDF_CHAR may be problematic.""""
%       ------------------------------------------------------------------------
%       https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/84
%       """"
%       TIME_MIN is ['2460091.5010676757']
%       TIME_MIN from metadata equates to 2000-01-01T11:58:55.818460091
%       Start time from filename: 20230527
%       Start time from epoch variable: 2023-05-27T00:01:32.247168768
%
%
%       TIME_MAX is ['2460092.5010688445']
%       TIME_MAX from metadata equates to 2000-01-01T11:58:55.818460092
%       End time from filename:
%       End time from epoch variable: 2023-05-28T00:01:32.348148608
%       """"
%       ------------------------------------------------------------------------
%       https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/85
%       has the same complaint as above.
%===============================================================================
OutGaSubset.TIME_MIN = bicas.utils.TT2000_to_UTC_str(OutputDataset.Zv.Epoch(1  ));
OutGaSubset.TIME_MAX = bicas.utils.TT2000_to_UTC_str(OutputDataset.Zv.Epoch(end));



OutGaSubset.MODS = bicas.const.GA_MODS_DB.get_MODS_strings_CA(outputDsi);



% ROC DFMD hints that value should not be set dynamically. (See meaning of
% non-italic black text for GA name in table.)
%DataObj.GlobalAttribute.CAVEATS = ?!!



% ~ASSERTION
if ~isscalar(OutGaSubset.Parents)
  [settingValue, settingKey] = Bso.get_fv(...
    'INPUT_CDF.GA_PARENTS_MISMATCH_POLICY');
  bicas.default_anomaly_handling(...
    L, settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
    ['The value of the input CDF files'' global attribute "Parents"', ...
    ' differ (and they should not, or?).'], ...
    'BICAS:DatasetFormat')
  % NOTE: Maybe the wrong choice of error ID, "DatasetFormat".
end

% ASSERTION: Required subset for every dataset
% --------------------------------------------
% NOTE: GAs can be conditional. Ex: MODS, Provider
% TODO-DEC: Is this assertion sensible? There are many permanent GAs not
%           mentioned here. Remove assertion?
irf.assert.struct(OutGaSubset, ...
  {'Parents', ...
  'Datetime', 'OBS_ID', 'SOOP_TYPE'}, 'all')
end







function ga_Parents = get_GA_Parents(InputDatasetsMap)

ga_Parents = {};

keysCa = InputDatasetsMap.keys;
for i = 1:numel(keysCa)

  InputDataset = InputDatasetsMap(keysCa{i});

  % NOTE: Parsing INPUT dataset filename to set some GAs.
  [logicalFileId, ~, ~, ~] = parse_dataset_filename(...
    irf.fs.get_name(InputDataset.filePath));
  ga_Parents{end+1} = ['CDF>', logicalFileId];

end    % for
end







function ga_SPICE_KERNELS = get_GA_SPICE_KERNELS(InputDatasetsMap)

keysCa = InputDatasetsMap.keys;
ga_SPICE_KERNELS = cell(0, 1);

for i = 1:numel(keysCa)

  InputDataset = InputDatasetsMap(keysCa{i});

  % Read GA, but convert to format which represents zero kernels as empty cell
  % array, since that is convenient for algorithm.
  if isfield(InputDataset.Ga, 'SPICE_KERNELS')
    parent_SPICE_KERNELS = InputDataset.Ga.SPICE_KERNELS;

    if isscalar(parent_SPICE_KERNELS) && any(ismember(parent_SPICE_KERNELS{1}, {'none', ' '}))
      parent_SPICE_KERNELS = cell(0, 1);
    end
  else
    parent_SPICE_KERNELS = cell(0, 1);
  end

  ga_SPICE_KERNELS = unique([ga_SPICE_KERNELS; parent_SPICE_KERNELS]);
end    % for

% Normalize to the data format used in datasets.
if isempty(ga_SPICE_KERNELS)
  ga_SPICE_KERNELS = {'none'};
end
end







% NOTE: Only works correctly for files that follow the official filenaming
% scheme. logicalFileId does not work for e.g. IRFU-internal filenaming
% extension, as e.g. for test files that might be sent to ROC as part of
% official RCS test package.
%
% NOTE: Does not change case.
% NOTE: Is effectively a wrapper around solo.adm.parse_dataset_filename().
%
function [logicalFileId, logicalSource, dataVersionStr, timeIntervalStr] ...
  = parse_dataset_filename(filename)
[~, basename, ~] = fileparts(filename);

% NOTE: Will include IRFU-internal filenaming extension.
logicalFileId = basename;

% Actually parse the dataset filename.
R = solo.adm.parse_dataset_filename(filename);
assert(~isempty(R), 'BICAS:Assertion', ...
  ['Can not parse dataset filename "%s" and therefore not', ...
  ' derive values for global attributes', ...
  ' (Logical_source, Data_version, Datetime). Filename does not appear', ...
  ' to follow filenaming conventions.'], filename)

logicalSource   = R.fnDatasetIdCdag;
dataVersionStr  = R.versionStr;
timeIntervalStr = R.timeIntervalStr;
end
