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
%
%
% ARGUMENTS
% =========
% InputDatasetsMap
%       NOTE: This function does not really need all of InputDatasetsMap as
%       input (contains ZVs) but the function uses that input argument since
%       it is easily accessible where this function is called.
% OutputDataset
%       Struct from processing with fields
%           .Ga.(globAttrName)
%               Subset of GA values that should be used.
%           .Zv.(zvName)
%               zVariables.
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
% PROPOSAL: Create class for GAs.
%   PRO: Can detect accidental overwriting/reuse of keys.

% ASSERTIONS
irf.assert.struct(OutputDataset.Ga, ...
  {'OBS_ID', 'SOOP_TYPE'}, {'Misc_calibration_versions'})



OutGaSubset = OutputDataset.Ga;

OutGaSubset.Parent_version = {};
OutGaSubset.Parents        = {};
OutGaSubset.Provider       = {};

%=============================
% Iterate over INPUT datasets
%=============================
keysCa = InputDatasetsMap.keys;
for i = 1:numel(keysCa)

  InputDatasetInfo = InputDatasetsMap(keysCa{i});
  InputGa          = InputDatasetInfo.Ga;

  % ASSERTION
  % NOTE: ROC DFMD is not completely clear on which version number should
  % be used.
  % NOTE: Stores all values to be safe.
  assert(isscalar(InputGa.Data_version), ...
    'BICAS:DatasetFormat', ...
    ['Global attribute "Data_version" for input dataset', ...
    ' with key=%s is not a MATLAB scalar (i.e. the global attribute is', ...
    ' not exactly ONE string).'], ...
    keysCa{i})

  % 2020-12-16, EJ: Has found input datasets to have global
  % attribute "Data_version" values which are either NUMERIC or STRINGS
  % (e.g. "02"). Varies.
  %
  %-----------------------------------------------------------------------
  % Ex: solo_L2_rpw-lfr-surv-bp1-cdag_20200625_V10.cdf:
  % NOTE: Seems to take last two digits from basename as Parent_version
  % values, even for RCTs (i.e. a BUG)!!!
  % Parent_version (9 entries):
  %       0 (CDF_CHAR/2):         "09"
  %       1 (CDF_CHAR/2):         "08"
  %       2 (CDF_CHAR/2):         "08"
  %       3 (CDF_CHAR/2):         "20"
  %       4 (CDF_CHAR/2):         "20"
  %       5 (CDF_CHAR/2):         "43"
  %       6 (CDF_CHAR/2):         "00"
  %       7 (CDF_CHAR/2):         "07"
  %       8 (CDF_CHAR/2):         "00"
  % Parents (9 entries):
  %       0 (CDF_CHAR/46):        "solo_L1_rpw-lfr-surv-bp1-cdag_20200625_V09.cdf"
  %       1 (CDF_CHAR/32):        "solo_HK_rpw-lfr_20200625_V08.cdf"
  %       2 (CDF_CHAR/32):        "solo_HK_rpw-bia_20200625_V08.cdf"
  %       3 (CDF_CHAR/40):        "SOLO_CAL_RCT-LFR-SCM_V20190123171020.cdf"
  %       4 (CDF_CHAR/41):        "SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf"
  %       5 (CDF_CHAR/40):        "SOLO_CAL_RCT-LFR-VHF_V20200720165743.cdf"
  %       6 (CDF_CHAR/55):        "SOLO_CAL_RCT-SCM_RPW_SCM-FM-MEB-PFM_V20190519120000.cdf"
  %       7 (CDF_CHAR/35):        "SOLO_CAL_RPW_BIAS_V202003101607.cdf"
  %       8 (CDF_CHAR/42):        "SOLO_CAL_RPW-HF-PREAMP_V20200624000000.cdf"
  %-----------------------------------------------------------------------
  % Ex: solo_L1_rpw-lfr-surv-bp2-cdag_20200625_V09.cdf
  % Parent_version (1 entry):
  %   0 (CDF_INT8/1):         8
  % Parents (1 entry):
  %   0 (CDF_CHAR/33):        "CDF>solo_L0_rpw-cdag_20200625_V08"
  %-----------------------------------------------------------------------
  % Ex: solo_L1_rpw-lfr-surv-bp2-cdag_20201225_V01.cdf, at internal
  % reprocessing (2021-01), not part of regular versioning.
  % Parent_version (1 entry):
  %      0 (CDF_CHAR/2):         "05"
  % Parents (1 entry):
  %      0 (CDF_CHAR/33):        "CDF>solo_L0_rpw-cdag_20201225_V05"
  %-----------------------------------------------------------------------
  % NOTE: Skeletons imply that Parent_version should be strings
  % (CDF_CHAR), though that setting should be inherited from some other
  % dataset, which might be a good or bad source.
  % Ex: SOLO_L2_RPW-LFR-SBM1-CWF-E_V11.skt:
  %   "Parent_version"      1:    CDF_CHAR     { " " }
  %-----------------------------------------------------------------------

  if isfield(InputGa, 'Provider')
    OutGaSubset.Provider = union(OutGaSubset.Provider, InputGa.Provider);
  else
    % IMPLEMENTATION NOTE: MAG datasets have been observed to not have
    % glob.attr. "Provider". VHT datasets have MAG datasets as parents.
    % /2021-05-05
    % Ex: solo_L2_mag-srf-normal_20200701_V02.cdf
    L.logf('warning', ...
      'Input dataset "%s"\ndoes not have CDF global attribute "Provider".\n', ...
      InputDatasetInfo.filePath)
  end

  % NOTE: Parsing INPUT dataset filename to set some GAs.
  [logicalFileId, ~, dataVersionStr, ~] = parse_dataset_filename(...
    irf.fs.get_name(InputDatasetInfo.filePath));
  % Sets string, not number. Correct?
  OutGaSubset.Parent_version{end+1} = dataVersionStr;
  OutGaSubset.Parents       {end+1} = ['CDF>', logicalFileId];

end



OutGaSubset.Software_name    = bicas.const.SWD_METADATA('SWD.identification.name');
OutGaSubset.Software_version = bicas.const.SWD_METADATA('SWD.release.version');

% BUG? Assigns local time, not UTC!!! ROC DFMD does not mention time zone.
OutGaSubset.Generation_date  = char(datetime("now","Format","uuuu-MM-dd'T'HH:mm:ss"));

% NOTE: Parsing OUTPUT dataset filename to set some GAs.
[logicalFileId, logicalSource, dataVersionStr, timeIntervalStr] = parse_dataset_filename(outputFilename);

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
% IMPLEMENTATION NOTE: Unclear if it is wise to overwrite GA Dataset_ID and
% Logical_source. In principle, the skeletons should contain the correct
% values. In principle, the ideal solution is to assert that GA Dataset_ID
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



MODS = bicas.const.GA_MODS_DB.get_MODS_strings_CA(outputDsi);
OutGaSubset.MODS = MODS;



% ROC DFMD hints that value should not be set dynamically. (See meaning of
% non-italic black text for GA name in table.)
%DataObj.GlobalAttribute.CAVEATS = ?!!



% ~ASSERTION
if ~isscalar(OutGaSubset.Parents)
  [settingValue, settingKey] = Bso.get_fv(...
    'INPUT_CDF.GA_PROVIDER_MISMATCH_POLICY');
  bicas.default_anomaly_handling(...
    L, settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
    ['The value of the input CDF files'' global attribute "Provider"', ...
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
  {'Parents', 'Parent_version', 'Provider', ...
  'Datetime', 'OBS_ID', 'SOOP_TYPE'}, 'all')
end



% NOTE: Only works correctly for files that follow the official filenaming
% scheme. logicalFileId does not work for e.g. IRFU-internal filenaming
% extension, as e.g. for test files that might be sent to ROC as part of
% official RCS test package.
%
% NOTE: Does not change case.
% NOTE: Wrapper around solo.adm.parse_dataset_filename().
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
