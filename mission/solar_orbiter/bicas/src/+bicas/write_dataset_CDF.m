%
% Function that writes one dataset CDF file.
%
%
% NOTE: Replaces NaN-->fill value (according to master CDF) for floats
% (write_dataobj() does not do this).
%
%
% BUG: write_nominal_dataset_CDF seems to only be able to overwrite pre-existing
% (proper) CDF files, but not empty pre-existing files.
%
%
% ARGUMENTS
% =========
% ZvsSubset
%       Struct with those ZVs which should be written to CDF. "Subset" since
%       it excludes those ZVs in the master file, and which should NOT be
%       overwritten.
% GaSubset
%       Struct with fields representing a subset of the CDF global attributes. A
%       specific set of fields is required.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-06-24 as a separate file, by moving the function out from
% other file.
%
function write_dataset_CDF(...
  ZvsSubset, GaSubset, outputFile, masterCdfPath, Bso, L)

% PROPOSAL: Automatic tests for sub-functions. ==> Make subfunctions public.
%   NOTE: Uses irf.cdf.write_dataobj() for generic functionality for actual
%         writing of CDF.
%   PROPOSAL: Convert to class.
%       TODO-DEC: Naming of class & main function?
%           cdf
%           out
%           write
%           dataset
%           --
%           cdfout.write_dataset()
%   PROPOSAL: Convert bicas.write_dataset_CDF and bicas.read_dataset_CDF into
%             combined package (class?!).
%       CON: Read & write make up two distinct categories.
%       TODO-NI: Is there any overlap in functionality?
%===========================================================================
% This function needs GlobalAttributes values from the input files:
%    Data_version ??!!
%
% PROPOSAL: Accept GlobalAttributes for all input datasets?!
% PROBLEM: Too many arguments.
% TODO-DEC: Should function find the master CDF file itself?
%===========================================================================

%============
% ASSERTIONS
%============
% UI ASSERTION: Check for directory collision. Always error.
if exist(outputFile, 'dir')     % Checks for directory.
  error(...
    'BICAS:PathNotAvailable', ...
    'Intended output dataset file path matches a pre-existing directory.')
end
% UI ASSERTION: Check for output file path collision with pre-existing file
%               or directory.
% Command checks for file and directory (should not do just file).
if exist(outputFile, 'file')
  [settingValue, settingKey] = Bso.get_fv(...
    'OUTPUT_CDF.PREEXISTING_OUTPUT_FILE_POLICY');

  anomalyDescrMsg = sprintf(...
    'Intended output dataset file path "%s" matches a pre-existing file.', ...
    outputFile);
  bicas.default_anomaly_handling(L, settingValue, settingKey, 'ERROR_WARNING_ILLEGAL_SETTING', ...
    anomalyDescrMsg, 'BICAS:PathNotAvailable')
end



%===========================
% Create (modified) dataobj
%===========================
DataObj = init_modify_dataobj(...
  ZvsSubset, GaSubset, masterCdfPath, outputFile, Bso, L);



%=====================================
% CASE: ACTUALLY WRITE OUTPUT DATASET
%=====================================
write_nominal_dataset_CDF(DataObj, outputFile, Bso, L)



%====================================
% Assert expected CDF format version
%====================================
% IMPLEMENTATION NOTE: There does not seem to be a way of checking this without
% reading a CDF file. Also, this check would make a bad automated test since
% BICAS is developed in irfu-matlab, which may have one CDF format version,
% while the officially compliant may use another CDF format version.
CdfMetadata = spdfcdfinfo(outputFile);

% spdfcdfinfo help text:
%   FormatVersion        A string containing the version of the CDF
%                        library used to create the file
CdfVersion       = CdfMetadata.FormatVersion;
CdfVersionRegexp = Bso.get_fv('OUTPUT_CDF.FORMAT_VERSION_REGEXP');

L.logf('info', 'Output dataset CDF format version: %s', CdfVersion)

if isempty(regexp(CdfVersion, ['^', CdfVersionRegexp, '$'], 'once'))
  error( ...
    'BICAS:IllegalOutputCdfFormatVersion', ...
    ['Generated file "%s" has CDF format version "%s" which incompatible', ...
    ' with the permitted CDF format version in regular expression "%s"', ...
    ' (setting OUTPUT_CDF.FORMAT_VERSION_REGEXP).'], ...
    outputFile, CdfVersion, CdfVersionRegexp)
end

end







% Create a modified dataobj, based on a master CDF, that can be written to file.
%
% NOTE: Only uses global attribute values from
%   (1) GaSubset, and
%   (2) master CDF.
% bicas.execute_SWM: get_output_dataset_GAs() which sets
% global attributes dynamically.
%
% NOTE: Assertions require that ZvsSubset contains records of data. Can not
% easily submit "no data" for debugging purposes (deactivate processing but
% still write file).
%
function DataObj = init_modify_dataobj(...
  ZvsSubset, GaSubset, ...
  masterCdfPath, outputFile, Bso, L)

%============
% ASSERTIONS
%============
if ~isfield(ZvsSubset, 'Epoch')
  error('BICAS:SWMProcessing', ...
    'Data for output dataset "%s" has no zVariable Epoch.', ...
    outputFile)
end
if isempty(ZvsSubset.Epoch)
  error('BICAS:SWMProcessing', ...
    'Data for output dataset "%s" contains an empty zVariable Epoch.', ...
    outputFile)
end
if ~issorted(ZvsSubset.Epoch, 'strictascend')
  error('BICAS:SWMProcessing', ...
    ['Data for output dataset "%s" contains a zVariable Epoch', ...
    ' that does not increase monotonically.'], ...
    outputFile)
end



%======================
% Read master CDF file
%======================
L.logf('info', 'Reading master CDF file: "%s"', masterCdfPath)
DataObj = dataobj(masterCdfPath);



% IMPLEMENTATION NOTE: VHT datasets do not have a zVar QUALITY_FLAG.
% /2023-08-10
if isfield(ZvsSubset, 'QUALITY_FLAG')
  ZvsSubset.QUALITY_FLAG = modify_QUALITY_FLAG(ZvsSubset.QUALITY_FLAG, Bso, L);
end



%==================================================================
% Iterate over all OUTPUT PD field names
% Set corresponding dataobj zVariables
% --------------------------------------
% NOTE: ~zVariables from processing, i.e. not ZVs in master CDF.
%==================================================================
% NOTE: Only sets a SUBSET of the zVariables in master CDF.
L.log('info', 'Converting PDV to dataobj (CDF data structure)')
ZvsLog  = struct();   % ZVs for logging.
pdFieldNameList = fieldnames(ZvsSubset);
for iPdFieldName = 1:length(pdFieldNameList)
  zvName    = pdFieldNameList{iPdFieldName};
  zvValuePd = ZvsSubset.(zvName);

  % HACK: Normalize for the purpose of (old) ZV logging
  % FPA --> Array.
  % Use Master CDF FVs, except floats use NaN.
  % ---------------------------------------------------
  % IMPLEMENTATION NOTE: This is to avoid having to modify other old
  % non-FPA-aware code for now. Long term, should probably adapt logging
  % to only work on FPAs.
  if isa(zvValuePd, 'bicas.utils.FPArray')
    %------------------------------
    % CASE: ZV is stored as an FPA
    %------------------------------

    % Normalize FPA --> array with CDF FV.
    [fv, ~, masterMc] = bicas.get_dataobj_FV_pad_value_MC(DataObj, zvName);
    assert(...
      strcmp(masterMc, zvValuePd.mc), ...
      'Mismatching MATLAB classes for ZV "%s" between master CDF ("%s") and MATLAB (FPA) variable ("%s"). Master CDF: %s', ...
      zvName, masterMc, zvValuePd.mc, masterCdfPath)

    if ismember(zvValuePd.mc, {'single', 'double'})
      fv = cast(NaN, zvValuePd.mc);
    end
    assert(all(~ismember(fv, zvValuePd.NFP_1D_array())))

    zvValueLog = zvValuePd.array(fv);
  else
    %--------------------------------
    % CASE: ZV is stored as an array
    %--------------------------------
    zvValueLog = zvValuePd;
  end

  ZvsLog.(zvName) = zvValueLog;

  DataObj = overwrite_dataobj_ZV(DataObj, zvName, zvValuePd, L);
end



% Log data to be written to CDF file.
bicas.utils.log_ZVs(ZvsLog, Bso, L)



%======================================================================
% Use GaSubset to overwrite pre-existing (assertion) global attributes
%======================================================================
fnList = fieldnames(GaSubset);
for iFn = 1:numel(fnList)
  fn = fnList{iFn};

  assert(isfield(DataObj.GlobalAttributes, fn), ...
    'BICAS:DatasetFormat', ...
    ['Master CDF does not appear to contain global attribute', ...
    ' "%s" which the BICAS processing has produced/set.'], fn)

  DataObj.GlobalAttributes.(fn) = GaSubset.(fn);
end



%================================================================
% Handle still-empty zVariables (zero records; always anomalies)
%================================================================
for fn = fieldnames(DataObj.data)'
  zvName = fn{1};

  if isempty(DataObj.data.(zvName).data)

    DataObj = handle_empty_ZV_anomaly(DataObj, zvName, ...
      masterCdfPath, Bso, L);

  end    % if isempty(DataObj.data.(zvName).data)
end    % for

end    % init_modify_dataobj







% Enforce global max value for zVar QUALITY_FLAG specified in BSO.
%
% NOTE: Ignore fill positions/values.
% NOTE: Applies to both L2 and L3.
% NOTE: Using shortened zVariable name QF = QUALITY_FLAG to make algorithm
%       clearer.
function QfFpa = modify_QUALITY_FLAG(QfFpa, Bso, L)

assert(isa(QfFpa, 'bicas.utils.FPArray'))

% ============================================
% Obtain global QUALITY_FLAG cap from settings
% ============================================
[qfMax, key] = Bso.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX');
% ASSERT: Valid setting.
assert(isfinite(qfMax))
qfMax = uint8(qfMax);
assert(bicas.utils.validate_ZV_QUALITY_FLAG(qfMax), ...
  'BICAS:Assertion:ConfigurationBug', ...
  'Illegal BICAS setting "%s"=%i.', key, qfMax)

% Log if setting is lower than allowed ZV maximum.
if qfMax < bicas.const.QUALITY_FLAG_MAX
  L.logf('warning', ...
    'Using setting %s = %i to set global max value for zVariable QUALITY_FLAG.', ...
    key, qfMax);
end

% Enforce QUALITY_FLAG cap.
bTooHighQfFpa                     = (QfFpa >= qfMax);
QfFpa(bTooHighQfFpa.array(false)) = bicas.utils.FPArray(qfMax);

end    % modify_QUALITY_FLAG







% Function used by init_modify_dataobj() for using ZVs from processing, to
% overwrite ZVs in dataobj (from master CDF).
%
% ARGUMENTS
% =========
% zvValuePd
%       zVar value from processing (unaltered). PD = Processing Data.
%
%
% SHORTENINGS
% ===========
% ZVA = zVariable Attribute
%
function DataObj = overwrite_dataobj_ZV(DataObj, zvName, zvValuePd, L)

% ASSERTION: Master CDF already contains the zVariable.
if ~isfield(DataObj.data, zvName)
  error('BICAS:Assertion:SWMProcessing', ...
    ['Trying to write to zVariable "%s" that does not exist', ...
    ' in the master CDF file.'], zvName)
end



%======================================================================
% Prepare PDV zVariable value to save to CDF:
% (1a) FPAs --> array
% (1b) floats: Replace NaN-->fill value
% (3)  Convert zVar variable to the corresponding MATLAB class specified in
%      the master CDF.
%
% NOTE: If both fill values and pad values have been replaced with NaN
% (when reading CDF), then the code can not distinguish between fill
% values and pad values writing the CDF.
%======================================================================
[fv, ~, mc] = bicas.get_dataobj_FV_pad_value_MC(DataObj, zvName);
if isa(zvValuePd, 'bicas.utils.FPArray')
  % Normalize FPA --> array with CDF FV.
  assert(strcmp(mc, zvValuePd.mc))
  assert(all(~ismember(fv, zvValuePd.NFP_1D_array())))
  zvValueTemp = zvValuePd.array(fv);
elseif isfloat(zvValuePd)
  % Normalize array --> array with CDF FV.
  zvValueTemp = irf.utils.replace_value(zvValuePd, NaN, fv);
else
  zvValueTemp = zvValuePd;
end
cdfMc = irf.cdf.convert_CDF_type_to_MATLAB_class(...
  DataObj.data.(zvName).type, 'Permit MATLAB classes');
zvValueCdf = cast(zvValueTemp, cdfMc);



%===========================================================================
% Set zVar attrs. SCALEMIN, SCALEMAX according to zVar data
% ---------------------------------------------------------
% NOTE: Must handle fill values when deriving min & max.
%   NOTE: For floats: Use zvValuePd (not cvValueCdf).
%   NOTE: Must handle ZVs with ONLY fill values.
%
% TODO-DEC: How handle ZVs with only fill value/NaN?
%       Ex: EAC when no AC diff data
%       Ex: BW=0
%           Ex: SBM1, SBM2 test data
%   PROPOSAL: Set SCALEMIN=SCALEMAX=0
%   PROPOSAL: Keep SCALEMIN, SCALEMAX from skeleton.
%   PROPOSAL: Set to fill value.
%       CON: Legal?
%
% PROPOSAL: Move code into separate function.
%
% NOTE: Could be simplified with FPAs?
%
% NOTE: For some ZVs, SCALEMIN & SCALEMAX will not be wrong, but may also
% not be very meaningful.
% NOTE: Some ZVs might not change, i.e. min=max.
%   Ex: IBIAS1/2/3, SYNCHRO_FLAG
%
% NOTE: Implementation seems to work for floats.
%===========================================================================
if isnumeric(zvValueCdf)
  iZv = find(strcmp(DataObj.VariableAttributes.SCALEMAX(:,1), zvName));
  assert(isscalar(iZv))

  % Derive 1D array without fill values for inspection and derivation of
  % ZVAs.
  zvValueCdfLin = zvValueCdf(:);
  zvValueCdfLin(zvValueCdfLin == fv) = [];
  assert(all(~isnan(zvValueCdfLin)), ...
    'BICAS:Assertion', ...
    'zvValueCdfLin for zvName="%s" contains NaN despite being expected not to.', ...
    zvName)

  % SCALEMIN/-MAX from master CDFs.
  SCALEMIN_zvaMaster = DataObj.VariableAttributes.SCALEMIN{iZv, 2};
  SCALEMAX_zvaMaster = DataObj.VariableAttributes.SCALEMAX{iZv, 2};

  % NOTE: Epoch SCALEMAX is string (dataset bug?) /2021-02-05
  % ==> SCALEMAX_master not numeric when zvValue is.
  if isnumeric(SCALEMIN_zvaMaster) && isnumeric(SCALEMAX_zvaMaster)

    % SCALEMIN/-MAX from processed data.
    SCALEMIN_zvaPd = min(zvValueCdfLin, [], 'all');
    SCALEMAX_zvaPd = max(zvValueCdfLin, [], 'all');

    %L.logf('debug', 'zvName             = %s', zvName)
    %L.logf('debug', 'SCALEMIN_zvaMaster = %g', SCALEMIN_zvaMaster)
    %L.logf('debug', 'SCALEMIN_zvaPd     = %g', SCALEMIN_zvaPd)
    %L.logf('debug', 'SCALEMAX_zvaMaster = %g', SCALEMAX_zvaMaster)
    %L.logf('debug', 'SCALEMAX_zvaPd     = %g', SCALEMAX_zvaPd)

    %===============================================================
    % DECISION POINT: Set/update SCALEMIN & SCALEMAX used in actual
    % output CDF.
    %===============================================================
    if ~isempty(zvValueCdfLin)
      % CASE: There is a min & max.
      assert(~isempty(SCALEMIN_zvaPd))
      assert(~isempty(SCALEMAX_zvaPd))

      SCALEMIN_zvaCdf = SCALEMIN_zvaPd;
      SCALEMAX_zvaCdf = SCALEMAX_zvaPd;
    else
      % CASE: There is no min/max due to absence of non-fill value
      % data.
      SCALEMIN_zvaCdf = 0;    % NOTE: Must later be typecast.
      SCALEMAX_zvaCdf = 0;
    end

    % NOTE: zvValue has already been typecast to CDF type, but any
    % (future?) algorithm (above) for setting values may cancel that
    % (e.g. for constants). Must therefore typecast again, just to be
    % sure.
    SCALEMIN_zvaCdf = cast(SCALEMIN_zvaCdf, cdfMc);
    SCALEMAX_zvaCdf = cast(SCALEMAX_zvaCdf, cdfMc);

    DataObj.VariableAttributes.SCALEMIN{iZv,2} = SCALEMIN_zvaCdf;
    DataObj.VariableAttributes.SCALEMAX{iZv,2} = SCALEMAX_zvaCdf;
  end
end

% Set zVariable.
DataObj.data.(zvName).data = zvValueCdf;
end







% Function used by init_modify_dataobj() when finding empty zVar.
%
% ARGUMENTS
% =========
% masterCdfPath
%       NOTE: Only needed for anomaly description message.
%
function DataObj = handle_empty_ZV_anomaly(...
  DataObj, zvName, masterCdfPath, Bso, L)

%==============================================================
% CASE: zVariable has zero records.
% This indicates that it should have been set using PDV field.
%==============================================================

% NOTE: Useful to specify master CDF path in the case of having
% multiple output datasets. Will otherwise not know which output
% dataset is referred to. Note: Can still read master CDF from
% preceding log messages.
anomalyDescrMsg = sprintf(...
  ['Master CDF "%s" contains zVariable "%s" which has not been', ...
  ' set (i.e. it has zero records) after adding ', ...
  'processing data. This should only happen for incomplete processing.'], ...
  masterCdfPath, zvName);

mc = irf.cdf.convert_CDF_type_to_MATLAB_class(...
  DataObj.data.(zvName).type, 'Permit MATLAB classes');
isNumericZVar = isnumeric(cast(0.0, mc));

if isNumericZVar
  %====================
  % CASE: Numeric zVar
  %====================
  [settingValue, settingKey] = Bso.get_fv(...
    'OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY');
  switch(settingValue)

    case 'USE_FILLVAL'
      %=========================================================
      % Create correctly-sized zVariable using only fill values
      %=========================================================
      % NOTE: Assumes that
      % (1) there is a PD fields/zVariable Epoch, and
      % (2) this zVariable should have as many records as Epoch.
      L.logf('warning', ...
        ['Setting numeric master/output CDF zVariable', ...
        ' "%s" to presumed correct size using fill', ...
        ' values due to setting "%s" = "%s".'], ...
        zvName, settingKey, settingValue)

      nEpochRecords  = size(ZvsSubset.Epoch, 1);
      [fv, ~, ~] = bicas.get_dataobj_FV_pad_value_MC(DataObj, zvName);
      zvSize      = [nEpochRecords, DataObj.data.(fn{1}).dim];
      zvValueTemp = cast(zeros(zvSize), mc);
      zvValueCdf  = irf.utils.replace_value(zvValueTemp, 0, fv);

      DataObj.data.(zvName).data = zvValueCdf;

    otherwise
      bicas.default_anomaly_handling(L, ...
        settingValue, settingKey, ...
        'ERROR_WARNING_ILLEGAL_SETTING', anomalyDescrMsg, ...
        'BICAS:SWMProcessing:DatasetFormat')
  end

else
  %========================
  % CASE: Non-numeric zVar
  %========================
  [settingValue, settingKey] = Bso.get_fv(...
    'OUTPUT_CDF.EMPTY_NONNUMERIC_ZV_POLICY');
  bicas.default_anomaly_handling(L, ...
    settingValue, settingKey, ...
    'ERROR_WARNING_ILLEGAL_SETTING', anomalyDescrMsg, ...
    'BICAS:SWMProcessing:DatasetFormat')
end    % if isNumericZVar

end







% NOTE: Unclear if can overwrite CDFs. Varies depending on file.
%
function write_nominal_dataset_CDF(DataObj, outputFile, Bso, L)
%===========================================
% Write to CDF file using write_CDF_dataobj
%===========================================

[strictNumericZvSizePerRecord, settingKey] = Bso.get_fv(...
  'OUTPUT_CDF.write_dataobj.strictNumericZvSizePerRecord');
if ~strictNumericZvSizePerRecord
  L.logf('warning', [...
    '========================================================================================================\n', ...
    'Permitting master CDF zVariable size per record to differ from the output CDF zVariable size per record.\n', ...
    'This is due to setting %s = "%g"\n', ...
    '========================================================================================================\n'], ...
    settingKey, strictNumericZvSizePerRecord);
end

L.logf('info', 'Writing dataset CDF file: %s', outputFile)
irf.cdf.write_dataobj( ...
  outputFile, ...
  DataObj.GlobalAttributes, ...
  DataObj.data, ...
  DataObj.VariableAttributes, ...
  DataObj.Variables, ...
  'calculateMd5Checksum',              true, ...
  'strictEmptyZvClass',                Bso.get_fv('OUTPUT_CDF.write_dataobj.strictEmptyZvClass'), ...
  'strictEmptyNumericZvSizePerRecord', Bso.get_fv('OUTPUT_CDF.write_dataobj.strictEmptyNumericZvSizePerRecord'), ...
  'strictNumericZvSizePerRecord',      strictNumericZvSizePerRecord)
end
