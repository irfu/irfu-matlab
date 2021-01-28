%
% Execute a "S/W mode" as (indirectly) specified by the CLI arguments. This
% function should be agnostic of CLI syntax.
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% SwModeInfo
% InputFilePathMap  : containers.Map with
%    key   = prodFuncInputKey
%    value = Path to input file
% OutputFilePathMap : containers.Map with
%    key   = prodFuncOutputKey
%    value = Path to output file
%
%
% "BUGS"
% ======
% - Sets GlobalAttributes.Generation_date in local time (no fixed time zone,
%   e.g. UTC+0).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-06-09
%
function execute_sw_mode(...
        SwModeInfo, InputFilePathMap, OutputFilePathMap, ...
        masterCdfDir, rctDir, NsoTable, SETTINGS, L)

    % TODO-NI: How verify dataset ID and dataset version against constants?
    %    NOTE: Need to read CDF first.
    %    NOTE: Need S/W mode.
    %
    % PROPOSAL: Verify output zVariables against master CDF zVariable dimensions
    %           (accessible with dataobj, even for zero records).
    %   PROPOSAL: function matchesMaster(DataObj, MasterDataobj)
    %       PRO: Want to use dataobj to avoid reading file (input dataset) twice.
    %
    % NOTE: Things that need to be done when writing PDV-->CDF
    %       Read master CDF file.
    %       Compare PDV variables with master CDF variables (only write a subset).
    %       Check variable types, sizes against master CDF.
    %       Write GlobalAttributes: Calibration_version, Parents, Parent_version,
    %           Generation_date, Logical_file_id,
    %           Software_version, SPECTRAL_RANGE_MIN/-MAX (optional?), TIME_MIN/-MAX
    %       Write VariableAttributes: pad value? (if master CDF does not contain
    %           a correct value), SCALE_MIN/-MAX
    %
    % PROPOSAL: Print variable statistics also for zVariables which are created with fill values.
    %   NOTE: These do not use NaN, but fill values.



    % ASSERTION: Check that all input & output dataset paths (strings) are
    % unique.
    % NOTE: Manually entering CLI argument, or copy-pasting BICAS call, can
    % easily lead to reusing the same path by mistake, and e.g. overwriting an
    % input file.
    datasetFileList = [InputFilePathMap.values(), OutputFilePathMap.values()];
    assert(numel(unique(datasetFileList)) == numel(datasetFileList), ...
        'BICAS:execute_sw_mode:CLISyntax', ...
        ['Input and output dataset paths are not all unique.', ...
        ' This hints of a manual mistake', ...
        ' in the CLI arguments in call to BICAS.'])



    %=================================
    % READ CDFs
    % ---------
    % Iterate over all the INPUT CDFs
    %=================================
    InputDatasetsMap = containers.Map();
    for i = 1:length(SwModeInfo.inputsList)
        prodFuncInputKey = SwModeInfo.inputsList(i).prodFuncInputKey;
        inputFilePath    = InputFilePathMap(prodFuncInputKey);

        %=======================
        % Read dataset CDF file
        %=======================
        [Zv, GlobalAttributes]             = bicas.read_dataset_CDF(...
            inputFilePath, SETTINGS, L);
        InputDatasetsMap(prodFuncInputKey) = struct(...
            'Zv', Zv, ...
            'Ga', GlobalAttributes);



        %===========================================
        % ASSERTIONS: Check GlobalAttributes values
        %===========================================
        % NOTE: Can not use
        % bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(Zv) since
        % not all zVariables have same number of records.
        % Ex: Metadata such as ACQUISITION_TIME_UNITS.
        if ~isfield(GlobalAttributes, 'Dataset_ID')
            error('BICAS:execute_sw_mode:Assertion:DatasetFormat', ...
                ['Input dataset does not contain (any accepted variation', ...
                ' of) the global attribute Dataset_ID.\n    File: "%s"'], ...
                inputFilePath)
        end
        cdfDatasetId = GlobalAttributes.Dataset_ID{1};

        if ~strcmp(cdfDatasetId, SwModeInfo.inputsList(i).datasetId)
            [settingValue, settingKey] = SETTINGS.get_fv(...
                'INPUT_CDF.GA_DATASET_ID_MISMATCH_POLICY');
            anomalyDescrMsg = sprintf(...
                ['The input CDF dataset''s stated DATASET_ID does', ...
                ' not match value expected from the S/W mode.\n', ...
                '    File: %s\n', ...
                '    Global attribute GlobalAttributes.Dataset_ID{1} : "%s"\n', ...
                '    Expected value:                                 : "%s"\n'], ...
                inputFilePath, cdfDatasetId, SwModeInfo.inputsList(i).datasetId);
            bicas.default_anomaly_handling(L, ...
                settingValue, settingKey, ...
                'E+W+illegal', anomalyDescrMsg, 'BICAS:DatasetFormat')
        end
    end



    %==============
    % PROCESS DATA
    %==============
    [settingNpefValue, settingNpefKey] = SETTINGS.get_fv(...
        'OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE');
    if ~settingNpefValue
        %==========================
        % CALL PRODUCTION FUNCTION
        %==========================
        OutputDatasetsMap = SwModeInfo.prodFunc(InputDatasetsMap, rctDir, NsoTable);
    else
        OutputDatasetsMap = [];
        L.logf('warning', ...
            'Disabled processing due to setting %s.', settingNpefKey)
    end



    %==================================
    % WRITE CDFs
    % ----------
    % Iterate over all the OUTPUT CDFs
    %==================================
    % ASSERTION: The output datasets generated by processing equal the datasets
    % required by the s/w mode.
    EJ_library.assert.castring_sets_equal(...
        OutputDatasetsMap.keys, ...
        {SwModeInfo.outputsList.prodFuncOutputKey});
    % ITERATE
    for iOutputCdf = 1:length(SwModeInfo.outputsList)
        OutputInfo = SwModeInfo.outputsList(iOutputCdf);

        prodFuncOutputKey = OutputInfo.prodFuncOutputKey;
        outputFilePath    = OutputFilePathMap(prodFuncOutputKey);



        %========================
        % Write dataset CDF file
        %========================
        masterCdfPath = fullfile(...
            masterCdfDir, ...
            bicas.get_master_CDF_filename(...
                OutputInfo.datasetId, ...
                OutputInfo.skeletonVersion));

        if ~settingNpefValue
            % CASE: Nominal
            OutputDataset = OutputDatasetsMap(OutputInfo.prodFuncOutputKey);

            ZvsSubset = OutputDataset.Zv;

            GaSubset = derive_output_dataset_GlobalAttributes(...
                InputDatasetsMap, OutputDataset, ...
                EJ_library.fs.get_name(outputFilePath), SETTINGS, L);
        else
            % CASE: No processing.
            ZvsSubset = [];
        end
        bicas.write_dataset_CDF( ...
            ZvsSubset, GaSubset, outputFilePath, masterCdfPath, ...
            SETTINGS, L );
    end



end   % execute_sw_mode



% Function for dynamically deriving global attributes for a specific output dataset
% given the global attributes of multiple input datasets.
%
% NOTE: Some of the global attribute values determined here are
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
%       input (contains zVars) but the function uses that input argument since
%       it is easily accessible where this function is called.
% OutputDataset
%       Struct from processing with fields
%           Ga.(globAttrName)
%               Subset of global attribute values that should be used. Currently
%               only includes:
%                   .Datetime
%           .Zv.(zvName)
%               zVariables.
%
%
% RETURN VALUE
% ============
% OutGaSubset
%       Struct where each field name corresponds to a CDF global atttribute.
%       NOTE: Deviates from the usual variable naming conventions.
%       GlobalAttributesSubset field names have
%               the exact names of CDF global attributes.
% outputFilename
%       Output dataset filename. Could potentially be used for deriving
%       Glob.attrs. Datetime (time interval string), Data_version,
%       (DATASET_ID).
%       NOTE: Not yet used. 
%
function OutGaSubset = derive_output_dataset_GlobalAttributes(...
        InputDatasetsMap, OutputDataset, outputFilename, SETTINGS, L)

    % PGA = Parents' GlobalAttributes.

    if ~isscalar(OutputDataset.Ga.Datetime)
        [settingValue, settingKey] = SETTINGS.get_fv(...
            'OUTPUT_CDF.GLOBAL_ATTRIBUTES.Datetime_NOT_SCALAR_POLICY');
        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
            ['Global attribute "Datetime" for output dataset', ...
            ' is not a MATLAB scalar (i.e. the global attribute does not consist', ...
            ' of exactly ONE string). This may be due to the corresponding input', ...
            ' dataset value being similarily incorrect.'], ...
            'BICAS:execute_sw_mode:Datetime')
    end

    OutGaSubset = OutputDataset.Ga;



    OutGaSubset.Parent_version = {};
    OutGaSubset.Parents        = {};
    OutGaSubset.Provider       = {};
    keysCa = InputDatasetsMap.keys;
    for i = 1:numel(keysCa)
        InputGa = InputDatasetsMap(keysCa{i}).Ga;

        % ASSERTION
        % NOTE: ROC DFMD is not completely clear on which version number should
        % be used.
        % NOTE: Stores all values to be safe.
        assert(isscalar(InputGa.Data_version), ...
            'BICAS:execute_sw_mode:DatasetFormat', ...
            ['Global attribute "Data_version" for input dataset', ...
            ' with key=%s is not a MATLAB scalar (i.e. the global attribute is', ...
            ' not exactly ONE string).'], ...
            keysCa{i})

        % 2020-12-16, EJ: Has found input datasets to have global
        % attribute "Data_version" values which are either NUMERIC or STRINGS
        % (e.g. "02"). Varies.

        % NOTE: Using Data_version to set Parent_version.
        OutGaSubset.Parent_version{end+1} = InputGa.Data_version{1};
        OutGaSubset.Parents       {end+1} = ['CDF>', InputGa.Logical_file_id{1}];
        OutGaSubset.Provider              = union(OutGaSubset.Provider, InputGa.Provider);
    end



    OutGaSubset.Software_name       = bicas.constants.SWD_METADATA('SWD.identification.name');
    OutGaSubset.Software_version    = bicas.constants.SWD_METADATA('SWD.release.version');
    % Static value?!!
    OutGaSubset.Calibration_version = SETTINGS.get_fv('OUTPUT_CDF.GLOBAL_ATTRIBUTES.Calibration_version');
    % BUG? Assigns local time, not UTC!!! ROC DFMD does not mention time zone.
    OutGaSubset.Generation_date     = datestr(now, 'yyyy-mm-ddTHH:MM:SS');         
    OutGaSubset.Logical_file_id     = get_logical_file_id(outputFilename);
    %DataObj.GlobalAttributes.SPECTRAL_RANGE_MIN
    %DataObj.GlobalAttributes.SPECTRAL_RANGE_MAX

    %---------------------------------------------------------------------------
    % "Metadata Definition for Solar Orbiter Science Data", SOL-SGS-TN-0009:
    %   "TIME_MIN   The date and time of the beginning of the first acquisition
    %               for the data contained in the file"
    %   "TIME_MAX   The date and time of the end of the last acquisition for the
    %               data contained in the file"
    %   States that TIME_MIN, TIME_MAX should be "Julian day" (not "modified
    %   Julian day", which e.g. OVT uses internally).
    %
    % NOTE: Implementation does not consider the integration time of each
    % sample.
    % NOTE: juliandate() is consistent with Julian date converter at
    % https://www.onlineconversion.com/julian_date.htm
    % NOTE: ZvsSubset.Epoch already asserted to be monotonically increasing.
    %
    % NOTE: Exact format unclear from documentation, autochecks.
    % NOTE: Issue for autochecks on L3:
    %       https://gitlab.obspm.fr/ROC/DataPool/-/issues/16
    % check_cdf_istp.solo_L3_rpw-bia.txt:
    %   """"
    % 	Global attribute TIME_MAX is of type CDF_DOUBLE.
    % 	    Datatypes other than CDF_CHAR may be problematic.
    % 	Global attribute TIME_MIN is of type CDF_DOUBLE.
    % 	    Datatypes other than CDF_CHAR may be problematic.""""
    % NOTE: ROC data reprocessed ~2021-01-25,
    % solo_L1_rpw-bia-current-cdag_20201201-20201231_V01.cdf (version number
    % probably not part of official versioning) uses
    %     TIME_MIN (1 entry):
    %         0 (CDF_CHAR/17):        "2459184.982450046"
    %     TIME_MAX (1 entry):
    %         0 (CDF_CHAR/17):        "2459215.007218565"
    % Note the number of decimals. No exponent. Other files with ten decimals.
    %
    % PROPOSAL: Copy values from the corresponding values from the relevant input dataset.
    %   CON: Does not work for downsampled.
    %   CON: There has historically been problems with copying bad values from
    %        not-up-to-date input datasets.
    %---------------------------------------------------------------------------
    % NOTE: Choosing 10 decimals (instead of 9) so that time resolution is
    % higher than highest LFR sampling frequency (not sure of highest for
    % TDS-LFM).
    TIME_MINMAX_FORMAT = '%.10f';
    gaTimeMinNbr = juliandate(EJ_library.cdf.TT2000_to_datevec(OutputDataset.Zv.Epoch(1  )));
    gaTimeMaxNbr = juliandate(EJ_library.cdf.TT2000_to_datevec(OutputDataset.Zv.Epoch(end)));
    OutGaSubset.TIME_MIN = sprintf(TIME_MINMAX_FORMAT, gaTimeMinNbr);
    OutGaSubset.TIME_MAX = sprintf(TIME_MINMAX_FORMAT, gaTimeMaxNbr);
    
    % ROC DFMD hints that value should not be set dynamically. (See meaning of
    % non-italic black text for global attribute name in table.)
    %DataObj.GlobalAttribute.CAVEATS = ?!!

    
    
    % ~ASSERTION
    if ~isscalar(OutGaSubset.Parents)
        [settingValue, settingKey] = SETTINGS.get_fv(...
            'INPUT_CDF.GA_PROVIDER_MISMATCH_POLICY');
        bicas.default_anomaly_handling(...
            L, settingValue, settingKey, 'E+W+illegal', ...
            ['The value of the input CDF files'' global attribute "Provider"', ...
            ' differ (and they should not, or?).'], ...
            'BICAS:execute_sw_mode:DatasetFormat')
        % NOTE: Maybe wrong choice of error ID "DatasetFormat".
    end

    % ASSERTION: Required subset for every dataset.
    EJ_library.assert.struct(OutGaSubset, ...
        {'Parents', 'Parent_version', 'Provider', ...
        'Datetime', 'OBS_ID', 'SOOP_TYPE'}, 'all')
end



% NOTE: Only works correctly for files that follow the official filenaming scheme.
%
% NOTE: Does not change case.
function logicalFileId = get_logical_file_id(filePath)
    % Use the filename without suffix.
    [~, basename, ~] = fileparts(filePath);
    logicalFileId = basename;
end
