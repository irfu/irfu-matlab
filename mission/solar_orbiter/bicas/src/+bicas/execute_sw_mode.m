% Execute a "S/W mode" as (indirectly) specified by the CLI arguments.
% This function should be agnostic of CLI syntax.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
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
% - Sets GlobalAttributes.Generation_date in local time (no fixed time zone).
% - Calls derive_output_dataset_GlobalAttributes for ALL input dataset and uses the result for ALL output datasets.
%   ==> If a S/W mode has multiple output datasets based on different sets of input datasets, then the GlobalAttributes
%   might be wrong. Should ideally be run on the exact input datasets (~EIn PDs) used to produce a specific output
%   dataset.
%
function execute_sw_mode(SwModeInfo, InputFilePathMap, OutputFilePathMap, masterCdfDir, calibrationDir, SETTINGS, L)
    %
    % QUESTION: How verify dataset ID and dataset version against constants?
    %    NOTE: Need to read CDF first.
    %    NOTE: Need S/W mode.
    %
    % PROPOSAL: Verify output zVariables against master CDF zVariable dimensions (accessible with dataobj, even for zero records).
    %   PROPOSAL: function matchesMaster(DataObj, MasterDataobj)
    %       PRO: Want to use dataobj to avoid reading file (input dataset) twice.
    %
    % QUESTION: What should be the relationship between data manager and S/W modes really?
    %           Should data manager check anything?
    %
    % NOTE: Things that need to be done when writing PDV-->CDF
    %       Read master CDF file.
    %       Compare PDV variables with master CDF variables (only write a subset).
    %       Check variable types, sizes against master CDF.
    %       Write GlobalAttributes: Calibration_version, Parents, Parent_version, Generation_date, Logical_file_id,
    %           Software_version, SPECTRAL_RANGE_MIN/-MAX (optional?), TIME_MIN/-MAX
    %       Write VariableAttributes: pad value? (if master CDF does not contain a correct value), SCALE_MIN/-MAX
    %
    % PROPOSAL: BUG FIX: Move global attributes into PDs somehow to let the processing functions collect the values during processing?
    %   PROPOSAL: Have PDs include global attributes in new struct structure.
    %             EIn PD:            EInPD(GlobalAttributes,          zVariables)   // All input dataset GAs.
    %             Intermediary PDs:     PD(GlobalAttributesCellArray, data)         // All input datasets GAs (multiple datasets).
    %             EOut PDs:         EOutPD(GlobalAttributesSubset,    data)         // Only those GAs that should be set. Should have been "collected" at this stage.
    %       PROBLEM: When collecting lists of GAs, must handle any overlap of input datasets when merging lists.
    %           Ex: (EIn1+EIn2-->Interm1; EIn1+EIn2-->Interm2; Interm1+Interm2-->EOut)
    %
    % PROPOSAL: Print variable statistics also for zVariables which are created with fill values.
    %   NOTE: These do not use NaN, but fill values.
    %
    % PROPOSAL: read_dataset_CDF, write_dataset_CDF as separate function files.
    %
    % PROPOSAL: Abolish get_fill_pad_values and using EJ_library.cdf.get_zvs_metadata_struct instead.
    
    
    
    GlobalAttributesCellArray = {};   % Use cell array since CDF global attributes may in principle contain different sets of attributes (field names).
    
    Cal = bicas.calib(calibrationDir, SETTINGS, L);
    
    
    
    % ASSERTION: Check that all input & output dataset paths (strings) are unique.
    % NOTE: Manually entering CLI argument, or copy-pasting BICAS call, can easily lead to reusing the same path by mistake,
    % and e.g. overwriting an input file.
    datasetFileList = [InputFilePathMap.values(), OutputFilePathMap.values()];
    assert(numel(unique(datasetFileList)) == numel(datasetFileList), 'BICAS:execute_sw_mode:CLISyntax', ...
        'Input and output dataset paths are not all unique. This hints of a manual mistake in the CLI arguments in call to BICAS.')
    
    
    
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
        [Zv, GlobalAttributes]             = read_dataset_CDF(inputFilePath, SETTINGS, L);
        InputDatasetsMap(prodFuncInputKey) = struct('Zv', Zv, 'Ga', GlobalAttributes);
        
        
        
        %===========================================
        % ASSERTIONS: Check GlobalAttributes values
        %===========================================
        % NOTE: Can not use bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(Zv) since not all zVariables have same number of
        % records. Ex: Metadata such as ACQUISITION_TIME_UNITS.
        if ~isfield(GlobalAttributes, 'Dataset_ID')
            error('BICAS:execute_sw_mode:Assertion:DatasetFormat', ...
                'Input dataset does not contain (any accepted variation of) the global attribute Dataset_ID.\n    File: "%s"', ...
                inputFilePath)
        end
        cdfDatasetId = GlobalAttributes.Dataset_ID{1};
        
        if ~strcmp(cdfDatasetId, SwModeInfo.inputsList(i).datasetId)
            [settingValue, settingKey] = SETTINGS.get_fv('INPUT_CDF.GA_DATASET_ID_MISMATCH_POLICY');
            anomalyDescrMsg = sprintf(...
                ['The input CDF dataset''s stated DATASET_ID does not match value expected from the S/W mode.\n', ...
                '    File: %s\n', ...
                '    Global attribute GlobalAttributes.Dataset_ID{1} : "%s"\n', ...
                '    Expected value:                                 : "%s"\n'], ...
                inputFilePath, cdfDatasetId, SwModeInfo.inputsList(i).datasetId);
            bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', anomalyDescrMsg, 'BICAS:DatasetFormat')
        end
        
        GlobalAttributesCellArray{end+1} = GlobalAttributes;
    end
    
    
    
    globalAttributesSubset = derive_output_dataset_GlobalAttributes(GlobalAttributesCellArray, SETTINGS, L);
    
    
    
    %==============
    % PROCESS DATA
    %==============
    OutputDatasetsMap = SwModeInfo.prodFunc(InputDatasetsMap, Cal);
    
    
    
    %==================================
    % WRITE CDFs
    % ----------
    % Iterate over all the OUTPUT CDFs
    %==================================
    for iOutputCdf = 1:length(SwModeInfo.outputsList)
        OutputInfo = SwModeInfo.outputsList(iOutputCdf);
        
        prodFuncOutputKey = OutputInfo.prodFuncOutputKey;
        outputFilePath    = OutputFilePathMap(prodFuncOutputKey);
        
        %========================
        % Write dataset CDF file
        %========================
        masterCdfPath = fullfile(...
            masterCdfDir, ...
            bicas.get_master_CDF_filename(OutputInfo.datasetId, OutputInfo.skeletonVersion));
        write_dataset_CDF ( ...
            OutputDatasetsMap(OutputInfo.prodFuncOutputKey), globalAttributesSubset, outputFilePath, masterCdfPath, ...
            OutputInfo.datasetId, SETTINGS, L );
    end
    
    
    
end   % execute_sw_mode







function GlobalAttributesSubset = derive_output_dataset_GlobalAttributes(GlobalAttributesCellArray, SETTINGS, L)
    % Function for global attributes for an output dataset from the global attributes of multiple input datasets (if there
    % are several).
    %
    % PGA = Parents' GlobalAttributes.
    %
    % RETURN VALUE
    % ============
    % GlobalAttributesSubset : Struct where each field name corresponds to a CDF global atttribute.
    %                          NOTE: Deviates from the usual variable naming conventions. GlobalAttributesSubset field names
    %                          have the exact names of CDF global attributes.
    %
    
    %ASSERT_MATCHING_TEST_ID = SETTINGS.get_fv('INPUT_CDF.GA_TEST_IDS_MISMATCH_POLICY');
    
    GlobalAttributesSubset.Parents        = {};            % Array in which to collect value for this file's GlobalAttributes (array-sized GlobalAttribute).
    GlobalAttributesSubset.Parent_version = {};
    %pgaTestIdList   = {};   % List = List with one value per parent.
    pgaProviderList = {};
    for i = 1:length(GlobalAttributesCellArray)
        GlobalAttributesSubset.Parents       {end+1} = ['CDF>', GlobalAttributesCellArray{i}.Logical_file_id{1}];
        
        % NOTE: ROC DFMD is not completely clear on which version number should be used.
        GlobalAttributesSubset.Parent_version{end+1} = GlobalAttributesCellArray{i}.Data_version{1};
        
        %pgaTestIdList                        {end+1} = GlobalAttributesCellArray{i}.Test_id{1};
        pgaProviderList                      {end+1} = GlobalAttributesCellArray{i}.Provider{1};
    end
    
    % NOTE: Test_id values can legitimately differ. E.g. "eeabc1edba9d76b08870510f87a0be6193c39051" and "eeabc1e".
    bicas.utils.assert_strings_equal(L, 0,                       pgaProviderList, 'The input CDF files'' GlobalAttribute "Provider" values differ.')
    %bicas.utils.assert_strings_equal(L, ASSERT_MATCHING_TEST_ID, pgaTestIdList,   'The input CDF files'' GlobalAttribute "Test_id" values differ.')
    
    % IMPLEMENTATION NOTE: Uses shortened "Test id" value in case it is a long one, e.g. "eeabc1edba9d76b08870510f87a0be6193c39051". Uncertain
    % how "legal" that is but it seems to be at least what people use in the filenames.
    % IMPLEMENTATION NOTE: Does not assume a minimum length for TestId since empty Test_id strings have been observed in
    % datasets. /2020-01-07
    GlobalAttributesSubset.Provider = pgaProviderList{1};
    %GlobalAttributesSubset.Test_Id  = pgaTestIdList{1}(1:min(7, length(pgaTestIdList{1})));
    
end







function [Zvs, GlobalAttributes] = read_dataset_CDF(filePath, SETTINGS, L)
    % Read elementary input process data from a CDF file. Copies all zVariables into fields of a regular structure.
    %
    %
    % RETURN VALUES
    % =============
    % Zvs              : ZVS = zVariables Struct. One field per zVariable (using the same name). The content of every such field equals the
    %                    content of the corresponding zVar.
    % GlobalAttributes : Struct returned from "dataobj".
    %
    %
    % NOTE: Fill & pad values are replaced with NaN for numeric data types.
    %       Other CDF data (attributes) are ignored.
    % NOTE: Uses irfu-matlab's dataobj for reading the CDF file.
    
    % NOTE: HK TIME_SYNCHRO_FLAG can be empty.
    
    
    disableReplacePadValue        = SETTINGS.get_fv('INPUT_CDF.REPLACE_PAD_VALUE_DISABLED');
    [ofvZvList,  ovfZvSettingKey] = SETTINGS.get_fv('INPUT_CDF.OVERRIDE_FILL_VALUE.ZV_NAMES');
    [ofvFillVal, ovfFvSettingKey] = SETTINGS.get_fv('INPUT_CDF.OVERRIDE_FILL_VALUE.FILL_VALUE');
    
    %===========
    % Read file
    %===========
    L.logf('info', 'Reading CDF file: "%s"', filePath)
    DataObj = dataobj(filePath);
    
    
    
    %=========================================================================
    % Copy zVariables (only the data) into analogous fields in smaller struct
    %=========================================================================
    L.log('info', 'Converting dataobj (CDF data structure) to PDV.')
    Zvs               = struct();
    ZvsLog            = struct();   % zVariables for logging.
    zVariableNameList = fieldnames(DataObj.data);
    for iZv = 1:length(zVariableNameList)
        zvName  = zVariableNameList{iZv};
        zvValue = DataObj.data.(zvName).data;
        
        ZvsLog.(zvName) = zvValue;
        
        %=================================================
        % Replace fill/pad values with NaN for FLOAT data
        %=================================================
        % QUESTION: How does/should this work with integer fields that should also be stored as integers internally?!!!
        %    Ex: ACQUISITION_TIME, Epoch.
        % QUESTION: How distinguish integer zVariables that could be converted to floats (and therefore use NaN)?
        if isfloat(zvValue)
            [fillValue, padValue] = get_fill_pad_values(DataObj, zvName);
            if ~isempty(fillValue)
                % CASE: There is a fill value.
                
                if any(ismember(zvName, ofvZvList))
                    L.logf('warning', ...
                        'Overriding input CDF fill value with %d due to settings "%s" and "%s".', ...
                        ofvFillVal, ovfZvSettingKey, ovfFvSettingKey)
                    fillValue = ofvFillVal;
                end
                
                zvValue = EJ_library.utils.replace_value(zvValue, fillValue, NaN);
            end
            if ~disableReplacePadValue
                zvValue = EJ_library.utils.replace_value(zvValue, padValue,  NaN);
            end
        else
            % Disable?! Only print warning if actually finds fill value which is not replaced?
            %L.logf('warning', 'Can not handle replace fill/pad values for zVariable "%s" when reading "%s".', zVariableName, filePath))
        end
        
        Zvs.(zvName) = zvValue;
    end
    
    
    
    % Log data read from CDF file
    bicas.proc_utils.log_zVars(ZvsLog, L)
    
    
    
    %=================================================================================
    % Normalize the field/zVar names
    % ------------------------------
    % NOTE: At least the test files
    % solo_L1R_rpw-tds-lfm-cwf-e_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
    % solo_L1R_rpw-tds-lfm-rswf-e_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
    % do not contain "DATASET_ID", only "Dataset_ID".
    %
    % NOTE: Has not found document that specifies the global attribute. /2020-01-16
    % https://gitlab.obspm.fr/ROC/RCS/BICAS/issues/7#note_11016
    % states that the correct string is "Dataset_ID".
    %=================================================================================
    [GlobalAttributes, fnChangeList] = EJ_library.utils.normalize_struct_fieldnames(DataObj.GlobalAttributes, ...
        {{{'DATASET_ID', 'Dataset_ID'}, 'Dataset_ID'}}, 'Assert one matching candidate');
    msgFunc = @(oldFn, newFn) (sprintf(...
        'Global attribute in input dataset\n    "%s"\nuses illegal alternative "%s" instead of "%s".\n', ...
        filePath, oldFn, newFn));
    bicas.handle_struct_name_change(fnChangeList, SETTINGS, L, msgFunc, 'Dataset_ID', 'INPUT_CDF.USING_GA_NAME_VARIANT_POLICY')
    
    
    
    %=================
    % Checks on Epoch
    %=================
    if ~isfield(Zvs, 'Epoch')
        error('BICAS:execute_sw_mode:DatasetFormat', 'Input dataset "%s" has no zVariable Epoch.', filePath)
    end
    if isempty(Zvs.Epoch)
        error('BICAS:execute_sw_mode:DatasetFormat', 'Input dataset "%s" contains an empty zVariable Epoch.', filePath)
    end
    
    
    
    %===============================================================================================
    % Check for increasing Epoch values
    % ---------------------------------
    % Examples:
    % solo_L1_rpw-lfr-surv-cwf-cdag_20200212_V01.cdf   (decrements 504 times)
    % solo_L1_rpw-lfr-surv-swf-cdag_20200212_V01.cdf   (1458 identical consecutive pairs of values)
    % solo_HK_rpw-bia_20200212_V01.cdf                 (decrements once)
    %===============================================================================================
    % IMPLEMENTATION NOTE: SOLO_L1_RPW-BIA-CURRENT have increasing Epoch, but not always MONOTONICALLY increasing Epoch.
    if ~issorted(Zvs.Epoch)   % Check for increasing values, but NOT monotonically increasing.
        
        anomalyDescrMsg = sprintf('Input dataset "%s"\ncontains an Epoch zVariable which values do not monotonically increment.\n', filePath);
        
        [settingValue, settingKey] = SETTINGS.get_fv('INPUT_CDF.NON-INCREMENTING_ZV_EPOCH_POLICY');
        switch(settingValue)
            case 'SORT'
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'other', ...
                    anomalyDescrMsg)
                
                % Sort (data) zVariables according to Epoch.
                [~, iSort] = sort(Zvs.Epoch);
                Zvs = select_ZVS_indices(Zvs, iSort);
                
                %             % NOTE: Sorting Epoch does not remove identical values. Must therefore check again.
                %             if ~issorted(Zvs.Epoch, 'strictascend')
                %                 error('BICAS:execute_sw_mode:DatasetFormat', ...
                %                     ['zVariable Epoch in input dataset "%s"\n does not increase non-monotonically even after sorting.', ...
                %                     ' It must contain multiple identical values (or the sorting algorithm does not work).'], ...
                %                     filePath)
                %             end
                
            otherwise
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:execute_sw_mode:DatasetFormat')
        end
    end
    
    
    
    L.logf('info', 'File''s Global attribute: Dataset_ID       = "%s"', GlobalAttributes.Dataset_ID{1})
    L.logf('info', 'File''s Global attribute: Skeleton_version = "%s"', GlobalAttributes.Skeleton_version{1})
    
end



function Zvs = select_ZVS_indices(Zvs, iArray)
% Function that modifies ZVS to only contain specified records in specified order.
%
% Can be used for
% ** Re-ordering records (sorting Epoch).
% ** Filtering records (only keeping some).
%
% NOTE: Only want to modify the zVariables that contain data, i.e. for which CDF variable attribute DEPEND_0=Epoch, not
% metadata e.g. ACQUISITION_TIME_UNITS. Code does not use rigorous condition. Should ideally use zVariable attribute
% DEPEND_0. Is therefore not a generic function.

    % NOTE: Can not use bicas.proc_utils.assert_struct_num_fields_have_same_N_rows(S); since want to ignore but permit
    % fields/zVars with other number of records.
    
    fnList = fieldnames(Zvs);
    
    for iZv = 1:numel(fnList)
        fn = fnList{iZv};
        Zv = Zvs.(fn);
        
        % IMPLEMENTATION NOTE: Using size to distinguish data & metadata zVariables.
        if size(Zv, 1) == size(Zvs.Epoch, 1)
            Zv = Zv(iArray, :,:,:,:,:,:,:);
        end
        Zvs.(fn) = Zv;
    end    
    
end



function write_dataset_CDF(...
        ZvsSubset, GlobalAttributesSubset, outputFile, masterCdfPath, datasetId, SETTINGS, L)
    %
    % Function that writes one ___DATASET___ CDF file.
    %
    
    %=======================================================================================================================
    % This function needs GlobalAttributes values from the input files:
    %    One value per file:      Data_version (for setting Parent_version).
    %    Data_version ??!!
    %
    % PROPOSAL: Accept GlobalAttributes for all input datasets?!
    % PROBLEM: Too many arguments.
    % QUESTION: Should function find the master CDF file itself?
    %   Function needs the dataset ID for it anyway.
    %   Function should check the master file anyway: Assert existence, GlobalAttributes (dataset ID, SkeletonVersion, ...)
    %=======================================================================================================================
    
    
    
    %===================
    % ASSERTIONS: Epoch
    %===================
    if ~isfield(ZvsSubset, 'Epoch')
        error('BICAS:execute_sw_mode', 'Data for output dataset "%s" has no zVariable Epoch.', outputFile)
    end
    if isempty(ZvsSubset.Epoch)
        error('BICAS:execute_sw_mode', 'Data for output dataset "%s" contains an empty zVariable Epoch.', outputFile)
    end
    if ~issorted(ZvsSubset.Epoch, 'strictascend')
        error('BICAS:execute_sw_mode', 'Data for output dataset "%s" contains a zVariable Epoch that does not increase monotonically.', outputFile)
    end
    
    
    
    %======================
    % Read master CDF file
    %======================
    L.logf('info', 'Reading master CDF file: "%s"', masterCdfPath)
    DataObj = dataobj(masterCdfPath);
    ZvsLog  = struct();   % zVars for logging.
    
    %=============================================================================================
    % Iterate over all OUTPUT PD field names (~zVariables) - Set corresponding dataobj zVariables
    %=============================================================================================
    % NOTE: Only sets a SUBSET of the zVariables in master CDF.
    pdFieldNameList = fieldnames(ZvsSubset);
    L.log('info', 'Converting PDV to dataobj (CDF data structure)')
    for iPdFieldName = 1:length(pdFieldNameList)
        zvName = pdFieldNameList{iPdFieldName};
        
        % ASSERTION: Master CDF already contains the zVariable.
        if ~isfield(DataObj.data, zvName)
            error('BICAS:execute_sw_mode:Assertion:SWModeProcessing', ...
                'Trying to write to zVariable "%s" that does not exist in the master CDF file.', zvName)
        end
        
        zvValue = ZvsSubset.(zvName);
        ZvsLog.(zvName) = zvValue;
        
        %================================================================================================================
        % Prepare PDV zVariable value:
        % (1) Replace NaN-->fill value
        % (2) Convert to the right MATLAB class
        %
        % NOTE: If both fill values and pad values have been replaced with NaN (when reading CDF), then the code can not
        % distinguish between fill values and pad values.
        %================================================================================================================
        if isfloat(zvValue)
            [fillValue, ~] = get_fill_pad_values(DataObj, zvName);
            zvValue        = EJ_library.utils.replace_value(zvValue, NaN, fillValue);
        end
        matlabClass = EJ_library.cdf.convert_CDF_type_to_MATLAB_class(DataObj.data.(zvName).type, 'Permit MATLAB classes');
        zvValue     = cast(zvValue, matlabClass);
        
        % Set zVariable.
        DataObj.data.(zvName).data = zvValue;
    end
    
    
    
    % Log data to be written to CDF file.
    bicas.proc_utils.log_zVars(ZvsLog, L)
    
    
    
    %==========================
    % Set CDF GlobalAttributes
    %==========================
    DataObj.GlobalAttributes.Software_name       = SETTINGS.get_fv('SWD.identification.name');
    DataObj.GlobalAttributes.Software_version    = SETTINGS.get_fv('SWD.release.version');
    DataObj.GlobalAttributes.Calibration_version = SETTINGS.get_fv('OUTPUT_CDF.GLOBAL_ATTRIBUTES.Calibration_version');         % "Static"?!!
    DataObj.GlobalAttributes.Generation_date     = datestr(now, 'yyyy-mm-ddTHH:MM:SS');         % BUG? Assigns local time, not UTC!!! ROC DFMD does not mention time zone.
    DataObj.GlobalAttributes.Logical_file_id     = get_logical_file_id(outputFile);
    DataObj.GlobalAttributes.Parents             = GlobalAttributesSubset.Parents;
    DataObj.GlobalAttributes.Parent_version      = GlobalAttributesSubset.Parent_version;
    %DataObj.GlobalAttributes.Data_version        = SETTINGS.get_fv('OUTPUT_CDF.DATA_VERSION');     % ROC DFMD says it should be updated in a way which can not be automatized?!!!
    DataObj.GlobalAttributes.Provider            = GlobalAttributesSubset.Provider;             % ROC DFMD contradictive if it should be set.
    %if SETTINGS.get_fv('OUTPUT_CDF.GLOBAL_ATTRIBUTES.SET_TEST_ID')
    %    DataObj.GlobalAttributes.Test_id         = GlobalAttributesSubset.Test_Id;              % ROC DFMD says that it should really be set by ROC.
    %end
    %DataObj.GlobalAttributes.SPECTRAL_RANGE_MIN
    %DataObj.GlobalAttributes.SPECTRAL_RANGE_MAX
    %DataObj.GlobalAttributes.TIME_MIN
    %DataObj.GlobalAttributes.TIME_MAX
    %DataObj.GlobalAttribute.CAVEATS ?!! ROC DFMD hints that value should not be set dynamically. (See meaning of non-italic black text for global attribute name in table.)
    
    
    
    %==============================================
    % Handle still-empty zVariables (zero records)
    %==============================================
    for fn = fieldnames(DataObj.data)'
        zvName = fn{1};
        
        if isempty(DataObj.data.(zvName).data)
            %===========================================================================================
            % CASE: zVariable has zero records, indicating that it should have been set using PDV field
            %===========================================================================================
            
            anomalyDescrMsg = sprintf(['Master CDF contains zVariable "%s" which has not been set (i.e. it has zero records) after adding ', ...
                'processing data. This should only happen for incomplete processing.'], ...
                zvName);
            
            matlabClass   = EJ_library.cdf.convert_CDF_type_to_MATLAB_class(DataObj.data.(zvName).type, 'Permit MATLAB classes');
            isNumericZVar = isnumeric(cast(0.000, matlabClass));
            
            if isNumericZVar
                %====================
                % CASE: Numeric zVar
                %====================
                [settingValue, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY');
                switch(settingValue)
                    case 'USE_FILLVAL'
                        %========================================================
                        % Create correctly-sized zVariable data with fill values
                        %========================================================
                        % NOTE: Assumes that
                        % (1) there is a PD fields/zVariable Epoch, and
                        % (2) this zVariable should have as many records as Epoch.
                        L.logf('warning', ...
                            ['Setting numeric master/output CDF zVariable "%s" to presumed correct size using fill', ...
                            ' values due to setting "%s" = "%s".'], ...
                            zvName, settingKey, settingValue)
                        
                        nEpochRecords = size(ZvsSubset.Epoch, 1);
                        [fillValue, ~] = get_fill_pad_values(DataObj, zvName);
                        zvSize  = [nEpochRecords, DataObj.data.(fn{1}).dim];
                        zvValue = cast(zeros(zvSize), matlabClass);
                        zvValue = EJ_library.utils.replace_value(zvValue, 0, fillValue);
                        
                        DataObj.data.(zvName).data = zvValue;
                        
                    otherwise
                        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', anomalyDescrMsg, ...
                            'BICAS:execute_sw_mode:SWModeProcessing:DatasetFormat')
                end
                
            else
                %========================
                % CASE: Non-numeric zVar
                %========================
                [settingValue, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.EMPTY_NONNUMERIC_ZV_POLICY');
                bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', anomalyDescrMsg, ...
                    'BICAS:execute_sw_mode:SWModeProcessing:DatasetFormat')
            end
        end
    end
    
    settingWriteFileDisabled = SETTINGS.get_fv('OUTPUT_CDF.WRITE_FILE_DISABLED');
    
    
    
    %==============================
    % Checks before writing to CDF
    %==============================
    % UI ASSERTION: Check for directory collision. Always error.
    if exist(outputFile, 'dir')     % Checks for directory.
        error('BICAS:execute_sw_mode', 'Intended output dataset file path matches a pre-existing directory.')
    end
    
    % Check if file writing is deliberately disabled.
    if settingWriteFileDisabled
        L.logf('warning', 'Writing output CDF file is disabled via setting OUTPUT_CDF.WRITE_FILE_DISABLED.')
        return
    end
    
    % Behaviour w.r.t. output file path collision with pre-existing file.
    if exist(outputFile, 'file')    % Checks for file and directory.
        [settingValue, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.PREEXISTING_OUTPUT_FILE_POLICY');
        
        anomalyDescrMsg = sprintf('Intended output dataset file path "%s" matches a pre-existing file.', outputFile);
        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
            anomalyDescrMsg, 'BICAS:execute_sw_mode')
        
    end
    
    %===========================================
    % Write to CDF file using write_CDF_dataobj
    %===========================================
    
    [strictNumericZvSizePerRecord, settingName] = SETTINGS.get_fv('OUTPUT_CDF.write_CDF_dataobj.strictNumericZvSizePerRecord');
    if strictNumericZvSizePerRecord
        logger.logf('warning', [...
            '=========================================================================================================', ...
            'Permitting master CDF zVariable size per record to differ from the output CDF zVariable size per record.\n'
            'This is due to setting %s = "%s"\n', ...
            '========================================================================================================='], ...
            strictNumericZvSizePerRecord, settingName);
    end
    
    L.logf('info', 'Writing dataset CDF file: %s', outputFile)
    EJ_library.cdf.write_dataobj( ...
        outputFile, ...
        DataObj.GlobalAttributes, ...
        DataObj.data, ...
        DataObj.VariableAttributes, ...
        DataObj.Variables, ...
        'calculateMd5Checksum',              true, ...
        'strictEmptyZvClass',                SETTINGS.get_fv('OUTPUT_CDF.write_CDF_dataobj.strictEmptyZvClass'), ...
        'strictEmptyNumericZvSizePerRecord', SETTINGS.get_fv('OUTPUT_CDF.write_CDF_dataobj.strictEmptyNumericZvSizePerRecord'), ...
        'strictNumericZvSizePerRecord',      strictNumericZvSizePerRecord)
    
end



function logicalFileId = get_logical_file_id(filePath)
    % Use the filename without suffix.
    [~, basename, ~] = fileparts(filePath);
   logicalFileId = basename;
end
% function logicalFileId = get_logical_file_id(datasetId, testId, provider, dataVersion)
%     % Construct a "Logical_file_id" as defined in the ROC DFMD.
%     %   NOTE 2020-03-19: Can not find in ROC DFMD 02/02,
%     % "The name of the CDF file without the ‘.cdf’ extension, using the file naming convention."
%     
%     bicas.assert_DATASET_ID(datasetId)
%     
%     if ~ischar(dataVersion ) || length(dataVersion)~=2
%         error('BICAS:execute_sw_mode:Assertion:IllegalArgument', 'Illegal dataVersion')
%     end
%     
%     providerParts = strsplit(provider, '>');
%     logicalFileId = [datasetId, '_', testId, '_', providerParts{1}, '_V', dataVersion];
% end



% fillValue : Empty if there is no fill value.
%
function [fillValue, padValue] = get_fill_pad_values(do, zVariableName)
    % NOTE: Uncertain how it handles the absence of a fill value. (Or is fill value mandatory?)
    % PROPOSAL: Remake into general-purpose function.
    % PROPOSAL: Remake into just using the do.Variables array?
    %    NOTE: Has to derive CDF variable type from do.Variables too.
    % PROPOSAL: Incorporate into dataobj?! Isn't there a buggy function/method there already?
    
    fillValue = getfillval(do, zVariableName);        % NOTE: Special function for dataobj.
    % NOTE: For unknown reasons, the fill value for tt2000 zVariables (or at least "Epoch") is stored as a UTC(?) string.
    if strcmp(do.data.(zVariableName).type, 'tt2000')
        fillValue = spdfparsett2000(fillValue);   % NOTE: Uncertain if this is the correct conversion function.
    end
    
    iZVariable = strcmp(do.Variables(:,1), zVariableName);
    padValue   = do.Variables{iZVariable, 9};
    % Comments in "spdfcdfinfo.m" should indirectly imply that column 9 is pad values since the structure/array
    % commented on should be identical.
end
