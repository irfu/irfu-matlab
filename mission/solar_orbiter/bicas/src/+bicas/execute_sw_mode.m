% Execute a "S/W mode" as (indirectly) specified by the CLI arguments.
% This function should be agnostic of CLI syntax.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
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
    % PROPOSAL: read_dataset_CDF as separate function file.
    
    
    
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
    [settingNpefValue, settingNpefKey] = SETTINGS.get_fv('OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE');
    if ~settingNpefValue
        OutputDatasetsMap = SwModeInfo.prodFunc(InputDatasetsMap, Cal);
    else
        OutputDatasetsMap = [];
        L.logf('warning', 'Disabled processing due to setting %s.', settingNpefKey)
    end
    
    
    
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
        
        if ~settingNpefValue        
            ZvsSubset = OutputDatasetsMap(OutputInfo.prodFuncOutputKey);
        else
            ZvsSubset = [];
        end
        bicas.write_dataset_CDF( ...
            ZvsSubset, globalAttributesSubset, outputFilePath, masterCdfPath, ...
            SETTINGS, L );
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
    pgaProviderList = {};
    for i = 1:length(GlobalAttributesCellArray)
        GlobalAttributesSubset.Parents       {end+1} = ['CDF>', GlobalAttributesCellArray{i}.Logical_file_id{1}];
        
        % NOTE: ROC DFMD is not completely clear on which version number should be used.
        GlobalAttributesSubset.Parent_version{end+1} = GlobalAttributesCellArray{i}.Data_version{1};
        
        pgaProviderList                      {end+1} = GlobalAttributesCellArray{i}.Provider{1};
    end
    
    pgaUniqueProviderList = unique(pgaProviderList);
    if ~isscalar(pgaUniqueProviderList)
        [settingValue, settingKey] = SETTINGS.get_fv('INPUT_CDF.GA_PROVIDER_MISMATCH_POLICY');
        bicas.default_anomaly_handling(...
            L, settingValue, settingKey, 'E+W+illegal', ...
            'The value of the input CDF files'' global attribute "Provider" differ (and they should not, or?).', ...
            'BICAS:execute_sw_mode:DatasetFormat')
        % NOTE: Maybe wrong choice of error ID "DatasetFormat".
    end
    
    GlobalAttributesSubset.Provider = pgaProviderList{1};
    
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
            [fillValue, padValue] = bicas.get_fill_pad_values(DataObj, zvName);
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
    bicas.proc_utils.log_zVars(ZvsLog, SETTINGS, L)
    
    
    
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
    % ASSERTION: Increasing Epoch values
    % ----------------------------------
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
