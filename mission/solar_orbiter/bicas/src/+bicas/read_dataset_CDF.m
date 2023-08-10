%
% Function that reads one dataset CDF file. Read elementary input process data
% from a CDF file.
%
%
% RETURN VALUES
% =============
% Dataset
%       Instance of bicas.InputDataset
%
%
% NOTE: Fill & pad values are replaced with NaN for float ZVs.
%       Other CDF data (attributes) are ignored.
% NOTE: Uses irfu-matlab's dataobj for reading the CDF file.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-09-17 as a separate file, by moving out the function from
% other file.
%
function Dataset = read_dataset_CDF(filePath, SETTINGS, L)
    
    tTicToc = tic();
    
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
    Zvs    = struct();
    ZvFvs  = struct();
    ZvsLog = struct();   % zVariables for logging.
    zVariableNameList = fieldnames(DataObj.data);
    for iZv = 1:length(zVariableNameList)
        zvName  = zVariableNameList{iZv};
        zvValue = DataObj.data.(zvName).data;
        
        ZvsLog.(zvName) = zvValue;
        
        [fillValue, padValue] = bicas.get_fill_pad_values(DataObj, zvName);
        ZvFvs.(zvName) = fillValue;
            
        %=================================================
        % Replace fill/pad values with NaN for FLOAT data
        %=================================================
        % TODO-NI: How does/should this work with integer fields that MUST
        %          also be stored as integers internally?!!!
        %    Ex: Epoch, ACQUISITION_TIME.
        if isfloat(zvValue)
            
            if ~isempty(fillValue)
                % CASE: There is a fill value.
                
                if any(ismember(zvName, ofvZvList))
                    L.logf('warning', ...
                        ['Overriding input CDF fill value with %d', ...
                        ' due to settings "%s" and "%s".'], ...
                        ofvFillVal, ovfZvSettingKey, ovfFvSettingKey)
                    fillValue = ofvFillVal;
                end
                
                zvValue = irf.utils.replace_value(zvValue, fillValue, NaN);
            end
            if ~disableReplacePadValue
                zvValue = irf.utils.replace_value(zvValue, padValue,  NaN);
            end
        end
        
        Zvs.(zvName) = zvValue;
    end

    
    
    % Log data read from CDF file
    bicas.utils.log_ZVs(ZvsLog, SETTINGS, L)
    
    
    
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
    [GlobalAttributes, fnChangeList] = irf.ds.normalize_struct_fieldnames(...
        DataObj.GlobalAttributes, ...
        {{{'DATASET_ID', 'Dataset_ID'}, 'Dataset_ID'}}, ...
        'Assert one matching candidate');
    msgFunc = @(oldFn, newFn) (sprintf(...
        ['Global attribute in input dataset', ...
        '\n    "%s"\nuses illegal alternative "%s" instead of "%s".\n'], ...
        filePath, oldFn, newFn));
    bicas.handle_struct_name_change(fnChangeList, SETTINGS, L, ...
        msgFunc, 'Dataset_ID', 'INPUT_CDF.USING_GA_NAME_VARIANT_POLICY')
    
    
    
    %===================
    % ASSERTIONS: Epoch
    %===================
    if ~isfield(Zvs, 'Epoch')
        error('BICAS:DatasetFormat', ...
            'Input dataset "%s" has no zVariable Epoch.', filePath)
    end
    if isempty(Zvs.Epoch)
        error('BICAS:DatasetFormat', ...
            'Input dataset "%s" contains an empty zVariable Epoch.', filePath)
    end
    
    
    
    %=========================================================================
    % ASSERTION: Increasing Epoch values
    % ----------------------------------
    % Examples:
    % solo_L1_rpw-lfr-surv-cwf-cdag_20200212_V01.cdf   (decrements 504 times)
    % solo_L1_rpw-lfr-surv-swf-cdag_20200212_V01.cdf   (1458 identical
    %                                                   consecutive pairs of
    %                                                   values)
    % solo_HK_rpw-bia_20200212_V01.cdf                 (decrements once)
    %=========================================================================
    % IMPLEMENTATION NOTE: SOLO_L1_RPW-BIA-CURRENT have increasing Epoch, but
    % not always MONOTONICALLY increasing Epoch.
    %
    % Check for increasing values, but NOT monotonically increasing.
    if ~issorted(Zvs.Epoch)
        
        anomalyDescrMsg = sprintf(...
            ['Input dataset "%s"\ncontains an Epoch zVariable', ...
            ' which values do not monotonically increment.\n'], ...
            filePath);
        
        [settingValue, settingKey] = SETTINGS.get_fv(...
            'INPUT_CDF.NON-INCREMENTING_ZV_EPOCH_POLICY');
        switch(settingValue)
            
            case 'SORT'
                bicas.default_anomaly_handling(...
                    L, settingValue, settingKey, 'other', ...
                    anomalyDescrMsg)
                
                % Sort (data) zVariables according to Epoch.
                [~, iSort] = sort(Zvs.Epoch);
                Zvs = select_ZVS_indices(Zvs, iSort);
                
                % % NOTE: Sorting Epoch does not remove identical values. Must therefore check again.
                % if ~issorted(Zvs.Epoch, 'strictascend')
                %     error('BICAS:DatasetFormat', ...
                %         ['zVariable Epoch in input dataset "%s"\n does not increase non-monotonically even after sorting.', ...
                %         ' It must contain multiple identical values (or the sorting algorithm does not work).'], ...
                %         filePath)
                % end

            otherwise
                bicas.default_anomaly_handling(...
                    L, settingValue, settingKey, ...
                    'E+W+illegal', ...
                    anomalyDescrMsg, 'BICAS:DatasetFormat')
        end
    end
    
    
    
    L.logf('info', 'File''s Global attribute: Dataset_ID       = "%s"', ...
        GlobalAttributes.Dataset_ID{1})
    L.logf('info', 'File''s Global attribute: Skeleton_version = "%s"', ...
        GlobalAttributes.Skeleton_version{1})



    % Create return value.
    Dataset = bicas.InputDataset(Zvs, ZvFvs, GlobalAttributes, filePath);
    
    
    
    % Just printing filename to avoid that actual time information is too far
    % right.
    bicas.log_speed_profiling(L, ...
        sprintf('%s: %s', mfilename, irf.fs.get_name(filePath)), ...
        tTicToc)
end



% Function that modifies ZVS to only contain specified records in specified
% order.
%
% Can be used for
% ** Re-ordering records (sorting Epoch).
% ** Filtering records (only keeping some).
%
% NOTE: Only want to modify the zVariables that contain data, i.e. for which CDF
% variable attribute DEPEND_0=Epoch, not metadata e.g. ACQUISITION_TIME_UNITS.
% Code does not use rigorous condition. Should ideally use zVariable attribute
% DEPEND_0. Is therefore not a generic function.
%
function Zvs = select_ZVS_indices(Zvs, iArray)
    % NOTE: Can not use
    % bicas.proc.utils.assert_struct_num_fields_have_same_N_rows(S); since want
    % to ignore but permit fields/ZVs with other number of records.
    
    fnList = fieldnames(Zvs);
    
    for iZv = 1:numel(fnList)
        fn = fnList{iZv};
        Zv = Zvs.(fn);
        
        % IMPLEMENTATION NOTE: Using size to distinguish data & metadata
        % zVariables.
        if size(Zv, 1) == size(Zvs.Epoch, 1)
            Zv = Zv(iArray, :,:,:,:,:,:,:);
        end
        Zvs.(fn) = Zv;
    end    
    
end
