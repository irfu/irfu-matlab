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
    % IMPLEMENTATION NOTE: Generating FPAs for all ZVs in parallel with the old
    % representation of ZV values. The intention is that this will be a move
    % towards only using FPAs for all ZVs.
    %
    % PROBLEM: How convert some ZVs to FPAs at a time?
    %   PROPOSAL: Use list(s) for ZVs that should be FPA and/or not.
    %       PROPOSAL: List of ZVs that should not be FPAs.
    %           PRO: Will shrink with time. When it is empty then one know that the process is complete.
    %           CON: Difficult to begin with since needs to fill with long list.
    %
    % PROPOSAL: Test code:
    %   PROPOSAL: Separate ~inner function that accepts ~dataobj as input argument.
    %       NOTE: bicas.get_dataobj_FV_pad_value_MC() takes DO as
    %             argument.
    %   PROPOSAL: Write CDF file as part of test(!).
    
    % List of ZVs that should be represented as FPAs (and not as plain arrays).
    FPA_ZV_NAME_BIAS_HK_CA = {'HK_BIA_MODE_MUX_SET', 'HK_BIA_MODE_DIFF_PROBE', 'HK_BIA_DIFF_GAIN'};
    FPA_ZV_NAME_LFR_SCI_CA = {'BIAS_MODE_MUX_SET', 'QUALITY_FLAG', 'QUALITY_BITMASK'};
    
    FPA_ZV_NAME_CA = [FPA_ZV_NAME_BIAS_HK_CA, FPA_ZV_NAME_LFR_SCI_CA];

    tTicToc = tic();



    %===========
    % Read file
    %===========
    L.logf('info', 'Reading CDF file: "%s"', filePath)
    DataObj = dataobj(filePath);
    
    
    
    %=========================================================================
    % Copy zVariables (only the data) into analogous fields in smaller struct
    %=========================================================================
    L.log('info', 'Converting dataobj (CDF data structure) to PDV.')
    ZvFpa  = struct();
    Zvs    = struct();
    ZvFv   = struct();
    ZvsLog = struct();   % zVariables (name+value) for logging.

    zVariableNameList = fieldnames(DataObj.data);
    for iZv = 1:length(zVariableNameList)
        zvName    = zVariableNameList{iZv};
        zvValueDo = DataObj.data.(zvName).data;
        ZvsLog.(zvName) = zvValueDo;    % NOTE: Do = As found in dataobj (before typecasting & replacing FV-->NaN)).
        
        [fv, ~, mc] = bicas.get_dataobj_FV_pad_value_MC(DataObj, zvName);
        
        % =========================
        % Normalize ZV MATLAB class
        % =========================
        zvValueTyped = zvValueDo;    % ZV value after normalizing type/class
        if ~strcmp(mc, class(zvValueTyped))
            % Ex: /nonhome_data/SOLAR_ORBITER/bicas_test_input/LFR-SURV-CWF-SWF_test/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200213_V07.cdf
            %     ZV "COMMON_BIA_STATUS_FLAG".
            %     zvValue = [] (double), mc = 'int8' (?)
            L.logf('warning', ...
                'Loaded CDF ZV "%s" (size [%s]) has MATLAB class "%s" which is different from the stated MATLAB class "%s". Typecasting ZV to correct.', ...
                zvName, num2str(size(zvValueTyped)), class(zvValueTyped), mc)
            zvValueTyped = cast(zvValueTyped, mc);
        end


        %=================================================
        % Replace fill/pad values with NaN for FLOAT data
        %=================================================
        % TODO-NI: How does/should this work with integer fields that MUST
        %          also be stored as integers internally?!!!
        %    Ex: Epoch, ACQUISITION_TIME.
        zvValueTypedNan = zvValueTyped;
        if isfloat(zvValueTypedNan)
            
            if ~isempty(fv)
                % CASE: There is a fill value.
                
                zvValueTypedNan = irf.utils.replace_value(zvValueTypedNan, fv, NaN);
            end
        end
        
        if ismember(zvName, FPA_ZV_NAME_CA)
            % ===============================
            % Derive FPA representation of ZV
            % ===============================
            if isempty(fv)
                Fpa = bicas.utils.FillPositionsArray(zvValueTyped, 'NO_FILL_POSITIONS');
            else
                Fpa = bicas.utils.FillPositionsArray(zvValueTyped, 'FILL_VALUE', fv);
            end
            
            ZvFpa.(zvName) = Fpa;
        else
            Zvs.(zvName)   = zvValueTypedNan;
        end
        % IMPLEMENTATION NOTE: Should not really need to store explicit FV for
        % ZVs converted to FPAs, but temporary conversions FPA-->array require
        % this any way.
        ZvFv.(zvName) = fv;
    end

    
    
    % Log data read from CDF file
    bicas.utils.log_ZVs(ZvsLog, SETTINGS, L)
    
    
    
    %=================================================================================
    % Normalize the field/zVar names
    % ------------------------------
    % NOTE: At least the test files
    % solo_L1R_rpw-tds-lfm-cwf-e_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
    % solo_L1R_rpw-tds-lfm-rswf-e_20190523T080316-20190523T134337_V02_les-7ae6b5e.cdf
    % do not contain GA "DATASET_ID", only "Dataset_ID".
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
    Dataset = bicas.InputDataset(ZvFpa, Zvs, ZvFv, GlobalAttributes, filePath);
    
    
    
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
% DEPEND_0. Is therefore not a truly generic function.
%
% RETURN VALUE
% ============
% Zvs
%       Modified version of input argument. All fields with same number of rows
%       as field "Epoch" modified by only keeping the rows specified in iArray.
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
