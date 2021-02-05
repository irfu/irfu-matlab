%
% Function that writes one dataset CDF file.
%
%
% BUG: write_nominal_dataset_CDF seems to only be able to overwrite pre-existing
% (proper) CDF files, but not empty pre-existing files.
%
%
% ARGUMENTS
% =========
% ZvsSubset
%       Struct with those zVars which should be written to CDF. "Subset" since
%       it excludes those zVars in the master file, and which should NOT be
%       overwritten.
% GaSubset
%       Struct with fields representing a subset of the CDF global attributes. A
%       specific set of fields is required.
% 
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-06-24 as a separate file, by moving out the function from
% bicas.executed_sw_mode.
%
function write_dataset_CDF(...
        ZvsSubset, GaSubset, outputFile, masterCdfPath, SETTINGS, L)
    
    % PROPOSAL: Use setting PROCESSING.ZV_QUALITY_FLAG_MAX to cap zVar
    %   QUALITY_FLAG.
    %   PRO: Replaces other code in multiple locations.
    %   PRO/CON: Presumes that every output dataset has a zVar QUALITY_FLAG.
    
    %===========================================================================
    % This function needs GlobalAttributes values from the input files:
    %    One value per file:      Data_version (for setting Parent_version).
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
            'BICAS:write_dataset_CDF', ...
            'Intended output dataset file path matches a pre-existing directory.')
    end
    % UI ASSERTION: Check for output file path collision with pre-existing file
    %               or directory.
    % Command checks for file and directory (should not do just file).
    if exist(outputFile, 'file')
        [settingValue, settingKey] = SETTINGS.get_fv(...
            'OUTPUT_CDF.PREEXISTING_OUTPUT_FILE_POLICY');
        
        anomalyDescrMsg = sprintf(...
            'Intended output dataset file path "%s" matches a pre-existing file.', ...
            outputFile);
        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
            anomalyDescrMsg, 'BICAS:write_dataset_CDF')
    end
    
    
    
    %===========================
    % Create (modified) dataobj
    %===========================
    % NPEF = No Processing Empty File
    [settingNpefValue, settingNpefKey] = SETTINGS.get_fv(...
        'OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE');
    if ~settingNpefValue
        
        %===================================================================
        % Set global max value for zVar QUALITY_FLAG
        % ------------------------------------------
        % NOTE: min(... 'includeNaN') implies that NaN always counts as the
        % lowest value.
        %===================================================================
        % PROPOSAL: Turn into generic function for capping QUALITY_FLAG based on
        % arbitrary setting.
        [value, key] = SETTINGS.get_fv('PROCESSING.ZV_QUALITY_FLAG_MAX');
        assert((0 < value) && (value <= 3), 'Illegal setting "%s"=%i.', key, value)
        if value < 3
            L.logf('warning', ...
                'Using setting %s = %i to set a zVar QUALITY_FLAG global max value.', ...
                key, value);
        end
        ZvsSubset.QUALITY_FLAG = min(...
            ZvsSubset.QUALITY_FLAG, ...
            value, 'includeNaN');
        
        
        DataObj = init_modif_dataobj(...
            ZvsSubset, GaSubset, masterCdfPath, outputFile, SETTINGS, L);
        % NOTE: This call will fail if setting
        % OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE=1 since processing is disabled and
        % therefore ZvsSubset=[] (can not be generated).
    end
    
    
    
    %===========================================
    % ASSERTIONS / Checks before writing to CDF
    %===========================================
    
    % Check if file writing is deliberately disabled.
    % NOTE: Do this as late as possible, in order to be able to test as much
    % code as possible without writing file.
    [settingValue, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.WRITE_FILE_DISABLED');
    if settingValue
        L.logf('warning', ...
            'Writing output CDF file is disabled via setting %s.', settingKey)
        return
    end
    
    if ~settingNpefValue
        %=====================================
        % CASE: ACTUALLY WRITE OUTPUT DATASET
        %=====================================
        write_nominal_dataset_CDF(DataObj, outputFile, SETTINGS, L)
    else
        L.logf('warning', ...
            'Writing empty output file due to setting %s.', settingNpefKey)
        write_empty_file(outputFile)
    end

end



% Create a modified dataobj that can be written to file. The dataobj is based on
% the master CDF.
%
% NOTE: Only uses global attribute values from
%   (1) GaSubset, and
%   (2) master CDF.
% bicas.execute_sw_mode: derive_output_dataset_GlobalAttributes() which sets
% global attributes dynamically.
%
%
% NOTE: Assertions require that ZvsSubset contains records of data. Can not
% easily submit "no data" for debugging purposes (deactivate processing but
% still write file).
%
function DataObj = init_modif_dataobj(...
        ZvsSubset, GaSubset, ...
        masterCdfPath, outputFile, SETTINGS, L)
    
    %============
    % ASSERTIONS
    %============
    if ~isfield(ZvsSubset, 'Epoch')
        error('BICAS:write_dataset_CDF', ...
            'Data for output dataset "%s" has no zVariable Epoch.', ...
            outputFile)
    end
    if isempty(ZvsSubset.Epoch)
        error('BICAS:write_dataset_CDF', ...
            'Data for output dataset "%s" contains an empty zVariable Epoch.', ...
            outputFile)
    end
    if ~issorted(ZvsSubset.Epoch, 'strictascend')
        error('BICAS:write_dataset_CDF', ...
            ['Data for output dataset "%s" contains a zVariable Epoch', ...
            ' that does not increase monotonically.'], ...
            outputFile)
    end
    
    
    
    %======================
    % Read master CDF file
    %======================
    L.logf('info', 'Reading master CDF file: "%s"', masterCdfPath)
    DataObj = dataobj(masterCdfPath);
    ZvsLog  = struct();   % zVars for logging.
    
    %==================================================================
    % Iterate over all OUTPUT PD field names
    % Set corresponding dataobj zVariables
    % --------------------------------------
    % NOTE: ~zVariables from processing, i.e. not zVars in master CDF.
    %==================================================================
    % NOTE: Only sets a SUBSET of the zVariables in master CDF.
    pdFieldNameList = fieldnames(ZvsSubset);
    L.log('info', 'Converting PDV to dataobj (CDF data structure)')
    for iPdFieldName = 1:length(pdFieldNameList)
        zvName = pdFieldNameList{iPdFieldName};
        
        % ASSERTION: Master CDF already contains the zVariable.
        if ~isfield(DataObj.data, zvName)
            error('BICAS:write_dataset_CDF:Assertion:SWModeProcessing', ...
                ['Trying to write to zVariable "%s" that does not exist', ...
                ' in the master CDF file.'], zvName)
        end
        
        % PD = Processing Data. Indicates that this is the unaltered value from
        %      processing.
        zvValuePd       = ZvsSubset.(zvName);
        ZvsLog.(zvName) = zvValuePd;
        
        %======================================================================
        % Prepare PDV zVariable value to save to CDF:
        % (1) Replace NaN-->fill value
        % (2) Convert to the right MATLAB class
        %
        % NOTE: If both fill values and pad values have been replaced with NaN
        % (when reading CDF), then the code can not distinguish between fill
        % values and pad values.
        %======================================================================
        if isfloat(zvValuePd)
            [fillValue, ~] = bicas.get_fill_pad_values(DataObj, zvName);
            zvValueTemp    = EJ_library.utils.replace_value(zvValuePd, NaN, fillValue);
        else
            zvValueTemp    = zvValuePd;
        end
        matlabClass = EJ_library.cdf.convert_CDF_type_to_MATLAB_class(...
            DataObj.data.(zvName).type, 'Permit MATLAB classes');
        zvValueCdf  = cast(zvValueTemp, matlabClass);
        
        
        
        %==============================================================
        % Modify zVar attrs. SCALEMIN, SCALEMAX according to zVar data
        %
        % EXPERIMENTAL - INCOMPLETE IMPLEMENTATION
        % ------------------------------------------------------------
        % NOTE: Must handle fill values when deriving min & max.
        %   NOTE: For floats: Use zvValuePd (not cvValueCdf).
        %   NOTE: Must handle zVars with ONLY fill values.
        %
        % NOTE: For some zVars, SCALEMIN & SCALEMAX will not be wrong, but may
        % also not be very meaningful.
        % NOTE: Some zVars might not change, i.e. min=max.
        %   Ex: IBIAS1/2/3, SYNCHRO_FLAG
        %
        % PROPOSAL: Assertion for absence of NaN.
        %==============================================================        
        if isfloat(zvValuePd) && 0    % DISABLED
            i = find(strcmp(DataObj.VariableAttributes.SCALEMAX(:,1), zvName));
            assert(isscalar(i))
            
            zva_SCALEMIN_master = DataObj.VariableAttributes.SCALEMIN{i,2};
            zva_SCALEMAX_master = DataObj.VariableAttributes.SCALEMAX{i,2};
            
            % NOTE: Epoch SCALEMAX is string (dataset bug?) /2021-02-05
            % ==> SCALEMAX_master not numeric when zvValue is.
            if isnumeric(zva_SCALEMIN_master) && isnumeric(zva_SCALEMAX_master)
            
                zva_SCALEMIN_Pd = min(zvValuePd, [], 'all');
                zva_SCALEMAX_Pd = max(zvValuePd, [], 'all');
                assert(~isnan(zva_SCALEMIN_Pd))
                assert(~isnan(zva_SCALEMAX_Pd))
            
                L.logf('debug', 'zvName              = %s', zvName)
                L.logf('debug', 'zva_SCALEMIN_master = %g', zva_SCALEMIN_master)
                L.logf('debug', 'zva_SCALEMIN_Pd     = %g', zva_SCALEMIN_Pd)
                L.logf('debug', 'zva_SCALEMAX_master = %g', zva_SCALEMAX_master)
                L.logf('debug', 'zva_SCALEMAX_Pd     = %g', zva_SCALEMAX_Pd)
                
                %=======================================================
                % DECISION POINT: How to set/update SCALEMIN & SCALEMAX
                %=======================================================
                zva_SCALEMIN_new = zva_SCALEMIN_Pd;
                zva_SCALEMAX_new = zva_SCALEMAX_Pd;
                
                % NOTE: zvValue has already been typecast to CDF type, but any
                % (future?) algorithm (above) for setting values may cancel
                % that. Must therefore typecast again, just to be sure.
                zva_SCALEMIN_new = cast(zva_SCALEMIN_new, matlabClass);
                zva_SCALEMAX_new = cast(zva_SCALEMAX_new, matlabClass);
            
                DataObj.VariableAttributes.SCALEMIN{i,2} = zva_SCALEMIN_new;
                DataObj.VariableAttributes.SCALEMAX{i,2} = zva_SCALEMAX_new;
            end
        end
        
        % Set zVariable.
        DataObj.data.(zvName).data = zvValueCdf;
    end
    
    
    
    % Log data to be written to CDF file.
    bicas.proc_utils.log_zVars(ZvsLog, SETTINGS, L)


    
    %======================================================================
    % Use GaSubset to overwrite pre-existing (assertion) global attributes
    %======================================================================
    fnList = fieldnames(GaSubset);
    for iFn = 1:numel(fnList)
        fn = fnList{iFn};
        
        assert(isfield(DataObj.GlobalAttributes, fn), ...
            'BICAS:write_dataset_CDF:init_modif_dataobj:DatasetFormat', ...
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
            
            matlabClass   = EJ_library.cdf.convert_CDF_type_to_MATLAB_class(...
                DataObj.data.(zvName).type, 'Permit MATLAB classes');
            isNumericZVar = isnumeric(cast(0.000, matlabClass));
            
            if isNumericZVar
                %====================
                % CASE: Numeric zVar
                %====================
                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'OUTPUT_CDF.EMPTY_NUMERIC_ZV_POLICY');
                switch(settingValue)
                    
                    case 'USE_FILLVAL'
                        %========================================================
                        % Create correctly-sized zVariable data with fill values
                        %========================================================
                        % NOTE: Assumes that
                        % (1) there is a PD fields/zVariable Epoch, and
                        % (2) this zVariable should have as many records as Epoch.
                        L.logf('warning', ...
                            ['Setting numeric master/output CDF zVariable', ...
                            ' "%s" to presumed correct size using fill', ...
                            ' values due to setting "%s" = "%s".'], ...
                            zvName, settingKey, settingValue)
                        
                        nEpochRecords  = size(ZvsSubset.Epoch, 1);
                        [fillValue, ~] = bicas.get_fill_pad_values(DataObj, zvName);
                        zvSize  = [nEpochRecords, DataObj.data.(fn{1}).dim];
                        zvValueTemp = cast(zeros(zvSize), matlabClass);
                        zvValueCdf  = EJ_library.utils.replace_value(zvValueTemp, 0, fillValue);
                        
                        DataObj.data.(zvName).data = zvValueCdf;
                        
                    otherwise
                        bicas.default_anomaly_handling(L, ...
                            settingValue, settingKey, ...
                            'E+W+illegal', anomalyDescrMsg, ...
                            'BICAS:write_dataset_CDF:SWModeProcessing:DatasetFormat')
                end
                
            else
                %========================
                % CASE: Non-numeric zVar
                %========================
                [settingValue, settingKey] = SETTINGS.get_fv(...
                    'OUTPUT_CDF.EMPTY_NONNUMERIC_ZV_POLICY');
                bicas.default_anomaly_handling(L, ...
                    settingValue, settingKey, ...
                    'E+W+illegal', anomalyDescrMsg, ...
                    'BICAS:write_dataset_CDF:SWModeProcessing:DatasetFormat')
            end
        end
    end
    
end    % init_modif_dataobj



% NOTE: Unclear if can overwrite CDFs. Varies depending on file.
%
function write_nominal_dataset_CDF(DataObj, outputFile, SETTINGS, L)
    %===========================================
    % Write to CDF file using write_CDF_dataobj
    %===========================================
    
    [strictNumericZvSizePerRecord, settingKey] = SETTINGS.get_fv(...
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
    EJ_library.cdf.write_dataobj( ...
        outputFile, ...
        DataObj.GlobalAttributes, ...
        DataObj.data, ...
        DataObj.VariableAttributes, ...
        DataObj.Variables, ...
        'calculateMd5Checksum',              true, ...
        'strictEmptyZvClass',                SETTINGS.get_fv('OUTPUT_CDF.write_dataobj.strictEmptyZvClass'), ...
        'strictEmptyNumericZvSizePerRecord', SETTINGS.get_fv('OUTPUT_CDF.write_dataobj.strictEmptyNumericZvSizePerRecord'), ...
        'strictNumericZvSizePerRecord',      strictNumericZvSizePerRecord)
end



% NOTE: Should always overwrite any pre-existing file.
%
function write_empty_file(filePath)
    fileId = fopen(filePath, 'w');    
    
    % ~ASSERTION
    if fileId == -1
        % NOTE: Technically non-BICAS error ID.
        error(...
            'BICAS:write_dataset_CDF:CanNotOpenFile', ...
            'Can not open file: "%s"', filePath)
    end
    
    % NOTE: Does not have to write any data to create empty file.
    fclose(fileId);
end
