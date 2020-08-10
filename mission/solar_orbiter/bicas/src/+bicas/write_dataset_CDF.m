%
% Function that writes one ___DATASET___ CDF file.
%
%
% BUG: write_nominal_dataset_CDF seems to only be able to overwrite pre-existing (proper) CDF files, but not empty
% files.
% 
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-06-24 as a separate file, moved from bicas.executed_sw_mode.
%
function write_dataset_CDF(...
        ZvsSubset, GlobalAttributesSubset, outputFile, masterCdfPath, SETTINGS, L)
% PROPOSAL: Be possible to overwrite files with write_nominal_dataset_CDF (when desired).
    
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
    
    
    
    [settingNpefValue, settingNpefKey] = SETTINGS.get_fv('OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE');
    if ~settingNpefValue
        DataObj = init_modif_dataobj(ZvsSubset, GlobalAttributesSubset, masterCdfPath, outputFile, SETTINGS, L);
        % NOTE: This call will fail if setting OUTPUT_CDF.NO_PROCESSING_EMPTY_FILE=1 since processing is disabled and
        % therefore ZvsSubset=[] (can not be generated).
    end
    
    
    
    %==============================
    % Checks before writing to CDF
    %==============================
    % UI ASSERTION: Check for directory collision. Always error.
    if exist(outputFile, 'dir')     % Checks for directory.
        error('BICAS:execute_sw_mode', 'Intended output dataset file path matches a pre-existing directory.')
    end
    
    % Check if file writing is deliberately disabled.
    [settingValue, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.WRITE_FILE_DISABLED');
    if settingValue
        L.logf('warning', 'Writing output CDF file is disabled via setting %s.', settingKey)
        return
    end
    
    % Check for output file path collision with pre-existing file.
    if exist(outputFile, 'file')    % Command checks for file and directory (can not do just file).
        [settingValue, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.PREEXISTING_OUTPUT_FILE_POLICY');
        
        anomalyDescrMsg = sprintf('Intended output dataset file path "%s" matches a pre-existing file.', outputFile);
        bicas.default_anomaly_handling(L, settingValue, settingKey, 'E+W+illegal', ...
            anomalyDescrMsg, 'BICAS:execute_sw_mode')
    end
    
    if ~settingNpefValue
        write_nominal_dataset_CDF(DataObj, outputFile, SETTINGS, L)
    else
        L.logf('warning', 'Writing empty output file due to setting %s.', settingNpefKey)
        write_empty_file(outputFile)
    end

end



% NOTE: Should always overwrite file.
function write_empty_file(filePath)
    fileId = fopen(filePath, 'w');    
    
    % ~ASSERTION
    if fileId == -1
        % NOTE: Technically non-BICAS error ID.
        error('BICAS:write_dataset_CDF:CanNotOpenFile', 'Can not open file: "%s"', filePath)
    end
    
    % NOTE: Does not have to write any data to create empty file.
    fclose(fileId);
end



% Create a modified dataobj that can be written to file.
%
%
% NOTE: Assertions require that ZvsSubset contains records of data. Can not easily submit "no data" for debugging
% purposes (deactivate processing but still write file).
%
function DataObj = init_modif_dataobj(ZvsSubset, GlobalAttributesSubset, masterCdfPath, outputFile, SETTINGS, L)
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
            [fillValue, ~] = bicas.get_fill_pad_values(DataObj, zvName);
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
    % Static value?!!
    DataObj.GlobalAttributes.Calibration_version = SETTINGS.get_fv('OUTPUT_CDF.GLOBAL_ATTRIBUTES.Calibration_version');
    % BUG? Assigns local time, not UTC!!! ROC DFMD does not mention time zone.
    DataObj.GlobalAttributes.Generation_date     = datestr(now, 'yyyy-mm-ddTHH:MM:SS');         
    DataObj.GlobalAttributes.Logical_file_id     = get_logical_file_id(outputFile);
    DataObj.GlobalAttributes.Parents             = GlobalAttributesSubset.Parents;
    DataObj.GlobalAttributes.Parent_version      = GlobalAttributesSubset.Parent_version;
    DataObj.GlobalAttributes.Provider            = GlobalAttributesSubset.Provider;                 
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
                        
                        nEpochRecords  = size(ZvsSubset.Epoch, 1);
                        [fillValue, ~] = bicas.get_fill_pad_values(DataObj, zvName);
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
    
end    % init_modif_dataobj



% NOTE: Unclear if can overwrite CDFs. Varies depending on file.
%
function write_nominal_dataset_CDF(DataObj, outputFile, SETTINGS, L)
    %===========================================
    % Write to CDF file using write_CDF_dataobj
    %===========================================
    
    [strictNumericZvSizePerRecord, settingKey] = SETTINGS.get_fv('OUTPUT_CDF.write_CDF_dataobj.strictNumericZvSizePerRecord');
    if strictNumericZvSizePerRecord
        L.logf('warning', [...
            '========================================================================================================\n', ...
            'Permitting master CDF zVariable size per record to differ from the output CDF zVariable size per record.\n'
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
