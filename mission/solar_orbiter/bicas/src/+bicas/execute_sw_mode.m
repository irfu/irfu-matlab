% execute_sw_mode   Execute a "S/W mode" as (implicitly) specified by the CLI arguments.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% swModeCliParameter : The CLI argument value that represents the S/W mode.
% InputFilePathMap   : containers.Map with
%    keys   = PDIDs
%    values = Paths to input files.
%
%
% IMPORTANT NOTE
% ==============
% As of 2016-10-26: This function temporarily uses both an old way and a new way of writing CDFs in parallel for the
% purpose of testing. The old way is meant to be phased out eventually.
%
%
% "BUGS:"
% =======
% - Sets GlobalAttributes.Generation_date in local time (no fixed time zone).
% - Calls derive_output_dataset_GlobalAttributes for ALL input dataset and uses the result for ALL output datasets.
%   ==> If a S/W mode has multiple output datasets based on different sets of input datasets, then the GlobalAttributes
%   might be wrong. Should ideally be run on the exact input datasets (~EIn PDs) used to produce a specific output
%   dataset.
%
function execute_sw_mode(DataManager, swModeCliParameter, InputFilePathMap, outputDir)
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
% Separate functions for converting CDF file<-->PDV?
%   NOTE: May be different procedures (GlobalAttributes etc.) for all post-calibration datasets.
%   NOTE: There are values which are identical for all output datasets (for a given set of input files).
%   --
%   NOTE: Things that need to be done when reading CDF-->PDV
%       PDV assertions: Contains the right zVariables, with the right types & sizes.
%           NOTE: Some could be checked against the (input!) master CDF.
%           Should be more PDV-specific.
%       Assign the data to a specific data manager PDV/PDID.
%       Convert fill/pad values-->NaN for numeric types, or at least for floating point zVariables.
%           QUESTION: How should code know which policy applies for integers? (Ex: ACQUSITION_TIME)
%               PROPOSAL: Data manager gets information on fill/pad values and does the selected substitutions itself.
%               PROPOSAL: read-CDF function get information on which conversions should be made.
%
%       PROPOSAL: Separate function for all datasets should do...
%           S/W mode assertions: GlobalAttributes: DATASET_ID, Skeleton_version
%           Read GlobalAttributes to be used later: Source_name, Data_version, Provider, Test_id, 
%   --
%   NOTE: Things that need to be done when writing PDV-->CDF
%       Read master CDF file.
%       Compare PDV variables with master CDF variables (only write a subset).
%       Check variable types, sizes against master CDF.
%       Write GlobalAttributes: Calibration_version, Parents, Parent_version, Generation_date, Logical_file_id,
%           Software_version, SPECTRAL_RANGE_MIN/-MAX (optional?), TIME_MIN/-MAX
%       Write VariableAttributes: pad value? (if master CDF does not contain a correct value), SCALE_MIN/-MAX
%       Determine the output filename. ==> 
%       Collect the output filenames into the JSON output filenames list.
%
% PROPOSAL: BUG FIX: Move global attributes into PDs somehow to let the data_manager collect the values during processing?
%   PROPOSAL: Have PDs include global attributes in new struct structure.
%             EIn PD:            EInPD(GlobalAttributes,          zVariables)   // All input dataset GAs.
%             Intermediary PDs:     PD(GlobalAttributesCellArray, data)         // All input datasets GAs (multiple datasets).
%             EOut PDs:         EOutPD(GlobalAttributesSubset,    data)         // Only those GAs that should be set. Should have been "collected" at this stage.
%       PROBLEM: When collecting lists of GAs, must handle any overlap of input datasets when merging lists.
%           Ex: (EIn1+EIn2-->Interm1; EIn1+EIn2-->Interm2; Interm1+Interm2-->EOut)
%       PROPOSAL: Have get_processing_info do the work if it is the same for every PDID.
%           CON: There is no such functionality in get_processing_info. All the processing is supposed to be done by the processing function.
%               PROPOSAL: Can use switch-case for the data processing function, then wrap another processing function
%                         for global attributes around it and return that function.
%   PROPOSAL: Have all PDs use struct PD(GlobalAttributes, data).
%   PROPOSAL: Separate DM dependancies to handle GAs directly EIn-EOut without intermediary PDs.
%       PRO: Can collect all GAs in one step.
%       CON: DM must explicitly state all the EIn PDs, rather than collect them automatically indirectly (recursively).x
%   PROPOSAL: Since functionality is identical for all (L2R-->L2S) processing, use data_manager.get_elementary_input_PDIDs.
%       CON: There is no proper link from EIn PDs to dataset global attributes unless these are in the EIn PDs, but they
%            are not really needed for the processing so it would be contradictive if they were.
%
% PROPOSAL: Print variable statistics also for zVariables which are created with fill values.
%   NOTE: These do not use NaN, but fill values.

global SETTINGS

irf.log('n', sprintf('Output directory = "%s"', outputDir));       

% ASSERTION
if ~exist(outputDir, 'dir')
    error('BICAS:execute_sw_mode:Assertion:PathNotFound', 'Output directory "%s" does not exist.', outputDir)
end



%===========================================================================
% Give all input CDF files (from the CLI argument list) to the data manager
%===========================================================================
inputPdidList = InputFilePathMap.keys;
GlobalAttributesCellArray = {};   % Use cell array since CDF global attributes may in principle contain different sets of attributes (field names).

% testIdList = {};
for i = 1:length(inputPdidList)
    eInPdid = inputPdidList{i};
    inputFilePath = InputFilePathMap(eInPdid);
    
    
    [processData, GlobalAttributes] = read_dataset_CDF(eInPdid, inputFilePath);
    DataManager.set_elementary_input_process_data(eInPdid, processData);
    
    GlobalAttributesCellArray{end+1} = GlobalAttributes;
end



globalAttributesSubset = derive_output_dataset_GlobalAttributes(GlobalAttributesCellArray);
SwModeInfo = DataManager.get_extended_sw_mode_info(swModeCliParameter);



%==================================
% Iterate over all the output CDFs
%==================================
JsonOutputCdfFilenameListStruct = struct;    % Struct(!) representing a JSON object.
for iOutputCdf = 1:length(SwModeInfo.outputs)
    OutputInfo = SwModeInfo.outputs{iOutputCdf};
    
    eOutPdid = OutputInfo.PDID;
    
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    % PROCESS DATA
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%
    ProcessData = DataManager.get_process_data_recursively(eOutPdid);

    % Write dataset CDF file.
    masterCdfPath = bicas.get_master_CDF_path(OutputInfo.DATASET_ID, OutputInfo.SKELETON_VERSION_STR);
    [outputFilename] = write_dataset_CDF ( ...
        ProcessData, globalAttributesSubset, outputDir, @get_output_filename, masterCdfPath, OutputInfo.DATASET_ID );
    
    % Collect list (struct) of output files.
    JsonOutputCdfFilenameListStruct.( OutputInfo.SWD_OUTPUT_FILE_IDENTIFIER ) = outputFilename;
end



%============================================================
% Print JSON object describing the created file(s) to stdout
% ----------------------------------------------------------
% Required by the RCS ICD iss2rev2, section 3.3.
%============================================================
str = bicas.utils.JSON_object_str(...
    JsonOutputCdfFilenameListStruct, ...
    SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'), ...
    SETTINGS.get_fv('JSON_OBJECT_STR.VALUE_POSITION'));
bicas.stdout_disp(str);

end   % execute_sw_mode







function GlobalAttributesSubset = derive_output_dataset_GlobalAttributes(GlobalAttributesCellArray)
% Function for global attributes for an output dataset from the global attributes of multiple input datasets (if there
% are several).
%
% RETURN VALUE:
% GlobalAttributesSubset : Struct where each field name corresponds to a CDF global atttribute.
%                          NOTE: Deviates from the usual variable naming conventions. GlobalAttributesSubset field names
%                          have the exact names of CDF global attributes.

global SETTINGS
ASSERT_MATCHING_TEST_ID = SETTINGS.get_fv('INPUT_CDF_ASSERTIONS.MATCHING_TEST_ID');

GlobalAttributesSubset.Parents        = {};            % Array in which to collect value for this file's GlobalAttributes (array-sized GlobalAttribute).
GlobalAttributesSubset.Parent_version = {};
pgaTestIdList   = {};   % PGA = parents' GlobalAttributes. List = List with one value per parent.
pgaProviderList = {};
for i = 1:length(GlobalAttributesCellArray)
    GlobalAttributesSubset.Parents       {end+1} = ['CDF>', GlobalAttributesCellArray{i}.Logical_file_id{1}];
    
    % NOTE: ROC DFMD is not completely clear on which version number should be used.
    GlobalAttributesSubset.Parent_version{end+1} = GlobalAttributesCellArray{i}.Data_version{1};
    
    pgaTestIdList                        {end+1} = GlobalAttributesCellArray{i}.Test_id{1};
    pgaProviderList                      {end+1} = GlobalAttributesCellArray{i}.Provider{1};
end

% NOTE: Test_id values can legitimately differ. E.g. "eeabc1edba9d76b08870510f87a0be6193c39051" and "eeabc1e".
bicas.utils.assert_strings_equal(0,                       pgaProviderList, 'The input CDF files'' GlobalAttribute "Provider" values differ.')
bicas.utils.assert_strings_equal(ASSERT_MATCHING_TEST_ID, pgaTestIdList,   'The input CDF files'' GlobalAttribute "Test_id" values differ.')

% NOTE: Uses shortened "Test id" value in case it is a long one, e.g. "eeabc1edba9d76b08870510f87a0be6193c39051". Uncertain
% how "legal" that is but it seems to be at least what people use in the filenames.
GlobalAttributesSubset.Provider = pgaProviderList{1};
GlobalAttributesSubset.Test_Id  = pgaTestIdList{1}(1:7);

end







function [outputFilename] = write_dataset_CDF(...
    ProcessData, GlobalAttributesSubset, outputFileParentDir, FilenamingFunction, masterCdfPath, datasetId)
% Function that writes one ___dataset___ CDF file.
%

%==========================================================================
% This function needs GlobaAttributes values from the input files:
%    One value per file:      Data_version (for setting Parent_version).
%    One value for all files: Test_id
%    Data_version ??!!
%
% PROPOSAL: Accept GlobalAttributes for all input datasets?!
% PROBLEM: Too many arguments.
% QUESTION: Should function find the master CDF file itself?
%   Function needs the dataset ID for it anyway.
%   Function should check the master file anyway: Assert existence, GlobalAttributes (dataset ID, SkeletonVersion, ...)
%==========================================================================

global SETTINGS



%======================
% Read master CDF file
%======================
irf.log('n', sprintf('Reading master CDF file: "%s"', masterCdfPath))
DataObj = dataobj(masterCdfPath);

%=============================================================================================
% Iterate over all OUTPUT PD field names (~zVariables) - Set corresponding dataobj zVariables
%=============================================================================================
% NOTE: Only sets a SUBSET of the zVariables in master CDF.
pdFieldNameList = fieldnames(ProcessData);
irf.log('n', 'Converting PDV to dataobj (CDF data structure)')
%bicas.dm_utils.log_array('explanation')
for iPdFieldName = 1:length(pdFieldNameList)
    zVariableName = pdFieldNameList{iPdFieldName};
    
    % ASSERTION: Master CDF already contains the zVariable.
    if ~isfield(DataObj.data, zVariableName)
        error('BICAS:execute_sw_mode:Assertion:SWModeProcessing', ...
        'Trying to write to zVariable "%s" that does not exist in the master CDF file.', zVariableName)
    end
    
    zVariableData = ProcessData.(zVariableName);
    
    % Prepare PDV zVariable data:
    % (1) Replace NaN-->fill value
    % (2) Convert to the right MATLAB class.
    if isfloat(zVariableData)
        [fillValue, ~] = get_fill_pad_values(DataObj, zVariableName);
        zVariableData = bicas.utils.replace_value(zVariableData, NaN, fillValue);
    end
    matlabClass = bicas.utils.convert_CDF_type_to_MATLAB_class(DataObj.data.(zVariableName).type, 'Permit MATLAB classes');
    zVariableData = cast(zVariableData, matlabClass);
    
    %=================================================================================================
    % Log data to be written to CDF file
    % ----------------------------------
    % NOTE: Log messages should reflect the values WRITTEN to file.
    %       Process data variables are logged separately by data_manager.
    %=================================================================================================
    %bicas.dm_utils.log_struct_arrays(zVariableName, zVariableData);
    
    % Set zVariable.
    DataObj.data.(zVariableName).data = zVariableData;
end


    
%==========================
% Set CDF GlobalAttributes
%==========================
DataObj.GlobalAttributes.Software_name       = SETTINGS.get_fv('SWD_IDENTIFICATION.name');
DataObj.GlobalAttributes.Software_version    = SETTINGS.get_fv('SWD_RELEASE.version');
DataObj.GlobalAttributes.Calibration_version = SETTINGS.get_fv('CALIBRATION_VERSION');         % "Static"?!!
DataObj.GlobalAttributes.Generation_date     = datestr(now, 'yyyy-mm-ddTHH:MM:SS');         % BUG? Assigns local time, not UTC!!! ROC DFMD does not mention time zone.
DataObj.GlobalAttributes.Logical_file_id     = logical_file_id(...
    datasetId, GlobalAttributesSubset.Test_Id, ...
    GlobalAttributesSubset.Provider, ...
    SETTINGS.get_fv('OUTPUT_CDF.DATA_VERSION'));
DataObj.GlobalAttributes.Parents             = GlobalAttributesSubset.Parents;
DataObj.GlobalAttributes.Parent_version      = GlobalAttributesSubset.Parent_version;
DataObj.GlobalAttributes.Data_version        = SETTINGS.get_fv('OUTPUT_CDF.DATA_VERSION');     % ROC DFMD says it should be updated in a way which can not be automatized?!!!
DataObj.GlobalAttributes.Provider            = GlobalAttributesSubset.Provider;             % ROC DFMD contradictive if it should be set.
if SETTINGS.get_fv('OUTPUT_CDF.SET_TEST_ID')
    DataObj.GlobalAttributes.Test_id         = GlobalAttributesSubset.Test_Id;              % ROC DFMD says that it should really be set by ROC.
end
%DataObj.GlobalAttributes.SPECTRAL_RANGE_MIN
%DataObj.GlobalAttributes.SPECTRAL_RANGE_MAX
%DataObj.GlobalAttributes.TIME_MIN
%DataObj.GlobalAttributes.TIME_MAX
%DataObj.GlobalAttribute CAVEATS ?!! ROC DFMD hints that value should not be set dynamically. (See meaning of non-italic black text for global attribute name in table.)



%==============================================
% Handle still-empty zVariables (zero records)
%==============================================
for fn = fieldnames(DataObj.data)'
    zVariableName = fn{1};
    
    if isempty(DataObj.data.(zVariableName).data)
        % CASE: zVariable has zero records, indicating that should have been set using PDV field.
        
        logMsg = sprintf(['CDF contains zVariable "%s" which has not been set (i.e. it has zero records) after adding ', ...
            'processing data to the master CDF. This should only happen for incomplete processing.'], ...
            zVariableName);
        
        if SETTINGS.get_fv('OUTPUT_CDF.EMPTY_ZVARIABLES_SET_TO_FILL')
            irf.log('w', logMsg)
            matlabClass = bicas.utils.convert_CDF_type_to_MATLAB_class(DataObj.data.(zVariableName).type, 'Permit MATLAB classes');
            
            % ASSERTION: Require numeric type.
            if ~isnumeric(cast(0.000, matlabClass))
                error('BICAS:sw_execute_sw_mode:SWModeProcessing', ...
                    'zVariable "%s" is non-numeric. Can not set it to correctly-sized data with fill values (not implemented).', zVariableName)
            end
            
            %========================================================
            % Create correctly-sized zVariable data with fill values
            %========================================================
            % NOTE: Assumes that
            % (1) there is a PD fields/zVariable Epoch, and
            % (2) this zVariable should have as many records as Epoch.
            irf.log('w', sprintf('Setting zVariable "%s" to correctly-sized data with fill values.', zVariableName))
            nEpochRecords = size(ProcessData.Epoch, 1);
            [fillValue, ~] = get_fill_pad_values(DataObj, zVariableName);
            zVariableSize = [nEpochRecords, DataObj.data.(fn{1}).dim];
            zVariableData = cast(zeros(zVariableSize), matlabClass);
            zVariableData = bicas.utils.replace_value(zVariableData, 0, fillValue);
            
            DataObj.data.(zVariableName).data = zVariableData;
        else
            error('BICAS:execute_sw_mode:SWModeProcessing', logMsg)
        end
    end
end

%===========================================
% Write to CDF file using write_CDF_dataobj
%===========================================
outputFilename = FilenamingFunction(...
    datasetId, GlobalAttributesSubset.Test_Id, GlobalAttributesSubset.Provider, SETTINGS.get_fv('OUTPUT_CDF.DATA_VERSION'));
filePath = fullfile(outputFileParentDir, outputFilename);
irf.log('n', sprintf('Writing dataset CDF file: %s', filePath))
bicas.utils.write_CDF_dataobj( filePath, ...
    DataObj.GlobalAttributes, ...
    DataObj.data, ...
    DataObj.VariableAttributes, ...
    DataObj.Variables ...
    )

end







function [ProcessData, GlobalAttributes] = read_dataset_CDF(pdid, filePath)
% Read elementary input process data from a CDF file and convert it to a format suitable as a data_manager "process data".
% Copies all zVariables into fields of a regular structure.
%
%
% ARGUMENTS
% =========
% pdid
% filePath
%
%
% RETURN VALUES
% =============
% ProcessData      : Struct with one field per zVar (named after the zVar). The content of every such field equals the
%                    content of the corresponding zVar.
% GlobalAttributes : Struct returned from "dataobj".
%
%
% NOTE: Fill & pad values are replaced with NaN for numeric data types.
%       Other CDF data (attributes) are ignored.
% NOTE: Uses irfu-matlab's dataobj for reading the CDF file.

% NOTE: HK TIME_SYNCHRO_FLAG can be empty.

global SETTINGS CONSTANTS

irf.log('n', sprintf('pdid=%s', pdid))    % NOTE: irf.log adds the method name.

%===========
% Read file
%===========
irf.log('n', sprintf('Reading CDF file: "%s"', filePath))
do = dataobj(filePath);                 % do=dataobj, i.e. irfu-matlab's dataobj!!!



%=========================================================================
% Copy zVariables (only the data) into analogous fields in smaller struct
%=========================================================================
irf.log('n', 'Converting dataobj (CDF data structure) to PDV.')
ProcessData       = struct();
zVariableNameList = fieldnames(do.data);
%bicas.dm_utils.log_array('explanation')
for i = 1:length(zVariableNameList)
    zVariableName = zVariableNameList{i};
    zVariableData = do.data.(zVariableName).data;
    
    %=================================================================================================
    % Log data to be written to CDF file
    % ----------------------------------
    % NOTE: Log messages should reflect the values READ from file.
    %       Process data variables are logged separately by data_manager.
    %=================================================================================================
    %bicas.dm_utils.log_struct_arrays(zVariableName, zVariableData);
    
    %=================================================
    % Replace fill/pad values with NaN for FLOAT data
    %=================================================
    % QUESTION: How does/should this work with integer fields that should also be stored as integers internally?!!!
    %    Ex: ACQUISITION_TIME, Epoch.
    % QUESTION: How distinguish integer zVariables that could be converted to floats (and therefore use NaN)?
    if isfloat(zVariableData)
        [fillValue, padValue] = get_fill_pad_values(do, zVariableName);
        zVariableData = bicas.utils.replace_value(zVariableData, fillValue, NaN);
        zVariableData = bicas.utils.replace_value(zVariableData, padValue,  NaN);
    else
        % Disable?! Only print warning if finds fill value which is not replaced?
        %irf.log('w', sprintf('Can not handle replace fill/pad values for zVariable "%s" when reading "%s".', zVariableName, filePath))
    end
    
    ProcessData.(zVariableName) = zVariableData;
end



fileDatasetId          = do.GlobalAttributes.DATASET_ID{1};
fileSkeletonVersionStr = do.GlobalAttributes.Skeleton_version{1};
irf.log('n', sprintf('File: DATASET_ID       = "%s"', fileDatasetId))
irf.log('n', sprintf('File: Skeleton_version = "%s"', fileSkeletonVersionStr))



%===================================================
% ASSERTIONS: Check GlobalAttributes values
%===================================================
% NOTE: Does print file name since it has only been previously been logged as "notice".
%bicas.dm_utils.assert_unvaried_N_rows(processData);
InputInfo = bicas.utils.select_cell_array_structs(CONSTANTS.INPUTS_INFO_LIST, 'PDID', {pdid});
InputInfo = InputInfo{1};
bicas.utils.assert_strings_equal(...
    SETTINGS.get_fv('INPUT_CDF_ASSERTIONS.STRICT_DATASET_ID'), ...
    {fileDatasetId, InputInfo.DATASET_ID}, ...
    sprintf('The input CDF file''s stated DATASET_ID does not match the value expected for the S/W mode.\n    File: %s\n    ', filePath))
bicas.utils.assert_strings_equal(...
    SETTINGS.get_fv('INPUT_CDF_ASSERTIONS.STRICT_SKELETON_VERSION'), ...
    {fileSkeletonVersionStr, InputInfo.SKELETON_VERSION_STR}, ...
    sprintf('The input CDF file''s stated Skeleton_version does not match the value expected for the S/W mode.\n    File: (%s)\n    ', filePath))



GlobalAttributes = do.GlobalAttributes;   % Assign return value.

end







function filename = get_output_filename(datasetId, testId, provider, dataVersion)
% Function that decides the filename to use for any given output dataset CDF file
%
% NOTE: ROC-TST-GSE-NTT-00017, "Data format and metadata definition for the ROC-SGSE data", iss2,rev1, Section 3.4
% specifies a file naming convention. Note: The version number should be the "data version"!
%
% dataVersion : two-digit string

% PROPOSAL: Include date and time?
%
% NOTE: May need to know the global attributes of the master CDF! See file naming convention, Test_id.
% Some may be identical to dataset ID, but still..



%current_time = datestr(now, 'yyyy-mm-dd_hhMMss.FFF');
%filename = [datasetId, '_V', skeletonVersionStr, '_', current_time, '.cdf'];
%filename = [datasetId, '_V', skeletonVersionStr, '___OUTPUT.cdf'];

filename = [logical_file_id(datasetId, testId, provider, dataVersion), '.cdf'];
end







function logicalFileId = logical_file_id(datasetId, testId, provider, dataVersion)
% Construct a "Logical_file_id" as defined in the ROC DFMD, global attribute+file name convention.

global CONSTANTS

CONSTANTS.assert_dataset_ID(datasetId)
if ~ischar(dataVersion ) || length(dataVersion)~=2
    error('BICAS:execute_sw_mode:Assertion:IllegalArgument', 'Illegal dataVersion')
end

providerParts = strsplit(provider, '>');
logicalFileId = [datasetId, '_', testId, '_', providerParts{1}, '_V', dataVersion];
end







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

iZVariable=find(strcmp(do.Variables(:,1), zVariableName));
padValue = do.Variables{iZVariable, 9};
% Comments in "spdfcdfinfo.m" should indirectly imply that column 9 is pad values since the structure/array
% commented on should be identical.
end

