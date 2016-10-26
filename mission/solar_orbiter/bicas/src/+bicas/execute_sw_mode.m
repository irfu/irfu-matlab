% execute_sw_mode   Execute a "S/W mode" as (implicitly) specified by the CLI arguments.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% swModeCliParameter : String
% InputFileList      : containers.Map with
%    keys   = Dataset IDs
%    values = Paths to input files.
%
%
% IMPORTANT NOTE: As of 2016-10-26:
% ================================
% This function temporarily uses both an old way and a new way of writing CDFs in
% parallel for testing. The old way is meant to be phased out eventually.
%
function execute_sw_mode(swModeCliParameter, InputFileList, outputDir)
%
% QUESTION: How verify dataset ID and dataset version against constants?
%    NOTE: Need to read CDF first.
%    NOTE: Need S/W mode.
%
% PROPOSAL: Verify output zVariables against master CDF zVariable dimensions (accessible with dataobj, even for zero records).
%
% TODO: Function for determining output filename.
%
% QUESTION: What should be the relationship between data manager and S/W modes really?
% Should data manager check anything?
%

global CONSTANTS

irf.log('n', sprintf('Output directory = "%s"', outputDir));       
if ~exist(outputDir, 'dir')
    error('BICAS:execute_sw_mode:Assertion:PathNotFound', 'Output directory "%s" does not exist.', outputDir)
end



DM = bicas.data_manager();



%===========================================================================
% Give all input CDF files (from the CLI argument list) to the data manager
%===========================================================================
input_PDIDs = InputFileList.keys;
for i = 1:length(input_PDIDs)
    PDID = input_PDIDs{i};
    inputFile = InputFileList(PDID);
    
    DM.set_elementary_input_CDF(PDID, inputFile);
end



C_sw_mode = bicas.data_manager.get_C_sw_mode_full(swModeCliParameter);



%==================================
% Iterate over all the output CDFs
%==================================
JSON_output_CDF_filenames = [];
for i = 1:length(C_sw_mode.outputs)
    C_output = C_sw_mode.outputs{i};
    
    outputFilename = get_output_filename(C_output.dataset_ID, C_output.skeleton_version_str);
    outputFilePath = fullfile(outputDir, outputFilename);
    [parentDir, basename, ext] = fileparts(outputFilePath);
    outputFilePathOld = fullfile(parentDir, [basename, '_old', ext]);
    outputFilePathNew = outputFilePath;
    
    PDID = C_output.PDID;
    processData = DM.get_process_data_recursively(PDID, C_sw_mode.ID);
    
    % Read master CDF file via spdfcdfread.
    masterCdfPath = bicas.get_master_CDF_path(C_output.dataset_ID, C_output.skeleton_version_str);
    MasterSpdf = [];
    [MasterSpdf.data, MasterSpdf.info] = spdfcdfread(masterCdfPath, 'Structure', 1, 'KeepEpochAsIs', 1);
    % Read master CDF file via dataobj.
    MasterDataobj = dataobj(masterCdfPath);
    
    
    
    %==================================================
    % Iterate over all output process data field names
    %==================================================
    fieldNameList = fieldnames(processData);
    for i = 1:length(fieldNameList)
        zVariableName = fieldNameList{i};
        irf.log('n', sprintf('Converting "%s" to CDF zVariable.', zVariableName))
        
        %===========================================
        % Add process field to CDF as a zVariable
        % (assuming there is one in the master CDF)
        %===========================================
        % NOTE: THIS CODE IS STILL INCOMPLETE.
        % 1) The data supplied to write_CDF/spdfcdfwrite is not complete.
        % 1a) The data from the master CDF does not contain attributes for zero-record zVariables in MasterSpdf.data, but it
        % does in MasterSpdf.info and that seems to be enough for write_CDF/spdfcdfwrite but it is not obvious that it should be so.
        % 1b) Not all zVariables are explicitly set (QUALITY_FLAG, QUALITY_BITMASK, DELTA_PLUS_MINUS).
        % 2) There should be a check/assertion that all zVariables expected to be set are also set.
        % 3) Using an index from MasterSpdf.info.Variables to set a variable in MasterSpdf.data(i_zVar).Data.
        % It is not obvious that this is correct although it seems to work.
        
        
        
        % Modify MasterSpdf.data, MasterSpdf.info
        % -----------------------------------------
        masterZVariableNameList = MasterSpdf.info.Variables(:,1);
        i_zVar = find(strcmp(zVariableName, masterZVariableNameList));
        if isempty(i_zVar)
            error('BICAS:execute_sw_mode:Assertion:DatasetFormat', 'Can not find zVariable "%s" in master file.', zVariableName)
        elseif ~isempty(MasterSpdf.data(i_zVar).Data) || ~isempty(MasterSpdf.data(i_zVar).VariableName)
            error('BICAS:execute_sw_mode:SWModeProcessing', 'Trying to overwrite zero-record CDF zVariable "%s" copied from master CDF file. Expected empty zVariable in master CDF file.', zVariableName)
        end
        MATLAB_class = bicas.utils.convert_CDF_type_to_MATLAB_class(MasterSpdf.info.Variables(i_zVar,4), 'Permit MATLAB classes');
        zVariableData = processData.(zVariableName);
        % Compensating for odd behaviour in spdfcdfread and hence in write_CDF.
        % This 
        if ~ischar(zVariableData) && ndims(zVariableData) == 3
            zVariableData = permute(zVariableData, [2,3,1]);
        end
        MasterSpdf.data(i_zVar).Data = cast(zVariableData, MATLAB_class);
        MasterSpdf.data(i_zVar).VariableName = zVariableName;
        
        
        
        % Modify MasterDataobj.data
        % --------------------------
        MATLAB_class = bicas.utils.convert_CDF_type_to_MATLAB_class(MasterDataobj.data.(zVariableName).type, 'Permit MATLAB classes');
        MasterDataobj.data.(zVariableName).data = cast(processData.(zVariableName), MATLAB_class);
    end



    % ASSERTION: MasterSpdf.data, MasterSpdf.info
    for i = 1:length(MasterSpdf.data)
        if isempty(MasterSpdf.data(i).Data)
            % Ugly error message. Uses two zVariable names since not sure if either is going to be available or correct.
            zVariableName1 = MasterSpdf.data(i).VariableName;
            zVariableName2 = MasterSpdf.info.Variables{i,1};
            error('BICAS:execute_sw_mode:SWModeProcessing', ...
                'Master CDF contains zVariable which has not been set (i.e. it has zero records). "%s", "%s".', ...
                zVariableName1, zVariableName2)
        end
    end
    % ASSERTION: MasterDataobj
    for fn = fieldnames(MasterDataobj.data)'
        if isempty(MasterDataobj.data.(fn{1}).data)
            error('BICAS:execute_sw_mode:SWModeProcessing', ...
                'Master CDF contains zVariable which has not been set (i.e. it has zero records). "%s".', ...
                MasterDataobj.data.(fn{1}).Name)
        end
    end



    bicas.utils.write_CDF_spdfcdfread(outputFilePathOld, MasterSpdf.data, MasterSpdf.info)    % NOTE: No "fill_empty" option.
    bicas.utils.write_CDF_dataobj(outputFilePathNew, ...
        MasterDataobj.GlobalAttributes, ...
        MasterDataobj.data, ...
        MasterDataobj.VariableAttributes, ...
        MasterDataobj.Variables ...
    )    % NOTE: No "fill_empty" option.
    JSON_output_CDF_filenames.(C_output.JSON_output_file_identifier) = outputFilename;
    
end



%============================================================
% Print JSON object describing the created file(s) to stdout
% ----------------------------------------------------------
% Required by the RCS ICD iss2rev2, section 3.3.
%============================================================
str = bicas.utils.JSON_object_str(JSON_output_CDF_filenames, CONSTANTS.C.JSON_object_str);
bicas.stdout_disp(str);

end



%=======================================================================
% This function decides what filename to use for any given output file.
%=======================================================================
function filename = get_output_filename(datasetId, skeleton_version_str)
    % PROPOSAL: Include date and time?
    % PROPOSAL: Move into constants.
    % QUESTION: Use the "Data_version" instead? ROC-TST-GSE-SPC-00017-LES, i1r3, Section 3.4/3.1 seems to imply this.
    
    %current_time = datestr(now, 'yyyy-mm-dd_hhMMss.FFF');
    %filename = [datasetId, '_V', skeleton_version_str, '_', current_time, '.cdf'];
    filename = [datasetId, '_V', skeleton_version_str, '___OUTPUT.cdf'];
end
