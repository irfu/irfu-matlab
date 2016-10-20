% execute_sw_mode   Execute a "S/W mode" as (implicitly) specified by the CLI arguments.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% sw_mode_CLI_parameter : String
% input_files : containers.Map with
%    keys   = Dataset IDs
%    values = Paths to input files.
%
function execute_sw_mode(sw_mode_CLI_parameter, input_files, output_dir)
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

irf.log('n', sprintf('Output directory = "%s"', output_dir));       
if ~exist(output_dir, 'dir')
    error('BICAS:execute_sw_mode:Assertion:PathNotFound', 'Output directory "%s" does not exist.', output_dir)
end



DM = bicas.data_manager();



%===========================================================================
% Give all input CDF files (from the CLI argument list) to the data manager
%===========================================================================
input_PDTs = input_files.keys;
for i = 1:length(input_PDTs)
    PDT = input_PDTs{i};
    input_file = input_files(PDT);
    
    DM.set_elementary_input_CDF(PDT, input_file);
end



C_sw_mode = bicas.data_manager.get_C_sw_mode_full(sw_mode_CLI_parameter);



%==================================
% Iterate over all the output CDFs
%==================================
JSON_output_CDF_filenames = [];
for i = 1:length(C_sw_mode.outputs)
    C_output = C_sw_mode.outputs{i};
    
    output_filename = get_output_filename(C_output.dataset_ID, C_output.skeleton_version_str);
    output_file_path = fullfile(output_dir, output_filename);
    [parent_dir, basename, ext] = fileparts(output_file_path);
    output_file_path_old = fullfile(parent_dir, [basename, '_old', ext]);
    output_file_path_new = output_file_path;
    
    PDT = C_output.PDT;
    process_data = DM.get_process_data_recursively(PDT, C_sw_mode.ID);
    
    % Read master CDF file via spdfcdfread.
    master_CDF_path = bicas.get_master_CDF_path(C_output.dataset_ID, C_output.skeleton_version_str);
    master_spdf = [];
    [master_spdf.data, master_spdf.info] = spdfcdfread(master_CDF_path, 'Structure', 1, 'KeepEpochAsIs', 1);
    % Read master CDF file via dataobj.
    master_dataobj = dataobj(master_CDF_path);
    
    
    
    %==================================================
    % Iterate over all output process data field names
    %==================================================
    fn_list = fieldnames(process_data);
    for i = 1:length(fn_list)
        zVar_name = fn_list{i};
        irf.log('n', sprintf('Converting "%s" to CDF zVariable.', zVar_name))
        
        %===========================================
        % Add process field to CDF as a zVariable
        % (assuming there is one in the master CDF)
        %===========================================
        % NOTE: THIS CODE IS STILL INCOMPLETE.
        % 1) The data supplied to write_CDF/spdfcdfwrite is not complete.
        % 1a) The data from the master CDF does not contain attributes for zero-record zVariables in master_spdf.data, but it
        % does in master_spdf.info and that seems to be enough for write_CDF/spdfcdfwrite but it is not obvious that it should be so.
        % 1b) Not all zVariables are explicitly set (QUALITY_FLAG, QUALITY_BITMASK, DELTA_PLUS_MINUS).
        % 2) There should be a check/assertion that all zVariables expected to be set are also set.
        % 3) Using an index from master_spdf.info.Variables to set a variable in master_spdf.data(i_zVar).Data.
        % It is not obvious that this is correct although it seems to work.
        
        
        
        % Modify master_spdf.data, master_spdf.info
        % ----------------------------------
        master_zVar_names = master_spdf.info.Variables(:,1);
        i_zVar = find(strcmp(zVar_name, master_zVar_names));
        if isempty(i_zVar)
            error('BICAS:execute_sw_mode:Assertion:DatasetFormat', 'Can not find zVariable "%s" in master file.', zVar_name)
        elseif ~isempty(master_spdf.data(i_zVar).Data) || ~isempty(master_spdf.data(i_zVar).VariableName)
            error('BICAS:execute_sw_mode:SWModeProcessing', 'Trying to overwrite zero-record CDF zVariable "%s" copied from master CDF file. Expected empty zVariable in master CDF file.', zVar_name)
        end
        MATLAB_class = bicas.utils.convert_CDF_type_to_MATLAB_class(master_spdf.info.Variables(i_zVar,4), 'Permit MATLAB classes');
        zVar_data = process_data.(zVar_name);
        % Compensating for odd behaviour in spdfcdfread and hence in write_CDF.
        % This 
        if ~ischar(zVar_data) && ndims(zVar_data) == 3
            zVar_data = permute(zVar_data, [2,3,1]);
        end
        master_spdf.data(i_zVar).Data = cast(zVar_data, MATLAB_class);
        master_spdf.data(i_zVar).VariableName = zVar_name;
        
        
        
        % Modify master_dataobj.data
        % --------------------------
        MATLAB_class = bicas.utils.convert_CDF_type_to_MATLAB_class(master_dataobj.data.(zVar_name).type, 'Permit MATLAB classes');
        master_dataobj.data.(zVar_name).data = cast(process_data.(zVar_name), MATLAB_class);
    end



    % ASSERTION: master_spdf.data, master_spdf.info
    for i = 1:length(master_spdf.data)
        if isempty(master_spdf.data(i).Data)
            % Ugly error message. Uses two zVariable names since not sure if either is going to be available or correct.
            zVar_name1 = master_spdf.data(i).VariableName;
            zVar_name2 = master_spdf.info.Variables{i,1};
            error('BICAS:execute_sw_mode:SWModeProcessing', ...
                'Master CDF contains zVariable which has not been set (i.e. it has zero records). "%s", "%s".', ...
                zVar_name1, zVar_name2)
        end
    end
    % ASSERTION: master_dataobj
    for fn = fieldnames(master_dataobj.data)'
        if isempty(master_dataobj.data.(fn{1}).data)
            error('BICAS:execute_sw_mode:SWModeProcessing', ...
                'Master CDF contains zVariable which has not been set (i.e. it has zero records). "%s".', ...
                master_dataobj.data.(fn{1}).Name)
        end
    end



    bicas.utils.write_CDF_spdfcdfread(output_file_path_old, master_spdf.data, master_spdf.info)    % NOTE: No "fill_empty" option.
    bicas.utils.write_CDF_dataobj(output_file_path_new, ...
        master_dataobj.GlobalAttributes, ...
        master_dataobj.data, ...
        master_dataobj.VariableAttributes, ...
        master_dataobj.Variables ...
    )    % NOTE: No "fill_empty" option.
    JSON_output_CDF_filenames.(C_output.JSON_output_file_identifier) = output_filename;
    
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
function filename = get_output_filename(dataset_ID, skeleton_version_str)
    % PROPOSAL: Include date and time?
    % PROPOSAL: Move into constants.
    % QUESTION: Use the "Data_version" instead? ROC-TST-GSE-SPC-00017-LES, i1r3, Section 3.4/3.1 seems to imply this.
    
    %current_time = datestr(now, 'yyyy-mm-dd_hhMMss.FFF');
    %filename = [dataset_ID, '_V', skeleton_version_str, '_', current_time, '.cdf'];
    filename = [dataset_ID, '_V', skeleton_version_str, '___OUTPUT.cdf'];
end
