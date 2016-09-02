% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-09
%
% Execute a "S/W mode" as (implicitly) specified by the CLI arguments.
%
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
% PROPOSAL: In principle, if using the MMS data manager concept, then the code does not need
%       to explicitly check for input files, only ask for a set of datasets.
%    PROBLEM: How check input dataset version? Need to look up S/W mode and constants.
%       Could in principle be different for different modes(?!)
%       QUESTION: Check where? In data_manager.set_input? In execute_sw_mode.
%          data_manager does logically need some knowledge of dataset versions since its algorithms
%          are tailored for different dataset versions.
%       --
%       PROPOSAL: Set input data not one CDF at a time, but one S/W mode at a time.
%
% TODO: Function for determining output filename.
%
% QUESTION: What should be the relationship between data manager and S/W modes really?
% Should data manager check anything?
%

global ERROR_CODES CONSTANTS

irf.log('n', sprintf('Output directory = "%s"', output_dir));       
if ~exist(output_dir, 'dir')
    errorp(ERROR_CODES.PATH_NOT_FOUND, 'Output directory "%s" does not exist.', output_dir)
end



DM = data_manager();



%===========================================================================
% Give all input CDF files (from the CLI argument list) to the data manager
%===========================================================================
input_process_data_types = input_files.keys;
for i = 1:length(input_process_data_types)
    process_data_type = input_process_data_types{i};
    input_file = input_files(process_data_type);
    
    DM.set_elementary_input_CDF(process_data_type, input_file);
end



C_sw_mode = CONSTANTS.get_C_sw_mode_full(sw_mode_CLI_parameter);



%==================================
% Iterate over all the output CDFs
%==================================
JSON_filenames_obj = [];
for i = 1:length(C_sw_mode.outputs)
    C_output = C_sw_mode.outputs{i};
    
    output_filename = get_output_filename(C_output.dataset_ID, C_output.skeleton_version_str);
    output_file_path = fullfile(output_dir, output_filename);
    
    JSON_filenames_obj.(C_output.JSON_output_file_identifier) = output_filename;
    
    process_data_type = C_output.process_data_type;
    process_data = DM.get_process_data_recursively(process_data_type, C_sw_mode.ID);
    
    % Read master CDF file.
    master_CDF_path = get_master_CDF_path(C_output.dataset_ID, C_output.skeleton_version_str);
    [master_data, master_info] = spdfcdfread(master_CDF_path, 'Structure', 1, 'KeepEpochAsIs', 1);
    
    
    
    %===============================================================
    % Iterate over all output process data field names (zVariables)
    %===============================================================
    fn_list = fieldnames(process_data);
    for i = 1:length(fn_list)
        zVar_name = fn_list{i};
        master_zVar_names = master_info.Variables(:,1);
        
        %===========================================
        % Add process field to CDF as zVariable
        % (assuming there is one in the master CDF)
        %===========================================
        % NOTE: THIS CODE IS STILL INCOMPLETE.
        % 1) The data supplied to write_CDF/spdfcdfwrite is not complete.
        % 1a) The data from the master CDF does not contain attributes for zero-record zVariables in master_data, but it
        % does in master_info and that seems to be enough for write_CDF/spdfcdfwrite but it is not obvious.
        % 1b) Not all zVariables are explicitly set (QUALITY_FLAG, QUALITY_BITMASK, DELTA_PLUS_MINUS) and up with one
        % record.
        % 2) There should be a check/assertion that all zVariables expected to be set are also set.
        % 3) Using an index from master_info.Variables to set a variable in master_data(i_zVar).Data.
        % It is not obvious that this is correct although it seems to work.
        i_zVar = find(strcmp(zVar_name, master_zVar_names));
        if isempty(i_zVar)
            errorp(ERROR_CODES.CDF_ERROR, 'Can not find zVariable "%s" in master file.', zVar_name)
        elseif ~isempty(master_data(i_zVar).Data) || ~isempty(master_data(i_zVar).VariableName)
            errorp(ERROR_CODES.CDF_ERROR, 'Trying to overwrite zero-record CDF zVariable "%s". Expected empty zVariable.', zVar_name)
        end
        MATLAB_class = convert_CDF_type_to_MATLAB_class(master_info.Variables(i_zVar,4), 'Permit MATLAB classes');
        master_data(i_zVar).Data = cast(process_data.(zVar_name), MATLAB_class);
        master_data(i_zVar).VariableName = zVar_name;
    end

    for i = 1:length(master_data)
        if isempty(master_data(i).Data)
            % Ugly error message. Uses two zVariable names since not sure if either is going to be available or correct.
            zVar_name1 = master_data(i).VariableName;
            zVar_name2 = master_info.Variables{i,1};
            errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Master CDF contains zVariable which has not been set. "%s", "%s".', zVar_name1, zVar_name2)
        end
    end



    write_CDF(output_file_path, master_data, master_info)    % NOTE: No "fill_empty" option.
    
end



%============================================================
% Print JSON object describing the created file(s) to stdout
% ----------------------------------------------------------
% Required by the RCS ICD iss2rev2, section 3.3.
%============================================================
str = JSON_object_str(JSON_filenames_obj, CONSTANTS.C.JSON_object_str);
stdout_disp(str);

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
