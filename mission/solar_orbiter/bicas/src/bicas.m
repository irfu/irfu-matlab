% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-03-xx
%
% BICAS = BIAS CAlibration Software (temporary name)
%
%
% IMPORTANT NOTE: The general interface that this software must comply with
% is described in the ROC-TST-GSE-ICD-00023-LES document (the "RCS ICD").
% Much documentation can thus be found there.
%
%
%
% INPUT / OUTPUT VARIABLES:
% -------------------------
% This function expects only string arguments corresponding to CLI arguments.
% This function expects the arguments defined in the RCS ICD and possibly additional inoffical
% arguments.
%
% - The official CLI parameter syntax is defined in
%   RCS ICD, Iss02 Rev02, Section 3.2.
% - NOTE: The official parameter syntax ("Specific input parameters", i.e.
%   "combinations of s/w mode and input parameters) used by this
%   particular code must be in agreement with "roc_sw_descriptor.js".
% - NOTE: The parameter syntax may contain additional inofficial parameters,
%   which are useful for development/debugging.
%
% SYNTAX/PARAMETERS:
% [--developer-bash-settings] ( --identification | --version | --help ) [--log <path>] | <S/W mode> <Common & specific input parameters>
%
% --developer-bash-settings : Used by the wrapper bash script. Permitted but ignored by this code.
%
% <S/W mode> : 
% Common input parameters
%    --output <absolute path to directory>
%    --log    <absolute path to directory>   (Note: Ignored by this code; The log file is not written to from here)
%    --config <absolute path to file>
% Specific input parameters
%    Specify input and output files and which are which.
%    Depends on the exact <S/W mode>.
%
% RETURN VALUE: error_code = The error code that is to be passed on to the OS/shell.
%
%
%
% NOTE: This function is the main MATLAB function, i.e. it is called by no other MATLAB code during
% regular use. It is intended to be wrapped in, and called from non-MATLAB code, e.g. a bash script.
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-06-02) but may very well work with other
% versions of MATLAB.
%
% NOTE: The code uses persistent variables. One may want to call "clear all" to clear these before
% calling the code again when working in MATLAB. (One can not put "clear all" in the code since that
% clears the function parameters.)
%
% IMPLEMENTATION NOTE: This code does not quit/exit using the MATLAB function "quit" since that
% always exits all of MATLAB which is undesirable when developing in the MATLAB IDE. The function
% returns the error code to make it possible for the bash wrapper to quit with an exit code instead.
%
% IMPLEMENTATION NOTE: The RCS ICD specifies tightly what should go to stdout. This code
% 1) prints all log messages to stdout, and
% 2) prints all messages intended for stdout to stdout but with a prefix so they can be filtered
% out by the calling bash script.
% Reasons: See the bash wrapper script.
%
function error_code = bicas( varargin )
%
% PROPOSAL: Set flag for MATLAB warnings. Disable?
%    NOTE: TN claims warnings are sent to stdout.
% PROPOSAL: Extra (inofficial) flag for setting the log level.
% QUESTION: Is the applicaton allowed to overwrite output files?
%


global ERROR_CODES
global REQUIRED_MATLAB_VERSION
global CONSTANTS                 % Initialized later
[ERROR_CODES, REQUIRED_MATLAB_VERSION] = error_safe_constants();

try
    % Among other things: Sets up paths to within irfu-matlab (excluding .git/).
    % NOTE: Prints to stdout. Can not deactivate this behaviour!
    irf('check_path');
    irf.log('debug')      % Set log level.
    
    % Check MATLAB version
    found_version = version('-release');
    if ~strcmp(found_version, REQUIRED_MATLAB_VERSION)
        errorp(ERROR_CODES.MISC_ERROR, ...
            'Wrong MATLAB version. Found %s. Requires %s.\n', found_version, REQUIRED_MATLAB_VERSION)
    end
    
    CONSTANTS = constants();   % Do not run initialization until the MATLAB version has been checked for. Could fail.
    
    % Log arguments
    for i = 1:length(varargin)
        irf.log('n', sprintf('CLI argument %2i: "%s"', i, varargin{i}))
    end
    irf.log('n', sprintf('Current working directory: "%s"', pwd));   % Useful for debugging the use of relative directory arguments.
    
    % Derive the root path of the software (BICAS directory structure root).
    [matlab_src_path, ~, ~] = fileparts(mfilename('fullpath'));
    sw_root_path = get_abs_path([matlab_src_path, filesep, '..']);
    CONSTANTS.SW_root_dir(sw_root_path)
    irf.log('n', sprintf('MATLAB source code path: "%s"', matlab_src_path))
    irf.log('n', sprintf('Software root path:      "%s"', sw_root_path))



    %================================================================
    % 1) Parse arguments, and
    % 2) select which of BICAS' different modes of operation.
    %================================================================
    arguments = varargin;
    
    % First argument: Flag to permit but ignore.
%     if (length(arguments) >= 1) && (strcmp(arguments{1}, '--developer-bash-settings'))
%         arguments = arguments(2:end);
%     end
    
    % Start configuring requirements on (remaining) arguments.
    flags = containers.Map;
    flags('log_path')         = struct('CLI_str', '--log',    'is_required', 0, 'expects_value', 1);   % Flag+value to permit but ignore.
    flags('config_file_path') = struct('CLI_str', '--config', 'is_required', 0, 'expects_value', 1);   % Flag+calue to permit but ignore.

    % Select mode of operations.
    if (length(arguments) < 1)

        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Not enough arguments found.')

    elseif (strcmp(arguments{1}, '--identification'))
        %============================
        % CASE: Print identification
        %============================
        [~] = parse_CLI_flags(arguments(2:end), flags);  % Check CLI syntax but ignore results.
        print_identification()

    elseif (strcmp(arguments{1}, '--version'))
        %============================
        % CASE: Print version
        %============================
        [~] = parse_CLI_flags(arguments(2:end), flags);  % Check CLI syntax but ignore results.
        print_version()

    elseif (strcmp(arguments{1}, '--help'))
        %============================
        % CASE: Print help
        %============================
        [~] = parse_CLI_flags(arguments(2:end), flags);  % Check CLI syntax but ignore results.
        print_help()

    else
        %==============================================
        % CASE: Should be a S/W mode (error otherwise)
        %==============================================

        try
            C_sw_mode = CONSTANTS.get_C_sw_mode_full(arguments{1});
        catch exception
            % NOTE: The message is slightly inaccurate. Argument "--version" etc. would have worked too.
            errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not interpret argument "%s" as a S/W mode.', arguments{1});
        end

        %-----------------------------------------------------------------
        % Configure CLI flags to expect, partly depending on the S/W mode
        %
        % NOTE: The flags are identified by strings (container.Map keys) which are a combination of the namespaces for
        % 1) identifiers for misc. flags e.g. "output_dir", "log_path".
        % 2) dataset IDs!
        %-----------------------------------------------------------------
        flags('output_dir') = struct('CLI_str', '--output', 'is_required', 1, 'expects_value', 1);
        C_inputs = C_sw_mode.inputs;    % C = Constants structure.
        input_process_data_types = {};   % List of keys used for input files.
        for i_input = 1:length(C_inputs)
            process_data_type = C_inputs{i_input}.process_data_type;
            
            % Configure one flag+value pair
            flag = [];
            flag.CLI_str = ['--', C_inputs{i_input}.CLI_parameter];
            flag.is_required = 1;
            flag.expects_value = 1;
            flags(process_data_type) = flag;
            
            input_process_data_types{end+1} = process_data_type;
        end

        %-----------------------------
        % Parse (remaining) arguments
        %-----------------------------
        parsed_flags = parse_CLI_flags(arguments(2:end), flags);
        
        
        
        input_files = containers.Map(input_process_data_types, parsed_flags.values(input_process_data_types));   % Extract subset of parsed arguments.
        
        output_dir = get_abs_path(parsed_flags('output_dir'));
        
        

        %===============================================================
        % EXECUTE S/W MODE
        %
        % CHOOSE IMPLEMENTATION TO USE.
        %===============================================================
        execute_sw_mode(C_sw_mode.CLI_parameter, input_files, output_dir)  % The intended real implementation
        %execute_sw_mode_TEST_IMPLEMENTATION(C_sw_mode.CLI_parameter, output_dir)   % IMPLEMENTATION FOR TESTING. OUTPUTS NONSENSE CDFs.
        %errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not completely implemented: Use S/W mode')
    end

        

    % EXIT
    error_code = ERROR_CODES.NO_ERROR;   % Default RETURN value.


    
catch exception
    
    try
        irf.log('critical', 'Main function caught an exception. Beginning error handling.');   % Print to stdout.
        
        %========================================================================================
        % Parse exception message - Extract and set error code
        % ----------------------------------------------------
        % ASSUMES: exception.message is on format: <Error code><Whitespace><Message string>.
        % errorp produces error messages on that format.
        % NOTE: Code does not seem to produce error for non-parsable strings, not even empty
        % strings. error_code = [], if can not parse.
        %========================================================================================
        message = exception.message;
        temp = strsplit(message);
        error_code = str2double(temp{1});                 % NOTE: Set error code for when exiting MATLAB.
        if isnan(error_code)
            error_code = ERROR_CODES.ERROR_IN_MATLAB_ERROR_HANDLING;       % NOTE: Set error code for when exiting MATLAB.
        else
            message = message(length(temp{1})+2:end);
        end
        
        %======================
        % Print the call stack
        %======================
        len = length(exception.stack);
        fprintf(2, 'MATLAB call stack:\n');    % Print to stderr.
        if (~isempty(len))
            for i=1:len
                sc = exception.stack(i);   % sc = stack call
                temp = strsplit(sc.file, filesep);
                filename = temp{end};
                
                fprintf(2, '    %-25s %-30s row %i,\n', [filename, ','], [sc.name, ','], sc.line);
            end
        end
        
        fprintf(2, [message, '\n']);    % Print to stderr.
        
        default_msg = sprintf('Exiting MATLAB application with error code %i.\n', error_code);
        fprintf(2, default_msg);        % Print to stderr.
        
        return
        
    catch exception
        %===================================================
        % CASE: There was an error in the error handling(!)
        %===================================================
        
        % NOTE: Only use very, very error safe code here.
        fprintf(2, 'Unknown error. Error in the MATLAB code''s error handling.\nException message: "%s"\n', ...
            exception.message');   % Print to stderr.
        
        error_code = ERROR_CODES.ERROR_IN_MATLAB_ERROR_HANDLING;   % Not even use hardcoded constant for this error?!!
        return
    end
end



end



%===================================================================================================
% TEST IMPLEMENTATION
% Implements a S/W mode by simply creating nonsense cdf output files where expected.
% Code should satisfy the RCS ICD with this implementation.
%
% NOTE: Will overwrite output file. Not necessary desirable in a real implementation but is
% practical for testing.
%
% ASSUMES
%===================================================================================================
function execute_sw_mode_TEST_IMPLEMENTATION(sw_mode_CLI_parameter, output_dir)

global ERROR_CODES CONSTANTS

irf.log('c', 'USING TEST IMPLEMENTATION FOR S/W MODES. ONLY CREATES NONSENSE CDF FILES.')

%C_mode = get_C_sw_mode(sw_mode_CLI_parameter);
temp = select_structs(C.sw_modes, 'CLI_parameter', {sw_mode_CLI_parameter});
C_mode = temp{1};
output_JSON = [];

% Iterate over OUTPUTS
for i = 1:length(C_mode.outputs)
    C_mode_output = C_mode.outputs{i};
    master_cdf_filename = C_mode_output.master_cdf_filename;
    output_filename = [C_mode_output.dataset_ID, '_', C_mode_output.dataset_version_str, '.cdf'];
    
    src_file  = fullfile(bias_constants.sw_root_dir(), C.master_cdfs_dir_rel, master_cdf_filename);
    dest_file = fullfile(output_dir, output_filename);
    
    irf.log('n', 'Trying to copy file')
    irf.log('n', sprintf('   from %s', src_file))
    irf.log('n', sprintf('   to   %s', dest_file))
    [success, copyfile_msg, ~] = copyfile(src_file, dest_file);   % Overwrites any pre-existing file.
    if ~success
        errorp(ERROR_CODES.MISC_ERROR, ...
            'Failed to copy file\n    from "%s"\n    to   "%s".\n"copyfile" error message: "%s"', ...
            src_file, dest_file, copyfile_msg)
    end
    
    output_JSON.(C_mode_output.JSON_output_file_identifier) = output_filename;
end

% Print list of files produced in the form of a JSON object.
str = JSON_object_str(output_JSON);
stdout_printf(str);

end


%===================================================================================================
function print_version()

% IMPLEMENTATION NOTE: Uses the software version in the S/W descriptor rather than the in the BICAS
% constants since the RCS ICD specifies that it should be that version.

swd = get_sw_descriptor();
stdout_printf('Version "%s"\n', swd.release.version)

end



%===================================================================================================
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07
%
% Print the JSON S/W descriptor.
%
function print_identification()
global CONSTANTS

D = get_sw_descriptor();
str = JSON_object_str(D, CONSTANTS.C.JSON_object_str);
stdout_printf(str);

end



%===================================================================================================
function print_help()
%
% PROPOSAL: Print error codes. Can use implementation to list them?
%    PROPOSAL: Define error codes with description strings?! Map?! Check for doubles?!
% PROPOSAL: Print CLI syntax incl. all modes?

global ERROR_CODES

D = get_sw_descriptor();
stdout_printf('%s\n%s\n', D.identification.name, D.identification.description)

stdout_printf('\nError codes (internal constants):\n')
for sfn = fieldnames(ERROR_CODES)'
    error_code = ERROR_CODES.(sfn{1});
    error_name = sfn{1};
    stdout_printf('   %3i = %s\n', error_code, error_name)
end

%errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not implemented: --help.')

end
