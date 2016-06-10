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
% --developer-bash-settings : Used by the wrapper bash script. Ignored by this code.
%
% <S/W mode> : 
% Common input parameters
%    --output <absolute path to directory>
%    --log    <absolute path to directory>   (Note: Ignored by this code; The log file not written from here)
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
% regular use. It is intended to be wrapped in and called from non-MATLAB code, e.g. a bash script.
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-06-02) but may very well work with other
% versions of MATLAB.
%
% NOTE: The code uses persistent variables. One may want to call "clear all" to clear these before
% calling the code again when working in MATLAB. (One can not put "clear all" in the code since that
% clears the function parameters.)
%
% IMPLEMENTATION NOTE: This code does not quit/exit using the MATLAB function "quit" since that
% always exits all of MATLAB which is undesirebly when developing in the MATLAB IDE. The function
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
% PROPOSAL: Have --help print the error codes.

init_global_constants
global ERROR_CODES
global REQUIRED_MATLAB_VERSION
error_code = ERROR_CODES.NO_ERROR;   % Default exit error code.
   
try
    % Among other things: Sets up paths to within irfu-matlab (excluding .git/).
    % NOTE: Prints to stdout. Can not deactivate this behaviour!
    irf('check_path');
    
    irf.log('debug')      % Set log level.
    
    arguments = varargin;
    for i = 1:length(arguments)
        irf.log('n', sprintf('CLI argument %2i: "%s"', i, arguments{i}))
    end
    irf.log('n', sprintf('Current working directory: "%s"', pwd));   % Useful for debugging the use of relative directory arguments.
    
    % Check MATLAB version
    found_version = version('-release');
    if ~strcmp(found_version, REQUIRED_MATLAB_VERSION)
        errorp(ERROR_CODES.MISC_ERROR, ...
            'Wrong MATLAB version. Found %s. Requires %s.\n', found_version, REQUIRED_MATLAB_VERSION)
    end

    % Derive the root path of the software (BICAS directory structure root).
    [matlab_src_path, ~, ~] = fileparts(mfilename('fullpath'));
    sw_root_path = get_abs_path([matlab_src_path, filesep, '..']);
    irf.log('n', sprintf('MATLAB source code path: "%s"', matlab_src_path))
    irf.log('n', sprintf('Software root path:      "%s"', sw_root_path))
    


    %=================
    % Parse arguments
    %=================
    
    % Flag to permit, but ignore.
    if (length(arguments) >= 1) && (strcmp(arguments{1}, '--developer-bash-settings'))
        arguments = arguments(2:end);
    end
    
    flags = containers.Map;
    flags('log_path') = struct('CLI_str', '--log', 'is_required', 0, 'expects_value', 1);    % Flag+value to permit but ignore.

    if (length(arguments) < 1)

        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Not enough arguments found.')

    elseif (strcmp(arguments{1}, '--identification'))

        [~] = parse_CLI_flags(arguments(2:end), flags);  % Check CLI syntax but ignore results.
        print_identification()

    elseif (strcmp(arguments{1}, '--version'))

        [~] = parse_CLI_flags(arguments(2:end), flags);  % Check CLI syntax but ignore results.
        print_version()

    elseif (strcmp(arguments{1}, '--help'))

        [~] = parse_CLI_flags(arguments(2:end), flags);  % Check CLI syntax but ignore results.
        print_help()

    else
        
        %=================================================
        % CASE: Must be a S/W mode (or invalid arguments)
        %=================================================

        flags('output_dir') = struct('CLI_str', '--output', 'is_required', 1, 'expects_value', 1);
        %flags('config_file_path') = struct('CLI_str', '--config', 'is_required', 0, 'expects_value', 1);

        % Figure out which S/W mode (which index among the constants).
        C = bicas_constants;
        temp = {};
        for i = 1:length(C.sw_modes)
            temp{end+1} = C.sw_modes{i}.CLI_parameter;
        end
        i_mode = find(strcmp(arguments{1}, temp));        
        if length(i_mode) ~= 1
            % NOTE: The message is slightly inaccurate. Argument "--version" etc. would have worked too.
            errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not interpret argument "%s" as a S/W mode.', arguments{1});
        end

        % Add CLI flags depending on the S/W mode.
        % ----------------------------------------
        C_inputs = C.sw_modes{i_mode}.inputs;    % C = Constants structure.
        input_keys = {};
        for i_input = 1:length(C_inputs)
            key = C_inputs{i_input}.CLI_parameter_name;    % Replace key with dataset_ID?!
            flag = [];
            flag.CLI_str = ['--', C_inputs{i_input}.CLI_parameter_name];
            flag.is_required = 1;
            flag.expects_value = 1;
            
            flags(key) = flag;
            input_keys{end+1} = key;
        end

        % Parse remaining arguments.
        parsed_flags = parse_CLI_flags(arguments(2:end), flags);
        
        input_files = containers.Map(input_keys, parsed_flags.values(input_keys));
        output_dir = get_abs_path(parsed_flags('output_dir'));
        
        irf.log('n', sprintf('Output directory = "%s"', output_dir));
        
        if ~exist(output_dir, 'dir')
            errorp(ERROR_CODES.PATH_NOT_FOUND, 'Output directory "%s" does not exist.', output_dir)
        end



        %===============================================================
        % TEST CODE
        % 
        % USING TEST CODE TO IMPLEMENT S/W MODES.
        % IGNORES INPUT FILES AND PRODUCES INVALID .CDF FILE AS OUTPUT.
        %===============================================================
        execute_sw_mode_TEST_IMPLEMENTATION(i_mode, output_dir, sw_root_path)
        %errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not completely implemented: Use S/W mode')
    end


    
catch exception
    
    try
        irf.log('critical', 'Main function caught an exception. Error handling. Exiting.');   % Print to stdout.
        
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
            error_code = ERROR_CODES.UNKNOWN_ERROR;       % NOTE: Set error code for when exiting MATLAB.
        else
            message = message(length(temp{1})+2:end);
        end
        
        %======================
        % Print the call stack
        %======================
        len = length(exception.stack);
        fprintf(2, 'MATLAB call stack:\n');
        if (~isempty(len))
            for i=1:len
                sc = exception.stack(i);   % sc = stack call
                temp = strsplit(sc.file, filesep);
                filename = temp{end};
                
                fprintf(2, '    %s (%s), row %i,\n', sc.name, filename, sc.line);
            end
        end
        
        fprintf(2, [message, '\n']);    % Print to stderr.
        
        default_msg = sprintf('Exiting MATLAB application with error code %i.\n', error_code);
        fprintf(2, default_msg);        % Print to stderr.
        
        return
        
    catch exception
        % CASE: There was an error in the error handling(!).
        
        % NOTE: Only use very, very error safe code here.
        msg = sprintf('Unknown error. Error in the MATLAB script''s error handling.\nException message: "%s"\n', ...
            exception.message');
        fprintf(2, msg);   % Print to stderr.
        
        error_code = ERROR_CODES.UNKNOWN_ERROR;   % Not even use hardcoded constant for error code?!!
        return
    end
end



end
%===================================================================================================
% 
% TEST IMPLEMENTATION
% Implements a S/W mode by simply creating nonsense cdf output files where expected.
% Code should satisfy the RCS ICD with this implementation.
%
% NOTE: Will overwrite output file. Not necessary desirable in a real implementation but is
% practical for testing.
%
function execute_sw_mode_TEST_IMPLEMENTATION(i_mode, output_dir, sw_root_path)

global ERROR_CODES

irf.log('c', 'USING TEST IMPLEMENTATION FOR S/W MODES. ONLY CREATES NONSENSE CDF FILES.')
        
C = bicas_constants;
C_mode = C.sw_modes{i_mode};
output_JSON = [];

% Iterate over OUTPUTS
for i = 1:length(C_mode.outputs)
    C_mode_output = C_mode.outputs{i};
    master_cdf_filename = C_mode_output.master_cdf_filename;
    output_filename = [C_mode_output.dataset_ID, '_', C_mode_output.dataset_version_str, '.cdf'];
    
    src_file  = fullfile(sw_root_path, C.master_cdfs_dir_rel, master_cdf_filename);
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

str = JSON_object_str(output_JSON);
stdout_printf(str);

end
%===================================================================================================
function print_version()

% IMPLEMENTATION NOTE: Uses the software version in the S/W descriptor rather than the in the BICAS
% constants since the RCS ICD specifies that it should be that version.

swd = get_sw_descriptor();
stdout_printf('Release version "%s"\n', swd.release.version)

end

%===================================================================================================
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07
%
% Print the JSON S/W descriptor.
%
function print_identification()

D = get_sw_descriptor();
str = JSON_object_str(D);
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
    stdout_printf('   %3i = %s\n', ERROR_CODES.(sfn{1}), sfn{1})
end

%errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not implemented: --help.')

end
