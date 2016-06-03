% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-03-xx
%
% BIASRCS = BIAS RCS (temporary name)
%
% (RCS = RPW Calibration Software; abbrev. from ROC-TST-GSE-ICD-00023-LES)
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
% NOTE: This function is the main MATLAB function, i.e. it is called by no other MATLAB code during
% regular use. It is intended to be wrapped in and called from non-MATLAB code, e.g. a bash script.
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-06-02) but may very well work in other
% versions.
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
function error_code = biasrcs( varargin )
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
    
    [matlab_src_path, ~, ~] = fileparts(mfilename('fullpath'));
    irf.log('n', sprintf('MATLAB source code path: "%s"', matlab_src_path))
    
    
    
    %=================
    % Parse arguments
    %=================
    
    % Argument to permit, but ignore.
    if (length(arguments) >= 1) && (strcmp(arguments{1}, '--developer-bash-settings'))
        arguments = arguments(2:end);
    end
    
    % Check MATLAB version
    found_version = version('-release');
    if ~strcmp(found_version, REQUIRED_MATLAB_VERSION)
        errorp(ERROR_CODES.MISC_ERROR, ...
            'Wrong MATLAB version. Found %s. Requires %s.\n', found_version, REQUIRED_MATLAB_VERSION)
    end


    
    flags = [];
    flags(1).return_field_name = 'log_path';
    flags(1).CLI_name = '--log';
    flags(1).is_required = 0;
    flags(1).expects_value = 1;
    
    if (length(arguments) < 1)

        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Not enough arguments found.')

    elseif (strcmp(arguments{1}, '--identification'))

        [~] = parse_CLI_flags(arguments(2:end), flags);
        print_identification()

    elseif (strcmp(arguments{1}, '--version'))

        [~] = parse_CLI_flags(arguments(2:end), flags);
        print_version()

    elseif (strcmp(arguments{1}, '--help'))
            
        errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not implemented: --help.')
            
    else
        
        %=================================================
        % CASE: Must be a S/W mode (or invalid arguments)
        %=================================================
        
        flag = [];
        flag.return_field_name = 'output_dir';
        flag.CLI_name = '--output';
        flag.is_required = 1;
        flag.expects_value = 1;
        flags(end+1) = flag;

%         flag = [];
%         flag.return_field_name = 'config_file_path';
%         flag.CLI_name = '--config';
%         flag.is_required = 0;
%         flag.expects_value = 1;
%         flags(end+1) = flag;

        % Figure out which S/W mode (which index among the constants).
        C = biasrcs_constants;
        temp = {};
        for i = 1:length(C.sw_modes)
            temp{end+1} = C.sw_modes{i}.CLI_parameter;
        end
        i_mode = find(strcmp(arguments{1}, temp));        
        if length(i_mode) ~= 1
            errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not interpret argument "%s" as S/W mode.', arguments{1});
        end
        
        % Add CLI flags depending on the S/W mode.
        C_inputs = C.sw_modes{i_mode}.inputs;    % C = Constants structure.
        for i_input = 1:length(C_inputs)
            % ADD flags to look for.
            flag = [];
            flag.return_field_name = C_inputs{i_input}.CLI_parameter_name;
            flag.CLI_name = ['--', flag.return_field_name];
            flag.is_required = 1;
            flag.expects_value = 1;
            flags(end+1) = flag;
        end
        
        % Parse remaining arguments.
        parsed_flags = parse_CLI_flags(arguments(2:end), flags);
        
        
        
        %===============================================================
        % TEST CODE
        % 
        % USING TEST CODE TO IMPLEMENT S/W MODES.
        % IGNORES INPUT FILES AND PRODUCES INVALID .CDF FILE AS OUTPUT.
        %===============================================================
        warning('S/W modes are not implemented yet. Using test implementation.')   % Use irf.log warning/critical?
        C_mode = C.sw_modes{i_mode};
        for i = 1:length(C_mode.outputs)
            C_mode_output = C_mode.outputs{i};
            master_cdf_filename = C_mode_output.master_cdf_filename;
            output_dir = parsed_flags.output_dir.value;
            irf.log('n', sprintf('Output directory = "%s"', output_dir));
            
            % NOTE: Good to check existence of directory in particular for relative paths, since the
            % bash wrapper might change the directory.
            if ~exist(output_dir, 'dir')
                errorp(ERROR_CODES.PATH_NOT_FOUND, 'Output directory "%s" does not exist.', output_dir)
            end
            
            output_filename = [C_mode_output.dataset_ID, '_', C_mode_output.dataset_version_str, '.cdf'];
            
            % TEST!!
            %src_file = ['~/work_files/SOLAR_ORBITER/irfu-matlab/mission/solar_orbiter/biasrcs/data/', master_cdf_filename];
            src_file  = fullfile(matlab_src_path, C.master_cdfs_dir_rel, master_cdf_filename);
            dest_file = fullfile(output_dir, output_filename);
            
            
            [success, copyfile_msg, ~] = copyfile(src_file, dest_file);
            if ~success
                errorp(ERROR_CODES.MISC_ERROR, ...
                    'Failed to copy file\n    from "%s"\n    to   "%s".\n"copyfile" error message: "%s"', ...
                    src_file, dest_file, copyfile_msg)
            end
        end
        
        
        
        %errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not completely implemented: Use S/W mode')
    end


    
catch exception
    
    % Parse error code from exception message.
    % ASSUMES: exception.message is on format: <error code><one whitespace><message string>.
    
    try
        
        % NOTE: Does not seem to produce error for non-parsable strings, not even
        % empty strings. error_code = [], if can not parse.
        message = exception.message;
        temp = strsplit(message);
        error_code = str2double(temp{1});
        if isnan(error_code)
            error_code = ERROR_CODES.UNKNOWN_ERROR;
        else
            message = message(length(temp{1})+2:end);
        end
        
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
        
        fprintf(2, [message, '\n']);
        irf.log('critical', message);
        
        default_msg = sprintf('Exiting MATLAB application with error code %i.\n', error_code);
        fprintf(2, default_msg);
        irf.log('critical', default_msg);
        
        return
        
    catch exception
        % CASE: There was an error in the error handling(!).
        
        % NOTE: Only use very, very error safe code here.
        msg = sprintf('Unknown error. Error in the MATLAB script''s error handling.\nException message: "%s"\n', ...
            exception.message');
        fprintf(1, msg);
        fprintf(2, msg);
        
        error_code = ERROR_CODES.UNKNOWN_ERROR;   % Not even use hardcoded constant for error code?!!
        return
    end
end



end
%===================================================================================================
function print_version()

swd = get_sw_descriptor();
stdout_printf('Release version "%s"\n', swd.release.version)

end
%===================================================================================================
function print_identification()

descr = get_sw_descriptor();
str = JSON_object_str(descr);
stdout_printf(str);

end
