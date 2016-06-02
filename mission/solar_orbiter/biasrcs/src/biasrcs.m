% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-03-xx
%
% BIASRCS = BIAS RCS (temporary name)
%
% (RCS = RPW Calibration Software; abbrev. from ROC-TST-GSE-ICD-00023-LES)
%
%
%
% IMPORTANT NOTE: The general interface that this software must comply with
% is described in the ROC-TST-GSE-ICD-00023-LES document (the "RCS ICD").
% Much documentation can thus be found there.
%
%
%
% PARAMETERS:
% -----------
% This function expects all the arguments defined in the RCS ICD and possibly additional inoffical
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
% SYNTAX:
% [--matlab] ( --identification | --version | --help ) [--log <path>] | <S/W mode> <Common & specific input parameters>
%
% --matlab : Flag for working (developing) from within the MATLAB IDE. Prevents the application
%            from exiting MATLAB itself. Useful for debugging.
%            Inofficial and not specified by the RCS ICD.
% <S/W mode> : 
% Common input parameters
%    --output <absolute path to directory>
%    --log    <absolute path to directory>   (Note: Ignored by this code; The log file not written from here)
%    --config <absolute path to file>
% Specific input parameters
%    Specify input and output files and which are which.
%    Depends on the exact <S/W mode>.
%
%
%
% NOTE: This function is the main MATLAB function, i.e. it is called by no other MATLAB code during
% regular use. It is intended to be wrapped in and called from non-MATLAB code, e.g. a bash script.
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-06-02) but may very well work in other
% versions.
%
% IMPLEMENTATION NOTE: The RCS ICD specifies tightly what should go to stdout.
% This code
% 1) prints all log messages to stdout, and
% 2) prints all messages intentended for stdout to stdout but with a prefix so they can be filtered
% out by the calling bash script.
% Reasons: See the bash wrapper script.
%
function biasrcs( varargin )
%
% PROPOSAL: Require extra argument when NOT working in MATLAB rather than the opposite?
% PROPOSAL: Set flag for MATLAB warnings. Disable?
%    NOTE: TN claims warnings are sent to stdout.
% PROPOSAL: Extra (inofficial) flag for setting the log level.
%

init_global_constants
global ERROR_CODES
global REQUIRED_MATLAB_VERSION
   
try
    % Among other things: Sets up paths to within irfu-matlab (excluding .git/).
    % NOTE: Prints to stdout. Can not deactivate!
    irf('check_path');
    
    irf.log('debug')      % Set log level.
    
    arguments = varargin;
    for i = 1:length(arguments)
        irf.log('n', sprintf('CLI argument %i: "%s"', i, arguments{i}))
    end
    
    
    
    %=================
    % Parse arguments
    %=================
    
    % Argument to permit, but ignore.
    if (length(arguments) >= 1) && (strcmp(arguments{1}, '--developer-bash-settings'))
        arguments = arguments(2:end);
    end
    
    % Set DO_NOT_EXIT_MATLAB flag based on CLI parameter.
    % ----------------------------------------------------
    % NOTE: We want to find this argument early in the code.
    % If an error occurs before this point, MATLAB could quit
    % while developing in the MATLAB IDE which is what the flag is supposed to prevent.
    DO_NOT_EXIT_MATLAB = 0;
    if (length(arguments) >= 1) && (strcmp(arguments{1}, '--matlab'))
        DO_NOT_EXIT_MATLAB = 1;
        arguments = arguments(2:end);
    end
    
    % Check MATLAB version
    % --------------------
    % IMPLEMENTATION NOTE: Might want to do this after DO_NOT_EXIT_MATLAB has been set.    
    found_version = version('-release');
    if ~strcmp(found_version, REQUIRED_MATLAB_VERSION)
        errorp(ERROR_CODES.MISC_ERROR, 'Wrong MATLAB version. Found %s. Requires %s.\n', found_version, REQUIRED_MATLAB_VERSION)
    end
   
    
    
    opt_single_flags = {};
    opt_double_flags = {'--log'};
    
    if (length(arguments) < 1)

        errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Not enough arguments found.')

    elseif (strcmp(arguments{1}, '--identification'))

        [~] = parse_single_double_arguments(arguments(2:end), opt_single_flags, opt_double_flags);
        print_identification()

    elseif (strcmp(arguments{1}, '--version'))

        [~] = parse_single_double_arguments(arguments(2:end), opt_single_flags, opt_double_flags);
        print_version()

    elseif (strcmp(arguments{1}, '--help'))
            
        errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not implemented.')
            
    else
        
        % CASE: Must be a S/W mode (or invalid arguments).

        % Figure out which S/W mode (index among the constants).
        C = biasrcs_constants;
        temp = {};
        for i = 1:length(C.sw_modes)
            temp{end+1} = C.sw_modes{i}.CLI_parameter;
        end
        i_mode = find(strcmp(arguments{1}, temp));        
        if length(i_mode) ~= 1
            errorp(ERROR_CODES.CLI_ARGUMENT_ERROR, 'Can not interpret argument "%s" as S/W mode.', arguments{1});
        end
        
        % Add (mandatory) flags depending on the S/W mode.
        C_inputs = C.sw_modes{i_mode}.inputs;    % C = Constants structure.
        for i_input = 1:length(C_inputs)
            % ADD flags to look for.
            opt_double_flags{end+1} = ['--', C_inputs{i_input}.CLI_parameter_name];
        end
        [~] = parse_single_double_arguments(arguments(2:end), opt_single_flags, opt_double_flags);
        
        errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Operation not implemented.')
    end


    
catch exception
    
    % Parse exit code from exception message.
    % Assumes exception.message on format: <exit code><one whitespace><message string>.
    %
    % QUESTION: How handle errors due to not finding "irf.log"?
    
    try
        
        % NOTE: Does not seem to produce error for non-parsable strings, not even
        % empty strings. exit_code = [], if can not parse.
        message = exception.message;
        temp = strsplit(message);
        exit_code = str2double(temp{1});
        if isnan(exit_code)
            exit_code = ERROR_CODES.UNKNOWN_ERROR;
        else
            message = message(length(temp{1})+2:end);
        end
        
        len = length(exception.stack);
        fprintf(2, 'Call stack:\n');
        if (~isempty(len))
            for i=1:len
                sc = exception.stack(i);   % sc = stack call
                temp = strsplit(sc.file, filesep);
                filename = temp{end};
                
                %fprintf(2, '    %s, row %i,\n', sc.name, sc.line);
                fprintf(2, '    %s (%s), row %i,\n', sc.name, filename, sc.line);
            end
        end
        
        fprintf(2, [message, '\n']);
        irf.log('critical', message);
        
        default_msg = sprintf('Exiting MATLAB application with exit code %i.\n', exit_code);
        fprintf(2, default_msg);
        irf.log('critical', default_msg);
        
        
        if ~DO_NOT_EXIT_MATLAB
            quit(exit_code)
        end
        
        
    catch exception
        fprintf(2, 'Unknown error. Error in the MATLAB script''s error handling.\n');
        quit(1)    % Not even use hardcoded constant for exit code?!!
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
