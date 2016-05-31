% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden, 2016
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
% This function expects all the arguments defined in the RCS ICD and
% possibly additional inoffical arguments.
%
% - The official CLI parameter syntax is defined in
%   ROC-TST-GSE-ICD-00023-LES, Iss02 Rev02, Section 3.2.
% - NOTE: The official parameter syntax ("Specific input parameters", i.e.
%   "combinations of s/w mode and input parameters) used by this
%   particular code must be in agreement with "roc_sw_descriptor.js".
% - NOTE: The parameter syntax may contain additional inofficial parameters,
%   which are useful for development/debugging.
%
% Parameters: [ --matlab ] <Syntax 1 or 2>
%    Syntax 1: --identification | --version | --help
%    Syntax 2: <S/W mode>
%              <Common and specific input parameters; arbitrary order>
%
% --matlab : Flag for working (developing) from within the MATLAB IDE. Prevents the application
%            from exiting MATLAB itself. Useful for debugging.
%            Inofficial and not specified by the RCS ICD.
%
% Common input parameters
%    --output <absolute path to directory>
%    --log    <absolute path to directory>   (Note: Ignored by this code; Log file not written here)
%    --config <absolute path to file>
% Specific input parameters
%    Specify input and output files and which are which.
%    Depends on the exact <S/W mode>. See RCS ICD and the "roc_sw_descriptor.js" for this software.
%
%
%
% NOTE: This function is the main MATLAB function, i.e. it is called by no other
% MATLAB code during regular use. It is intended to be wrapped in
% and called from non-MATLAB code, e.g. a bash script.
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-05-19).
%
% IMPLEMENTATION NOTE: The RCS ICD specifies tightly what should go to stdout.
%
% IMPLEMENTATION NOTE: Could be useful to print/log MATLABs regular printouts (e.g. when
%       omitting/forgetting semicolon at end of line). Could be useful to
%       temporarily merge that output with irf.log output and RCS-ICD specified output.
%



% PROBLEM: irf('check_path'); prints to stdout and that can not be
% redirected(?). There may be other irfu-matlab functions with the same problem.
%    NOTE: Should be same problem with matlab warnings.
%    PROPOSAL: Use evalc.
%       NOTE: How handle errors?
% PROPOSAL: Somehow explicitly separate RCS-ICD specified stdout output from other stdout output?
%    PROPOSAL: Print RCS-ICD output with prefix and have wrapper bash script separate it out?
%       PRO: Can catch log/warning/error messages from the wrapper script.
%          Ex: Not finding matlab, failed call to matlab, matlab not finding the matlab script etc.
%    PROPOSAL: Print RCS-ICD output to separate stream and have the wrapper bash script separate it out?
%       CON: There is no such thing as another stream (?) unless it is a file with fprintf.
%       CON: The wrapper bash script then has to parse the CLI arguments to find the log file and be able to send
%            filtered output to the log file.
%          PROPOSAL: Let the wrapper script "filter" out stdout and send NOTHING to the log file.
%
% PROPOSAL: Build in MATLAB version check (warning?). If not, at least document the intended version.
% PROPOSAL: Require extra argument when NOT working in MATLAB rather than the opposite?
% QUESTION: If log file already exists, overwrite, amend, or error?
% PROPOSAL: Set flag for MATLAB warnings. Disable?
%    NOTE: TN claims warnings are sent to stdout.
%
% PROPOSAL: Rewrite JSON writing code to generate a string which can then be passed on to code that "prints official"
% stdout which adds prefixes.
function biasrcs( varargin )

try

    UNKNOWN_ERROR_CODE = 1;
    ARGUMENT_ERROR_CODE = 100;
    OPERATION_NOT_IMPLEMENTED_ERROR_CODE = 101;
    
    % Among other things: Sets up paths to within irfu-matlab (excluding .git/).
    % NOTE: Prints to stdout. Can not deactivate!
    irf('check_path');
    
    irf.log('debug')      % Set log level.
    %irf.log('log_out', '~/temp/biasrcslog.txt')    % NOTE: Amends to any preexisting file.
    
    % TEST
    %irf.log('c', 'Critical-level log message') % Goes to stdOUT!
    %irf.log('w', 'Warning-level log message')  % Goes to stdOUT!
    %irf.log('n', 'Notice-level log message')   % Goes to stdOUT!
    %irf.log('d', 'Debug-level log message')    % Goes to stdOUT!
    
    
    
    arguments = varargin;
    for i = 1:length(arguments)
        irf.log('n', sprintf('Command line argument %i: %s', i, arguments{i}))
    end
    
    % Argument to permit, but ignore.
    if (length(arguments) >= 1) && (strcmp(arguments{1}, '--developer-bash-settings'))
        arguments = arguments(2:end);
    end
    
    % Set DO_NOT_EXIT_MATLAB flag based on CLI parameters.
    DO_NOT_EXIT_MATLAB = 0;
    if (length(arguments) >= 1) && (strcmp(arguments{1}, '--matlab'))
        DO_NOT_EXIT_MATLAB = 1;
        arguments = arguments(2:end);
    end
    
    
    
    if (length(arguments) < 1)
        
        error('%i No arguments found. At least one argument is required.', ARGUMENT_ERROR_CODE)
        
    elseif (length(arguments) == 1)
        if (strcmp(arguments{1}, '--identification'))
            
            print_identification()
            
        elseif (strcmp(arguments{1}, '--version'))
            
            print_version()
            
        elseif (strcmp(arguments{1}, '--help'))
            
            error('%i Operation not implemented.', OPERATION_NOT_IMPLEMENTED_ERROR_CODE)
            
        else
            
            error('%i Can not interpret CLI parameters.', OPERATION_NOT_IMPLEMENTED_ERROR_CODE)
        end
    else
    
        % CASE: length(arguments) >= 2)
        parse_calibration_arguments(arguments)
    end
    

    
catch exception
    
    % Parse exit code from exception message.
    % Assumes exception.message on format: <exit code><one whitespace><message string>.
    %
    % QUESTION: How handle errors due to not finding "irf.log"?
    
    try
        
        % NOTE: Does not seem to produce error for non-parsable strings, not even
        % empty strings. exit_code = [], if can not parse.
        temp = strsplit(exception.message);
        exit_code = str2num(temp{1});
        message = exception.message(length(temp{1}+2):end);
        if isempty(exit_code)
            exit_code = UNKNOWN_ERROR_CODE;
        end
        
        len = length(exception.stack);
        fprintf(2, 'Call stack:\n');
        if (~isempty(len))
            for i=1:len
                fprintf(2, '    %s, row %i,\n', exception.stack(i).name, exception.stack(i).line);
            end
        end
        
        fprintf(2, [message, '\n']);
        irf.log('critical', message);
        
        default_msg = sprintf('Exiting MATLAB application with exit code %i\n', exit_code);
        fprintf(2, default_msg);
        irf.log('critical', default_msg);
        
        
        if ~DO_NOT_EXIT_MATLAB
            quit(exit_code)
        end
        
        
    catch exception
        fprintf(2, 'Unknown error. Error in the MATLAB script''s error handling.\n');
        quit(1)
    end
end



%==========================================================================
% Parse arguments associated with a calibration run.
    function parameters = parse_calibration_arguments(arguments)
        if length(arguments) <= 2
            error('%i Found un-paired argument.', ARGUMENT_ERROR_CODE)
        end
        
        parameters = struct;
        
        common_parameter_names = {'output', 'config', 'log'};
        
        % Read first argument: Interpret as "mode" name.
        mode = arguments{1};
        arguments = arguments(2:end);
        mode_parameter_names = {};
        if (strcmp(mode, 'testmode1'))
            mode_parameter_names{end+1} = 'input';
        else
            error('%i Can not recognize software mode "%s".', ARGUMENT_ERROR_CODE, mode)
        end
        
        % Iterate over remaining arguments: Interpret as pairs of parameter name + parameter value.
        [parameter_names, parameter_values] = parse_argument_names_values(arguments);
        
        for i = 1:length(parameter_names)
            name = parameter_names{i};
            
            i_c = find(strcmp(name, common_parameter_names));   % c = common
            i_m = find(strcmp(name, mode_parameter_names));     % m = mode
            
            if ~isempty(i_c) || ~isempty(i_m)
                parameters.(name) = parameter_values{i};
            else
                error('%i Can not recognize parameter name "%s".', ARGUMENT_ERROR_CODE, name)
            end
        end
        
        common_parameter_names
        mode_parameter_names
    end
%==========================================================================
% Iterate over remaining arguments: Interpret as pairs of parameter name + parameter value.
% NOTE: Does not check for doubles.
    function [parameter_names, parameter_values] = parse_argument_names_values(arguments)
        parameter_names  = [];
        parameter_values = [];
        
        while (true)
            if length(arguments) == 0
                break
            elseif length(arguments) == 1
                error('%i Found flag "%s" but no corresponding flag value.', ARGUMENT_ERROR_CODE, arguments{1})
            end
            
            name  = arguments{1};
            value = arguments{2};
            arguments = arguments(3:end);
            
            % Check for flag name beginning with "--".
            if ~(length(name) >=2 && strcmp(name(1:2), '--'))
                error('%i Can not interpret argument as parameter name: "%s".', ARGUMENT_ERROR_CODE, flag_name)
            end
            name = name(3:end);
            
            parameter_names{end+1}  = name;
            parameter_values{end+1} = value;
        end
    end
%==========================================================================
%==========================================================================
    function print_version()
        descr = define_descriptor();
        stdoutprintf('Release version "%s"\n', descr.release.version)
    end
%==========================================================================
    function print_identification()
        descr = define_descriptor();
        str = JSON_object_str(descr);
        stdoutprintf(str);
    end
%==========================================================================
% Returns structure representing the JSON object represented in the file "roc_sw_descriptor.js"
% specified in the RCS ICD.
    function obj = define_descriptor()
        ERIK_P_G_JOHANSSON = 'Erik P G Johansson';

        obj = struct();
        
        identification.project     = 'ROC-SGSE';
        identification.name        = 'BIASRCS (temporary name)';
        identification.identifier  = 'ROC-SGSE-BIASRCS';    % Temporary
        identification.description = 'BIAS calibration software (temporary description)';
        obj.identification = identification;
        
        release.version =      '0.0.1';
        release.date =         '2016-05-19';
        release.author =       ERIK_P_G_JOHANSSON;
        release.contact =      'erik.johansson@irfu.se';
        release.institute =    'IRF-U';
        release.modification = 'None (Initial release)';
        obj.release = release;
        
        environment.executable = 'bin/biasrcs';
        obj.environment = environment;
        
        obj.modes = [];
        
        mode = [];
        mode.name = 'testmode1';
        mode.purpose =  'Mode 1 (temporary purpose description)';
        mode.inputs.input_SCI.identifier = 'ROC-SGSE_L2R_RPW-LFR-SBM1-CWF';
        mode.inputs.input_SCI.version    = '01';
        mode.inputs.input_HK.identifier  = 'ROC-SGSE_HK_RPW-BIA';
        mode.inputs.input_HK.version     = '01';
        
        mode.outputs.output_SCI.identifier  = 'ROC-SGSE_L2S_RPW-BIA-xxxxx';
        mode.outputs.output_SCI.name        = 'xxxxx (temporary name)';
        mode.outputs.output_SCI.description = 'Contains xxxxx (temporary description)';
        mode.outputs.output_SCI.level       = 'L2S';
        mode.outputs.output_SCI.release.author       = ERIK_P_G_JOHANSSON;
        mode.outputs.output_SCI.release.date         = '2016-05-19';
        mode.outputs.output_SCI.release.version      = '01';
        mode.outputs.output_SCI.release.contact      = 'erik.johansson@irfu.se';
        mode.outputs.output_SCI.release.institute    = 'IRF-U';
        mode.outputs.output_SCI.release.modification = 'None (initial release)';
        obj.modes{end+1} = mode;
        
        mode = [];
        mode.name = 'testmode2';
        mode.purpose =  'Mode 2 (temporary purpose description)';
        mode.inputs.input_SCI.identifier = 'ROC-SGSE_L2R_RPW-TDS-SURV-RSWF';
        mode.inputs.input_SCI.version    = '01';
        mode.outputs.output_SCI.identifier  = 'ROC-SGSE_L2S_RPW-BIA-xxxxx';
        mode.outputs.output_SCI.name        = 'xxxxx (temporary name)';
        mode.outputs.output_SCI.description = 'Contains xxxxx (temporary description)';
        mode.outputs.output_SCI.level       = 'L2S';
        mode.outputs.output_SCI.release.author       = ERIK_P_G_JOHANSSON;
        mode.outputs.output_SCI.release.date         = '2016-05-19';
        mode.outputs.output_SCI.release.version      = '01';
        mode.outputs.output_SCI.release.contact      = 'erik.johansson@irfu.se';
        mode.outputs.output_SCI.release.institute    = 'IRF-U';
        mode.outputs.output_SCI.release.modification = 'None (initial release)';
        obj.modes{end+1} = mode;
    end

end
