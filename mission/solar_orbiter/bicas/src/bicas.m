% errorCode = bicas( varargin )   Main function that launches BICAS.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-03-xx
%
% BICAS = BIAS CAlibration Software
%
%
%
% IMPORTANT NOTE: INFORMATION IMPORTANT FOR JUST USING THIS CODE CAN BE FOUND IN "readme.txt" AND OTHER DOCUMENTATION
% TEXT FILES (*.txt). To prevent the duplication of documentation, comments in this file tries to only cover subjects
% important for understanding the implementation and information not already present in other documentation text files.
%
%
%
% This function is the main MATLAB function, i.e. it is called by no other MATLAB code during regular use. It is
% intended to be wrapped in, and called from non-MATLAB code, e.g. a bash script.
%
%
% ASSUMES: The current file is in the <BICAS>/src directory.
%
%
% ARGUMENTS AND RETURN VALUE:
% ---------------------------
% This function expects exactly the CLI arguments submitted to the bash launcher script. This function therefore expects
% the arguments defined in the RCS ICD and possibly additional inoffical arguments.
%
% RETURN VALUE: errorCode = The error code that is to be passed on to the OS/shell.
%
% Notes:
% - The official parameter syntax for S/W modes must be in agreement with "roc_sw_descriptor.js" as specified by the RCS ICD.
% - The parameter syntax may contain additional inofficial parameters, which are useful for development/debugging, but which are still compatible with the RCS ICD.
% - The (MATLAB) code ignore but permits the CLI flags --log and --config.
%
%
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-06-02) but may very well work with other
% versions of MATLAB.
%
% IMPLEMENTATION NOTE: This code does not quit/exit using the MATLAB function "quit" since that
% always exits all of MATLAB which is undesirable when developing in the MATLAB IDE. The function
% returns the error code to make it possible for the bash wrapper to quit with an exit code instead.
%
% IMPLEMENTATION NOTE: The RCS ICD specifies tightly what should go to stdout. This code
% 1) prints all log messages to stdout, and
% 2) prints all messages intended for the final stdout to stdout but with a prefix so they can be filtered
% out by the calling wrapper bash script.
% Reasons: See the bash wrapper script.
%
function errorCode = bicas( varargin )
%
% PROPOSAL: Set flag for MATLAB warnings. Disable?
%    NOTE: TN claims warnings are sent to stdout.
% PROPOSAL: Extra (inofficial) flag for setting the log level.
% QUESTION: Is the applicaton allowed to overwrite output files?
%
% PROPOSAL: Check that all master cdf files are present/available.
%
% PROPOSAL: Rename to "bicas_main.m", or bicas.main (+bicas/main.m).
% PROPOSAL: Move prescribed MATLAB version to the config file.
% PROPOSAL: Put a summarized version of CLI syntax in "bicas --help" (somethinger easier that the S/W descriptor).
%    PRO: Useful when S/W descriptor becomes big and complex.
%
% PROPOSAL: Do not print just one error message based on msgID, pick several possible one.
%   CON: Does not match with picking exactly one error code to return.



global CONSTANTS                 % Gobal structure "CONSTANTS" is initialized later.
[ERROR_CODES, REQUIRED_MATLAB_VERSION] = bicas.error_safe_constants();



try
    
    startTimeTicSeconds = tic;
    
    % Among other things: Sets up paths to within irfu-matlab (excluding .git/).
    % NOTE: Prints to stdout. Can not deactivate this behaviour!
    irf('check_path');
    irf.log('notice')      % Set initial log level value until it is later overridden by the Config value.
    
    %======================
    % Check MATLAB version
    %======================
    matlabVersionString = version('-release');
    if ~strcmp(matlabVersionString, REQUIRED_MATLAB_VERSION)
        error('BICAS:BadMATLABVersion', 'Using bad MATLAB version. Found version "%s". BICAS requires version "%s".\n', ...
            matlabVersionString, REQUIRED_MATLAB_VERSION)
    end
    
    %=======================================================================
    % Derive the root path of the software (BICAS directory structure root)
    %=======================================================================
    % ASSUMES: The current file is in the <BICAS>/src directory.
    [matlabSrcPath, ~, ~] = fileparts(mfilename('fullpath'));
    bicasRootPath = bicas.utils.get_abs_path(fullfile(matlabSrcPath, '..'));
    
    %=============================
    % Initialize global constants
    %=============================
    % NOTE: Does not initialize CONSTANTS until the MATLAB version has been checked for. The init. code could otherwise
    % fail.
    CONSTANTS = bicas.constants(bicasRootPath);
    irf.log(CONSTANTS.C.IRF_LOG_LEVEL);
    
    %=======================================
    % Log misc. paths and all CLI arguments
    %=======================================
    irf.log('n', sprintf('BICAS software root path:      "%s"', bicasRootPath))
    irf.log('n', sprintf('BICAS MATLAB source code path: "%s"', matlabSrcPath))
    irf.log('n', sprintf('Current working directory:     "%s"', pwd));   % Useful for debugging the use of relative directory arguments.
    for i = 1:length(varargin)
        irf.log('n', sprintf('CLI argument %2i: "%s"', i, varargin{i}))    % PROPOSAL: Combine into a single multiline log message?
    end



    %================================================================
    % 1) Parse arguments, and
    % 2) select which of BICAS' different modes of operation.
    %================================================================
    cliArgumentsArray = varargin;
    
    % Start configuring requirements on (remaining) arguments.
    FlagsConfigMap = containers.Map;
    FlagsConfigMap('log_path')         = struct('cliString', '--log',    'isRequired', 0, 'expectsValue', 1);   % Flag+value to permit but ignore.
    FlagsConfigMap('config_file_path') = struct('cliString', '--config', 'isRequired', 0, 'expectsValue', 1);   % Flag+calue to permit but ignore.

    % Select mode of operations.
    if (length(cliArgumentsArray) < 1)

        error('BICAS:CLISyntax', 'Not enough arguments found.')

    elseif (strcmp(cliArgumentsArray{1}, '--identification'))
        %============================
        % CASE: Print identification
        %============================
        [~] = bicas.utils.parse_CLI_flags(cliArgumentsArray(2:end), FlagsConfigMap);  % Check CLI syntax but ignore results.
        print_identification()

    elseif (strcmp(cliArgumentsArray{1}, '--version'))
        %============================
        % CASE: Print version
        %============================
        [~] = bicas.utils.parse_CLI_flags(cliArgumentsArray(2:end), FlagsConfigMap);  % Check CLI syntax but ignore results.
        print_version()

    elseif (strcmp(cliArgumentsArray{1}, '--help'))
        %============================
        % CASE: Print help
        %============================
        [~] = bicas.utils.parse_CLI_flags(cliArgumentsArray(2:end), FlagsConfigMap);  % Check CLI syntax but ignore results.
        print_help(ERROR_CODES)

    else
        %==============================================
        % CASE: Should be a S/W mode (error otherwise)
        %==============================================
        DataManager = bicas.data_manager();
        %try
            C_sw_mode = DataManager.get_C_sw_mode_full(cliArgumentsArray{1});
        %catch exception
            % NOTE: The message is slightly inaccurate since it assumes that the first argument is a S/W mode.
            % Argument "--version" etc. would have worked too.
            %error('BICAS:CLISyntax', 'Can not interpret argument "%s" as a S/W mode.', cliArgumentsArray{1});
        %end

        %-----------------------------------------------------------------
        % Configure CLI flags to expect, partly depending on the S/W mode
        %
        % NOTE: The flags are identified by strings (container.Map keys) which are a combination of the namespaces for
        % 1) identifiers for misc. flags e.g. "output_dir", "log_path".
        % 2) dataset IDs!
        %-----------------------------------------------------------------
        FlagsConfigMap('output_dir') = struct('cliString', '--output', 'isRequired', 1, 'expectsValue', 1);
        C_inputs = C_sw_mode.inputs;      % C = Constants structure.
        inputPdids = {};                  % List of keys used for input files.
        for iInput = 1:length(C_inputs)
            pdid = C_inputs{iInput}.PDID;
            
            % Configure one flag+value pair
            FlagConfig = [];
            FlagConfig.cliString    = ['--', C_inputs{iInput}.CLI_parameter];
            FlagConfig.isRequired   = 1;
            FlagConfig.expectsValue = 1;
            FlagsConfigMap(pdid) = FlagConfig;
            
            inputPdids{end+1} = pdid;
        end

        %-----------------------------
        % Parse (remaining) arguments
        %-----------------------------
        ParsedCliArgumentsMap = bicas.utils.parse_CLI_flags(cliArgumentsArray(2:end), FlagsConfigMap);
        
        
        
        InputFilesMap = containers.Map(inputPdids, ParsedCliArgumentsMap.values(inputPdids));   % Extract subset of parsed arguments.
        
        outputDir = bicas.utils.get_abs_path(ParsedCliArgumentsMap('output_dir'));
        
        

        %===============================================================
        % EXECUTE S/W MODE
        %
        % CHOOSE IMPLEMENTATION TO USE.
        %===============================================================
        bicas.execute_sw_mode(DataManager, C_sw_mode.CLI_parameter, InputFilesMap, outputDir)   % The intended real implementation
        %execute_sw_mode_TEST_IMPLEMENTATION(C_sw_mode.CLI_parameter, outputDir)   % IMPLEMENTATION FOR TESTING. OUTPUTS NONSENSE CDFs.
        
    end

        

    executionWallTimeSeconds = toc(startTimeTicSeconds);
    irf.log('n', sprintf('Execution took %g s (wall time).', executionWallTimeSeconds));    % Always log (-->critical)?
        
        
        
    % EXIT
    errorCode = ERROR_CODES.NO_ERROR;   % Default RETURN value.



catch exception
    
    try
        irf.log('critical', 'Main function caught an exception. Beginning error handling.');   % Print to stdout.
        
        message = exception.message;
        
        %==================================================================
        % Convert MATLAB error message identifiers into return error codes
        %==================================================================
        % NOTE: The order in which tests occur matters since the same error message identifier may contain multiple
        % matching components. Therefore, the matching of msg IDs with error codes is also not an exact science.
        errorId = strsplit(exception.identifier, ':');
        if     any(strcmpi(errorId, 'OperationNotImplemented'));   errorCode = ERROR_CODES.OPERATION_NOT_IMPLEMENTED;
        elseif any(strcmpi(errorId, 'PathNotFound'));              errorCode = ERROR_CODES.PATH_NOT_FOUND;
        elseif any(strcmpi(errorId, 'CLISyntax'));                 errorCode = ERROR_CODES.CLI_SYNTAX_ERROR;
        elseif any(strcmpi(errorId, 'SWModeProcessing'));          errorCode = ERROR_CODES.SW_MODE_PROCESSING_ERROR;
        elseif any(strcmpi(errorId, 'DatasetFormat'));             errorCode = ERROR_CODES.DATASET_FORMAT_ERROR;
        elseif any(strcmpi(errorId, 'IllegalConfiguration'));      errorCode = ERROR_CODES.CONFIGURATION_ERROR;
        elseif any(strcmpi(errorId, 'Assertion')) ;                errorCode = ERROR_CODES.ASSERTION_ERROR;
        %elseif any(strcmpi(errorId, ''))
        %    errorCode = ERROR_CODES.;
        else
            errorCode = ERROR_CODES.MISC_ERROR;
            %errorCode = ERROR_CODES.ERROR_IN_MATLAB_ERROR_HANDLING;
        end
        
        %======================
        % Print the call stack
        %======================
        callStackLength = length(exception.stack);
        fprintf(2, 'MATLAB call stack:\n');    % Print to stderr.
        if (~isempty(callStackLength))
            for i=1:callStackLength
                stackCall = exception.stack(i);
                temp      = strsplit(stackCall.file, filesep);
                filename  = temp{end};
                
                fprintf(2, '    %-25s %-55s row %i,\n', [filename, ','], [stackCall.name, ','], stackCall.line);
            end
        end
        
        fprintf(2, [message, '\n']);    % Print to stderr.
        
        fprintf(2, 'Exiting MATLAB application with error code %i.\n', errorCode);        % Print to stderr.
        
        return
        
    catch exception
        %===================================================
        % CASE: There was an error in the error handling(!)
        %===================================================
        
        % NOTE: Only use very, very error safe code here.
        fprintf(2, 'Unknown error. Error in the MATLAB code''s error handling.\nException message: "%s"\n', ...
            exception.message');   % Print to stderr.
        
        errorCode = ERROR_CODES.ERROR_IN_MATLAB_ERROR_HANDLING;   % Not even use hardcoded constant for this error?!!
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
%===================================================================================================
% function execute_sw_mode_TEST_IMPLEMENTATION(sw_mode_CLI_parameter, outputDir)
% 
% global ERROR_CODES CONSTANTS
% 
% irf.log('c', 'USING TEST IMPLEMENTATION FOR S/W MODES. ONLY CREATES NONSENSE CDF FILES.')
% 
% %C_mode = get_C_sw_mode(sw_mode_CLI_parameter);
% temp = bicas.utils.select_structs(C.sw_modes, 'CLI_parameter', {sw_mode_CLI_parameter});
% C_mode = temp{1};
% output_JSON = [];
% 
% % Iterate over OUTPUTS
% for i = 1:length(C_mode.outputs)
%     C_mode_output = C_mode.outputs{i};
%     master_cdf_filename = C_mode_output.master_cdf_filename;
%     output_filename = [C_mode_output.dataset_ID, '_V', C_mode_output.skeleton_version_str, '.cdf'];
%     
%     src_file  = fullfile(bias_constants.sw_root_dir(), C.master_cdfs_dir_rel, master_cdf_filename);
%     dest_file = fullfile(outputDir, output_filename);
%     
%     irf.log('n', 'Trying to copy file')
%     irf.log('n', sprintf('   from %s', src_file))
%     irf.log('n', sprintf('   to   %s', dest_file))
%     [success, copyfile_msg, ~] = copyfile(src_file, dest_file);   % Overwrites any pre-existing file.
%     if ~success
%         error('BICAS:FailedToCopyFile',, ...
%             'Failed to copy file\n    from "%s"\n    to   "%s".\n"copyfile" error message: "%s"', ...
%             src_file, dest_file, copyfile_msg)
%     end
%     
%     output_JSON.(C_mode_output.JSON_output_file_identifier) = output_filename;
% end
% 
% % Print list of files produced in the form of a JSON object.
% str = bicas.utils.JSON_object_str(output_JSON);
% bicas.stdout_disp(str);
% 
% end



%===================================================================================================
function print_version()

% IMPLEMENTATION NOTE: Uses the software version in the S/W descriptor rather than the in the BICAS
% constants since the RCS ICD specifies that it should be that specific version.
% This in principle inefficient but precise.

global CONSTANTS

D = bicas.get_sw_descriptor();
bicas.stdout_printf('Version %s\n', D.release.version)

end



%===================================================================================================
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07
%
% Print the JSON S/W descriptor.
%
function print_identification()
global CONSTANTS

D = bicas.get_sw_descriptor();
str = bicas.utils.JSON_object_str(D, CONSTANTS.C.JSON_object_str);
bicas.stdout_disp(str);

end



%===================================================================================================
function print_help(ERROR_CODES)
%
% PROPOSAL: Print error codes. Can use implementation to list them?
%    PROPOSAL: Define error codes with description strings?! Map?! Check for doubles?!
% PROPOSAL: Print CLI syntax incl. for all modes? More easy to parse than the S/W descriptor.

%error('BICAS:OperationNotImplemented', 'Operation not implemented: --help.')

D = bicas.get_sw_descriptor();
bicas.stdout_printf('%s\n%s\n', D.identification.name, D.identification.description)
bicas.stdout_printf('\nError codes (internal constants):\n')
for sfn = fieldnames(ERROR_CODES)'
    errorCode = ERROR_CODES.(sfn{1});
    errorName = sfn{1};
    bicas.stdout_printf('   %3i = %s\n', errorCode, errorName)
end
bicas.stdout_printf('\nSee "readme.txt" for more help.\n')

end
