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
    irf.log('notice')      % Set initial log level value until it is later overridden by the config value.
    
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
    
    
    %===================================================================
    % Configure permitted flags COMMON for all BICAS modes of operation
    %===================================================================
    FlagsConfigMap = containers.Map;
    FlagsConfigMap('log_path')         = struct('cliFlagString', '--log',      'occurrenceRequirement', '0-1',   'nValues', 1);   % NOTE: Flag+value to permit but ignore since handled by bash launcher script.
    FlagsConfigMap('config_file_path') = struct('cliFlagString', '--config',   'occurrenceRequirement', '0-1',   'nValues', 1);   % NOTE: Flag+value to permit but ignore since handled by bash launcher script.
    FlagsConfigMap('modified_settings') = struct('cliFlagString', '--setting', 'occurrenceRequirement', '0-inf', 'nValues', 2);

    

    DataManager = bicas.data_manager();
    cliArgumentsList = varargin;
    
    
    
    %=====================================================================
    % Configure permitted CLI flags for SPECIFIC BICAS modes of operation
    %=====================================================================
    if (length(cliArgumentsList) < 1)
        error('BICAS:CLISyntax', 'Not enough arguments found.')

    elseif (strcmp(cliArgumentsList{1}, '--version'))
        bicasModeOfOperation = 'Print version';

    elseif (strcmp(cliArgumentsList{1}, '--identification'))
        bicasModeOfOperation = 'Print S/W descriptor';

    elseif (strcmp(cliArgumentsList{1}, '--help'))
        bicasModeOfOperation = 'Print help';

    else
        bicasModeOfOperation = 'Processing S/W mode';
        %==============================================
        % CASE: Should be a S/W mode (error otherwise)
        %==============================================
        try
            ExtendedSwModeInfo = DataManager.get_extended_sw_mode_info(cliArgumentsList{1});
        catch exception
            % NOTE: The message is slightly inaccurate since it assumes that the first argument is a S/W mode.
            % Argument "--veersson" (misspelled --version) etc. would have produced error here too.
            error('BICAS:CLISyntax', 'Can not interpret argument "%s" as a S/W mode (or any other legal argument).', cliArgumentsList{1});
        end

        
        
        %=======================================================================================
        % Configure requirements on (remaining) CLI arguments depending on the S/W mode
        % -----------------------------------------------------------------------------
        % NOTE: The flags are identified by strings (container.Map keys) which are a in reality a combination of the
        % namespaces for
        % 1) identifiers for misc. flags e.g. "output_dir", "log_path".
        % 2) dataset IDs!
        % This is not really appropriate.
        %=======================================================================================
        FlagsConfigMap('output_dir') = struct('cliFlagString', '--output', 'occurrenceRequirement', '1', 'nValues', 1);
        inputsInfoList = ExtendedSwModeInfo.inputs;      % C = Constants structure.
        inputPdidsList = {};                  % List of keys used for input files.
        for iInput = 1:length(inputsInfoList)
            pdid = inputsInfoList{iInput}.PDID;
            
            % Configure one flag+value pair
            FlagConfig = [];
            FlagConfig.cliFlagString         = ['--', inputsInfoList{iInput}.CLI_PARAMETER];
            FlagConfig.occurrenceRequirement = '1';
            FlagConfig.nValues               = 1;
            FlagsConfigMap(pdid) = FlagConfig;
            
            inputPdidsList{end+1} = pdid;
        end
    end



    %=======================================================================
    % Parse CLI arguments which are COMMON for all BICAS modes of operation
    %=======================================================================
    FlagValuesMap = bicas.utils.parse_CLI_flags(cliArgumentsList(2:end), FlagsConfigMap);
    
    % Extract the modified settings from the CLI arguments.
    % NOTE: ModifiedSettings is filled with values in the order of the CLI arguments.
    %       ==> A later flag (for the same setting) overwrites the value of a former flag.
    ModifiedSettings = containers.Map;
    valuesListsLists = FlagValuesMap('modified_settings');
    for iSetting = 1:length(valuesListsLists)
        ModifiedSettings(valuesListsLists{iSetting}{1}) = valuesListsLists{iSetting}{2};
    end

    
    
    %=============================
    % Initialize global constants
    %=============================
    % IMPLEMENTATION NOTE: Does not initialize CONSTANTS until here because:
    %    1) MATLAB version has been checked for. The initialization code could otherwise fail.
    %    2) Overriding settings from the CLI arguments are available here.
    % Should preferably not use irf.log before here so that the right logging level is used.
    CONSTANTS = bicas.constants(bicasRootPath, ModifiedSettings);
    irf.log(CONSTANTS.C.LOGGING.IRF_LOG_LEVEL);
    
    

    %=======================================
    % Log misc. paths and all CLI arguments
    %=======================================
    irf.log('n', sprintf('BICAS software root path:      "%s"', bicasRootPath))
    irf.log('n', sprintf('BICAS MATLAB source code path: "%s"', matlabSrcPath))
    irf.log('n', sprintf('Current working directory:     "%s"', pwd));   % Useful for debugging the use of relative directory arguments.
    for i = 1:length(varargin)
        irf.log('n', sprintf('CLI argument %2i: "%s"', i, varargin{i}))    % PROPOSAL: Combine into a single multiline log message?
    end



    %===========================================================================
    % Perform actions which are SPECIFIC for different BICAS modes of operation
    %===========================================================================
    if strcmp(bicasModeOfOperation, 'Print version')
        %============================
        % CASE: Print version
        %============================
        print_version(DataManager)

    elseif strcmp(bicasModeOfOperation, 'Print S/W descriptor')
        %============================
        % CASE: Print identification
        %============================
        print_identification(DataManager)

    elseif strcmp(bicasModeOfOperation, 'Print help')
        %============================
        % CASE: Print help
        %============================
        print_help(ERROR_CODES, DataManager)

    elseif strcmp(bicasModeOfOperation, 'Processing S/W mode')

        % Extract the input files (datasets) from CLI arguments.
        InputFilesMap = containers.Map;
        for iPdid = 1:length(inputPdidsList)
            valuesListsLists = FlagValuesMap(inputPdidsList{iPdid});
            InputFilesMap(inputPdidsList{iPdid}) = valuesListsLists{1}{1};   % Extract subset of parsed arguments.
        end
        
        % Extract the output directory from CLI arguments.
        valuesListsLists = FlagValuesMap('output_dir')
        outputDir = bicas.utils.get_abs_path(valuesListsLists{1}{1});
        
        

        %==================
        % EXECUTE S/W MODE
        %==================
        bicas.execute_sw_mode( DataManager, ExtendedSwModeInfo.CLI_PARAMETER, InputFilesMap, outputDir )
        
    else
        error('BICAS:Assertion', 'Can not interpret bicasModeOfOperation. This indicates a pure code bug.')
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
function print_version(DataManager)

% IMPLEMENTATION NOTE: Uses the software version in the S/W descriptor rather than the in the BICAS
% constants since the RCS ICD specifies that it should be that specific version.
% This is in principle inefficient but "precise".

swd = bicas.get_sw_descriptor(DataManager);
bicas.stdout_printf('Version %s\n', swd.release.version)

end



%===================================================================================================
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07
%
% Print the JSON S/W descriptor.
%
function print_identification(DataManager)
global CONSTANTS

swd = bicas.get_sw_descriptor(DataManager);
str = bicas.utils.JSON_object_str(swd, CONSTANTS.C.JSON_OBJECT_STR);
bicas.stdout_disp(str);

end



%===================================================================================================
function print_help(ERROR_CODES, DataManager)
%
% PROPOSAL: Print error codes. Can use implementation to list them?
%    PROPOSAL: Define error codes with description strings?! Map?! Check for doubles?!
% PROPOSAL: Print CLI syntax incl. for all modes? More easy to parse than the S/W descriptor.

%error('BICAS:OperationNotImplemented', 'Operation not implemented: --help.')

swd = bicas.get_sw_descriptor(DataManager);
bicas.stdout_printf('%s\n%s\n', swd.identification.name, swd.identification.description)
bicas.stdout_printf('\nError codes (internal constants):\n')
for sfn = fieldnames(ERROR_CODES)'
    errorCode = ERROR_CODES.(sfn{1});
    errorName = sfn{1};
    bicas.stdout_printf('   %3i = %s\n', errorCode, errorName)
end
bicas.stdout_printf('\nSee "readme.txt" for more help.\n')

end
