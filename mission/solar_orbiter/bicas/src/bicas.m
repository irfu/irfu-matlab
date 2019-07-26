%
% Main MATLAB function that launches BICAS.
%
% BICAS = BIAS CAlibration Software
%
% This function is BICAS' main MATLAB function, i.e. it is called by no other MATLAB code during regular use. It is
% intended to be wrapped in, and called from, non-MATLAB code, e.g. a bash script.
%
%
% IMPORTANT NOTE: DOCUMENTATION
% =============================
% Documentation important for using BICAS as a whole can be found in
% (1) "readme.txt" and other documentation text files (*.txt), and
% (2) the RCS SUM document.
% To prevent the duplication of documentation, comments in the source code tries to only cover subjects important for
% understanding the implementation and information not already present in other documentation text files.
%
%
% ARGUMENTS
% =========
% varargin: This function expects exactly the CLI arguments submitted to the bash launcher script as a sequence of
%           MATLAB strings. This function therefore expects the arguments defined in the RCS ICD and possibly additional
%           inoffical arguments.
% Notes:
% - The official parameter syntax for S/W modes must be in agreement with "roc_sw_descriptor.js" as specified by
%   the RCS ICD.
% - The parameter syntax may contain additional inofficial parameters, which are useful for development/debugging, but
%   which are still compatible with the RCS ICD.
% - The (MATLAB) code ignores but permits the CLI option --log.
%
%
% RETURN VALUE
% ============
% errorCode = The error code that is to be passed on to the OS/shell.
%
%
% NOTES
% =====
% ASSUMES: The current file is in the <BICAS>/src directory.
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-06-02) but may very well work with other
% versions of MATLAB.
%
% IMPLEMENTATION NOTE: This code does not quit/exit using the MATLAB function "quit" since that always exits all of
% MATLAB which is undesirable when developing code from in the MATLAB IDE. The function returns the error code to make
% it possible for the bash wrapper to quit with an exit code instead.
%
% IMPLEMENTATION NOTE: The RCS ICD specifies what should go to stdout. BICAS
% 1) prints all log messages to stdout, and
% 2) prints all messages intended for BICAS' final, actual stdout (as produced by the bash wrapper) to stdout but with a
% prefix so they can be filtered out by the calling wrapper bash script.
% RATIONALE: See the bash wrapper script.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-03-xx
%
function errorCode = bicas( varargin )
%
% PROPOSAL: Set option for MATLAB warnings. Disable?
%    NOTE: TN claims warnings are sent to stdout.
% TODO-NEED-INFO: Is the application allowed to overwrite output files?
%
% PROPOSAL: Check that all master cdf files are present/available.
%
% PROPOSAL: Rename to bicas.main (+bicas/main.m).
%   NOTE: Function already has an internal function named "main".
%   PROPOSAL: 
%
% PROPOSAL: Put a summarized version of CLI syntax in "bicas --help" (somethinger easier that the S/W descriptor).
%    PRO: Useful when S/W descriptor becomes big and complex.
%
% PROPOSAL: Split up the "main" function in several functions (outsource chunks of its code to smaller functions which are called).
%
% PROPOSAL: Better handling of errors in dataobj (reading CDF files).
%   PROPOSAL: Wrap dataobj in function and catch and rethrow errors with BICAS' error IDs.



% Clear any previous instance of global variables
% -----------------------------------------------
% This is useful to avoid mistakenly using a previously initialized version of CONSTANTS or SETTINGS when the
% initialization has failed and when developing in MATLAB. Must be done as early as possible in the execution.
clear -global CONSTANTS SETTINGS

[ERROR_TYPES_INFO, REQUIRED_MATLAB_VERSION, INOFFICIAL_ARGUMENTS_SEPARATOR] = bicas.error_safe_constants();



try

    errorCode = main(REQUIRED_MATLAB_VERSION, ERROR_TYPES_INFO, INOFFICIAL_ARGUMENTS_SEPARATOR, varargin);

catch Exception1
    %================================================================
    % CASE: Caught an error in the regular execution of the software
    %================================================================

    try
        irf.log('critical', 'Main function caught an exception. Beginning error handling.');   % Print to stdout.
        fprintf(2, 'exception1.identifier = "%s"\n', Exception1.identifier);    % Print to stderr.
        fprintf(2, 'exception1.message    = "%s"\n', Exception1.message);       % Print to stderr.

        %=================================================================================
        % Use MATLAB error message identifiers to identify one or multiple "error types".
        %=================================================================================
        msgIdentifierParts = strsplit(Exception1.identifier, ':');
        errorTypesList = msgIdentifierParts(ERROR_TYPES_INFO.isKey(msgIdentifierParts));    % Cell array of error types (strings) only.
        if isempty(errorTypesList)
            errorTypesList = {'UntranslatableErrorMsgId'};
        end

        %===================================
        % Print all identified error types.
        %===================================
        fprintf(2, 'Matching error types:\n');% Print to stderr.
        for i = 1:numel(errorTypesList)
            fprintf(2, '    %s\n', ERROR_TYPES_INFO(errorTypesList{i}).description);   % Print to stderr.
        end
        % NOTE: Choice - Uses the last part of the message ID for determining error code to return.
        errorCode = ERROR_TYPES_INFO(errorTypesList{end}).code;

        %======================
        % Print the call stack
        %======================
        callStackLength = length(Exception1.stack);
        fprintf(2, 'MATLAB call stack:\n');    % Print to stderr.
        if (~isempty(callStackLength))
            for i=1:callStackLength
                stackCall = Exception1.stack(i);
                temp      = strsplit(stackCall.file, filesep);
                filename  = temp{end};

                fprintf(2, '    %-25s %-55s row %i,\n', [filename, ','], [stackCall.name, ','], stackCall.line);
            end
        end



        fprintf(2, 'Exiting MATLAB application with error code %i.\n', errorCode);        % Print to stderr.
        return

    catch Exception2    % Deliberately use different variable name to distinguish the exception from the previous one.
        %========================================================
        % CASE: Caught an error in the regular error handling(!)
        %========================================================

        % NOTE: Only use very, very error-safe code here.
        fprintf(2, 'Error in the MATLAB code''s error handling.\n');   % Print to stderr.
        fprintf(2, 'exception2.identifier = "%s"\n', Exception2.identifier);          % Print to stderr.
        fprintf(2, 'exception2.message    = "%s"\n', Exception2.message);             % Print to stderr.

        errorCode = ERROR_TYPES_INFO('MatlabCodeErrorHandlingError').code;             % Use hardcoded constant for this error?!!

        fprintf(2, 'Exiting MATLAB application with error code %i.\n', errorCode);    % Print to stderr.
        return
    end
end



end    % bicas



% BICAS's de facto main function, without error handling.
function errorCode = main(REQUIRED_MATLAB_VERSION, ERROR_TYPES_INFO, INOFFICIAL_ARGUMENTS_SEPARATOR, cliArgumentsList)



startTimeTicSeconds = tic;



%==================================
% ~ASSERTION: Check MATLAB version
%==================================
matlabVersionString = version('-release');
if ~strcmp(matlabVersionString, REQUIRED_MATLAB_VERSION)
    error('BICAS:BadMatlabVersion', ...
        'Using bad MATLAB version. Found version "%s". BICAS requires version "%s".\n', ...
        matlabVersionString, REQUIRED_MATLAB_VERSION)
end
fprintf(1, 'Using MATLAB, version %s.\n', matlabVersionString);

%===================================================================================================================
% Initialize irfu-matlab "library"
%
% Among other things: Sets up paths to within irfu-matlab (excluding .git/).
% NOTE: Prints to stdout. Can not deactivate this behaviour!
% NOTE: Should not call irf('check') which looks for updates to irfu-matlab (can not distinguish between updates to
%       BICAS or the rest of irfu-matlab).
%===================================================================================================================
irf('check_path');
irf('check_os');              % Maybe not strictly needed.
irf('matlab');
irf('cdf_leapsecondstable');
irf('version')                % Print e.g. "irfu-matlab version: 2017-02-21,  v1.12.6".
irf.log('notice')             % Set initial log level value until it is later overridden by the config value.



%===============================
% Derive BICAS's directory root
%===============================
%ROC_RCS_PATH = getenv('ROC_RCS_PATH');     % Use environment variable.
%irf.log('n', sprintf('ROC_RCS_PATH = "%s"', ROC_RCS_PATH));
% ASSUMES: The current file is in the <BICAS>/src directory.
[matlabSrcPath, ~, ~] = fileparts(mfilename('fullpath'));   % Use path of the current MATLAB file.
bicasRootPath         = bicas.utils.get_abs_path(fullfile(matlabSrcPath, '..'));



%=======================================
% Log misc. paths and all CLI arguments
%=======================================
irf.log('n', sprintf('BICAS software root path:  "%s"', bicasRootPath))
irf.log('n', sprintf('Current working directory: "%s"', pwd));   % Useful for debugging the use of relative directory arguments.
for i = 1:length(cliArgumentsList)
    irf.log('n', sprintf('CLI argument %2i: "%s"', i, cliArgumentsList{i}))    % PROPOSAL: Combine into a single multiline log message?
end



%========================================
% Initialize global settings & constants
%========================================
% IMPLEMENTATION NOTE: Does not initialize CONSTANTS until here because:
%    1) MATLAB version should have been checked for first. The initialization code could otherwise fail.
%    2) Needs BICAS root path.
% NOTE: Constants will later be modified by the CLI arguments.
% Should preferably not use irf.log before here so that the right logging level is used.
global CONSTANTS
global SETTINGS
CONSTANTS = bicas.constants(bicasRootPath);
SETTINGS  = bicas.create_default_SETTINGS();



%=============================================
% First-round interpretation of CLI arguments
%=============================================
CliData = bicas.interpret_CLI_args(cliArgumentsList, INOFFICIAL_ARGUMENTS_SEPARATOR);



%=================================================
% Modify settings according to configuration file
%=================================================
if ~isempty(CliData.configFile)
    rowList = bicas.utils.read_text_file(CliData.configFile);
    ConfigFileSettingsVsMap = bicas.interpret_config_file(rowList);
    SETTINGS.set_preexisting_from_strings(ConfigFileSettingsVsMap);    % Modify SETTINGS
end



%=========================================================
% Modify settings according to (inofficial) CLI arguments
%=========================================================
SETTINGS.set_preexisting_from_strings(CliData.ModifiedSettingsMap);    % Modify SETTINGS
SETTINGS.make_read_only();
% CASE: SETTINGS has now been finalized and is read-only (by assertion) after this.



irf.log(SETTINGS.get_fv('LOGGING.IRF_LOG_LEVEL'));
irf.log('n', bicas.sprint_SETTINGS)                 % Prints/log the contents of SETTINGS.



%================================
% Set pipelineId, calibrationDir
%================================
% COMPLETE CODE, BUT NOT ALL NEEDED BY OTHER CODE YET.
%
%pipelineId     = read_env_variable('ROC_PIP_NAME',        'PROCESSING.ROC_PIP_NAME_OVERRIDE');
%calibrationDir = read_env_variable('ROC_RCS_CAL_PATH',    'PROCESSING.ROC_RCS_CAL_PATH_OVERRIDE');
masterCdfDir   = read_env_variable(SETTINGS, 'ROC_RCS_MASTER_PATH', 'PROCESSING.ROC_RCS_MASTER_PATH_OVERRIDE');
irf.log('n', sprintf('masterCdfDir = "%s"', masterCdfDir))



DataManager = bicas.data_manager();    % NOTE: Requires CONSTANTS (not necessarily SETTINGS) to be initialized.



switch(CliData.functionalityMode)
    case 'version'
        print_version(DataManager)
    case 'identification'
        print_identification(DataManager)
    case 'help'
        print_help(ERROR_TYPES_INFO, DataManager)
    case 'S/W mode'
        %==============================================================================
        % CASE: Should be a S/W mode (deduced from elimination of other possibilities)
        %==============================================================================
        try
            ExtendedSwModeInfo = DataManager.get_extended_sw_mode_info(CliData.swModeArg);    % NOTE: FIRST USE OF DataManager.
        catch Exception1
            % NOTE: Misspelled "--version" etc. would be interpreted as S/W mode and produce error here too.
            error('BICAS:CLISyntax', ...
                'Can not interpret first argument "%s" as a S/W mode (or any other legal first argument).', ...
                CliData.swModeArg);
        end
        
        
        
        %======================================================================
        % Parse CliData.SpecInputParametersMap arguments depending on S/W mode
        %======================================================================
        
        % Extract INPUT dataset files from arguments.
        inputsInfoList = ExtendedSwModeInfo.inputs;
        InputFilesMap  = containers.Map();
        for i = 1:numel(inputsInfoList)
            optionHeader = inputsInfoList{i}.CLI_OPTION_BODY;
            
            % UI ASSERTION
            if ~CliData.SpecInputParametersMap.isKey(optionHeader)
                error('BICAS:CLISyntax', 'Can not find CLI argument(s) for input "%s".', optionHeader)
            end
            
            inputFile = CliData.SpecInputParametersMap( optionHeader );
            InputFilesMap(inputsInfoList{i}.PDID) = inputFile;
        end
        
        % Extract OUTPUT dataset files from arguments.
        outputsInfoList = ExtendedSwModeInfo.outputs;
        OutputFilesMap  = containers.Map();
        for i = 1:numel(outputsInfoList)
            optionHeader = outputsInfoList{i}.CLI_OPTION_BODY;
            
            % UI ASSERTION
            if ~CliData.SpecInputParametersMap.isKey(optionHeader)
                error('BICAS:CLISyntax', 'Can not find CLI argument(s) for input "%s".', optionHeader)
            end
            
            outputFile = CliData.SpecInputParametersMap( optionHeader );
            OutputFilesMap(outputsInfoList{i}.PDID) = outputFile;
        end
        
        
        
        %==================
        % EXECUTE S/W MODE
        %==================
        bicas.execute_sw_mode( DataManager, ExtendedSwModeInfo.CLI_PARAMETER, InputFilesMap, OutputFilesMap, masterCdfDir )
        
    otherwise
        error('BICAS:Assertion', 'Illegal value functionalityMode="%s"', functionalityMode)
end    % if ... else ... / switch



executionWallTimeSeconds = toc(startTimeTicSeconds);
irf.log('n', sprintf('Execution took %g s (wall time).', executionWallTimeSeconds));    % Always log (-->critical)?



% EXIT
errorCode = ERROR_TYPES_INFO('NoError').code;   % Default RETURN value.
end



function print_version(DataManager)

% IMPLEMENTATION NOTE: Uses the software version in the S/W descriptor rather than the in the BICAS
% constants since the RCS ICD specifies that it should be that specific version.
% This is in principle inefficient but "precise".
% NOTE: Uses the s/w name from the s/w descriptor too (instead of SETTINGS) since available anyway.

swd = bicas.get_sw_descriptor(DataManager);
bicas.stdout_printf('%s version %s\n', swd.identification.name, swd.release.version)

end



% Print the JSON S/W descriptor.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07
%
function print_identification(DataManager)

global SETTINGS

swd = bicas.get_sw_descriptor(DataManager);
str = bicas.utils.JSON_object_str(swd, ...
    SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'), ...
    SETTINGS.get_fv('JSON_OBJECT_STR.VALUE_POSITION'));
bicas.stdout_disp(str);

end



% Print user-readable "help text".
%
% NOTE: Useful if the output printed by this function can be used for copy-pasting into RCS User Manual (RUM).
%
function print_help(ERROR_TYPES_INFO, DataManager)
%
% PROPOSAL: Print CLI syntax incl. for all modes? More easy to parse than the S/W descriptor.



% Print software name & description
Swd = bicas.get_sw_descriptor(DataManager);
print_version(DataManager)
bicas.stdout_printf('%s\n', Swd.identification.description)

%==========================
% Print error codes & types
%==========================
errorCodesList = cellfun(@(x) (x.code), ERROR_TYPES_INFO.values);   % Array of (unsorted) error codes.
[~, iSort] = sort(errorCodesList);
errorTypesInfoList = ERROR_TYPES_INFO.values;      % Cell array of structs (unsorted).
errorTypesInfoList = errorTypesInfoList(iSort);    % Cell array of structs sorted by error code.
bicas.stdout_printf('\nError codes:\n')
for i = 1:numel(errorTypesInfoList)
    bicas.stdout_printf('   %3i = %s\n', errorTypesInfoList{i}.code, errorTypesInfoList{i}.description)
end

% Print settings
bicas.stdout_disp(bicas.sprint_SETTINGS)   % Includes title

bicas.stdout_printf('\nSee "readme.txt" and user manual for more help.\n')
end



% Read environment variable, but allow the value to be overriden by a settings variable.
function v = read_env_variable(SETTINGS, envVarName, settingsOverrideName)
settingsOverrideValue = SETTINGS.get_fv(settingsOverrideName);

if isempty(settingsOverrideValue)
    v = getenv(envVarName);
else
    v = settingsOverrideValue;
end

% UI ASSERTION
if isempty(v)
    error('BICAS:Assertion', ...
        'Can not set internal variable corresponding to environment variable "%s" from either (1) the environment variable, or (2) settings key value "%s".', ...
        envVarName, settingsOverrideName)
end
end
