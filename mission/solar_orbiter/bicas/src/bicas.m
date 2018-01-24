% errorCode = bicas( varargin )   Main function that launches BICAS.
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
% ARGUMENTS AND RETURN VALUE
% ==========================
% ARGUMENTS: This function expects exactly the CLI arguments submitted to the bash launcher script as a sequence of
% MATLAB strings. This function therefore expects the arguments defined in the RCS ICD and possibly additional inoffical
% arguments.
%
% RETURN VALUE: errorCode = The error code that is to be passed on to the OS/shell.
%
% Notes:
% - The official parameter syntax for S/W modes must be in agreement with "roc_sw_descriptor.js" as specified by
%   the RCS ICD.
% - The parameter syntax may contain additional inofficial parameters, which are useful for development/debugging, but
%   which are still compatible with the RCS ICD.
% - The (MATLAB) code ignores but permits the CLI options --log and --config.
%
%
% NOTES
% =====
% ASSUMES: The current file is in the <BICAS>/src directory.
%
% NOTE: This code is designed for MATLAB 2016a (as of 2016-06-02) but may very well work with other
% versions of MATLAB.
%
% IMPLEMENTATION NOTE: This code does not quit/exit using the MATLAB function "quit" since that
% always exits all of MATLAB which is undesirable when developing in the MATLAB IDE. The function
% returns the error code to make it possible for the bash wrapper to quit with an exit code instead.
%
% IMPLEMENTATION NOTE: The RCS ICD specifies tightly what should go to stdout. BICAS
% 1) prints all log messages to stdout, and
% 2) prints all messages intended for BICAS' final, actual stdout (as produced by the bash wrapper) to stdout but with a
% prefix so they can be filtered out by the calling wrapper bash script.
% Reasons: See the bash wrapper script.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-03-xx
%
function errorCode = bicas( varargin )
%
% PROPOSAL: Set option for MATLAB warnings. Disable?
%    NOTE: TN claims warnings are sent to stdout.
% TODO-NEED-INFO: Is the applicaton allowed to overwrite output files?
%
% PROPOSAL: Check that all master cdf files are present/available.
%
% PROPOSAL: Rename to bicas.main (+bicas/main.m).
% PROPOSAL: Put a summarized version of CLI syntax in "bicas --help" (somethinger easier that the S/W descriptor).
%    PRO: Useful when S/W descriptor becomes big and complex.
%
% NOTE: Implementation of the parsing of CLI arguments is problematic for when processing s/w mode.
%   PROBLEM: Present, ICD & inofficial options mixed: The s/w mode argument influences which succeeding arguments are
%   permitted (depends on the s/w mode).
%   ==> Must invoke DataManager before SETTINGS is fully initialized (from the CLI arguments)
%   PROBLEM: Anticipated future changes: Arguments after s/w mode argument (choice of pipeline, test modes) influences
%   which s/w modes are allowed, and also which arguments those s/w modes in turn allow.
%   PROBLEM: --config (ICD option) could influence which ICD options that are permissible (e.g. pipeline).
%   --
%   PROPOSAL: Make it possible to separate sequence of inofficial arguments from other arguments before parsing individual options.
%       PRO: Makes it possible to first parse the inofficial arguments and modify SETTINGS, before parsing anything
%               s/w modes which might require an updated SETTINGS variable.
%       PROPOSAL: Inofficial arguments can only be added before s/w mode.
%           PRO: Does not need separator argument.
%           CON: Can not parse reliably since does not know where sequence of inofficial arguments ends.
%       PROPOSAL: Inofficial arguments can only be added after (inofficial, optional) separator argument, e.g. "---".
%           NOTE: Separator argument can be placed in SETTINGS itself, but not be read from CLI arguments.
%           PRO: Smoother for just adding/appending inofficial arguments after existing arguments.
%   PROPOSAL: Somehow try all possible interpretations of arguments to see if any one of them matches, e.g. try all s/w
%       modes.
%   PROPOSAL: Try to parse some options in argument sequence (e.g. --config) before parsing all.
%       CON: Can fail if searched-for options coincide with values of other options.
%           CON: Needs second parsing of all arguments. ==> Can give error for that case for syntactically correct argument lists.
%           PROPOSAL: Parse some options assuming that all other options have exactly N values argument.
%               NOTE: Does not need to know anything else about other options since this determines possible positions of
%                   the sought option: index = 1 + (N+1)*m, m=integer.
%               CON: Can still not mix inofficial and ICD arguments.
%
% PROPOSAL: Not declare SETTINGS as a global variable until it is certain that it has been updated/finalized.
%   PROPOSAL: Different names for global and local SETTINGS variable, even if temporary.


% Clear any previous instance of global variables (as early as possible).
% This is useful to avoid mistakenly using a previously initialized version of CONSTANTS or SETTINGS when the
% initialization has failed and when developing in MATLAB.
clear -global CONSTANTS SETTINGS

[ERROR_TYPES_INFO, REQUIRED_MATLAB_VERSION, INOFFICIAL_ARGUMENTS_SEPARATOR] = bicas.error_safe_constants();



try

    errorCode = main(REQUIRED_MATLAB_VERSION, ERROR_TYPES_INFO, INOFFICIAL_ARGUMENTS_SEPARATOR, varargin);

catch exception1
    
    try
        irf.log('critical', 'Main function caught an exception. Beginning error handling.');   % Print to stdout.
        fprintf(2, 'exception1.identifier = "%s"\n', exception1.identifier);    % Print to stderr.        
        fprintf(2, 'exception1.message    = "%s"\n', exception1.message);       % Print to stderr.
        
        %=================================================================================
        % Use MATLAB error message identifiers to identify one or multiple "error types".
        %=================================================================================
        msgIdentifierParts = strsplit(exception1.identifier, ':');
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
        callStackLength = length(exception1.stack);
        fprintf(2, 'MATLAB call stack:\n');    % Print to stderr.
        if (~isempty(callStackLength))
            for i=1:callStackLength
                stackCall = exception1.stack(i);
                temp      = strsplit(stackCall.file, filesep);
                filename  = temp{end};
                
                fprintf(2, '    %-25s %-55s row %i,\n', [filename, ','], [stackCall.name, ','], stackCall.line);
            end
        end
        

        
        fprintf(2, 'Exiting MATLAB application with error code %i.\n', errorCode);        % Print to stderr.        
        return
        
    catch exception2    % Deliberately use different variable name to distinguish the exception from the previous one.
        %===================================================
        % CASE: There was an error in the error handling(!)
        %===================================================
        
        % NOTE: Only use very, very error-safe code here.
        fprintf(2, 'Error in the MATLAB code''s error handling.\n');   % Print to stderr.
        fprintf(2, 'exception2.identifier = "%s"\n', exception2.identifier);          % Print to stderr.        
        fprintf(2, 'exception2.message    = "%s"\n', exception2.message);             % Print to stderr.
        
        errorCode = ERROR_TYPES_INFO('MatlabCodeErrorHandlingError').code;             % Use hardcoded constant for this error?!!
        
        fprintf(2, 'Exiting MATLAB application with error code %i.\n', errorCode);    % Print to stderr.
        return
    end
end



end    % bicas



% BICAS's de facto main function, without error handling.
function errorCode = main(REQUIRED_MATLAB_VERSION, ERROR_TYPES_INFO, INOFFICIAL_ARGUMENTS_SEPARATOR, cliArgumentsList)



startTimeTicSeconds = tic;

%=================================
% ASSERTION: Check MATLAB version
%=================================
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
irf('matlab');                % Maybe not strictly needed.
irf('cdf_leapsecondstable');
irf('version')                % Print e.g. "irfu-matlab version: 2017-02-21,  v1.12.6".
irf.log('notice')             % Set initial log level value until it is later overridden by the config value.

%=======================================================================
% Derive the root path of the software (BICAS directory structure root)
%=======================================================================
%ROC_RCS_PATH = getenv('ROC_RCS_PATH');     % Use environment variable.
%irf.log('n', sprintf('ROC_RCS_PATH = "%s"', ROC_RCS_PATH));
% ASSUMES: The current file is in the <BICAS>/src directory.
[matlabSrcPath, ~, ~] = fileparts(mfilename('fullpath'));   % Use path of the current MATLAB file.
bicasRootPathFromCode = bicas.utils.get_abs_path(fullfile(matlabSrcPath, '..'));
bicasRootPath = bicasRootPathFromCode;    % Select which path to use as BICAS root path.



%=======================================
% Log misc. paths and all CLI arguments
%=======================================
irf.log('n', sprintf('BICAS software root path:                    "%s"', bicasRootPath))
irf.log('n', sprintf('Current working directory:                   "%s"', pwd));   % Useful for debugging the use of relative directory arguments.
for i = 1:length(cliArgumentsList)
    irf.log('n', sprintf('CLI argument %2i: "%s"', i, cliArgumentsList{i}))    % PROPOSAL: Combine into a single multiline log message?
end



%=============================
% Initialize global constants
%=============================
% IMPLEMENTATION NOTE: Does not initialize CONSTANTS until here because:
%    1) MATLAB version should have been checked for first. The initialization code could otherwise fail.
%    2) Needs BICAS root path.
% NOTE: Constants will later be modified by the CLI arguments.
% Should preferably not use irf.log before here so that the right logging level is used.
global CONSTANTS
global SETTINGS
CONSTANTS = bicas.constants(bicasRootPath);
%SETTINGS  = bicas.settings;
SETTINGS = bicas.create_default_SETTINGS();
irf.log(SETTINGS.get_tv('LOGGING.IRF_LOG_LEVEL'));   % NOTE: May set the logging level to the same level as before. Remove?



%===============================================================================
% Separate CLI arguments into two different sequences/lists:
% (1) icdCliArgumentsList   = List of official arguments, as defined in RCS ICD.
% (2) inoffCliArgumentsList = List of inofficial arguments (may be empty).
%===============================================================================
iArgSeparator = find(strcmp(cliArgumentsList, INOFFICIAL_ARGUMENTS_SEPARATOR));
if numel(iArgSeparator) == 0
    icdCliArgumentsList   = cliArgumentsList;
    inoffCliArgumentsList = {};
elseif numel(iArgSeparator) == 1    % NOTE: Permit argument separator to be the very last argument.
    if (iArgSeparator <= 1)
        error('BICAS:CLISyntax', 'CLI argument separator at wrong position.')
    end
    icdCliArgumentsList   = cliArgumentsList( 1 : (iArgSeparator-1) );
    inoffCliArgumentsList = cliArgumentsList( (iArgSeparator+1) : end );
else
    error('BICAS:CLISyntax', 'Found more than one CLI argument separator.')
end



%=============================================================================
% Configure permitted ICD CLI options COMMON for all BICAS modes of operation
%=============================================================================
IcdOptionsConfigMap = containers.Map;
CONFIG_OPTION_HEADER = '--config';
% NOTE: log_path and config_file_path are both options to permit but ignore since they are handled by bash launcher script.
IcdOptionsConfigMap('log_path')            = struct('optionHeader', '--log',               'occurrenceRequirement', '0-1',   'nValues', 1);
IcdOptionsConfigMap('config_file_path')    = struct('optionHeader', CONFIG_OPTION_HEADER,  'occurrenceRequirement', '0-1',   'nValues', 1);



%==============================================================================================================
% Find the --config option among the ICD arguments in order to load the config file before parsing
% the remaining ICD CLI arguments.
% 
% NOTE: Implementation assumes that --config option is optional.
% ASSUMES: ICD CLI syntax implies that the option header can only be found at even index (2,4, ...) positions.
%==============================================================================================================
iArg = find(strcmp('CONFIG_OPTION_HEADER', icdCliArgumentsList));
if numel(iArg) > 1
    error('BICAS:bicas:CLISyntax', 'Found multiple instances of %s option.', CONFIG_OPTION_HEADER)
elseif numel(iArg) == 1
    % ASSERTION
    if ~((mod(iArg, 2) == 0) && (iArg < numel(icdCliArgumentsList)))
        error('BICAS:bicas:CLISyntax', 'Found %s option at illegal position.', CONFIG_OPTION_HEADER)
    end
    % CASE: There is one --config option.
    configFilePath = icdCliArgumentsList{iArg+1};
else
    % CASE: There is no --config option.
    configFilePath = fullfile(bicasRootPath, SETTINGS.get_tv('DEFAULT_CONFIG_FILE_RELATIVE_PATH'));
end



%=========================================================================================================
% Read environment variables specified by RCS ICD and store them in SETTINGS.
%
% NOTE: The placement of this code determines the precedence of other ways of setting the SETTINGS value.
% The placement of this code is chosen in anticipation of also reading SETTINGS from a config file.
%=========================================================================================================
env_var_2_SETTINGS('ROC_PIP_NAME',     'PROCESSING.ROC_PIP_NAME');
env_var_2_SETTINGS('ROC_RCS_CAL_PATH', 'PROCESSING.ROC_RCS_CAL_PATH');



InoffOptionsConfigMap = containers.Map;
InoffOptionsConfigMap('modified_settings') = struct('optionHeader', '--setting', 'occurrenceRequirement', '0-inf', 'nValues', 2);
% Parse inofficial CLI arguments.
InoffOptionValuesMap = bicas.utils.parse_CLI_options(inoffCliArgumentsList, InoffOptionsConfigMap);

%=======================================================================================================================
% Extract the modified settings from the inofficial CLI arguments.
%
% NOTE: ModifiedSettingsMap corresponds to one definition of ONE option (in the meaning of parse_CLI_options) and is
% filled with the corresponding option values in the order of the CLI arguments.
%       ==> A later occurrence of an option with the same first option value, overwrites previous occurrences of the
%       option with the same first option value.
%       E.g. --setting LOGGING.IRF_LOG_LEVEL w --setting LOGGING.IRF_LOG_LEVEL n
%=======================================================================================================================
ModifiedSettingsMap = containers.Map;
valuesListsLists = InoffOptionValuesMap('modified_settings');   % Should contain multiple occurrences of the same option. Is hence a "list of lists".
for iSetting = 1:length(valuesListsLists)
    ModifiedSettingsMap(valuesListsLists{iSetting}{1}) = valuesListsLists{iSetting}{2};
end
% Should preferably not use irf.log before here so that the right logging level is used.
SETTINGS.set_preexisting_from_strings(ModifiedSettingsMap);    % Modify SETTINGS
SETTINGS.make_read_only();
% CASE: SETTINGS has now been finalized and will never change after this.



irf.log(SETTINGS.get_fv('LOGGING.IRF_LOG_LEVEL'));
irf.log('n', bicas.sprint_SETTINGS)                 % Prints the contents of SETTINGS.



DataManager = bicas.data_manager();    % NOTE: Requires CONSTANTS (not necessarily SETTINGS) to be initialized.



%=====================================================================
% Read the first CLI argument -- Determine BICAS modes of operation
% -----------------------------------------------------------------
% ==> Configure permitted CLI options
%=====================================================================
if (length(icdCliArgumentsList) < 1)
    error('BICAS:CLISyntax', 'Not enough arguments found.')
    
elseif (strcmp(icdCliArgumentsList{1}, '--version'))
    %============================
    % CASE: Print version
    %============================
    %IcdOptionValuesMap = bicas.utils.parse_CLI_options(icdCliArgumentsList(2:end), IcdOptionsConfigMap);
    print_version(DataManager)
    
elseif (strcmp(icdCliArgumentsList{1}, '--identification'))
    %============================
    % CASE: Print s/w descriptor
    %============================
    %IcdOptionValuesMap = bicas.utils.parse_CLI_options(icdCliArgumentsList(2:end), IcdOptionsConfigMap);
    print_identification(DataManager)
    
elseif (strcmp(icdCliArgumentsList{1}, '--help'))
    %============================
    % CASE: Print help
    %============================
    %IcdOptionValuesMap = bicas.utils.parse_CLI_options(icdCliArgumentsList(2:end), IcdOptionsConfigMap);
    print_help(ERROR_TYPES_INFO, DataManager)

else
    %==============================================================================
    % CASE: Should be a S/W mode (deduced from elimination of other possibilities)
    %==============================================================================
    try
        ExtendedSwModeInfo = DataManager.get_extended_sw_mode_info(icdCliArgumentsList{1});    % NOTE: FIRST USE OF DataManager.
    catch exception1
        % NOTE: Argument "--verson" (misspelled "--version") etc. would have produced error here too.
        error('BICAS:CLISyntax', ...
            'Can not interpret first argument "%s" as a S/W mode (or any other legal first argument).', ...
            icdCliArgumentsList{1});
    end
    
    
    
    %==============================================================================================================
    % Configure requirements on (remaining) (ICD) CLI arguments depending on the S/W mode
    % -----------------------------------------------------------------------------------
    % NOTE/BUG RISK: The options are identified by option ID strings (container.Map keys). Here the code uses
    % (1) identifiers for misc. options (e.g. "output_dir", "log_path"), and
    % (2) dataset IDs
    % as option IDs. This is not really appropriate but works as long as there is no overlap between the two sets of
    % strings. Assertion checks for illegal option IDs.
    %==============================================================================================================
    IcdOptionsConfigMap('output_dir') = struct('optionHeader', '--output', 'occurrenceRequirement', '1', 'nValues', 1);
    inputsInfoList = ExtendedSwModeInfo.inputs;
    inputPdidsList = {};                  % List of keys used for input files.
    
    for iInput = 1:length(inputsInfoList)    % For every input dataset...
        pdid = inputsInfoList{iInput}.PDID;
        
        % Configure one option.
        OptionConfig = [];
        OptionConfig.optionHeader          = ['--', inputsInfoList{iInput}.OPTION_HEADER_SH];
        OptionConfig.occurrenceRequirement = '1';
        OptionConfig.nValues               = 1;
        
        % ASSERTION
        if IcdOptionsConfigMap.isKey(pdid)
            error('BICAS:Assertion:IllegalConfiguration', 'Dataset ID used as option identifier conflicts with other option identifier. Bad hardcoding.')
        end
        IcdOptionsConfigMap(pdid) = OptionConfig;
        
        inputPdidsList{end+1} = pdid;
    end
    
    %================================
    % Parse ICD CLI arguments (bulk)
    %================================
    IcdOptionValuesMap = bicas.utils.parse_CLI_options(icdCliArgumentsList(2:end), IcdOptionsConfigMap);
    
    
    
    % Extract the input files (datasets) from CLI arguments.
    InputFilesMap = containers.Map;
    for iPdid = 1:length(inputPdidsList)
        valuesListsLists = IcdOptionValuesMap(inputPdidsList{iPdid});
        InputFilesMap(inputPdidsList{iPdid}) = valuesListsLists{1}{1};   % Extract subset of parsed arguments.
    end
    
    % Extract the output directory from CLI arguments.
    valuesListsLists = IcdOptionValuesMap('output_dir');
    outputDir = bicas.utils.get_abs_path(valuesListsLists{1}{1});
    
    %==================
    % EXECUTE S/W MODE
    %==================
    bicas.execute_sw_mode( DataManager, ExtendedSwModeInfo.CLI_PARAMETER, InputFilesMap, outputDir )
    
end    % main



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



% Print "help text".
% 
% NOTE: Useful if this can be used for copy-pasting into RCS User Manual (RUM).
%
function print_help(ERROR_TYPES_INFO, DataManager)
%
% PROPOSAL: Print CLI syntax incl. for all modes? More easy to parse than the S/W descriptor.



% Print software name & description
swd = bicas.get_sw_descriptor(DataManager);
print_version(DataManager)
bicas.stdout_printf('%s\n', swd.identification.description)

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



% Read value of environment variable and use its value to set a corresponding value in SETTINGS.
%
% NOTE: Will only modify SETTINGS if the environment variable exists/is non-empty.
% NOTE: Will only set the SETTINGS value as a string value (not numeric).
%
% IMPLEMENTATION NOTE: Function may seem superfluous but it clarifies the code, and removes some mistyping risks.
function env_var_2_SETTINGS(envVarName, settingsKey)
global SETTINGS

envVarValue = getenv(envVarName);
if ~isempty(envVarValue)
    SETTINGS.set_prexisting(settingsKey, envVarValue);
end
end
