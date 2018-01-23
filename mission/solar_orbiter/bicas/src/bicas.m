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
% - The (MATLAB) code ignores but permits the CLI flags --log and --config.
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
% PROPOSAL: Set flag for MATLAB warnings. Disable?
%    NOTE: TN claims warnings are sent to stdout.
% PROPOSAL: Extra (inofficial) flag for setting the log level.
% TODO-NEED-INFO: Is the applicaton allowed to overwrite output files?
%
% PROPOSAL: Check that all master cdf files are present/available.
%
% PROPOSAL: Rename to bicas.main (+bicas/main.m).
% PROPOSAL: Put a summarized version of CLI syntax in "bicas --help" (somethinger easier that the S/W descriptor).
%    PRO: Useful when S/W descriptor becomes big and complex.
%
% NOTE: Implementation of the parsing of CLI arguments is problematic for when processing s/w mode.
%   Present: The s/w mode argument influences which succeeding arguments are permitted (depends on the s/w mode).
%   ==> Must invoke DataManager before SETTINGS is fully initialized (from the CLI arguments)
%   Anticipated future changes: Arguments after s/w mode argument (choice of pipeline, test modes) influences which s/w
%   modes are allowed, and also which arguments those s/w modes in turn allow.
%   PROPOSAL: Make it possible to separate sequence of inofficial arguments from other arguments before parsing individual flags.
%       PRO: Makes it possible to first parse the inofficial arguments and modofy SETTINGS, before parsing anything
%               s/w modes which might require an updated SETTINGS variable.
%       PROPOSAL: Inofficial arguments can only be added before s/w mode.
%           PRO: Does not need separator argument.
%           CON: Can not parse reliably since does not know where sequence of inofficial arguments ends.
%       PROPOSAL: Inofficial arguments can only be added after (inofficial, optional) separator argument, e.g. "---".
%           NOTE: Separator argument can be placed in SETTINGS itself, but not be read from CLI arguments.
%           PRO: Smoother for just adding/appending inofficial arguments after existing arguments.
%   PROPOSAL: Somehow try all possible interpretation of arguments to see if any one of them matches, e.g. try all s/w
%       modes.
%
% PROPOSAL: Not declare SETTINGS as a global variable until it is certain that it has been updated/finalized.
%   PROPOSAL: Different names for global and local SETTINGS variable, even if temporary.



% Clear any previous instance of global variables.
% This is useful to avoid mistakenly using a previously initialized version of CONSTANTS or SETTINGS when the
% initialization has failed and when developing in MATLAB.
clear -global CONSTANTS SETTINGS

[ERROR_TYPES_INFO, REQUIRED_MATLAB_VERSION] = bicas.error_safe_constants();



try

    errorCode = main(REQUIRED_MATLAB_VERSION, ERROR_TYPES_INFO, varargin);

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



end



% BICAS's de facto main function, without error handling.
function errorCode = main(REQUIRED_MATLAB_VERSION, ERROR_TYPES_INFO, cliArgumentsList)

global CONSTANTS                  % Global structure "CONSTANTS" is initialized later.
global SETTINGS                   % Global structure "SETTINGS"  is initialized later.

startTimeTicSeconds = tic;

% Among other things: Sets up paths to within irfu-matlab (excluding .git/).
% NOTE: Prints to stdout. Can not deactivate this behaviour!
% NOTE: Should not call irf('check') which looks for updates to irfu-matlab (can not distinguish between updates to
%       BICAS or the rest of irfu-matlab).
irf('check_path');
irf('check_os');              % Maybe not strictly needed.
irf('matlab');                % Maybe not strictly needed.
irf('cdf_leapsecondstable');
irf.log('notice')             % Set initial log level value until it is later overridden by the config value.
irf('version')                % Print e.g. "irfu-matlab version: 2017-02-21,  v1.12.6".

%======================
% Check MATLAB version
%======================
matlabVersionString = version('-release');
if ~strcmp(matlabVersionString, REQUIRED_MATLAB_VERSION)
    error('BICAS:BadMatlabVersion', ...
        'Using bad MATLAB version. Found version "%s". BICAS requires version "%s".\n', ...
        matlabVersionString, REQUIRED_MATLAB_VERSION)
end
fprintf(1, 'Using MATLAB, version %s.\n', matlabVersionString);

%=======================================================================
% Derive the root path of the software (BICAS directory structure root)
%=======================================================================
% ASSUMES: The current file is in the <BICAS>/src directory.
[matlabSrcPath, ~, ~] = fileparts(mfilename('fullpath'));
bicasRootPath = bicas.utils.get_abs_path(fullfile(matlabSrcPath, '..'));



%=============================
% Initialize global constants
%=============================
% IMPLEMENTATION NOTE: Does not initialize CONSTANTS until here because:
%    1) MATLAB version should have been checked for first. The initialization code could otherwise fail.
%    2) Needs BICAS root path.
% NOTE: Constants will later be modified by the CLI arguments.
% Should preferably not use irf.log before here so that the right logging level is used.
CONSTANTS = bicas.constants(bicasRootPath);
SETTINGS  = bicas.settings;
irf.log(SETTINGS.get('LOGGING.IRF_LOG_LEVEL'));   % NOTE: May set the logging level to the same level as before.



%===================================================================
% Configure permitted flags COMMON for all BICAS modes of operation
%===================================================================
FlagsConfigMap = containers.Map;
% NOTE: log_path and config_file_path are both flag+value to permit but ignore since they are handled by bash launcher script.
FlagsConfigMap('log_path')          = struct('cliFlagString', '--log',     'occurrenceRequirement', '0-1',   'nValues', 1);
FlagsConfigMap('config_file_path')  = struct('cliFlagString', '--config',  'occurrenceRequirement', '0-1',   'nValues', 1);
FlagsConfigMap('modified_settings') = struct('cliFlagString', '--setting', 'occurrenceRequirement', '0-inf', 'nValues', 2);



DataManager = bicas.data_manager();    % Requires CONSTANTS (not necessarily SETTINGS) to be initialized.



%=====================================================================
% Read the first CLI argument -- Determine BICAS modes of operation
% -----------------------------------------------------------------
% ==> Configure permitted CLI flags
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
        ExtendedSwModeInfo = DataManager.get_extended_sw_mode_info(cliArgumentsList{1});    % NOTE: FIRST USE OF DataManager.
    catch exception1
        % NOTE: Argument "--verson" (misspelled "--version") etc. would have produced error here too.
        error('BICAS:CLISyntax', 'Can not interpret first argument "%s" as a S/W mode (or any other legal first argument).', cliArgumentsList{1});
    end
    
    
    
    %==============================================================================================================
    % Configure requirements on (remaining) CLI arguments depending on the S/W mode
    % -----------------------------------------------------------------------------
    % NOTE/BUG RISK: The flags are identified by strings (container.Map keys) which are a in reality a combination
    % of the namespaces for:
    % (1) identifiers for misc. flags e.g. "output_dir", "log_path".
    % (2) dataset IDs!
    % This is not really appropriate but works as long as there is no overlap between the two sets of strings.
    %
    % PROPOSAL: Assertion for checking whether the map key has previously used.
    %==============================================================================================================
    FlagsConfigMap('output_dir') = struct('cliFlagString', '--output', 'occurrenceRequirement', '1', 'nValues', 1);
    inputsInfoList = ExtendedSwModeInfo.inputs;      % C = Constants structure.
    inputPdidsList = {};                  % List of keys used for input files.
    
    for iInput = 1:length(inputsInfoList)    % For every input dataset...
        pdid = inputsInfoList{iInput}.PDID;
        
        % Configure one flag+value pair.
        FlagConfig = [];
        FlagConfig.cliFlagString         = ['--', inputsInfoList{iInput}.CLI_PARAMETER];
        FlagConfig.occurrenceRequirement = '1';
        FlagConfig.nValues               = 1;
        
        % ASSERTION
        if FlagsConfigMap.isKey(pdid)
            error('BICAS:Assertion:IllegalConfiguration', 'Dataset ID used as flag identifier conflicts with other flag identifier. Bad hardcoding.')
        end
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
ModifiedSettingsAsStrings = containers.Map;
valuesListsLists = FlagValuesMap('modified_settings');
for iSetting = 1:length(valuesListsLists)
    ModifiedSettingsAsStrings(valuesListsLists{iSetting}{1}) = valuesListsLists{iSetting}{2};
end



% Modify settings
% Should preferably not use irf.log before here so that the right logging level is used.
SETTINGS.modify_settings(ModifiedSettingsAsStrings)

% CASE: SETTINGS has now been finalized and will never change after this.
irf.log(SETTINGS.get('LOGGING.IRF_LOG_LEVEL'));
irf.log('n', bicas.sprint_SETTINGS)                 % Prints the contents of SETTINGS.



%=======================================
% Log misc. paths and all CLI arguments
%=======================================
irf.log('n', sprintf('BICAS software root path:      "%s"', bicasRootPath))
irf.log('n', sprintf('BICAS MATLAB source code path: "%s"', matlabSrcPath))
irf.log('n', sprintf('Current working directory:     "%s"', pwd));   % Useful for debugging the use of relative directory arguments.
for i = 1:length(cliArgumentsList)
    irf.log('n', sprintf('CLI argument %2i: "%s"', i, cliArgumentsList{i}))    % PROPOSAL: Combine into a single multiline log message?
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
    print_help(ERROR_TYPES_INFO, DataManager)
    
elseif strcmp(bicasModeOfOperation, 'Processing S/W mode')
    
    % Extract the input files (datasets) from CLI arguments.
    InputFilesMap = containers.Map;
    for iPdid = 1:length(inputPdidsList)
        valuesListsLists = FlagValuesMap(inputPdidsList{iPdid});
        InputFilesMap(inputPdidsList{iPdid}) = valuesListsLists{1}{1};   % Extract subset of parsed arguments.
    end
    
    % Extract the output directory from CLI arguments.
    valuesListsLists = FlagValuesMap('output_dir');
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



% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07
%
% Print the JSON S/W descriptor.
%
function print_identification(DataManager)

global SETTINGS

swd = bicas.get_sw_descriptor(DataManager);
str = bicas.utils.JSON_object_str(swd, ...
    SETTINGS.get('JSON_OBJECT_STR.INDENT_SIZE'), ...
    SETTINGS.get('JSON_OBJECT_STR.VALUE_POSITION'));
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
