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
% ASSUMES: The current file is located in the <BICAS>/src directory.
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
function errorCode = main( varargin )
%
% PROPOSAL: Set option for MATLAB warnings. Disable?
%   NOTE: TN claims warnings are sent to stdout.
% TODO-NEED-INFO: Is the application allowed to overwrite output files?
%
% PROPOSAL: Put a summarized version of CLI syntax in "bicas --help" (somethinger easier that the S/W descriptor).
%   PRO: Useful when S/W descriptor becomes big and complex.
%
% PROPOSAL: Split up the "main_without_error_handling" function in several functions (outsource chunks of its code to smaller functions which are called).
%   PROPOSAL: Printing CLI arguments.
%
% PROPOSAL: Better handling of errors in dataobj (reading CDF files).
%   PROPOSAL: Wrap dataobj in function and catch and rethrow errors with BICAS' error IDs.
%
% PROPOSAL: When printing settings, print how each key value has been set: default, config file, CLI argument.
%
% PROPOSAL: Print Exception causes recursively.
%   NOTE: Technically a tree structure, not a chain/sequence.


% Clear any previous instance of global variables
% -----------------------------------------------
% This is useful to avoid mistakenly using a previously initialized version of SETTINGS when the
% initialization has failed and when developing in MATLAB. Must be done as early as possible in the execution.
clear -global CONSTANTS SETTINGS    % Clearing obsoleted variable CONSTANTS for safety.

C = bicas.error_safe_constants();



try
    errorCode = C.EMIDP_2_INFO('NoError').errorCode;
    main_without_error_handling(varargin);

catch Exception1
    %================================================================
    % CASE: Caught an error in the regular execution of the software
    %================================================================
    try
        % IMPLEMENTATION NOTE: The error handling collects one long string with log/error messages for one bicas.log
        % call, instead of making multiple bicas.log calls. This avoids having stdout and stderr messages mixed
        % (alternating rows with stdout and stderr) in the MATLAB GUI, making it easier to read.
        msg = '';
        msg = [msg, sprintf('Main function caught an exception. Starting error handling.\n')];
        msg = [msg, sprintf('Exception1.identifier = "%s"\n', Exception1.identifier)];
        msg = [msg, sprintf('Exception1.message    = "%s"\n', Exception1.message)];

        %=================================================================================
        % Use MATLAB error message identifiers to identify one or multiple "error types".
        %=================================================================================
        msgIdentifierParts = strsplit(Exception1.identifier, ':');
        emidpList = msgIdentifierParts(C.EMIDP_2_INFO.isKey(msgIdentifierParts));    % Cell array of message identifier parts (strings) only.
        if isempty(emidpList)
            emidpList = {'UntranslatableErrorMsgId'};
        end

        %===================================
        % Print all identified error types.
        %===================================
        msg = [msg, sprintf('Matching MATLAB error message identifier parts (error types derived from Exception1.identifier):\n')];
        for i = 1:numel(emidpList)
            emidp = emidpList{i};
            msg  = [msg, sprintf('    %-23s : %s\n', emidp, C.EMIDP_2_INFO(emidp).description)];
        end
        % NOTE: Choice - Uses the last part of the message ID for determining error code to return.
        errorCode = C.EMIDP_2_INFO(emidpList{end}).errorCode;

        %======================
        % Print the call stack
        %======================
        callStackLength = length(Exception1.stack);
        msg = [msg, sprintf('MATLAB call stack:\n')];
        if (~isempty(callStackLength))
            for i=1:callStackLength
                stackCall = Exception1.stack(i);
                temp      = strsplit(stackCall.file, filesep);
                filename  = temp{end};

                %msg = [msg, sprintf('    %-27s %-52s row %i,\n', [filename, ','], [stackCall.name, ','], stackCall.line)];
                msg = [msg, sprintf('    row %3i, %-27s %-52s\n', stackCall.line, [filename, ','], stackCall.name)];
            end
        end

        msg = [msg, sprintf('Exiting MATLAB application with error code %i.\n', errorCode)];
        bicas.log('error', msg)
        return

    catch Exception2    % Deliberately use different variable name to distinguish the exception from the previous one.
        %========================================================
        % CASE: Caught an error in the regular error handling(!)
        %========================================================

        % NOTE: Only use very, very error-safe code here.
        %       Does not use bicas.log() or similar.
        fprintf(2, 'Error in the MATLAB code''s error handling.\n');   % Print to stderr.
        fprintf(2, 'exception2.identifier = "%s"\n', Exception2.identifier);          % Print to stderr.
        fprintf(2, 'exception2.message    = "%s"\n', Exception2.message);             % Print to stderr.

        % NOTE: The RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 3.4.3 specifies
        %   error code 0 : No error
        %   error code 1 : Every kind of error (!)
        errorCode = 1;

        fprintf(2, 'Exiting MATLAB application with error code %i.\n', errorCode);    % Print to stderr.
        return
    end
end



end    % main



% BICAS's de facto main function, without error handling.
function main_without_error_handling(cliArgumentsList)



startTimeTicSeconds = tic;

C = bicas.error_safe_constants();



%==================================
% ~ASSERTION: Check MATLAB version
%==================================
matlabVersionString = version('-release');
if ~ismember(matlabVersionString, C.PERMITTED_MATLAB_VERSIONS)
    error('BICAS:main:BadMatlabVersion', ...
        'Using bad MATLAB version. Found version "%s". BICAS requires any of the following MATLAB versions: %s.\n', ...
        matlabVersionString, strjoin(C.PERMITTED_MATLAB_VERSIONS, ', '))
end
bicas.logf('info', 'Using MATLAB, version %s.\n\n', matlabVersionString);



% Log that BICAS (the MATLAB code) has started running.
% RATIONALE: This is useful when one manually looks through the log file and tries to identify the beginning of a
% particular run. The BICAS log is always amended to and may therefore contain log messages from multiple runs.
bicas.logf('info', [...
    '###########################################\n', ...
    '###########################################\n', ...
    '#### BICAS'' MATLAB CODE STARTS RUNNING ####\n', ...
    '###########################################\n', ...
    '###########################################\n'])



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



%===============================
% Derive BICAS's directory root
%===============================
% ASSUMES: The current file is in the <BICAS>/src/+bicas/ directory.
[matlabSrcPath, ~, ~] = fileparts(mfilename('fullpath'));   % Use path of the current MATLAB file.
bicasRootPath         = EJ_library.utils.get_abs_path(fullfile(matlabSrcPath, '..', '..'));



%=======================================
% Log misc. paths and all CLI arguments
%=======================================
% IMPLEMENTATION NOTE: Want this as early as possible, before interpreting arguments. Should therefore not merge this
% with printing settings. This might help debug why settings were not set.
bicas.logf('info', 'BICAS software root path:  bicasRootPath = "%s"', bicasRootPath)
bicas.logf('info', 'Current working directory: pwd           = "%s"', pwd);   % Useful for debugging the use of relative directory arguments.
bicas.logf('info', '\nCOMMAND-LINE INTERFACE (CLI) ARGUMENTS TO BICAS\n')
bicas.logf('info',   '===============================================')
cliArgumentsQuotedList = {};
for i = 1:length(cliArgumentsList)
    % UI ASSERTION
    % IMPLEMENTATION NOTE: This check useful when calling BICAS from MATLAB (not bash).
    if ~ischar(cliArgumentsList{i})
        error('BICAS:main', 'Argument %i is not a string.', i)
    end
    
    bicas.logf('info', '    CLI argument %2i: "%s"', i, cliArgumentsList{i})
    cliArgumentsQuotedList{i} = ['''', cliArgumentsList{i}, ''''];
end
cliArgStrWhSpaceSep = strjoin(cliArgumentsQuotedList, ' ');
cliArgStrCommaSep   = strjoin(cliArgumentsQuotedList, ', ');
% IMPLEMENTATION NOTE: Printing the entire sequence of arguments, quoted with apostophe, is useful for copy-pasting to
% both MATLAB command prompt and bash.
bicas.logf('info', '    CLI arguments for copy-pasting:\n')
bicas.logf('info', '        Quoted, whitespace-separated: %s\n\n', cliArgStrWhSpaceSep)
bicas.logf('info', '        Quoted, comma-separated:      %s\n\n', cliArgStrCommaSep)
bicas.logf('info', '\n\n')



%========================================
% Initialize global settings & constants
%========================================
global SETTINGS
SETTINGS  = bicas.create_default_SETTINGS();



%=============================================
% First-round interpretation of CLI arguments
%=============================================
CliData = bicas.interpret_CLI_args(cliArgumentsList);



%=================================================
% Modify settings according to configuration file
%=================================================
if ~isempty(CliData.configFile)
    configFile = CliData.configFile;
else
    configFile = fullfile(bicasRootPath, C.DEFAULT_CONFIG_FILE_RELATIVE_PATH);
end
bicas.logf('info', 'configFile = "%s"', configFile)
rowList                 = EJ_library.utils.read_text_file(configFile, '(\r\n|\r|\n)');
ConfigFileSettingsVsMap = bicas.interpret_config_file(rowList);
bicas.log('info', 'Overriding subset of in-memory settings using config file.')
SETTINGS = overwrite_settings_from_strings(SETTINGS, ConfigFileSettingsVsMap, 'configuration file');    % Modify SETTINGS



%=========================================================
% Modify settings according to (inofficial) CLI arguments
%=========================================================
bicas.log('info', 'Overriding subset of in-memory settings using (optional, inofficial) CLI arguments, if any.')
SETTINGS = overwrite_settings_from_strings(SETTINGS, CliData.ModifiedSettingsMap, 'CLI arguments');    % Modify SETTINGS
SETTINGS.make_read_only();
% CASE: SETTINGS has now been finalized and is read-only (by assertion) after this.



%======================
% ASSERTIONS: SETTINGS
%======================
EJ_library.utils.assert.castring_regexp(SETTINGS.get_fv('SWD.release.version'), '[0-9]+\.[0-9]+\.[0-9]+')
EJ_library.utils.assert.castring_regexp(SETTINGS.get_fv('SWD.release.date'),    '20[1-3][0-9]-[01][0-9]-[0-3][0-9]')
% Validate S/W release version
% ----------------------------
% RCS ICD 00037, iss1rev2, Section 5.3 S/W descriptor file validation scheme implies this regex.
% NOTE: It is hard to thoroughly follow the description, but the end result should be under
% release-->version-->pattern (not to be confused with release_dataset-->version--pattern).
EJ_library.utils.assert.castring_regexp(SETTINGS.get_fv('SWD.release.version'), '(\d+\.)?(\d+\.)?(\d+)')



bicas.log('info', bicas.sprint_SETTINGS(SETTINGS))    % Prints/log the contents of SETTINGS.



SwModeDefs = bicas.swmode_defs(SETTINGS);



switch(CliData.functionalityMode)
    case 'version'
        print_version(SwModeDefs.List, SETTINGS)
        
    case 'identification'
        print_identification(SwModeDefs.List, SETTINGS)
        
    case 'S/W descriptor'
        print_sw_descriptor(SwModeDefs.List, SETTINGS)
        
    case 'help'
        print_help(SETTINGS)
        
    case 'S/W mode'
        %==============================================================================
        % CASE: Should be a S/W mode (deduced from elimination of other possibilities)
        %==============================================================================
        try
            SwModeInfo = SwModeDefs.get_sw_mode_info(CliData.swModeArg);
        catch Exception1
            % NOTE: Misspelled "--version" etc. would be interpreted as S/W mode and produce error here too.
            error('BICAS:main:CLISyntax', ...
                'Can not interpret first argument "%s" as a S/W mode (or any other legal first argument).', ...
                CliData.swModeArg);
        end



        %======================================================================
        % Parse CliData.SpecInputParametersMap arguments depending on S/W mode
        %======================================================================
        
        % Extract INPUT dataset files from SIP arguments.
        InputFilesMap = extract_rename_Map_keys(...
            CliData.SpecInputParametersMap, ...
            {SwModeInfo.inputsList(:).cliOptionHeaderBody}, ...
            {SwModeInfo.inputsList(:).prodFuncInputKey});
        
        % Extract OUTPUT dataset files from SIP arguments.
        OutputFilesMap = extract_rename_Map_keys(...
            CliData.SpecInputParametersMap, ...
            {SwModeInfo.outputsList(:).cliOptionHeaderBody}, ...
            {SwModeInfo.outputsList(:).prodFuncOutputKey});
        
        % ASSERTION: Assume correct number of arguments (the only thing not implicitly checked by extract_rename_Map_keys above).
        nSipExpected = numel(SwModeInfo.inputsList) + numel(SwModeInfo.outputsList);
        nSipActual   = numel(CliData.SpecInputParametersMap.keys);
        if nSipExpected ~= nSipActual
            error('BICAS:main:CLISyntax', 'Illegal number of "specific input parameters" (input & output datasets). Expected %i, but got %i.', nSipExpected, nSipActual)
        end        

        
        
        %================================
        % Set pipelineId, calibrationDir
        %================================
        % NOTE: Reading environment variables first here, where they are needed.
        calibrationDir = read_env_variable(SETTINGS, 'ROC_RCS_CAL_PATH',    'ENV_VAR_OVERRIDE.ROC_RCS_CAL_PATH');
        %pipelineId     = read_env_variable(SETTINGS, 'ROC_PIP_NAME',        'ENV_VAR_OVERRIDE.ROC_PIP_NAME');   % RGTS or RODP
        masterCdfDir   = read_env_variable(SETTINGS, 'ROC_RCS_MASTER_PATH', 'ENV_VAR_OVERRIDE.ROC_RCS_MASTER_PATH');
        bicas.logf('info', 'calibrationDir = "%s"', calibrationDir)
        %bicas.logf('info', 'pipelineId     = "%s"', pipelineId)
        bicas.logf('info', 'masterCdfDir   = "%s"', masterCdfDir)



        %==================
        % EXECUTE S/W MODE
        %==================
        bicas.execute_sw_mode( SwModeInfo, InputFilesMap, OutputFilesMap, masterCdfDir, calibrationDir, SETTINGS )

    otherwise
        error('BICAS:main:Assertion', 'Illegal value functionalityMode="%s"', functionalityMode)
end    % if ... else ... / switch



executionWallTimeSeconds = toc(startTimeTicSeconds);
bicas.logf('info', 'Time used for execution (wall time): %g [s]', executionWallTimeSeconds);    % Always log (-->critical)?
end



function NewMap = extract_rename_Map_keys(SrcMap, srcKeysList, newKeysList)
assert(numel(srcKeysList) == numel(newKeysList))
NewMap = containers.Map();

for i = 1:numel(srcKeysList)
    srcKey = srcKeysList{i};
    
    if ~SrcMap.isKey(srcKey)
        error('BICAS:main:Assertion', 'Can not find source key "%s"', srcKey)
    end
    NewMap(newKeysList{i}) = SrcMap(srcKey);
end
end



% Print software version on JSON format.
%
% RCS ICD 00037 iss1/rev2, draft 2019-07-11, Section 5.2:
% The example (but not the text) makes it clear that code should print JSON object.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created <<2019-08-05
%
function print_version(SwModeDefsList, SETTINGS)

% IMPLEMENTATION NOTE: Uses the software version in the S/W descriptor rather than the in the BICAS
% constants since the RCS ICD specifies that it should be that specific version.
% This is in principle inefficient but also "precise".

JsonSwd = bicas.get_sw_descriptor(SwModeDefsList, SETTINGS);

JsonVersion = [];
JsonVersion.version = JsonSwd.release.version;

strVersion = bicas.utils.JSON_object_str(JsonVersion, ...
    SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'));
bicas.stdout_print(strVersion);
end



% Print the JSON S/W descriptor identification section.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07
%
function print_identification(SwModesDefsList, SETTINGS)

JsonSwd = bicas.get_sw_descriptor(SwModesDefsList, SETTINGS);
strSwd = bicas.utils.JSON_object_str(JsonSwd.identification, ...
    SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'));
bicas.stdout_print(strSwd);

end



% Print the JSON S/W descriptor.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-07/2019-09-24
%
function print_sw_descriptor(SwModesDefsList, SETTINGS)

JsonSwd = bicas.get_sw_descriptor(SwModesDefsList, SETTINGS);
strSwd = bicas.utils.JSON_object_str(JsonSwd, ...
    SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'));
bicas.stdout_print(strSwd);

end



% Print user-readable "help text".
%
% NOTE: Useful if the output printed by this function can be used for copy-pasting into RCS User Manual (RUM).
%
function print_help(SETTINGS)
%
% PROPOSAL: Print CLI syntax incl. for all modes? More easy to parse than the S/W descriptor.

C = bicas.error_safe_constants();



% Print software name & description
bicas.stdout_printf('\n%s version %s\n', SETTINGS.get_fv('SWD.identification.name'), SETTINGS.get_fv('SWD.release.version') )
bicas.stdout_print(SETTINGS.get_fv('SWD.identification.description'))

%==========================
% Print error codes & types
%==========================
errorCodesList = cellfun(@(x) (x.errorCode), C.EMIDP_2_INFO.values);   % Array of (unsorted) error codes.
[~, iSort] = sort(errorCodesList);
empidList          = C.EMIDP_2_INFO.keys;
errorTypesInfoList = C.EMIDP_2_INFO.values;        % Cell array of structs (unsorted).
empidList          = empidList(iSort);
errorTypesInfoList = errorTypesInfoList(iSort);    % Cell array of structs sorted by error code.
bicas.stdout_printf('\nERROR CODES, ERROR MESSAGE IDENTIFIERS, HUMAN-READABLE DESCRIPTIONS\n')
bicas.stdout_printf(  '===================================================================')
for i = 1:numel(errorTypesInfoList)
    errorType = errorTypesInfoList{i};
    bicas.stdout_printf(['    %1i : %s\n', ...
                         '        %s\n'], errorType.errorCode, empidList{i}, errorType.description)
end

% Print settings
bicas.stdout_print(bicas.sprint_SETTINGS(SETTINGS))   % Includes title

bicas.stdout_printf('See "readme.txt" and user manual for more help.\n')
end



% Read environment variable, but allow the value to be overriden by a settings variable.
function v = read_env_variable(SETTINGS, envVarName, settingsOverrideName)
settingsOverrideValue = SETTINGS.get_fv(settingsOverrideName);

if isempty(settingsOverrideValue)
    v = getenv(envVarName);
else
    bicas.logf('info', 'Environment variable "%s" overridden by setting\n    %s = "%s"\n', envVarName, settingsOverrideName, settingsOverrideValue)
    v = settingsOverrideValue;
end

% UI ASSERTION
if isempty(v)
    error('BICAS:main:Assertion', ...
        'Can not set internal variable corresponding to environment variable "%s" from either (1) the environment variable, or (2) settings key value "%s".', ...
        envVarName, settingsOverrideValue)
end
end




% Modify multiple settings, where the values are strings but converted to numerics as needed. Primarily intended
% for updating settings with values from CLI arguments (which by their nature are initially strings).
%
%
% ARGUMENTS
% =========
% ModifiedSettingsAsStrings : containers.Map
%   <keys>   = Settings keys (strings). Must pre-exist as a SETTINGS key.
%   <values> = Settings values AS STRINGS.
%              Preserves the type of settings value for strings and numerics. If the pre-existing value is
%              numeric, then the argument value will be converted to a number.
%              Numeric row vectors are represented as a comma separated-list (no brackets), e.g. "1,2,3".
%              Empty numeric vectors can not be represented.
%
%
% NOTE/BUG: No good checking (assertion) of whether the string format of a vector makes sense.
%
function SETTINGS = overwrite_settings_from_strings(SETTINGS, ModifiedSettingsMap, valueSource)
    
    keysList = ModifiedSettingsMap.keys;
    for iModifSetting = 1:numel(keysList)
        key              = keysList{iModifSetting};
        newValueAsString = ModifiedSettingsMap(key);
        
        % ASSERTION
        if ~isa(newValueAsString, 'char')
            error('BICAS:settings:Assertion:IllegalArgument', 'Map value is not a string.')
        end
        
        %==================================================
        % Convert string value to appropriate MATLAB class.
        %==================================================
        switch(SETTINGS.get_setting_value_type(key))
            case 'numeric'
                newValue = textscan(newValueAsString, '%f', 'Delimiter', ',');
                newValue = newValue{1}';    % Row vector.
            case 'string'
                newValue = newValueAsString;
            otherwise
                error('BICAS:settings:Assertion:ConfigurationBug', 'Can not handle the MATLAB class=%s of internal setting "%s".', class(oldValue), key)
        end
        
        % Overwrite old setting.
        SETTINGS.update_value(key, newValue, valueSource);
    end
    
end
