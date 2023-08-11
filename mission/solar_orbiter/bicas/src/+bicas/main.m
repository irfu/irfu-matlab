%
% Main MATLAB function that launches BICAS.
%
% BICAS = BIAS CAlibration Software
%
% This function is BICAS' main MATLAB function, i.e. it is called by no other
% MATLAB code during regular use. It is intended to be wrapped in, and called
% from, non-MATLAB code, e.g. a bash script.
%
%
% IMPORTANT NOTE: DOCUMENTATION
% =============================
% Documentation important for using BICAS as a whole can be found in
% (1) "readme.txt" and other documentation text files (*.txt), and
% (2) the RCS SUM document.
% To prevent the duplication of documentation, comments in the source code tries
% to only cover subjects important for understanding the implementation and
% information not already present in other documentation text files.
%
%
% ARGUMENTS
% =========
% varargin
%       This function expects exactly the CLI arguments submitted to the bash
%       launcher script as a sequence of MATLAB strings. This function therefore
%       expects the arguments defined in the RCS ICD and possibly additional
%       inoffical arguments.
% Notes:
% - The official parameter syntax for S/W modes must be in agreement with
%   "roc_sw_descriptor.js" as specified by the RCS ICD.
% - The parameter syntax may contain additional unofficial parameters, which are
%   useful for development/debugging, but which are still compatible with the
%   RCS ICD.
% - The (MATLAB) code ignores but permits the CLI option --log.
%
%
% RETURN VALUE
% ============
% errorCode : The error code that is to be passed on to the OS/shell.
%
%
% NOTES
% =====
% ASSUMES: The current file is located in the <BICAS>/src directory.
%
% NOTE: This code is designed for MATLAB 2019b (as of 2020-01-20) but may very
% well work with other versions of MATLAB.
%
% IMPLEMENTATION NOTE: This code does not quit/exit using the MATLAB function
% "quit" since that always exits all of MATLAB which is undesirable when
% developing code from in the MATLAB IDE. The function returns the error code to
% make it possible for the bash wrapper to quit with an exit code instead.
%
% IMPLEMENTATION NOTE: The RCS ICD specifies what should go to stdout. BICAS
% 1) prints all log messages to stdout, and
% 2) prints all messages intended for BICAS' final, actual stdout (as produced
%    by the bash wrapper) to stdout but with a prefix so they can be filtered
%    out by the calling wrapper bash script.
% RATIONALE: See the bash wrapper script.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-03-xx.
%
function errorCode = main( varargin )
    %
    % PROPOSAL: Set option for MATLAB warnings. Disable?
    %   NOTE: TN claims warnings are sent to stdout.
    % TODO-NI: Is the application allowed to overwrite output files?
    %
    % PROPOSAL: Put a summarized version of CLI syntax in "bicas --help" (somethinger easier that the S/W descriptor).
    %   PRO: Useful when S/W descriptor becomes big and complex.
    %
    % PROPOSAL: Split up the "main_without_error_handling" function into several functions
    %           (outsource chunks of its code to smaller functions which are called).
    %   PROPOSAL: Printing CLI arguments.
    %
    % PROPOSAL: Better handling of errors in dataobj (reading CDF files).
    %   PROPOSAL: Wrap dataobj in function and catch and rethrow errors with BICAS' error IDs.
    %
    % PROPOSAL: Print MATLAB path (return value from path()).
    %   CON: Too many rows.
    %
    % PROPOSAL: Option for overriding settings via CLI argument in MATLAB using a containers.Map value.
    %   PROBLEM: Not obvious which order of precedence makes sense. Complicated to use order among settings arguments.
    %       PROPOSAL: Applies setting BEFORE CLI settings args.
    %       PROPOSAL: Applies setting AFTER CLI settings args.
    % 
    % PROPOSAL: Log some kind of indicator of de facto code version.
    %   PROBLEM: How handle if runs outside git repo?
    %   PROPOSAL: git commit (latest)
    %       NOTE: Different on irfu-matlab and bicas_ROC.
    %   PROPOSAL: git branch
    %
    % PROPOSAL: Use irf.str.assist_print_table more.
    %   Ex: Logging settings, CLI arguments(?), error codes & messages(?)
    
    
    try
        
        % NOTE: Permitting logging to file from MATLAB instead of bash wrapper
        % in case of using unofficial option.
        L = bicas.Logger('bash wrapper', true);
        
        
        
        %========================================================================
        % Initialize irfu-matlab "library"
        % --------------------------------
        % Among other things: Sets up paths to within irfu-matlab (excluding
        % .git/).
        %
        % NOTE: Prints to stdout. Can not deactivate this behaviour!
        % NOTE: Should not call irf('check') which looks for updates to
        %       irfu-matlab (can not distinguish between updates to BICAS or the
        %       rest of irfu-matlab).
        %
        % IMPLEMENTATION NOTE: bicas.Logger.ICD_log_msg uses
        % irf.str.add_prefix_on_every_row().
        % ==> Must initialize paths for irf BEFORE using
        %     bicas.Logger.log/logf.
        %========================================================================
        irf('check_path');
        
        %=======================================================================
        % IMPLEMENTATION NOTE: Disabling irf(...) commands that produce file
        %   ~/.matlab_datastore_<hostname>
        % since this causes some kind of problems for ROC (Quynh Nhu NGUYEN
        % when running BICAS in parallel).
        % See https://gitlab.obspm.fr/ROC/RCS/BICAS/-/issues/71 .
        %
        % Not clear why this is a problem since the commands should be able to
        % read the file if it has already been created. Since the function calls
        % are not truly necessary, they have simply been disabled. Of these
        % commands, irf('cdf_leapsecondstable') is the most useful one, but the
        % only functionality it adds is still only to check which is the most
        % recent version leap second table of (1) the built-in one, and (2)
        % CDF_LEAPSECONDSTABLE, and give a warning if CDF_LEAPSECONDSTABLE is
        % not the most recent one.
        % /Erik P G Johansson 2022-09-16
        %=======================================================================
        % Not strictly needed.
        % NOTE: Creates ~/.matlab_datastore_<hostname>
        % irf('check_os');
        % NOTE: Creates ~/.matlab_datastore_<hostname>
        % irf('matlab');
        % NOTE: Creates ~/.matlab_datastore_<hostname>
        % NOTE: Sets environment variable CDF_LEAPSECONDSTABLE if it has not
        %       already been set.
        % irf('cdf_leapsecondstable');
        % Print e.g. "irfu-matlab version: 2017-02-21,  v1.12.6".
        % irf('version')



        % Default error code (i.e. no error).
        errorCode = bicas.constants.EMIDP_2_INFO('NoError').errorCode;
        main_without_error_handling(varargin, L);
        
    catch Exception1
        %================================================================
        % CASE: Caught an error in the regular execution of the software
        %================================================================
        try
            msg = sprintf('Main function caught an exception. Starting error handling.\n');
        
            [msgRecursive, errorCode] = recursive_exception_msg(Exception1);
            msg = [msg, msgRecursive];
            
            msg = [msg, sprintf('Exiting MATLAB application with error code %i.\n', errorCode)];
            L.log('error', msg)

            return
            
        % Deliberately use different variable name to distinguish the exception
        % from the previous one.
        catch Exception2
            %========================================================
            % CASE: Caught an error in the regular error handling(!)
            %========================================================
            % NOTE: Only use very, very error-safe code here.
            %       Does not use bicas.Logger or similar.
            % NOTE: Prints to stderr (not stdout).
            
            fprintf(2, 'Error in the MATLAB code''s error handling.\n');
            fprintf(2, 'exception2.identifier = "%s"\n', Exception2.identifier);
            fprintf(2, 'exception2.message    = "%s"\n', Exception2.message);
            
            % NOTE: The RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section
            % 3.4.3 specifies
            %   error code 0 : No error
            %   error code 1 : Every kind of error (!)
            errorCode = 1;
            
            fprintf(2, ...
                'Exiting MATLAB application with error code %i.\n', errorCode);
            return
        end
    end
    
    
    
end    % main



% Create logging/error message for a given exception, and which is recursive in
% Exception.cause.
%
function [msg, errorCode] = recursive_exception_msg(Exception)
    
    CAUSES_RECURSIVE_INDENTATION_LENGTH = 8;

    % IMPLEMENTATION NOTE: The error handling collects one long string with
    % log/error messages for one bicas.Logger.log call, instead of making
    % multiple bicas.Logger.log calls. This avoids having stdout and stderr
    % messages mixed (alternating rows with stdout and stderr) in the MATLAB
    % GUI, making it easier to read.
    msg = '';
    msg = [msg, sprintf('Exception.identifier = "%s"\n', Exception.identifier)];
    msg = [msg, sprintf('Exception.message    = "%s"\n', Exception.message)];
    
    %=========================================================================
    % Use MATLAB error message identifiers to identify one or multiple "error
    % types".
    %=========================================================================
    msgIdentifierParts = strsplit(Exception.identifier, ':');
    % Cell array of message identifier parts (strings) only.
    emidpList = msgIdentifierParts(bicas.constants.EMIDP_2_INFO.isKey(msgIdentifierParts));
    if isempty(emidpList)
        emidpList = {'UntranslatableErrorMsgId'};
    end
    
    %===================================
    % Print all identified error types.
    %===================================
    msg = [msg, sprintf(...
        ['Matching MATLAB error message identifier parts (error types', ...
        ' derived from Exception1.identifier):\n'])];
    for i = 1:numel(emidpList)
        emidp = emidpList{i};
        msg  = [msg, sprintf('    %-23s : %s\n', ...
            emidp, bicas.constants.EMIDP_2_INFO(emidp).description)];
    end
    % NOTE: Choice - Uses the last part of the message ID for determining error
    % code to return.
    errorCode = bicas.constants.EMIDP_2_INFO(emidpList{end}).errorCode;
    
    %======================
    % Print the call stack
    %======================
    callStackLength = length(Exception.stack);
    msg = [msg, sprintf('MATLAB call stack:\n')];
    if (~isempty(callStackLength))
        for i=1:callStackLength
            stackCall = Exception.stack(i);
            temp      = strsplit(stackCall.file, filesep);
            filename  = temp{end};
            
            msg = [msg, sprintf('    row %4i, %-27s %-52s\n', ...
                stackCall.line, [filename, ','], stackCall.name)];
        end
    end
    
    for iCause = 1:numel(Exception.cause)
        msg = [msg, sprintf('Logging Exception.cause{%i}:\n', iCause)];
        
        %================
        % RECURSIVE CALL
        %================
        % NOTE: Does not capture return value errorCode.
        recursiveMsg = recursive_exception_msg(Exception.cause{iCause});
        
        recursiveMsg = irf.str.indent(recursiveMsg, ...
            CAUSES_RECURSIVE_INDENTATION_LENGTH);
        msg = [msg, recursiveMsg];
    end
    
end



% BICAS's de facto main function, without error handling.
%
function main_without_error_handling(cliArgumentsList, L)
    
    
    
    tTicToc = tic();
    
    
    
    %==================================
    % ~ASSERTION: Check MATLAB version
    %==================================
    matlabVersionString = version('-release');
    if ~ismember(matlabVersionString, bicas.constants.PERMITTED_MATLAB_VERSIONS)
        error('BICAS:BadMatlabVersion', ...
            ['Using bad MATLAB version. Found version "%s".', ...
            ' BICAS requires any of the following MATLAB versions: %s.\n'], ...
            matlabVersionString, ...
        strjoin(bicas.constants.PERMITTED_MATLAB_VERSIONS, ', '))
    end
    L.logf('info', 'Using MATLAB, version %s.\n\n', matlabVersionString);
    
    
    
    % Log that BICAS (the MATLAB code) has started running.
    % RATIONALE: This is useful when one manually looks through the log file and
    % tries to identify the beginning of a particular run. The BICAS log is
    % always amended to and may therefore contain log messages from multiple
    % runs.
    L.logf('info', [...
        '############################################\n', ...
        '############################################\n', ...
        '#### BICAS'' MATLAB CODE STARTS RUNNING ####\n', ...
        '############################################\n', ...
        '############################################\n'])



    % IMPLEMENTATION NOTE: Runs before irf(...) commands. Added after a problem
    % of calling irf('check_os') which indirectly calls system('hostname') at
    % ROC:roc2-dev.
    L.logf('debug', 'OS environment variable PATH                 = "%s"', ...
        getenv('PATH'));
    
    % NOTE: Useful for seeing which leap second table was actually used, e.g. at
    % ROC.
    L.logf('debug', 'OS environment variable CDF_LEAPSECONDSTABLE = "%s"', ...
        getenv('CDF_LEAPSECONDSTABLE'));
    

    
    %===============================
    % Derive BICAS's directory root
    %===============================
    % ASSUMES: The current file is in the <BICAS>/src/+bicas/ directory.
    % Use path of the current MATLAB file.
    [matlabSrcPath, ~, ~] = fileparts(mfilename('fullpath'));
    bicasRootPath         = irf.fs.get_abs_path(...
        fullfile(matlabSrcPath, '..', '..'));
    
    
    
    %=======================================
    % Log misc. paths and all CLI arguments
    %=======================================
    % IMPLEMENTATION NOTE: Want this as early as possible, before interpreting
    % arguments. Should therefore not merge this with printing settings. This
    % might help debug why settings were not set.
    L.logf('info', 'BICAS software root path:  bicasRootPath = "%s"', bicasRootPath)
    % Working directory useful for debugging the use of relative directory
    % arguments.
    L.logf('info', 'Current working directory: pwd           = "%s"', pwd);
    L.logf('info', '\nCOMMAND-LINE INTERFACE (CLI) ARGUMENTS TO BICAS\n')
    L.logf('info',   '===============================================')
    cliArgumentsQuotedList = {};
    for i = 1:length(cliArgumentsList)
        % UI ASSERTION
        % IMPLEMENTATION NOTE: This check useful when calling BICAS from MATLAB
        % (not bash).
        if ~ischar(cliArgumentsList{i})
            error('BICAS:CLISyntax', 'Argument %i is not a string.', i)
        end
        
        L.logf('info', '%2i: "%s"', i, cliArgumentsList{i})
        cliArgumentsQuotedList{i} = ['''', cliArgumentsList{i}, ''''];
    end
    cliArgStrWhSpaceSep = strjoin(cliArgumentsQuotedList, ' ');
    cliArgStrCommaSep   = strjoin(cliArgumentsQuotedList, ', ');
    % IMPLEMENTATION NOTE: Printing the entire sequence of arguments, quoted
    % with apostophe, is useful for copy-pasting to both MATLAB command prompt
    % and bash.
    L.logf('info', '\n')
    L.logf('info', 'CLI arguments for copy-pasting\n')
    L.logf('info', '------------------------------\n')
    L.logf('info', 'Single-quoted, whitespace-separated: %s\n\n', cliArgStrWhSpaceSep)
    L.logf('info', 'Single-quoted, comma-separated:      %s\n\n', cliArgStrCommaSep)
    L.logf('info', '\n\n')
    
    
    
    %========================================
    % Initialize global settings & constants
    %========================================
    SETTINGS  = bicas.create_default_SETTINGS();
    
    
    
    %=============================================
    % First-round interpretation of CLI arguments
    %=============================================
    CliData = bicas.interpret_CLI_args(cliArgumentsList);
    
    
    
    %==============================================================
    % Configure unofficial log file, written to from within MATLAB
    %==============================================================
    if ~isempty(CliData.matlabLogFile)
        % NOTE: Requires that bicas.Logger has been initialized to permit
        % writing to log file.
        L.set_log_file(CliData.matlabLogFile);
    else
        % There should be no log file, (generated from within MATLAB).
        L.set_log_file([]);
    end
    
    
    
    %=================================================
    % Modify settings according to configuration file
    %=================================================
    if ~isempty(CliData.configFile)
        configFile = CliData.configFile;
    else
        configFile = fullfile(...
            bicasRootPath, ...
            bicas.constants.DEFAULT_CONFIG_FILE_RELATIVE_PATH);
    end
    L.logf('info', 'configFile = "%s"', configFile)
    L.log('info', 'Overriding subset of in-memory settings using config file.')
    bicas.override_settings_from_config_file(configFile, SETTINGS, L)
    
    
    
    %=========================================================
    % Modify settings according to (unofficial) CLI arguments
    %=========================================================
    L.log('info', ...
        ['Overriding subset of in-memory settings using', ...
        ' (optional, unofficial) CLI arguments, if any.'])
    SETTINGS.override_values_from_strings(...
        CliData.ModifiedSettingsMap, 'CLI arguments');
    
    
    
    SETTINGS.make_read_only();
    % CASE: SETTINGS has now been finalized and is read-only (by assertion)
    % after this.
    
    
    
    % Print/log the content of SETTINGS.
    L.log('info', bicas.sprint_SETTINGS(SETTINGS))
    
    % Print/log selected parts of bicas.constants.
    L.log('info', sprint_constants())
    
    
    
    SwmDefs = bicas.swm.get_SWML(SETTINGS, L);
    
    
    
    switch(CliData.functionalityMode)
        case 'version'
            print_version(SwmDefs.List, SETTINGS)
            
        case 'identification'
            print_identification(SwmDefs.List, SETTINGS)
            
        case 'S/W descriptor'
            print_SWD(SwmDefs.List, SETTINGS)
            
        case 'help'
            print_help(SETTINGS)
            
        case 'S/W mode'
            %============================
            % CASE: Should be a S/W mode
            %============================
            try
                SwmInfo = SwmDefs.get_SWM(CliData.swmArg);
            catch Exception1
                % NOTE: Misspelled "--version" etc. would be interpreted as S/W
                % mode and produce error here too.
                error('BICAS:CLISyntax', ...
                    ['Can not interpret first argument "%s" as a S/W mode', ...
                    ' (or any other legal first argument).'], ...
                    CliData.swmArg);
            end
            
            
            
            %=================================================================
            % Parse CliData.SpecInputParametersMap arguments depending on S/W
            % mode
            %=================================================================
            
            % Extract INPUT dataset files from SIP arguments.
            InputFilesMap = extract_rename_Map_keys(...
                CliData.SpecInputParametersMap, ...
                {SwmInfo.inputsList(:).cliOptionHeaderBody}, ...
                {SwmInfo.inputsList(:).prodFuncInputKey});
            
            % Extract OUTPUT dataset files from SIP arguments.
            OutputFilesMap = extract_rename_Map_keys(...
                CliData.SpecInputParametersMap, ...
                {SwmInfo.outputsList(:).cliOptionHeaderBody}, ...
                {SwmInfo.outputsList(:).prodFuncOutputKey});
            
            % ASSERTION: Assume correct number of arguments (the only thing not
            % implicitly checked by extract_rename_Map_keys above).
            nSipExpected = numel(SwmInfo.inputsList) + numel(SwmInfo.outputsList);
            nSipActual   = numel(CliData.SpecInputParametersMap.keys);
            if nSipExpected ~= nSipActual
                error('BICAS:CLISyntax', ...
                    ['Illegal number of "specific input parameters"', ...
                    ' (input & output datasets). Expected %i, but got %i.'], ...
                    nSipExpected, nSipActual)
            end
            
            
            
            %==========================
            % Set rctDir, masterCdfDir
            %==========================
            % NOTE: Reading environment variables first here, where they are
            % needed.
            rctDir       = read_env_variable(SETTINGS, L, ...
                'ROC_RCS_CAL_PATH',    'ENV_VAR_OVERRIDE.ROC_RCS_CAL_PATH');
            masterCdfDir = read_env_variable(SETTINGS, L, ...
                'ROC_RCS_MASTER_PATH', 'ENV_VAR_OVERRIDE.ROC_RCS_MASTER_PATH');
            L.logf('info', 'rctDir       = "%s"', rctDir)
            L.logf('info', 'masterCdfDir = "%s"', masterCdfDir)

            irf.assert.dir_exists(rctDir)
            irf.assert.dir_exists(masterCdfDir)
            
            
            
            %===================
            % Read RCS NSO file
            %===================
            rcsNsoRelativePath = SETTINGS.get_fv('PROCESSING.RCS_NSO.FILE.RELATIVE_PATH');
            rcsNsoOverridePath = SETTINGS.get_fv('PROCESSING.RCS_NSO.FILE.OVERRIDE_PATH');
            if isempty(rcsNsoOverridePath)
                rcsNsoPath = fullfile(bicasRootPath, rcsNsoRelativePath);
            else
                rcsNsoPath = rcsNsoOverridePath;
            end
            
            %L.logf('info', 'rcsNsoPath = "%s"', rcsNsoPath);
            L.logf('info', 'Loading RCS NSO table XML file "%s"', rcsNsoPath)
            NsoTable = bicas.NSO_table(rcsNsoPath);



            %==================
            % EXECUTE S/W MODE
            %==================
            bicas.execute_SWM(...
                SwmInfo, InputFilesMap, OutputFilesMap, ...
                masterCdfDir, rctDir, NsoTable, SETTINGS, L )
            
        otherwise
            error('BICAS:Assertion', ...
                'Illegal value functionalityMode="%s"', functionalityMode)
    end    % if ... else ... / switch
    
    
    
    bicas.log_speed_profiling(L, 'main_without_error_handling', tTicToc);
end    % main_without_error_handling



function NewMap = extract_rename_Map_keys(SrcMap, srcKeysList, newKeysList)
    assert(numel(srcKeysList) == numel(newKeysList))
    NewMap = containers.Map();
    
    for i = 1:numel(srcKeysList)
        srcKey = srcKeysList{i};
        
        if ~SrcMap.isKey(srcKey)
            error('BICAS:Assertion', 'Can not find source key "%s"', srcKey)
        end
        NewMap(newKeysList{i}) = SrcMap(srcKey);
    end
end



% Print software version on JSON format.
%
% RCS ICD 00037 iss1/rev2, draft 2019-07-11, Section 5.2:
% The example (but not the text) makes it clear that code should print JSON
% object.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created <<2019-08-05
%
function print_version(SwmDefsList, SETTINGS)
    
    % IMPLEMENTATION NOTE: Uses the software version in the S/W descriptor
    % rather than the in the BICAS constants since the RCS ICD specifies that it
    % should be that specific version. This is in principle inefficient but also
    % "precise".
    
    JsonSwd = bicas.get_SWD(SwmDefsList);
    
    JsonVersion = [];
    JsonVersion.version = JsonSwd.release.version;
    
    strVersion = bicas.utils.JSON_object_str(JsonVersion, ...
        SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'));
    bicas.stdout_print(strVersion);
end



% Print the JSON S/W descriptor identification section.
%
% NOTE: Argument is *not* an instance of bicas.swm.SWML.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-06-07
%
function print_identification(SwmDefsList, SETTINGS)
    
    JsonSwd = bicas.get_SWD(SwmDefsList);
    strSwd = bicas.utils.JSON_object_str(JsonSwd.identification, ...
        SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'));
    bicas.stdout_print(strSwd);
    
end



% Print the JSON S/W descriptor.
%
% NOTE: Argument is *not* an instance of bicas.swm.SWML.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-06-07/2019-09-24
%
function print_SWD(SwmDefsList, SETTINGS)
    
    JsonSwd = bicas.get_SWD(SwmDefsList);
    strSwd = bicas.utils.JSON_object_str(JsonSwd, ...
        SETTINGS.get_fv('JSON_OBJECT_STR.INDENT_SIZE'));
    bicas.stdout_print(strSwd);
    
end



% Print user-readable "help text".
%
% NOTE: Useful if the output printed by this function can be used for
% copy-pasting into RCS User Manual (RUM).
% NOTE: No logging.
%
function print_help(SETTINGS)
    %
    % PROPOSAL: Print CLI syntax incl. for all modes? More easy to parse than the S/W descriptor.
    
    
    
    % Print software name & description
    bicas.stdout_print(sprint_constants());
    
    %==========================
    % Print error codes & types
    %==========================
    % Array of (unsorted) error codes.
    errorCodesList = cellfun(@(x) (x.errorCode), bicas.constants.EMIDP_2_INFO.values);
    [~, iSort] = sort(errorCodesList);
    empidList          = bicas.constants.EMIDP_2_INFO.keys;
    errorTypesInfoList = bicas.constants.EMIDP_2_INFO.values;   % Cell array of structs (unsorted).
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



% Read environment variable, but allow the value to be overriden by a settings
% variable.
function v = read_env_variable(SETTINGS, L, envVarName, overrideSettingKey)
    settingsOverrideValue = SETTINGS.get_fv(overrideSettingKey);
    
    if isempty(settingsOverrideValue)
        v = getenv(envVarName);
    else
        L.logf('info', ...
            'Environment variable "%s" overridden by setting\n    %s = "%s"\n', ...
            envVarName, overrideSettingKey, settingsOverrideValue)
        v = settingsOverrideValue;
    end
    
    % UI ASSERTION
    if isempty(v)
        error('BICAS:Assertion', ...
            ['Can not set internal variable corresponding to', ...
            ' environment variable "%s" from either', ...
            ' (1) the environment variable, or (2) settings key value %s="%s".'], ...
            envVarName, overrideSettingKey, settingsOverrideValue)
    end
end



% Create string for logging. Summarizes relevant constants in
% bicas.constantants, but not all.
%
% In practice, only prints s/w descriptor values. Primarily want to log the
% version information.
%
function s = sprint_constants()
    %
    % NOTE: Does not print error codes (bicas.constants), but print_help() does.
    % PROPOSAL: PERMITTED_MATLAB_VERSIONS
    
    s = sprintf([...
        '\n', ...
        'SELECTED (HARD-CODED) INTERNAL BICAS CONSTANTS\n', ...
        '==============================================\n']);
    
    keysCa = bicas.constants.SWD_METADATA.keys;   % Always row vector.
    keysCa = sort(keysCa)';   % Column vector.
    nKeys  = numel(keysCa);
    
    valuesCa = cell(nKeys, 1);
    for i = 1:nKeys
        valuesCa{i, 1} = bicas.constants.SWD_METADATA(keysCa{i});
    end
    [~, dataCa, columnWidths] = irf.str.assist_print_table(...
        {'Constant', 'Value'}, [keysCa, valuesCa], {'left', 'left'});
    
    for iRow = 1:size(dataCa, 1)
        s = [s, sprintf('%s = %s\n', dataCa{iRow, 1}, dataCa{iRow, 2})];
    end
    s = [s, newline];
end
