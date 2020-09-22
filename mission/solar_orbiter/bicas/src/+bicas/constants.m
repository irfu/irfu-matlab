%
% Hard-coded constants which are needed
% -- for error handling
% -- early, before regular settings are initialized,
% and which  thus need to be initialized independent of settings and in a way which is unlikely to trigger errors.
%
% NOTE: This file contains the authoritative definitions of the meaning of error codes that should (maybe) be used in
% documentation.
%
%
% VARIABLE NAMING CONVENTIONS
% ===========================
% TBW   = To Bash Wrapper.
% EMIDP = (MATLAB) Error Message Identifier Part. One of the colon-separated parts of the MException .identifier
%         string field (error message ID).
%         NOTE: "Component" has a special meaning in the context of error message IDs. Therefore uses the term
%         "part" instead.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-07-09, as a replacement for the FUNCTION error_safe_constant created 2016-06-02.
%
classdef constants   % < handle
    % PROPOSAL: Error category for bad input datasets (both science and HK).
    %   PRO: Has similar for RCTs.
    
    
    
    properties(Constant)

        % Permissible string values returned by "version('-release')" when using
        % the correct MATLAB version.
        %
        % NOTE: BICAS originally required MATLAB R2016a but no longer does.
        % NOTE: ROC only needs MATLAB R2019b. Source:
        % https://gitlab.obspm.fr/ROC/RCS/BICAS/issues/2#note_10804
        PERMITTED_MATLAB_VERSIONS         = {'2019b'};        

        % Path to default config file relative to BICAS's directory root. Note
        % that this is also implicitly the constant for the default config file
        % filename.
        DEFAULT_CONFIG_FILE_RELATIVE_PATH = fullfile('config', 'bicas.conf');
        
        % MATLAB stdout prefix to signal to bash wrapper that the log message should be passed on to STDOUT (without the
        % prefix).
        STDOUT_PREFIX_TBW = 'STDOUT: ';
        
        % MATLAB stdout prefix to signal to bash wrapper that the log message should be passed on to LOG FILE (without
        % the prefix).
        LOG_FILE_PREFIX_TBW = 'LOG FILE: ';

        % Information to "interpret" and "translate" captured exceptions
        % --------------------------------------------------------------
        % containers.Map with
        %   key   = Any one of the colon-separated parts of a MATLAB error message identifier string (see
        %           "error" function).
        %   value = Struct with fields representing a type of error:
        %       .errorCode   = The error code/number to be returned from BICAS' main function.
        %                      IMPORTANT NOTE: A MATLAB error message identifier may match multiple "error types"
        %                      (keys). The error-handling code (try-catch) should decide whether every message
        %                      identifier should be used to identify only one error type if there are multiple ones to
        %                      choose from.
        %       .description = English human-readable text describing the error. Implicitly defines what
        %                      kinds of errors this error code should cover.
        %
        %
        EMIDP_2_INFO = bicas.constants.init_EMIDP_2_INFO;
        
        
        
        % Regular expression that the CLI name of a s/w mode must satisfy.
        %
        % The RCS ICD 00037, iss1rev2, draft 2019-07-11, section 5.3 seems to imply this regex for S/W mode
        % CLI parameters: ^[A-Za-z][\\w-]+$
        % NOTE: Only one backslash in MATLAB regex as opposed to in the RCS ICD.
        %
        % NOTE: Must not begin with "--" to be confused with CLI options, but the above constraint ensures this.
        %
        % NOTE: help regexp: "\w    A word character [a-z_A-Z0-9]"
        %
        SW_MODE_CLI_OPTION_REGEX = '[A-Za-z][\w-]+';
        
    end    % properties(Constant)
    
    

    methods(Static, Access=private)
        
        

        function MAP = init_EMIDP_2_INFO()
            % NOTE: The RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 3.4.3 specifies
            %   error code 0 : No error
            %   error code 1 : Every kind of error (!)
            
            MAP = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            MAP('NoError')                      = init_struct(0, 'No error');
            MAP('BadMatlabVersion')             = init_struct(1, 'Using the wrong MATLAB version.');
            MAP('UntranslatableErrorMsgId')     = init_struct(1, 'Error occurred, but code can not translate the error''s MATLAB message identifier into any of BICAS''s internal standard error codes.');
            MAP('MatlabCodeErrorHandlingError') = init_struct(1, 'The MATLAB code''s own error handling failed.');
            MAP('CLISyntax')                    = init_struct(1, 'Can not interpret command-line interface (CLI) arguments syntax.');
            MAP('PathNotFound')                 = init_struct(1, 'A specified directory or file does not exist.');
            MAP('OperationNotImplemented')      = init_struct(1, 'Execution has reached a portion of the code that has not been implemented yet.');
            MAP('Assertion')                    = init_struct(1, 'Detected an internal state that should never be possible in a bug-free code that receives correct inputs.');
            MAP('IllegalArgument')              = init_struct(1, 'An argument to an internal function had an illegal value.');
            MAP('SWModeProcessing')             = init_struct(1, 'Error in s/w mode processing (processing datasets).');
            MAP('DatasetFormat')                = init_struct(1, 'Error when interpreting (official CDF) input datasets, including master CDF files.');
            MAP('IllegalCodeConfiguration')     = init_struct(1, 'Bad hard-coded configuration (or possibly configurable setting but should not be), e.g. constants, S/W descriptor. This should ideally indicate a pure code bug, i.e. it is not triggered by certain user-controlled input.');
            MAP('CannotInterpretConfigFile')    = init_struct(1, 'Can not interpret the content of the configuration file. This implies a problem with the syntax.');
            MAP('ConfigurationBug')             = init_struct(1, 'Trying to configure BICAS in an illegal way.');
            MAP('FailedToReadInterpretRCT')     = init_struct(1, 'Can not interpret the content of the calibration file (RCT) file, e.g. because the RCT contains invalid calibration values.');
            MAP('FailedToReadInterpretNsOps')   = init_struct(1, 'Can not read or interpret the content of the non-standard operations file.');
            MAP('CannotFindRegexMatchingRCT')   = init_struct(1, 'Can not find any matching calibration file to read. No file matches regular expression.');
            
            % IMPLEMENTATION NOTE: Using a nested function merely to keep the function call short.
            function ErrorTypeInfo = init_struct(errorCode, errorDescription)
                ErrorTypeInfo = struct(...
                    'errorCode',   errorCode, ...
                    'description', errorDescription);
            end
            
        end
        
        
        
    end    % methods(Static)

end
