%
% Class used for logging. Create an instance which is then passed on to other
% functions for them to be able to print log messages.
%
% Print log message to MATLAB's stdout in standardized way.
%
% NOTE: This class provides some functionality that is meant to be useful when
%       calling BICAS code from outside BICAS in order to control the logging.
% Ex: Disable logging.
% Ex: Log to file from within MATLAB.
%
% NOTE: If logging to file from MATLAB, then one can not (in principle) log
%       MATLAB's own startup messages.
%
%
% RATIONALE: INSTANTIATED CLASS INSTEAD OF LOGGING FUNCTIONS
% ==========================================================
% Improve flexibility. No matter how complex the configuration of the logging
% is, it can be done once when configuring the class instance, and all the
% parameters are then passed along with the object and be used when generating
% log messages.
% --
% Ex: Switch between any combination of logging to
%   (1a) file and/or (1b) stdout, or (2) don't log at all.
% Ex: Switch between log prefixes or not.
% Ex: Non-BICAS code that uses BICAS code (e.g. bicas.proc.L1L2.cal.Cal) can have other
% logging, or none.
% Ex: Can implement accepting log messages before specifying the log file, by
% temporarily storing the messages.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-23
%
classdef Logger < handle
%
% PROBLEM: If logging to file from MATLAB, then can not log errors from not
%          being able to parse the CLI arguments, since log file is one of them.
%   CON: Can use "hack" to catch log file before doing the actual parsing.
%   CON: The types of CLI argument syntax errors not captured are very basic.
%        Does not include not finding paths.
%
% PROPOSAL: Log from bash as before for nominal execution, AND use separate
%           CLI log arguments to log to file from MATLAB.
%   Ex: --log-to-file-from-matlab <log file>
%
% PROPOSAL: irf.assert.trailing_LF
%   PROPOSAL: Simultaneously assert not trailing CR+LF.
%
% TODO-DEC: Should special, extra logging functionality be in this class or outside of it?
%   Ex: bicas.proc.L1L2.cal.utils.log_TF_function_handle
%   Ex: Logging for speed tests.
%   PRO: Shorter call. Can use L.method(...) instead of bicas.logfunc(L, ...)
%
% PROPOSAL: Way of emphasizing selected log messages.
%   PROPOSAL: Message is framed with "=".
%   PROPOSAL: All text is uppercase.
%   PROPOSAL: Unofficial logLevel='WARNING' (capitals) which is formatted and
%             translated into a regular logLevel='warning'.
%   PROPOSAL: Separate log method.



    properties(Access=private)
        % IMPLEMENTATION NOTE: Constant defined here and not centrally (e.g.
        % SETTINGS) to make sure that it is error-safe and always initialized.
        % Needed for early initialization and error handling (try-catch).
        %LOG_PREFIX = 'LOG FILE: ';
        
        LINE_FEED      = newline;
        
        stdoutOption   = 'none';
        
        logFileEnabled = false;
        logFileId      = [];
        logFileBuffer  = {};
    end



    methods(Access=public)

        % Constructor
        %
        %
        % ARGUMENTS
        % =========
        % stdoutOption
        %       String constant. Whether and how to log to stdout.
        %           'none'
        %               Do not log to stdout.
        %           'human-readable'
        %               Log to stdout as is most convenient for a human reader.
        %           'bash wrapper'
        %               Log to stdout as required by BICAS bash wrapper script.
        % logFileEnabled
        %       Logical/numerical. Whether to write to log file.
        %       NOTE: Log messages are stored in a buffer until the log file is
        %       specified later, using the designated method.
        %
        function obj = Logger(stdoutOption, logFileEnabled)
            % PROPOSAL: Separate arguments for stdout and log file behaviour
            %   CON: Want short call for no logging.
            
            % ASSERTIONS
            % IMPLEMENTATION NOTE: Assertion for number of arguments, since this
            % used to be variable.
            assert(nargin == 2, 'Not two arguments.')
            assert(isscalar(logFileEnabled) && islogical(logFileEnabled), ...
                'Argument logFileEnabled is not a scalar logical.')
            logFileEnabled = logical(logFileEnabled);
            
            switch(stdoutOption)
                case 'none'
                    obj.stdoutOption = 'none';
                    
                case 'human-readable'
                    obj.stdoutOption = 'human-readable';
                    
                case 'bash wrapper'
                    obj.stdoutOption = 'bash wrapper';
                    
                otherwise
                    error('BICAS:Assertion', 'Illegal argument "%s".', stdoutOption)
            end
            
            obj.logFileEnabled = logFileEnabled;
        end



        % Specify the log file to use, if logging to file has been enabled
        % (assertion). Previous log messages have been buffered.
        %
        % 
        % RATIONALE: NOT USING CONSTRUCTOR
        % ================================
        % It is useful to be able to specify log file AFTER that some logging
        % has been done.
        %
        %
        % ARGUMENTS
        % =========
        % logFile
        %       Path to log file.
        %       NOTE: If empty, then do not (ever) use log file. This option is
        %       useful if it is not known at the time of calling the constructor
        %       whether a log file should be used or not. This way the log
        %       message buffer can be cleared to potentially conserve RAM.
        % 
        function obj = set_log_file(obj, logFile)
            assert(obj.logFileEnabled, ...
                ['Trying to specify log file without having enabled', ...
                ' log file in constructor.'])
            assert(isempty(obj.logFileId), ...
                ['Trying to specify log file a second time by calling', ...
                ' this function a second time.'])

            if ~isempty(logFile)
                % CASE: Set log file.
                
                %===============
                % Open log file
                %===============
                % NOTE: Overwrite any pre-existing file.
                [fileId, fopenErrorMsg] = fopen(logFile, 'w');
                if fileId == -1
                    error(...
                        'BICAS:CanNotOpenFile', ...
                        ['Can not open log file "%s" for writing.', ...
                        ' fopen error message: "%s"'], ...
                        logFile, fopenErrorMsg)
                    % NOTE: Does not alter the object properties.
                end
                obj.logFileId = fileId;
                
                %=============================
                % Write buffered log messages
                %=============================
                for i = 1:numel(obj.logFileBuffer)
                    obj.write_to_log_file(obj.logFileBuffer{i});
                end
                obj.logFileBuffer = {};
                
            else
                % CASE: There should be no log file (despite constructor saying
                % there should/could be one).
                
                obj.logFileEnabled = false;
                obj.logFileBuffer  = {};
            end
        end



        % Fundemental method for logging. Other methods may wrap this method to
        % provide addition functionality.
        %
        %
        % ARGUMENTS
        % =========
        % logLevel
        %       String constant.
        %       NOTE: Value 'error' WILL NOT THROW ERROR. This is so that error
        %       handling code can log using this alternative. To actually THROW
        %       an error, use function error(...) or throw an exception
        %       directly.
        % msgStr
        %       Potentially multi-row string to be printed.
        %       NOTE: Multi-row strings must end with line feed (after last
        %       row).
        %
        %
        % Author: Erik P G Johansson, IRF, Uppsala, Sweden
        % First created 2019-07-26
        %
        function log(obj, logLevel, msg)
            % PROPOSAL: Be able to read "debug mode" flag so can choose whether to print or not.
            %   NOTE: Apropos RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.4 table.
            
            % RCS ICD compliant string.
            rcsIcdLogMsg = obj.ICD_log_msg(logLevel, msg);            

            %=================
            % Print to stdout
            %=================
            switch(obj.stdoutOption)
                case 'none'
                    % Do nothing
                    
                case 'human-readable'
                    obj.write_to_stdout(rcsIcdLogMsg)
                    
                case 'bash wrapper'
                    % String that is intended to be read by BICAS bash wrapper
                    % as stdout.
                    bashWrapperRecipientStr = irf.str.add_prefix_on_every_row(...
                        rcsIcdLogMsg, ...
                        bicas.const.LOG_FILE_PREFIX_TBW);
                    
                    obj.write_to_stdout(bashWrapperRecipientStr)
                    
                otherwise
                    error('BICAS:Assertion', ...
                        'Illegal property value obj.stdoutOption="%s".', ...
                        obj.stdoutOption)
            end
            
            %===================
            % Write to log file
            %===================
            if obj.logFileEnabled
                
                if isempty(obj.logFileId)
                    % CASE: Has not yet specified log file
                    obj.logFileBuffer{end+1} = rcsIcdLogMsg;
                else
                    % CASE: Has specified log file.
                    obj.write_to_log_file(rcsIcdLogMsg);
                end
            end
            
            %=========================================================
            % Print error messages to stderr (regardless of settings)
            %=========================================================
            if strcmp(logLevel, 'error')
                
                % Make sure string ends with line feed
                % ------------------------------------
                % IMPLEMENTATION NOTE: Necessary for stderr messages to end up
                % on separate lines. Not doing so produces some output rows (at
                % least inside the MATLAB GUI) with mixed stderr and std out
                % content which is hard to read.
                % NOTE: irf.str.add_prefix_on_every_row already does this
                % for the log messages.
                stderrStr = msg;
                if stderrStr(end) ~= obj.LINE_FEED
                    stderrStr = [stderrStr, obj.LINE_FEED];
                end

                % IMPLEMENTATION NOTE: Should not not use
                % bashWrapperRecipientStr, since it has the wrong prefix.
                obj.write_to_stderr(stderrStr)
            end

        end



        % Wrapper around bicas.Logger.log. Prints with pattern + parameters.
        %
        %
        % ARGUMENTS
        % =========
        % msgStr   : Log message, in the form of "format" argument to sprintf.
        % varargin : Variables used in msgStr
        %
        %
        % Author: Erik P G Johansson, IRF, Uppsala, Sweden
        % First created 2019-07-26
        %
        function logf(obj, logLevel, msgStr, varargin)
            
            obj.log( logLevel, sprintf(msgStr, varargin{:}) );
            
        end



    end    % methods(Access=public)
    
    
    
    methods(Access=private)
        
        
        
        function write_to_stdout(obj, str)
            % NOTE: Must print using function that reacts to trailing line feed.
            fwrite(1, str);
        end
        
        function write_to_stderr(obj, str)
            % NOTE: Must print using function that reacts to trailing line feed.
            fwrite(2, str);
        end
        
        function write_to_log_file(obj, str)
            % NOTE: Must print using function that reacts to trailing line feed.
            fwrite(obj.logFileId, str);
        end
        
        
        
        % NOTE: Partly defined by RCS ICD 00037, iss1/rev2, draft 2019-07-11,
        %       Section 4.2.3.
        % NOTE: RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.3 speaks
        %       of a "debug mode" not implemented here. Function always prints
        %       debug-level messages.
        % NOTE: Does NOT add bicas.const.LOG_FILE_PREFIX_TBW required for
        %       wrapper script to recognize log file messages in stdout. This is
        %       intentional since one may want both log message version with and
        %       without bicas.const.LOG_FILE_PREFIX_TBW.
        % NOTE: Current implementation could be a static method, but should
        %       probably remain an instance method as future changes might
        %       require it.
        %
        % RETURN VALUE
        % ============
        % rcsIcdLogMsg
        %       String that (1) end with line feed, (2) conforms to format for
        %       logs specified by the RCS ICD.
        %
        function rcsIcdLogMsg = ICD_log_msg(obj, logLevel, logMsg)
            
            switch(logLevel)
                case 'debug'   ; logLevelStr = 'DEBUG';
                case 'info'    ; logLevelStr = 'INFO';
                case 'warning' ; logLevelStr = 'WARNING';
                case 'error'   ; logLevelStr = 'ERROR';
                otherwise
                    error('BICAS:Assertion:IllegalArgument', ...
                        'Illegal logLevel="%s"', logLevel)
            end
            
            rcsIcdRowTimestamp = char(datetime("now","Format","uuuu-MM-dd'T'HH:mm:ss"));
            rcsIcdRowPrefix    = sprintf('%s -- %s -- ', ...
                rcsIcdRowTimestamp, logLevelStr);
            
            rcsIcdLogMsg = irf.str.add_prefix_on_every_row(...
                logMsg, rcsIcdRowPrefix);
        end



    end    % methods(Access=private)
    
    
    
end
