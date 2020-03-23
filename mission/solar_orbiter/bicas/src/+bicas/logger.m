%
% Class used for logging. Create an instance which is then passed on to other functions for them to be able
% to print log messages.
%
% Print log message to MATLAB's stdout in standardized way.
%
% NOTE: This class provides some functionality that is meant to be useful when calling BICAS code from outside BICAS in
% order to control the logging.
% Ex: Disable logging.
% Ex: Log to file from within MATLAB.
%
%
% RATIONALE: INSTANTIATED CLASS INSTEAD OF LOGGING FUNCTIONS
% ==========================================================
% Improve flexibility. No matter how complex the configuration of the logging is, it can be done once when configuring
% the class instance, and all the variables are then passed along with the object, or by using persistant/global
% variables (ugly!).
% Ex: Switch between any combination of logging to (1) file and/or (2) stdout, or (3) don't log at all.
% Ex: Switch between log prefixes or not.
% Ex: Non-BICAS code that uses BICAS code (e.g. bicas.calib) can have other logging, or none.
% Ex: Can implement accepting log messages before specifying the log file, by temporarily storing the messages.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-23
%
classdef logger < handle
%
% PROBLEM: If logging to file from MATLAB, then can not (in principle) log MATLAB's own startup messages.
% PROBLEM: If logging to file from MATLAB, then can not log errors from not being able to parse the CLI arguments, since
%          log file is one of them.
%   CON: Can use "hack" to catch log file before doing the actual parsing.
%   CON: The types of CLI argument syntax errors not captured are very basic. Does not include not finding paths.
%
% PROPOSAL: Log from bash as before for nominal execution, AND use separate CLI log arguments to log to file from MATLAB.
%   Ex: --log-to-file-from-matlab <log file>
%
% PROPOSAL: Move LOG_PREFIX to error_safe_constants.
%
% PROPOSAL: EJ_library.assert.trailing_LF
%   PROPOSAL: Simultaneously assert not trailing CR+LF.
%
% PROPOSAL: Initialization options
%   No logging.
%   Log to file, but specify it after instantiation.
%   Log to stdout.



    properties(Access=private)
        % IMPLEMENTATION NOTE: Constant defined here and not centrally (e.g. SETTINGS) to make sure that it is error-safe and
        % always initialized. Needed for early initialization and error handling (try-catch).
        LOG_PREFIX = 'LOG: ';
        
        LINE_FEED  = char(10);
        
        stdoutOption   = 'None';
        
        logFileEnabled = false;
        logFileId      = [];
        logFileBuffer  = {};
    end



    methods(Access=public)

        % Constructor
        %
        % ARGUMENTS
        % =========
        % No arguments   : No logging at all.
        % --
        % stdoutOption   : String constant.
        %                   'None'
        %                   'human-readable' : Log to stdout as is most convenient for a human reader.
        %                   'bash wrapper'   : Log to stdout as required by BICAS bash wrapper script.
        % logFileEnabled : Logical/numerical. Whether to write to log file.
        %                  NOTE: Log messages are stored in a buffer until the log file is specified later, using
        %                  special method.
        %
        function obj = logger(stdoutOption, logFileEnabled)
            % PROPOSAL: Separate arguments for stdout and log file behaviour
            %   CON: Want short call for no logging.
            
            if nargin == 0
                obj.stdoutOption   = 'None';
                obj.logFileEnabled = false;
                
            elseif nargin == 2
                assert(islogical(logFileEnabled) || isnumeric(logFileEnabled), 'Illegal argument logFileEnabled.')
                logFileEnabled = logical(logFileEnabled);
                
                switch(stdoutOption)
                    case 'None'
                        obj.stdoutOption = 'None';
                        
                    case 'human-readable'
                        obj.stdoutOption = 'human-readable';
                        
                    case 'bash wrapper'
                        obj.stdoutOption = 'bash wrapper';
                        
                    otherwise
                        error('Illegal argument "%s".', stdoutOption)
                end
                
                obj.logFileEnabled = logFileEnabled;
            else
                error('Illegal number of arguments.')
            end
        end



        % Specify the log file to use, if logging to file has been enabled. Previous log messages have been buffered.
        % 
        % RATIONALE: NOT USING CONSTRUCTOR
        % ================================
        % It is useful to be able to specify log file after that some logging has been done.
        % 
        function obj = set_log_file(obj, logFile)
            assert(obj.logFileEnabled,     'Trying to specify log file without having enabled log file in constructor.')
            assert(isempty(obj.logFileId), 'Trying to specify log file twice.')

            %===============
            % Open log file
            %===============
            % NOTE: Overwrite any pre-existing file.
            [fileId, fopenErrorMsg] = fopen(logFile, 'w');
            if fileId == -1
                error('BICAS:logger:Assertion', 'Can not open log file "%s". fopen error message: "%s"', logFile, fopenErrorMsg)
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
        end



        % ARGUMENTS
        % =========
        % logLevel : String constant.
        % msgStr   : Potentially multi-row string to be printed. NOTE: Multi-row strings must end with line feed.
        %
        %
        % Author: Erik P G Johansson, IRF, Uppsala, Sweden
        % First created 2019-07-26
        %
        function log(obj, logLevel, msg)
            % PROPOSAL: Be able to read "debug mode" flag so can choose whether to print or not.
            %   NOTE: Apropos RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.4 table.
            
            % RCS ICD compliant string.
            rcsIcdLogMsg            = obj.ICD_log_msg(logLevel, msg);            
            % String that is intended to be read by BICAS bash wrapper as stdout.
            bashWrapperRecipientStr = EJ_library.utils.add_prefix_on_every_row(rcsIcdLogMsg, obj.LOG_PREFIX);

            %=================
            % Print to stdout
            %=================
            switch(obj.stdoutOption)
                case 'None'
                case 'human-readable'
                    obj.write_to_stdout(rcsIcdLogMsg)
                case 'bash wrapper'
                    obj.write_to_stdout(bashWrapperRecipientStr)
                otherwise
                    error('BICAS:logger:Assertion', 'Illegal property value obj.stdoutOption="%s".', obj.stdoutOption)
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
                % IMPLEMENTATION NOTE: Necessary for stderr messages to end up on separate lines. Not doing so produces
                % some output rows (at least inside the MATLAB GUI) with mixed stderr and std out content which is hard
                % to read.
                % NOTE: EJ_library.utils.add_prefix_on_every_row already does this for the log messages.
                stderrStr = msg;
                if stderrStr(end) ~= obj.LINE_FEED
                    stderrStr = [stderrStr, obj.LINE_FEED];
                end

                % IMPLEMENTATION NOTE: Should not not use bashWrapperRecipientStr, since it has the wrong prefix.
                obj.write_to_stderr(stderrStr)
            end

        end



        % Wrapper around bicas.logger.log but prints with pattern + parameters.
        %
        %
        % ARGUMENTS
        % =========
        % msgStr   : Log message, in the form of "format" argument to sprintf.
        % varargin : Variables used in msgStr
        %
        %
        % Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
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
        
        
        
        % NOTE: Partly defined by RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.3.
        % NOTE: RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.3 speaks of a "debug mode" not implemented here.
        %       Function always prints debug-level messages.
        % NOTE: Does not add LOG_PREFIX required for wrapper script to recognize log messages in stdout. This is
        %       intentional since one may want both log message version with and without LOG_PREFIX.
        % NOTE: Could be a static method.
        %
        % RETURN VALUE
        % ============
        % rcsIcdLogMsg : String that (1) end with line feed, (2) conforms to format for logs specified by the RCS ICD.
        %
        function rcsIcdLogMsg = ICD_log_msg(obj, logLevel, logMsg)
            
            switch(logLevel)
                case 'debug'   ; logLevelStr = 'DEBUG';
                case 'info'    ; logLevelStr = 'INFO';
                case 'warning' ; logLevelStr = 'WARNING';
                case 'error'   ; logLevelStr = 'ERROR';
                otherwise
                    error('BICAS:logger:ICD_log_msg:Assertion:IllegalArgument', 'Illegal logLevel="%s"', logLevel)
            end
            
            rcsIcdRowTimestamp = datestr(clock, 'yyyy-mm-ddTHH:MM:SS');
            rcsIcdRowPrefix    = sprintf('%s -- %s -- ', rcsIcdRowTimestamp, logLevelStr);
            
            rcsIcdLogMsg = EJ_library.utils.add_prefix_on_every_row(logMsg, rcsIcdRowPrefix);
        end



    end    % methods(Access=private)
    
    
    
end    % classdef logger
