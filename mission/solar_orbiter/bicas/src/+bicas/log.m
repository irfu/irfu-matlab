%
% Print log message to MATLAB's stdout in standardized way.
% 
%
% NOTE: Can handle multi-row strings.
% NOTE: Adds its own constant prefix to row so that wrapper bash script can filter out the rows.
% NOTE: Partly defined by RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.3.
% NOTE: RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.3 speaks of a "debug mode" not implemented here.
%       Function always prints debug-level messages.
%
%
% ARGUMENTS
% =========
% logLevel : String constant.
% msgStr   : Potentially multi-row string to be printed. NOTE: Multi-row strings must end with line feed.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-07-26
%
function log(logLevel, msgStr)
% PROPOSAL: Be able to read "debug mode" flag so can choose whether to print or not.
%   NOTE: Apropos RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 4.2.4 table.
% PROPOSAL: Move LOG_PREFIX to error_safe_constants.

% IMPLEMENTATION NOTE: Constant defined here and not centrally (e.g. SETTINGS) to make sure that it is error-safe and
% always initialized. Needed for early initialization and error handling (try-catch).
LOG_PREFIX = 'LOG: ';



switch(logLevel)
    case 'debug'
        logLevelStr = 'DEBUG';
    case 'info'
        logLevelStr = 'INFO';
    case 'warning'
        logLevelStr = 'WARNING';
    case 'error'
        logLevelStr = 'ERROR';
    otherwise
        error('BICAS:log:Assertion:IllegalArgument', 'Illegal logLevel="%s"', logLevel)
end

timestamp = datestr(clock, 'yyyy-mm-ddTHH:MM:SS');
rowPrefix = sprintf('%s%s -- %s -- ', LOG_PREFIX, timestamp, logLevelStr);

printStr = bicas.utils.add_prefix_on_every_row(msgStr, rowPrefix);



%=============================
% Print log message to stdout
%=============================
fwrite(1, printStr);    % NOTE: Must print using function that reacts to trailing line feed.

%=============================================
% Additionally print error messages to stderr
%=============================================
if strcmp(logLevel, 'error')
    
    % Make sure string ends with line feed.
    % IMPLEMENTATION NOTE: Necessary for stderr messages to end up on separate lines. Not doing so produces some output
    % rows (at least inside the MATLAB GUI) with mixed stderr and std out content which is hard to read.
    % NOTE: bicas.utils.add_prefix_on_every_row already does this for the log messages.
    LINE_FEED = char(10); 
    if msgStr(end) ~= LINE_FEED
        msgStr = [msgStr, LINE_FEED];
    end
    

    fwrite(2, msgStr);    % NOTE: Must print using function that reacts to trailing line feed.
    % IMPLEMENTATION NOTE: Can not print printStr, since it has the wrong prefix, interpreted by the bash wrapper script.
end

end