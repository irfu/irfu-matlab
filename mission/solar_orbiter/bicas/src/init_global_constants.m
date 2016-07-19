% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
% Initialize (some) global constants.
%
% IMPORTANT: Only intented for initialization which is so trivial that no error handling (try-catch) is
% necessary so that its constants can be used in the error handling and be called outside try-catch.
%
% IMPLEMENTATION NOTE: It is useful to have this code separate so that it can be called separately
% before calling other functions separately (for testing; without launching the main function).
%
function  init_global_constants
%
% PROPOSAL: Return structure (lazy eval.) rather than initializing global constants.



% NOTE: These constants are used by the error handling (the main function's catch section) and
% should therefore be available in that code.
%
% NOTE: These constants are MATLAB exit error codes which are passed to the wrapper bash script
% which use them as exit codes.
global ERROR_CODES
ERROR_CODES.NO_ERROR = 0;
ERROR_CODES.MISC_ERROR = 1;
ERROR_CODES.UNKNOWN_ERROR = 2;          % Only use for error in error handling?
ERROR_CODES.CLI_ARGUMENT_ERROR = 100;   % Can not interpret command-line arguments.
ERROR_CODES.OPERATION_NOT_IMPLEMENTED = 101;   % Execution has reached a portion of the code that has not been implemented yet.
ERROR_CODES.ASSERTION_ERROR = 102;      % Detected an internal state that never be possible. This should indicate a pure code bug.
ERROR_CODES.PATH_NOT_FOUND = 103;       % Directory or file does not exist.
ERROR_CODES.SW_MODE_PROCESSING_ERROR = 104;

global REQUIRED_MATLAB_VERSION
REQUIRED_MATLAB_VERSION = '2016a';  % Value returned from "version('-release')".

end

