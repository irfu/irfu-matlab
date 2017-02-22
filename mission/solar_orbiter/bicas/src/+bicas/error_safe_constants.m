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
function [ERROR_CODES, REQUIRED_MATLAB_VERSION] = error_safe_constants
%
% PROPOSAL: Return structure used for initializing global variables. Rename ~"safe_constants", "error_safe_constants".
%    PROPOSAL: Return one "SCONSTANT", "SCONST" containing ".ERROR_CODES" and ".REQUIRED_MATLAB_VERSION".
%       PRO: "global SCONST" shorter than "global ERROR_CODES".
%       CON: Longer error commands.
%       --
%       PROPOSAL: ".ERR_CODES", ".ERRORS".
%
% QUESTION: How distinguish between assertion errors and other errors?
%
% PROPOSAL: Redefine CDF_ERROR as CDF_READ_VALIDATION_ERROR - Something is wrong with the data in a read CDF (input CDF, master CDF).
% PROPOSAL: I/O error.
%

% NOTE: These constants are used by the error handling (the main function's catch section) and should therefore be
% available in that code.
%
% NOTE: These constants are MATLAB exit error codes which are passed to the wrapper bash script which uses them as exit
% codes.
ERROR_CODES = [];

ERROR_CODES.NO_ERROR                       = 0;          % NOTE: The RCS ICD specifies error code==0 <==> no error.

ERROR_CODES.MISC_ERROR                     = 1;
ERROR_CODES.ERROR_IN_MATLAB_ERROR_HANDLING = 2;          % QUESTION: Not for the launch scripts error handling?
%ERROR_CODES.UNKNOWN_ERROR                  = 3;          % Error, and does not know how to translate error identifer to error code. Rename? Abolish?

ERROR_CODES.CLI_SYNTAX_ERROR          = 100;             % Can not interpret command-line arguments syntax.
ERROR_CODES.OPERATION_NOT_IMPLEMENTED = 101;             % Execution has reached a portion of the code that has not been implemented yet.
ERROR_CODES.ASSERTION_ERROR           = 102;             % Detected an internal state that should never be possible in a bug-free code, ideally even with any possibly input.
                                                         % This should ideally indicate a pure code bug.
ERROR_CODES.PATH_NOT_FOUND            = 103;             % Directory or file does not exist.
ERROR_CODES.SW_MODE_PROCESSING_ERROR  = 104;
ERROR_CODES.DATASET_FORMAT_ERROR      = 105;             % Error when interpreting (official CDF) datasets, including master CDF files.
ERROR_CODES.CONFIGURATION_ERROR       = 106;             % Bad configuration (in particular hard-coded), e.g. constants, S/W descriptor.



REQUIRED_MATLAB_VERSION = '2016a';   % Value returned from "version('-release')".

end
