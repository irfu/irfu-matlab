%
% Initialize (some) hardcoded constants which are needed
% -- error handling
% -- early, before regular settings are initialized,
% and which  thus need to be initialized independent of settings and in a way which is unlikely to trigger errors.
%
%
% NOTE: This file contains the authoritative definitions of the meaning of error codes that should (Maybe) be used in
% documentation.
%
%
% RETURN VALUES
% =============
% EMIDP_2_INFO : containers.Map. 
%                      key = Any one of the colon-separated parts of a MATLAB error message identifier string (see
%                            "error" function).
%                      value = Struct with fields representing a type of error:
%                           .code        = The error code/number to be returned from BICAS' main function.
%                           .description = English human-readable text describing the error. Implicitly defines what
%                                          kinds of errors this error code should cover.
%                      IMPORTANT NOTE: A MATLAB error message identifier may match multiple "error types" (keys). The
%                      error-handling code (try-catch) should decide whether every message identifier should be used to
%                      identify only one error type if there are multiple ones to choose from.
% REQUIRED_MATLAB_VERSION        : String value that should be identical to the value returned by "version('-release')"
%                                  when using the correct MATLAB version.
% INOFFICIAL_ARGUMENTS_SEPARATOR : String. Argument that separates the official (RCS ICD) arguments from inofficial
%                                  arguments.
%
%
% DEFINITIONS
% ===========
% EMIDP : (MATLAB) Error Message Identifier Part. One of the colon-separated parts of the exception .identifier string
%         field.
%         NOTE: "Component" has a special meaning in that context. Therefore uses term "part" instead.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
function C = error_safe_constants
% PROPOSAL: Change name. Something not just "error safe".
%   PROPOSAL: ~constants. (Change name of bicas.constants).

% NOTE: The RCS ICD 00037, iss1/rev2, draft 2019-07-11, Section 3.4.3 specifies
%   error code 0 : No error
%   error code 1 : Every kind of error (!)
MAP = containers.Map('KeyType', 'char', 'ValueType', 'any');
MAP('NoError')                      = info_struct(0, 'No error');
MAP('BadMatlabVersion')             = info_struct(1, 'Using the wrong MATLAB version.');
MAP('UntranslatableErrorMsgId')     = info_struct(1, 'Error occurred, but code can not translate the error''s MATLAB message identifier into any of BICAS''s internal standard error codes.');
MAP('MatlabCodeErrorHandlingError') = info_struct(1, 'The MATLAB code''s own error handling failed.');
MAP('CLISyntax')                    = info_struct(1, 'Can not interpret command-line interface (CLI) arguments syntax.');
MAP('PathNotFound')                 = info_struct(1, 'A specified directory or file does not exist.');
MAP('OperationNotImplemented')      = info_struct(1, 'Execution has reached a portion of the code that has not been implemented yet.');
MAP('Assertion')                    = info_struct(1, 'Detected an internal state that should never be possible in a bug-free code that receives correct inputs.');
MAP('IllegalArgument')              = info_struct(1, 'An argument to an internal function had an illegal value.');
MAP('SWModeProcessing')             = info_struct(1, 'Error in s/w mode processing (processing data sets).');
MAP('DatasetFormat')                = info_struct(1, 'Error when interpreting (official CDF) datasets, including master CDF files.');
MAP('IllegalCodeConfiguration')     = info_struct(1, 'Bad hard-coded configuration (or possibly configurable setting but should not be), e.g. constants, S/W descriptor. This should ideally indicate a pure code bug, i.e. it is not triggered by certain user-controlled input.');
MAP('CannotInterpretConfigFile')    = info_struct(1, 'Can not interpret the content of the configuration file. This implies a problem with the syntax.');
C.EMIDP_2_INFO = MAP;

C.REQUIRED_MATLAB_VERSION           = '2016a';
C.INOFFICIAL_ARGUMENTS_SEPARATOR    = '--';
C.DEFAULT_CONFIG_FILE_RELATIVE_PATH = fullfile('config', 'bicas.conf');    % Path (incl. filename) to default config file. Relative to BICAS's directory root.

end



function errorTypeInfo = info_struct(errorCode, errorDescription)
    errorTypeInfo = struct('errorCode', errorCode, 'description', errorDescription);
end
