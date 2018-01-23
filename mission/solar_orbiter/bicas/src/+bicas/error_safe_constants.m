%
% Initialize (some) hardcoded constants which are needed early, before regular settings constants are initialized, and
% thus need to be initialized in a way which is unlikely to trigger errors.
%
%
% NOTE: This file contains the authoritative definitions of the meaning of error codes that should be used in
% documentation.
%
%
% RETURN VALUES
% =============
% ERROR_TYPES_INFO : containers.Map
%                    key = Any one of the colon-separated parts of a error message identifier string (see "error" function).
%                    value = Struct with fields representing a type of error:
%                        .code        = The error code/number to be returned from BICAS' main function.
%                        .description = English human-readable text describing the error. Implicitly defines what kinds
%                                            of errors this error code should cover.
%                    IMPORTANT NOTE: A MATLAB error message identifier may match multiple "error types" (keys). The
%                    error-handling code (try-catch) should decide whether every message identifier should be used to
%                    identify only one error type if there are multiple ones to choose from.
% REQUIRED_MATLAB_VERSION        : String value that should be identical to the value returned by "version('-release')"
%                                  when using the correct MATLAB version.
% INOFFICIAL_ARGUMENTS_SEPARATOR : String. Argument that separates the official (ICD) arguments from inofficial
%                                  arguments.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
function [ERROR_TYPES_INFO, REQUIRED_MATLAB_VERSION, INOFFICIAL_ARGUMENTS_SEPARATOR] = error_safe_constants
%
% PROPOSAL: Redefine CDF_ERROR as CDF_READ_VALIDATION_ERROR - Something is wrong with the data in a read CDF (input CDF, master CDF).
% PROPOSAL: I/O error.
%
% PROPOSAL: Separate error code for hardcoded config error.
% PROPOSAL: Separate error code for not being able to identify error from error message identifier.
%
% NOTE: Uses capital letter for initialisms in message identifiers: CLI, SW.

MAP = containers.Map('KeyType', 'char', 'ValueType', 'any');

MAP('NoError')                      = info_struct(0, 'No error');  % NOTE: The RCS ICD specifies error code==0 <==> no error.

MAP('BadMatlabVersion')             = info_struct(100, 'Using the wrong MATLAB version.');
MAP('UntranslatableErrorMsgId')     = info_struct(101, 'Error occurred, but code can not translate the error''s MATLAB message identifier into any of BICAS'' internal standard error codes.');
MAP('MatlabCodeErrorHandlingError') = info_struct(102, 'The MATLAB code''s own error handling failed.');
MAP('CLISyntax')                    = info_struct(103, 'Can not interpret command-line interface (CLI) arguments syntax.');
MAP('PathNotFound')                 = info_struct(104, 'A specified directory or file does not exist.');
MAP('OperationNotImplemented')      = info_struct(105, 'Execution has reached a portion of the code that has not been implemented yet.');
MAP('Assertion')                    = info_struct(106, 'Detected an internal state that should never be possible in a bug-free code that receives correct inputs.');
MAP('IllegalArgument')              = info_struct(107, 'Argument passed to internal function had an illegal value.');
MAP('SWModeProcessing')             = info_struct(108, 'Error in s/w mode processing (processing data sets).');
MAP('DatasetFormat')                = info_struct(109, 'Error when interpreting (official CDF) datasets, including master CDF files.');
MAP('IllegalConfiguration')         = info_struct(110, 'Bad configuration (in particular hard-coded), e.g. constants, S/W descriptor. This should ideally indicate a pure code bug, i.e. it is not triggered by certain input.');

                     
ERROR_TYPES_INFO = MAP;

REQUIRED_MATLAB_VERSION = '2016a';

INOFFICIAL_ARGUMENTS_SEPARATOR = '--';

end



function errorTypeInfo = info_struct(errorCode, errorDescription)
    errorTypeInfo = struct('code', errorCode, 'description', errorDescription);
end
