% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
% Alternative error function.
% This function is intended to be used instead of MATLAB's own "error" function to standardize the
% behaviour. Should be able to handle multiline messages(?).
%
% errorp = error prime ("prime" as in the apostrophe after a variable to denote a different version
% of about the same functionality)
%
% NOTE: Prints to stderr, but NOT to stdout. The bash wrapper script sends both stdout and
% stderr to the log.
%
function errorp(error_code, message, varargin)
%
% PROPOSAL/TODO: Make it possible to submit a caught exception (thrown elsewhere).
%     Ex: get_abs_path
%

message = sprintf(['%i ', message], error_code, varargin{:});

error(message)   % Prints to stderr.

end
