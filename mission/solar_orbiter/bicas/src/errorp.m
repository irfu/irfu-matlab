% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-02
%
% Alternative error function.
% This function is intended to be used instead of MATLAB's own "error" function to standardize the
% behaviour. Should be able to handle multiline messages(?).
%
% errorp = error prime ("prime" as in the apostrophe after a variable)
%
function errorp(error_code, message, varargin)
%
% PROPOSAL/TODO: Make it possible to submit a caught exception (thrown elsewhere).
%     Ex: get_abs_path
%

message = sprintf(['%i ', message], error_code, varargin{:});

fprintf(1, ['ERROR: ', message, '\n']); % Use log.irf('c', message) instead? Can handle multiline messages?
error(message)

end
