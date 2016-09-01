% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-08-05
%
% Wrapper around stdout_disp but prints with pattern+parameters.
%
function stdout_printf(pattern, varargin)

    stdout_disp(sprintf(pattern, varargin{:}));
    
end