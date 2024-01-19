%
% Wrapper around bicas.stdout_print() but prints with pattern+parameters.
% See that function.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-08-05
%
function stdout_printf(pattern, varargin)

bicas.stdout_print( sprintf(pattern, varargin{:}) );

end
