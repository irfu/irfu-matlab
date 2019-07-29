%
% Wrapper around bicas.log but prints with pattern+parameters.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-07-26
%
function logf(logLevel, varargin)

bicas.log( logLevel, sprintf(varargin{:}) );

end