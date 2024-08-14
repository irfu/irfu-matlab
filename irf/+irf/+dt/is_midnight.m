%
% Check whether datetime array represents midnight.
%
% NOTE: Does not require any particular TimeZone.
%
%
% ARGUMENTS
% =========
% DT
%       datetime array.
%
%
% RETURN VALUES
% =============
% bIsMidnight
%       Logical array of the same size as argument.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function bIsMidnight = is_midnight(Dt)

bIsMidnight = dateshift(Dt, 'start', 'day') == Dt;

end
