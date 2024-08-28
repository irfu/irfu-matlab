%
% Convert string(s) YYYY-MM-DD to UTC datetime object with timestamp at
% midnight. (E.g. format 2024-01-01T00:00:00.000Z, with a trailing Z, does NOT
% work, deliberately).
%
% This function is mostly useful for test code and hence the short name.
%
% UM = UTC Midnight
%
%
% ARGUMENTS
% =========
% strCa
%       Either
%       (1) String (one timestamp)
%       (2) Cell array of strings (timestamps for corresponding elements).
%
%
% RETURN VALUES
% =============
% datetime with the same size as the argument.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function Dt = um(strCa)
% PROPOSAL: Rename UTC*.
%   PRO: More consistent with irf.dt.UTC().
%   CON: Longer.
%     PRO: Used many times in tests.

assert(iscell(strCa) || ischar(strCa))
% NOTE: datetime() also accepts other datetime objects, with any
% time-of-day. "Must" therefore forbid.

Dt          = datetime(strCa);
Dt.TimeZone = 'UTCLeapSeconds';

irf.dt.assert_UTC_midnight(Dt)
end
