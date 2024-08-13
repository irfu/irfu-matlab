%
% Assert that object (arbitrary size) is datetime, UTC, and only contains
% midnight timestamps.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function assert_UTC_midnight(Dt)
% PROPOSAL: Implement using additional separate function for midnight,
%           independent of TimeZone.

irf.dt.assert_UTC(Dt)
assert(all(Dt == dateshift(Dt, 'start', 'day'), 'all'), ...
  'datetime object does not exclusively contain timestamps representing midnight.')
end
