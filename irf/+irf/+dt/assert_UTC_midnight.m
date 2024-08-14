%
% Assert that object (arbitrary size) is datetime, UTC, and only contains
% midnight timestamps.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function assert_UTC_midnight(Dt)

irf.dt.assert_UTC(Dt)
assert(all(irf.dt.is_midnight(Dt), 'all'), ...
  'datetime object does not exclusively contain timestamps representing midnight.')

end
