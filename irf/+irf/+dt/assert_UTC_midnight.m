%
% Assert that object (arbitrary size) is datetime, UTC, and only contains
% midnight timestamps.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function assert_UTC_midnight(Dt)
  assert(isa(Dt, 'datetime'))
  assert(strcmp(Dt.TimeZone, 'UTCLeapSeconds'), ...
    'datetime object is not TimeZone=UTC.')
  assert(all(Dt == dateshift(Dt, 'start', 'day'), 'all'), ...
    'datetime object does not exclusively contain timestamps representing midnight.')
end
