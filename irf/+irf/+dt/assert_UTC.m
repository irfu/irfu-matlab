%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function assert_UTC(Dt)
% PROPOSAL: Automatic test code.

assert(isa(Dt, 'datetime'), ...
  'Argument is not a datetime objet.')
assert(strcmp(Dt.TimeZone, 'UTCLeapSeconds'), ...
  'datetime object is not TimeZone=UTC.')

end
