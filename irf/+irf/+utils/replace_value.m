%
% Replace every occurrence of a specific value (incl. NaN) in an array. Handles
% NaN as any other value.
%
%
% IMPLEMENTATION NOTE: Can not use "changem" since:
% (1) it does not handle replacement NaN-->x (needed when writing CDF files).
%     The reverse works though.
% (2) it does not check (assert) if the old or new values fit into the data
%     (given the data type).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-10-05
%
function x = replace_value(x, oldValue, newValue)
% TODO-DEC: Is it appropriate to use assert that x must be able to have
%                the value of oldValue (not newValue)? replace(x, NaN,
%                newValue) could be OK for integer x since it does not do
%                anything.

% ASSERTIONS
assert(isscalar(oldValue), 'Argument oldValue is not scalar.')
assert(isscalar(newValue), 'Argument newValue is not scalar.')
if ~isfloat(x) && (isnan(oldValue) || isnan(newValue))
  error('Using NaN for non-float data.')
end

% NOTE: Works for non-vectors.
if isnan(oldValue)
  b = isnan(x);
else
  b = (x==oldValue);
end
x(b) = newValue;

end
