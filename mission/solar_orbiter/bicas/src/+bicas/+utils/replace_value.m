% x = replace_value(x, oldValue, newValue)   Replace every occurrence of a specific value (incl. NaN) in an array.
%
% General-purpose routine, but created specifically for the purpose of converting fill/pad values<--->NaN when
% reading/writing CDF files.
%
% NOTE: Handles NaN as any other value.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-05
%
function x = replace_value(x, oldValue, newValue)
% QUESTION: Is it appropriate to use assert that x must be able to have the value of oldValue (not newValue)?
%     replace(x, NaN, newValue) could be OK for integer x since it does not do anything.

% ASSERTION
if ~isfloat(x) && (isnan(oldValue) || isnan(newValue))
    error('BICAS:replace_value:Assertion:IllegalArgument', 'Using NaN for non-float data.')
end

if isnan(oldValue)
    i = isnan(x);
else
    i = (x==oldValue);
end
x(i) = newValue;

end
