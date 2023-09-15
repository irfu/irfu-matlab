%
% Given a sequence of logical values (1D array), arbitrarily placed on a
% coordinate axis x (numeric 1D array), add an x margin to the x-ordered
% sequences of true values.
%
%
% ARGUMENTS
% =========
% x
%       Numeric 1D array. Does not need to be sorted. Finite.
% b1
%       Logical 1D array.
% xMargin
%       Scalar, positive, finite or  +inf.
%
%
% RETURN VALUES
% =============
% b2
%       Logical. Same size as b1. Every true value is within an x distance
%       xMargin of a true v1 value. Every false value is not.
%       --
%       (b2(i) == true)
%       <==>
%       There is at least one j such that b1(j)==true and
%       abs(x(i)-x(j)) <= xMargin.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-09.
%
function b2 = true_with_margin(x, b1, xMargin)
%
% PROPOSAL: Extend to use arbitrary x2<>x1, associated with y2.
%   CON: Less well defined what that means. How handle if x2(i) is in the
%        middle of an uninterrupted y1=true sequence, but more than xMargin
%        away from nearest x1(j)?
%
% PROPOSAL: Find "boundary" elements, b1=true, but false next to it (on at
%           least one side). Iterate over and set b2=true for all elements
%           within range.

% IMPLEMENTATION NOTE: It has proven hard to implement without some kind of
% loop over elements. Can probably not have a loop with fewer iterations
% than this.

irf.assert.vector(x)
irf.assert.vector(b1)
assert(numel(x) == numel(b1), 'Arguments x and b1 do not have the same number of columns.')

assert(all(isfinite(x)))
assert(islogical(b1))
assert(isscalar(xMargin) && (xMargin >= 0) && ~isnan(xMargin))

% Sort data in x order since using irf.utils.split_by_false requires
% it. Might also imply that it finds fewer longer same-value intervals =>
% faster execution.
[x, iSortArray] = sort(x);
b1 = b1(iSortArray);

[i1Array, i2Array] = irf.utils.split_by_false(b1);

b2   = b1;
for i = 1:numel(i1Array)
  xMin = (x(i1Array(i)) - xMargin);
  xMax = (x(i2Array(i)) + xMargin);
  b2 = b2 | ((xMin <= x) & (x <= xMax));
end

% Effectively use inverse permutation of iSortArray.
b2(iSortArray) = b2;
end
