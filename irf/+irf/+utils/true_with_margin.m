%
% Given a sequence of logical values (1D array), arbitrarily placed on a
% coordinate axis x (numeric 1D array), add an x margin to the x-ordered
% sequences of true values.
%
%
% ARGUMENTS
% =========
% x
%       Numeric column array. Does not need to be sorted. Finite.
% b1
%       Logical column array.
% xMargin1, xMargin2
%       Scalar number. Positive, finite or +inf.
%       Margins towards lower and higher x values respectively.
%
%
% RETURN VALUES
% =============
% b2
%       Logical column array. Same as b1, except that additional b2 elements
%       have been set to true iff the corresponding x values are within a
%       distance xMargin1 below or xMargin2 above a true b1 value.
%       --
%       (b2(i) == true)
%       <==>
%       There is at least one j such that b1(j)==true and
%       -xMargin1 <= x(i)-x(j) <= xMargin2.
%
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-09.
%
function b2 = true_with_margin(x, b1, xMargin1, xMargin2)
%
% PROPOSAL: Algorithm: Find "boundary" elements, b1=true, but false next to it
%           (on at least one side). Iterate over and set b2=true for all elements
%           within range.

% IMPLEMENTATION NOTE: It has proven hard to implement the functionality without
% some kind of loop over elements. One can probably not have a loop with fewer
% iterations than this.

assert(iscolumn(x))
assert(iscolumn(b1))
assert(numel(x) == numel(b1), 'Arguments x and b1 do not have the same number of elements.')

assert(all(isfinite(x)))
assert(islogical(b1))
assert(isscalar(xMargin1) && (xMargin1 >= 0) && ~isnan(xMargin1))
assert(isscalar(xMargin2) && (xMargin2 >= 0) && ~isnan(xMargin2))



% Sort data in x order since using irf.utils.split_by_false() requires
% it. Might also imply that it finds fewer longer same-value intervals =>
% faster execution.
[x, iSortArray] = sort(x);
b1 = b1(iSortArray);

[i1Array, i2Array] = irf.utils.split_by_false(b1);

b2 = b1;
for i = 1:numel(i1Array)
  % For the given x interval (indexed "i") of continuous true value, set all
  % values for that x interval plus margins.
  xMin = (x(i1Array(i)) - xMargin1);
  xMax = (x(i2Array(i)) + xMargin2);
  b2 = b2 | ((xMin <= x) & (x <= xMax));
end

% Effectively use inverse permutation of iSortArray.
b2(iSortArray) = b2;
end
