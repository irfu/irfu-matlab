%
% Given a logical 1D array, find uninterrupted sequences of true.
%
%
% ARGUMENTS
% =========
% bArray
%       1D array. Logical or numeric.
%
%
% RETURN VALUES
% =============
% i1Array, i2Array
%       Nx1 arrays, numeric, same size. Indices such that
%       bArray(i1Array(i) : i2Array(i)) == true,
%       and covers all true values in bArray exactly once.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-26.
%
function [i1Array, i2Array] = split_by_false(bArray)
% PROPOSAL: Function name?
%   PROPOSAL: split_by_*
%   PROPOSAL: group_by_*
%   PROPOSAL: ~boolean
%   PROPOSAL: ~false? ~true
%   PROPOSAL:
%   NOTE: Compare strsplit, irf.utils.split_by_jumps

irf.assert.vector(bArray)

% Add false components to
% (1) make it easy to handle empty (zero length) array,
% (2) make it easy to handle array's beginning and end.
b = [false; bArray(:); false];

% NOTE: Indices same as in "b".
bBegin = ~b(1:end-1) &  b(2:end);
bEnd   =  b(1:end-1) & ~b(2:end);

% NOTE: i1/2Array indices same as in "bArray" therefore -1.
i1Array = find(bBegin) + 1 - 1;
i2Array = find(bEnd)   + 0 - 1;

% Enforce 0x1 arrays for empty arrays.
if isempty(i1Array)
  i1Array = ones(0,1);
  i2Array = i1Array;
end
end
