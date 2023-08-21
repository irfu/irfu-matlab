%
% For a set of a set of N-D vectors with the same size in dimension 1, find all
% uninterrupted sequences of indices (in dimension 1) for which all variables
% are constant, i.e. var(i, :, :, ...) is identical for all values "i" in the
% sequence.
%
% NOTE: NaN counts as equal to itself.
% NOTE: Implementation uses irf.utils.find_equalities and therefore
% "isequaln", i.e.
%   -- IMPORTANT NOTE: Does not care about MATLAB class, not even recursively,
%      e.g. {'A'} == {65}.
%   -- Does not care about the order of fieldnames in structs
%   -- Counts NaN as equal to itself.
%
%
% ARGUMENTS
% =========
% varargin : At least one argument. All arguments must have the same size in
%            dimension 1.
%
%
% RETURN VALUES
% =============
% i1Array, i2Array : Numeric column arrays. First and last index in all
%                    sequences where arguments have constant values in the first
%                    dimension. Never contains zero-length sequences (e.g. end
%                    before beginning).
%                    NOTE: i1Array <= i2Array.
%                    NOTE: 0x1 if arguments have zero rows.
% n                : Number of subsequences. numel(i1Array).
%                    RATIONALE: This is provided as a pure convenience to the
%                    caller who would normally have to derive this anyway.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-08-13.
%
function [i1Array, i2Array, n] = split_by_change(varargin)
% NOTE: In principle, one might not want to split by every change.
%   Ex: BICAS UFV=1 ==> voltageNaN=true ==> One block, regardless of other
%       arguments.

% ASSERTION
% Require at least one argument, since size of return values is ~undefined
% (?!!) otherwise.
assert(numel(varargin) >= 1, 'Must have at least one argument.')
% NOTE: irf.utils.find_equalities checks that the arguments have the
% same number of rows.



fhArray = irf.utils.find_equalities(1, varargin{:});

if isempty(fhArray)
  i1Array = zeros(0,1);
  i2Array = zeros(0,1);
else
  bChange = ( fhArray(1:end-1) ~= fhArray(2:end) );
  iChange = find(bChange);
  i1Array = [1;       iChange+1     ];
  i2Array = [iChange; numel(fhArray)];
end

n = numel(i1Array);
end
