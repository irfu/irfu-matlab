%
% Given a set of same-sized arrays in one dimension, find all sets of indices
% for which the array components are all equal.
%
% NOTE: Counts NaN as equal to itself.
% NOTE: Is faster the more nearby indices (within searchDistance) are equal.
% NOTE: Implementation uses "isequaln", i.e.
%   -- IMPORTANT NOTE: Does not care about MATLAB class, not even recursively,
%      e.g. {'A'} == {65}.
%   -- Does not care about the order of fieldnames in structs
%   -- Counts NaN as equal to itself.
%
%
% NOTES ON PERFORMANCE
% ====================
% NOTE: The performance of the function should theoretically depend very much on
% (1) which rows are equal or not,
% (2) the number of input variables,
% (3) the data types of input variables,
% (4) the order if arguments,
% (5) argument "searchDistance".
% --
% Empirically, the speed varies widely depending on data type:
% (For vectors of unique values:)
%   time_datetime / time_double ~ 270x
%   time_string   / time_double ~ 15x-16x
%   time_int64    / time_double ~ 1x   (about same speed)
%
%
% ARGUMENTS
% =========
% searchDistance
%       Distance (in dimension 1) over which equalities will be searched for.
%       NOTE: Allowed to be inf.
% varargin
%       Arbitrary number of arrays with equal number of rows (size N in
%       dimension 1).
%
%
% RETURN VALUE
% ============
% fauxHashArray
%       Nx1 numeric array.
%       fauxHashArray(i) == fauxHashArray(j): Data are equal, always.
%       fauxHashArray(i) <> fauxHashArray(j): Data are inequal, if
%                                             abs(i-j) <= searchDistance
%       NOTE: The exact numerical values used may change as there is no
%       obviously natural way of setting them.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% Initially created 2020-04-10
%
function [fauxHashArray] = find_equalities(searchDistance, varargin)
%
% PROPOSAL: Better name
%   PROPOSAL: Imply equality between index values, not variables.
% PROPOSAL: Speed test.
%
% PROPOSAL: Terminology change: faux hash=fauxHash --> pseudohash
%
% PROPOSAL: Recursive call for each array separately. Then run algorithm to merge fhArrays.
%   PRO: More likely to find more equalities within each array separately. ==> Potentially faster.
%   PRO: Can merge multiple Nx1 fhArrays to one NxM array. ==> Potentially faster
%       PRO: Useful for hashing (unclear if makes sense though).
%   ~CON: The algorithm as it stands is the algorithm to merge fhArrays.
%
% PROPOSAL: For arrays which are not Nx1 numeric arrays (equal in format to fhArray), reduce them to an fhArray
%           first via a recursive call. Then merge the resulting fhArrays to one numeric NxM array which dimension
%           1-indices can be compared directly.
%   PRO: Does not need to iterate over arguments (derivatives of) for each comparison. ==> Potentially faster.
%
% PROPOSAL: Derive (simplistic) INTERNAL "hash" for data (used only by algorithm). Check for collisions.
%   NOTE: Does not to base hash on all data, just some that can be easily accessed.
%   PRO: Can avoid a lot of equality testing for data with few equalities.
%       PRO: Especially for high searchDistance.
%   PRO: For data with many equalities, can still sort by hash and only check equality between those.
%   CON: Hash needs to be calculated by reading every variable, can not just stop part-way as when calculating
%   equality (first variable inequality implies inequality for all variables as a whole).
%
% PROPOSAL: Return value for whether the algorithm managed to find all equalities (regardless of searchDistance).
%
% PROPOSAL: Have the algorithm predict the next equality given an already identified repeating sequence.
%   PRO: Could reduce the number of equality tests. ==> Speed up code
%
% PROPOSAL: Use unique() for numeric arrays.

nVariableArgs = numel(varargin);

% ASSERTIONS
assert(nVariableArgs >= 1)
assert(searchDistance >= 0)
for iArg = 1:nVariableArgs
  nRowsI = size(varargin{iArg}, 1);
  if iArg == 1
    nRows = nRowsI;
  end
  assert(nRows == nRowsI, ...
    'Data arguments do not have the same size in dimension 1.')
end

%===============================
% Iterate over all data indices
%===============================
%     guaranteesAllEqualities = true;
fauxHashArray = zeros(nRows,1);    % FH = Faux Hash
for i = 1:nRows

  % Assume unique hash until proven otherwise (tentative value).
  % NOTE: Using initial value "i" means that equal data indices refer to
  % the first occurrence. Not the only choice.
  fauxHashArray(i) = i;

  %=======================================================
  % Iterate over all preceding data indices within range
  %=======================================================
  jEarliest = max(1, i-searchDistance);
  for j = (i-1):-1:jEarliest

    %==========================================
    % Test equality between data indices i & j
    %==========================================
    ijEqual = true;   % Assume equality until proven otherwise.
    for iArg = 1:nVariableArgs
      % IMPLEMENTATION NOTE: Uses "isequaln", not isequal".
      if ~isequaln(varargin{iArg}(i, :), varargin{iArg}(j, :))
        ijEqual = false;
        break
      end
    end

    % Overwrite faux hash if equality.
    if ijEqual
      fauxHashArray(i) = fauxHashArray(j);
      break
    end

  end

  %         if ~ijEqual && (j == jEarliest)
  %             guaranteesAllEqualities = false;
  %         end
end

fauxHashArray = fauxHashArray(:);
end
