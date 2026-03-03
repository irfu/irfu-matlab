%
% Given a number of arrays with the same number of rows, find groups of rows
% (into all arrays) where the rows (for every array separately) are identical.
%
% Function is primarily intended for easily finding groups of CDF records where
% groups of settings (zVariables) are identical.
%
%
% PERFORMANCE NOTES
% =================
% The speed of this function should theoretically scale as t~O(n^2), if
% all rows are unequal. For real BICAS RPW data, there should not be so many
% unique rows that this creates an insurmountable problem.
%
% POTENTIAL PROBLEM: If a large initial block of rows are ~identical, and is
% followed by another block of rows unequal to all rows in the first block (e.g.
% at an RPW mode change), then the algorithm could potentially become slow since
% current rows in the second block are all compared with rows in the first block
% before it reaches the identical rows in the second block and the search can be
% stopped. This effect has not been tested though.
%
%
% ARGUMENTS
% =========
% varargin
%       All arguments must have the same size in dimension 1 (rows). May be zero
%       arguments
%
%
% RETURN VALUES
% =============
% iGroupArCa
%       Column cell array of row indices.
%       {iGroup} = Sorted column array of indices to rows in arguments.
%       NOTE: The sequence of values {iGroup}(1) is sorted.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
function iGroupArCa = group_unique_rows(varargin)
% PROPOSAL: Better name
%   ~find, group, rows
%   find_unique_rows
%   find_unique_row_groups
%   group_unique_rows
%
% NOTE: iFirstIdenticalRowAr should be the same (almost) as fauxHashArray
%       irf.utils.find_equalities() with searchDistance=Inf.
% PROPOSAL: Return iFirstIdenticalRowAr as pseudo-hashes (under some other
%           name) and use this function for replacing
%           irf.utils.find_equalities() with searchDistance=Inf.
%   NOTE: There are performance issues with the use of
%         irf.utils.find_equalities() in irf.utils.find_equalities() in
%         solo.adm.group_sort_DSMD_versions(). Should probably test for that
%         for both functions first.
%
% PROPOSAL: Search for identical rows, beginning at the current row and upwards
%           (decreasing row number).
%   PRO: Potentially faster for real data with mode changes.
%   CON: Effectively assumes that similar rows are close to each other.
% PROPOSAL: Divide rows into blocks. All blocks are separately searched for
%           groups first. Results are then combined.
%   CON: Effectively assumes that similar rows are close to each other.

nArrays = numel(varargin);

nRowsArray = zeros(nArrays, 1);
for iArray = 1:nArrays
  nRowsArray(iArray) = size(varargin{iArray}, 1);
end
nRows = unique(nRowsArray);
assert(length(nRows) <= 1, ...
  'Data arguments do not have the same size in dimension 1.')

% (iRow) = First identical row for iRow. Can be itself.
iFirstIdenticalRowAr = nan(nRows, 1);    % Preallocate.
% Preallocate. Needed by algorithm, because it is not filled in one index at a
% time.
iGroupArCa = cell(nRows, 1);

for iRowCurrent = 1:nRows
  % Given a row, compare it to all preceding rows. If the current row does not
  % match any preceding rows, then add a new group of rows. If the current row
  % does match any preceding row, then add it to the group of that preceding row
  % and stop searching for more matches.

  % Assume current row is unique, until proven otherwise.
  % IMPLEMENTATION NOTE: This design simplifies the implementation.
  iFirstIdenticalRowAr(iRowCurrent) =  iRowCurrent;
  iGroupArCa(          iRowCurrent) = {iRowCurrent};

  for iRowOther = 1:iRowCurrent-1
    if rows_equal(varargin, iRowOther, iRowCurrent)

      iFirstIdenticalRowAr(iRowCurrent) = iRowOther;
      iGroupArCa{iRowOther, 1} = [iGroupArCa{iRowOther}; iRowCurrent];
      break

    end
  end

end

% Only keep elements set to valid values.
iGroupArCa = iGroupArCa(unique(iFirstIdenticalRowAr), 1);
end



function b = rows_equal(arraysCa, iRow1, iRow2)
for iArray = 1:numel(arraysCa)
  array = arraysCa{iArray};

  % NOTE: Using isequaln() (not isequal()).
  % NOTE: Should be able to handle all array dimensionalities: All dimension
  %       starting from dimension 2 are linearized, which is fine for
  %       comparisons.
  if ~isequaln(array(iRow1, :), array(iRow2, :))
    b = false;
    return
  end

  b = true;
end
end
