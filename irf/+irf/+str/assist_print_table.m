%
% Help the caller create a simple text table for printing.
%
%
% NOTE: The caller has to add e.g. trailing commas in the table cells if the
%       caller wants them (the commas) to be left/center adjusted too.
% NOTE: Deliberately does NOT remove leading and trailing whitespace before
%       padding so that caller can use this for fine-tuning adjustment.
%   Ex: Right-adjust & pad "123  " and "3.5" ==> "123  " and "  3.5" so that
%       decimal point remains in the same place.
%
%
% ARGUMENTS
% =========
% headerStrs        : 1D cell array of strings. {iCol}.
% dataStrs          : 2D cell array of strings. {iRow, iCol}.
%                     Contains the strings to be displayed in the table.
%                     Strings may have independently arbitrary lengths.
% columnAdjustments : 1D cell array of string constants.
%                     'left', 'center', 'right'.
%
%
% RETURN VALUES
% =============
% headerStrs,
% dataStrs     : Same as corresponding arguments, but with padded whitespace so
%                that all cells in the same column (data+headers combined) have
%                the same column width.
% columnWidths : Numeric 1D array of column widths, i.e. length of the longest
%                string in respective columns.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-14.
%
function [headerStrs, dataStrs, columnWidths] = assist_print_table(headerStrs, dataStrs, columnAdjustments)
% ~PROBLEM: Will return last column that contains trailing whitespace that are unnecessary in log files etc.
%   NOTE: The caller can remove, add, or reorder columns. What is the last column is not self-evident.
%
% PROPOSAL: New name for headerStrs
%   headersCa, headerStrsCa
% PROPOSAL: New name for tableStrs.
%   ~str, ~table, ~values, ~Ca, ~data
%   strTable, valueTable, tableValuesCa, tableStrsCa
%   dataCa, dataStrsCa
%
% PROBLEM: Can not handle the case of "merged columns": One header for two columns
%   Ex: Numeric value + unit in two separate columns. Numeric value is right-adjusted, unit is left-adjusted, one
%       header for both.
%       Ex: "123 bytes"
%   PROPOSAL: headerStrs{iCol} == [] ==> Merged column headers.

HEADER_ADJUSTMENT = 'center';

% NORMALIZE
headerStrs        = headerStrs(:)';           % Row
columnAdjustments = columnAdjustments(:)';    % Row

% ASSERTIONS
[~, nCols] = irf.assert.sizes(...
  headerStrs,         [ 1, -2], ...
  dataStrs,           [-1, -2], ...
  columnAdjustments,  [ 1, -2]);

combinedStrs = [headerStrs; dataStrs];
for iCol = 1:nCols
  columnWidths(iCol)  = max(cellfun(@length, combinedStrs(:, iCol)));
  headerStrs(:, iCol) = pad_adjust(headerStrs(1, iCol), HEADER_ADJUSTMENT,       columnWidths(iCol));
  dataStrs(  :, iCol) = pad_adjust(dataStrs(  :, iCol), columnAdjustments{iCol}, columnWidths(iCol));
end
end



function strCa = pad_adjust(strCa, adjustment, newWidth)
% PROPOSAL: Make into generic function.
% NOTE: Deliberately does NOT remove leading and trailing whitespace.

switch(adjustment)
  case 'left'
    padArg = 'right';
  case 'right'
    padArg = 'left';
  case 'center'
    padArg = 'both';
  otherwise
    error('Illegal argument adjustment="%s".', adjustment)
end

strCa = pad(strCa, newWidth, padArg);
end
