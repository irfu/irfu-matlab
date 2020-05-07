%
% Help the caller create a simple text table for printing.
%
%
% NOTE: The caller has to add e.g. trailing commas in the table cells if the caller want them to left/center adjusted
% too.
% NOTE: Adjustment does not take pre-existing leading or trailing whitespace into account.
% 
%
% ARGUMENTS
% =========
% headerStrs        : 1D cell array of strings, one per column.
% tableCells        : 2D cell array of strings. Size nRow x nColums. Contains the strings to be displayed in the table.
%                     Strings may have independently arbitrary lengths.
% columnAdjustments : 1D cell array. Each component is a string constant with one of values 'left', 'center', 'right'.
%
%
% RETURN VALUES
% =============
% headerStrs   : Argument headerStrs, but with padded whitespace so that all cells in the same column have the same width.
% tableCells   : Argument tableCells, ---------------------------------------------"--------------------------------------
% columnWidths : 1D array of column widths, i.e. length of the longest string in respective columns.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2020-01-14.
%
function [headerStrs, tableStrs, columnWidths] = assist_print_table(headerStrs, tableStrs, columnAdjustments)
    % ~PROBLEM: Will return last column that contains trailing whitespace that are unnecessary in log files etc.
    %   NOTE: The caller can remove, add, or reorder columns. What is the last column is not self-evident.
    
    HEADER_ADJUSTMENT = 'center';
    
    nRows = size(tableStrs, 1);
    nCols = size(tableStrs, 2);
    
    assert(size(tableStrs, 2) == length(headerStrs))
    assert(size(tableStrs, 2) == length(columnAdjustments), 'Differing number of rows in tableStrs and columnAdjustments.')
    
    for iCol = 1:nCols
        columnWidths(iCol) = max(cellfun(@length, [headerStrs(iCol); tableStrs(:, iCol)]));
        %columnWidths(iCol) = min(minColWidths(iCol), columnWidths(iCol));
    end    
    
    for iCol = 1:nCols
        % Pad strings to equal width within the column.
        formatStr = sprintf('%%%is', columnWidths(iCol));
        for iRow = 1:nRows
            tableStrs{iRow, iCol} = sprintf(formatStr, tableStrs{iRow, iCol});
        end
        
        headerStr        = sprintf(formatStr, headerStrs{iCol});
        headerStrs{iCol} = strjust(headerStr, HEADER_ADJUSTMENT);
        
        % Left/center/right-adjust the column.
        tableStrs(:, iCol) = strjust(tableStrs(:, iCol), columnAdjustments{iCol});
    end 

end
