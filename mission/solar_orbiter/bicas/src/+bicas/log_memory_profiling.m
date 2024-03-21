%
% Log memory use in the caller workspace at the time of the call.
%
% IMPORTANT NOTE: SIMPLISTIC profiling, since it only uses "whos". It can
% therefore only access the memory use of variables in the caller workspace.
%
%
% ARGUMENTS
% =========
% locationName : Arbitrary string describing the location in the code/code
%                execution that is being profiled.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-27.
%
function log_memory_profiling(L, locationName)
% PROPOSAL: Automatic test code.
%
% PROPOSAL: Use generic code for creating tables.
%   NOTE: Already using irf.str.assist_print_table().
% PROPOSAL: Make more grep-friendly.
%   Ex: Sums, specific code locations.
%
% PROPOSAL: Column for largest variables.
%   PROPOSAL: Size ranking. 1=largest.

M = evalin('caller', 'whos');

firstRowStr = sprintf(...
  ['%s: Variable memory use in current workspace', ...
  ' (i.e. only a subset of memory use):\n'], ...
  locationName);



LOG_PREFIX         = 'MEMORY --';
HEADER_STRS        = {'Variable name', 'Memory use', 'Unit'};
COLUMN_ADJUSTMENTS = {'left', 'right', 'left'};
COLUMN_SEPARATOR   = ' ';
INDENT_SIZE        = 4;
DIVIDER_LINE_1_CHAR = '=';
DIVIDER_LINE_2_CHAR = '-';

dataStrs = cell(0,3);
for i = 1:numel(M)
  [valueStr, unit] = select_unit(M(i).bytes);
  dataStrs{i, 1} = sprintf('%s',   M(i).name);
  dataStrs{i, 2} = valueStr;
  dataStrs{i, 3} = sprintf('[%s]', unit);
end
[valueStr, unit] = select_unit(sum([M.bytes]));
dataStrs{i+1, 1} = 'SUM';
dataStrs{i+1, 2} = valueStr;
dataStrs{i+1, 3} = sprintf('[%s]', unit);
nDataRows = size(dataStrs, 1);



[headerStrs, dataStrs, ~] = irf.str.assist_print_table(...
  HEADER_STRS, dataStrs, COLUMN_ADJUSTMENTS);



headerDataCa = {strjoin(headerStrs, COLUMN_SEPARATOR)};
dividerStr1  = irf.str.repeat(DIVIDER_LINE_1_CHAR, numel(headerDataCa{1}));
dividerStr2  = irf.str.repeat(DIVIDER_LINE_2_CHAR, numel(headerDataCa{1}));
headerDataCa{end+1} = dividerStr1;
for iRow = 1:nDataRows
  if iRow == nDataRows
    headerDataCa{end+1} = dividerStr2;
  end
  headerDataCa{end+1} = strjoin(dataStrs(iRow, :), COLUMN_SEPARATOR);
end
tableStr = [strjoin(headerDataCa, newline), newline];
tableStr = irf.str.indent(tableStr, INDENT_SIZE);   % Indent table.

str = [firstRowStr, tableStr];
str = irf.str.add_prefix_on_every_row(str, LOG_PREFIX);
L.logf('debug', str)
end



% IMPLEMENTATION NOTE: Returns formatted value strings instead of numeric
% values, so that can add right-hand whitespace when not using decimals, thus
% placing the decimal point at the right location.
function [valueStr, unit] = select_unit(valueBytes)
if     valueBytes >= 2^20
  valueStr = sprintf('%.1f', valueBytes / 2^20);
  unit     = 'MiB';
elseif valueBytes >= 2^10
  valueStr = sprintf('%.1f', valueBytes / 2^10);
  unit     = 'kiB';
else
  % NOTE: Adds whitespace instead of decimals.
  valueStr = sprintf('%.0f  ', valueBytes);
  unit     = 'bytes';
end
end
