% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Print to stdout but add a standard prefix string to every row to make it possible to
% separate (grep) those particular rows from other stdout/log messages.
% 1) if varargin empty: Use string "pattern" directly.
% 2) if varargin not empty: Use sprintf with "pattern" and varargin as arguments.
% (This makes it possible to print "raw" strings directly without having sprintf interpreting them.)
%
% NOTE: Requires string that is finally printed to end with a line break.
% (The algorithm would be ambiguous without this criterion.)
% NOTE: Can handle strings with many line breaks.
%
function stdout_printf(pattern, varargin)

global ERROR_CODES
global CONSTANTS



% Line break string.
% Put into the "constants" structure?!!
% NOTE: sprintf('\n') is actually platform dependent?!!
% NOTE: Used for interpreting line breaks in the function argument.
% NOTE: Used for line breaks in the output.
LINE_BREAK = sprintf('\n');   


if isempty(varargin)
    str = pattern;
else
    str = sprintf(pattern, varargin{:});
end



if (length(str) < 1) || ~strcmp(str(end), LINE_BREAK)
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Not a legal string for printing. String must end with line feed.')
end

% strsplit options:
%       'CollapseDelimiters' - If true (default), consecutive delimiters in S
%         are treated as one. If false, consecutive delimiters are treated as
%         separate delimiters, resulting in empty string '' elements between
%         matched delimiters.
rows = strsplit(str, LINE_BREAK, 'CollapseDelimiters', false);

rows(end) = [];   % Remove last string since it corresponds to empty string after last line feed.

new_str = '';
for row = rows
    new_str = [new_str, CONSTANTS.C.stdout_prefix, row{1}, LINE_BREAK];
end

fprintf(1, new_str);



end
