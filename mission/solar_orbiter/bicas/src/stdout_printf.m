% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Print as with fprintf(1, ...) but add a standard prefix string to every row to make it possible to
% separate (grep) those particular rows from other stdout/log messages.
%
% NOTE: Requires string to end with line feed (either the LF character, or the string "\n").
% (The algorithm would be ambiguous without this criterion.)
% NOTE: Can handle strings with many line feeds.
%
function stdout_printf(pattern, varargin)
%
% IMPLEMENTATION NOTE: Interprets string with sprintf before making any actual string manipulations
% since this takes care of escape codes etc. After that, substitutes the actual line feed character,
% rather than the (sub)string '\n'.
%
% PROPOSAL: Modify to accept (only) string with LF directly (not string '\n')?

global ERROR_CODES
STDOUT_PREFIX = 'STDOUT: ';



LF = sprintf('\n');   % Line feed _character_.
str = sprintf(pattern, varargin{:});
if (length(str) < 1) || ~strcmp(str(end), LF)
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Not a legal string for printing. String must end with line feed.')
end

rows = strsplit(str, LF, 'CollapseDelimiters', false);
rows(end) = [];   % Remove last string since it corresponds to empty string after last line feed.

new_str = '';
for row = rows
    new_str = [new_str, STDOUT_PREFIX, row{1}, LF];
end

fprintf(1, new_str);



end
