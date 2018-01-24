% Print string to stdout but add a standard prefix string to every row to make it possible to
% separate (grep) those particular rows from other stdout/log messages.
% 
% NOTE: Can handle strings with many line feeds.
%
% ARGUMENTS
% =========
% str : String to be printed. NOTE: Assumes that string ends with a line feed.
%       (The intent of the algorithm would be ambiguous without this criterion.)
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31.
%
function stdout_disp(str)

global SETTINGS



% Character which when found in the argument string shall be interpreted as a line break.
% NOTE: Also used for line breaks in the string sent to "disp".
ARGUMENT_LINE_BREAK_CHAR = char(10);

STDOUT_PREFIX = SETTINGS.get_fv('STDOUT_PREFIX');


% ASSERTION: Require there to be a last character which represents a line break.
if length(str) < 1 || ~strcmp(str(end), ARGUMENT_LINE_BREAK_CHAR)
    error('BICAS:stdout_disp:Assertion:IllegalArgument', 'Not a legal string for printing. String must end with line feed.')
end



% NOTE: Removes the last character (line feed) since "disp" implicitly adds a line break at the end.
newStr = [...
    STDOUT_PREFIX, ...
    strrep(   str(1:end-1), ARGUMENT_LINE_BREAK_CHAR, [ARGUMENT_LINE_BREAK_CHAR, STDOUT_PREFIX]   )];

disp(newStr)

end
