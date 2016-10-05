% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Print string to stdout but add a standard prefix string to every row to make it possible to
% separate (grep) those particular rows from other stdout/log messages.
% 
% NOTE: Can handle strings with many line feeds.
%
% str : String to be printed. NOTE: Assumes that string ends with a line feed.
%       (The intent of the algorithm would be ambiguous without this criterion.)
%
function stdout_disp(str)

global CONSTANTS



% Line break string.
% Put into the "constants" structure?!!
% NOTE: Used for interpreting line breaks in the function argument.
% NOTE: Used for line breaks in the output.
LINE_FEED = char(10);



% ASSERTIONS
if length(str) < 1 || ~strcmp(str(end), LINE_FEED)
    error('BICAS:stdout_disp:Assertion:IllegalArgument', 'Not a legal string for printing. String must end with line feed.')
end



% NOTE: Removes the last character (line feed) since "disp" implicitly adds a line break at the end.
new_str = [CONSTANTS.C.stdout_prefix, strrep(str(1:end-1), LINE_FEED, [LINE_FEED, CONSTANTS.C.stdout_prefix])];

disp(new_str)



end
