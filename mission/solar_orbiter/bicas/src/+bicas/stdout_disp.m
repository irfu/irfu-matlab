% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-05-31
%
% Print to stdout but add a standard prefix string to every row to make it possible to
% separate (grep) those particular rows from other stdout/log messages.
% 
% str : String to be printed. NOTE: Assumes that sting ends with a line break.
%       (The intent of the algorithm would be ambiguous without this criterion.)
%
% NOTE: Can handle strings with many line breaks.
%
function stdout_disp(str)

global ERROR_CODES
global CONSTANTS



% Line break string.
% Put into the "constants" structure?!!
% NOTE: Seen some hint somewhere that sprintf('\n') is actually platform dependent?!! Verify.
% NOTE: Used for interpreting line breaks in the function argument.
% NOTE: Used for line breaks in the output.
LINE_BREAK = sprintf('\n');



if length(str) < 1 || ~strcmp(str(end), LINE_BREAK)
    errorp(ERROR_CODES.ASSERTION_ERROR, 'Not a legal string for printing. String must end with line break.')
end



% NOTE: Removes the last character (line break) since "disp" implicitly adds a line break at the end.
new_str = [CONSTANTS.C.stdout_prefix, strrep(str(1:end-1), LINE_BREAK, [LINE_BREAK, CONSTANTS.C.stdout_prefix])];

disp(new_str)



end
