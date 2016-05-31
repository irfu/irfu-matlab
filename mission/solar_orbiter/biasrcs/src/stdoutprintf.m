%
% Print as with fprintf(1, ...)
% but add a standard prefix string to every row to make it possible to separate (grep)
% those particular rows from other stdout/log messages.
%
% NOTE: Can handle strings with many line feeds.
% NOTE: Requires string to end with line feed (\n). The right algorithm would be ambiguous without
% this criterion.
%
function stdoutprintf(pattern, varargin)
%
% IMPLEMENTATION NOTE: Interprets string with sprintf before making any actual string manipulations
% since this takes care of escape codes etc. After that, substitutes the actual line feed character,
% rather than the (sub)string '\n'.
%
% PROPOSAL: Modify to accept (only) string with LF directly?
% PROPOSAL: Change name. soprintf? stdout_printf?

STDOUT_PREFIX = 'STDOUT: ';



LF = sprintf('\n');   % Line feed _character_.
str = sprintf(pattern, varargin{:});
if (length(str) < 1) || ~strcmp(str(end), LF)
    error('Not a legal string for printing. String must end with line feed.')
end

rows = strsplit(str, LF, 'CollapseDelimiters', false);
rows(end) = [];   % Remove last string since it corresponds to empty string after last line feed.

new_str = '';
for row = rows
    new_str = [new_str, STDOUT_PREFIX, row{1}, LF];
end

fprintf(1, new_str);



end
