%
% Add a prefix string to the beginning of every row in a multirow string.
%
% NOTE: Uses linefeed as linebreak.
%
%
% ARGUMENTS
% =========
% str    : Potentially multi-row string to be printed.
%          NOTE: If one-row   string: Argument may optionally end with line feed.
%                If multi-row string: Every row must end with line feed, including the last row.
%          NOTE: Not entirely unproblematic. Feeding unknown text string from an exception requires one to manually add
%                a line feed for this to work.
% prefix : One-row string without line feed.
%
%
% RETURN VALUE
% ============
% newStr : One or multi-row string
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-07-26
%
function newStr = add_prefix_on_every_row(str, prefix)
% PROPOSAL: Separate function for adding/ensuring trailing line feed for one-row strings.
%   PRO: bicas.Logger.log needs this separately for stderr printing.
%   PRO: "More clean".
%   CON: More hard-coded line feeds.

LINE_FEED = char(10);



% Convenience functionality: Add trailing line feed for strings without line feed (which have to be one-row strings).
% Multi-row strings must already have an ending linebreak.
if isempty(strfind(str, LINE_FEED))
  str = [str, LINE_FEED];
end



% ASSERTION: Require there to be a last character which represents a line break.
%
% IMPLEMENTATION NOTE: "Must" require input string to end with line feed to make the desired result more unambiguous.
% If the input string is not required to end with line feed, it is ambiguous whether an ending line feed (1) represents
% the end of the last row, or (2) precedes a last empty row consisting of zero characters (and itself not ending with a
% line feed).
if length(str) < 1 || ~strcmp(str(end), LINE_FEED)
  error('add_prefix_on_every_row:Assertion:IllegalArgument', 'Not a legal string for printing. String must end with line feed.')
end



% NOTE: Ignores the last character (line feed) when searching for line feed.
newStr = [...
  prefix, ...
  strrep(   str(1:end-1), LINE_FEED, [LINE_FEED, prefix]   ), ...
  LINE_FEED];

end
