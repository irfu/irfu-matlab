%
% Indent multirow string using irf.str.add_prefix_on_every_row.
% See irf.str.add_prefix_on_every_row for details.
%
% NOTE: Uses linefeed as linebreak.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-09.
%
function s = indent(s, nWhitespace)
% NOTE: "repmat" accepts negative sizes. Therefore wants to assert
% nonnegative size.
assert(nWhitespace >= 0, 'nWhitespace must be nonnegative.')

indentationStr = repmat(' ', 1, nWhitespace);
s = irf.str.add_prefix_on_every_row(s, indentationStr);
end
