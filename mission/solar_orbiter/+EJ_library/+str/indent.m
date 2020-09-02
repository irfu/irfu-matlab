%
% Indent multirow string using EJ_library.str.add_prefix_on_every_row.
% See EJ_library.str.add_prefix_on_every_row for details.
%
% NOTE: Uses linefeed as linebreak.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-09.
%
function s = indent(s, nWhitespace)
    assert(nWhitespace >= 0, 'nWhitespace must be nonnegative.')
    
    % NOTE: "repmat" accepts negative sizes. May therefore want to assert nonnegative size.
    indentationStr = repmat(' ', 1, nWhitespace);
    s = EJ_library.str.add_prefix_on_every_row(s, indentationStr);
end
