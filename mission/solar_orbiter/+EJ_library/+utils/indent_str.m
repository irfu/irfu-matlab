%
% Indent multirow string using EJ_library.utils.add_prefix_on_every_row.
% See EJ_library.utils.add_prefix_on_every_row for details.
%
% NOTE: Uses linefeed as linebreak.
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-09.
%
function s = indent_str(s, nWhitespace)
    assert(nWhitespace >= 0, 'nWhitespace mus be nonegative.')
    
    % NOTE: "repmat" accepts negative sizes.
    indentationStr = repmat(' ', 1, nWhitespace);
    s = EJ_library.utils.add_prefix_on_every_row(s, indentationStr);
end
