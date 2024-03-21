%
% Escape special characters in string before submitting it as title, label etc
% in plots in order to avoid special characters from being interpreted.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% s : String, or cell arrays of strings.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created <=2020-04-01.
%
function s = escape_str(s)
% PROPOSAL: Use regexprep.

s = strrep(s, '\', '\\');   % NOTE: Must come first!
s = strrep(s, '_', '\_');
s = strrep(s, '^', '\^');
s = strrep(s, '{', '\{');
s = strrep(s, '}', '\}');
end
