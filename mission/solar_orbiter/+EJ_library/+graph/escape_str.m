%
% Escape special characters in string before submitting it as title, label etc in plots in order to avoid special
% characters from being interpreted.
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% s : String, or cell arrays of strings.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created <=2020-04-01.
%
function s = escape_str(s)
    s = strrep(s, '_', '\_');
end
