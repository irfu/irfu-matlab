%
% Escape string before submitting it as title, label etc in plots.
%
% Also works for cell arrays of strings.
function s = escape_str(s)
    s = strrep(s, '_', '\_');
end
