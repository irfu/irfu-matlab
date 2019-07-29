%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created <=2016-10-28
%
function react_to_false_assertion(giveError, msg)
% Function for either giving a warning, or an error depending on a setting (presumably a global setting).

if giveError
    error('BICAS:Assertion', msg)
else
    LINE_FEED = char(10);
    bicas.log('warning', ['FALSE ASSERTION: ', msg, LINE_FEED])
end

end
