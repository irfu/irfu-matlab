%
% Create a string which is another string repeated n times.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created <=2013-06-25.
%
function str = repeat(s, n)
assert(isnumeric(n) && n>=0)

str = repmat(s, [1, n]);
end
