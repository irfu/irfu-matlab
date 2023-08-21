%
% Remove trailing forward slash from path (string).
% Pure string function.
%
%
% ARGUMENTS
% =========
% path : Path.
%        NOTE: "/" is a permitted, despite that it is uncertain what the
%        desired result should be.
%
%
% RETURN VALUES
% =============
% path : Input argument without trailing slash.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-09-16.
%
function path = remove_trailing_slash(path)
% PROPOSAL: Automatic test code.
%
% PROBLEM: How handle "/"?
%   PROPOSAL: Permanently permit input path = "/" and remove trailing slash from it.
%       PRO: Well-defined string function without special cases.
%       PRO: Might be   desired if "/" is a  relative path.
%       CON: Might be undesired if "/" is an absolute path.
%           NOTE: Empty string can be seen as a representation of both
%           absolute path / and a relative path .. The caller know which in
%           both cases.
%   PROPOSAL: Assertion to exclude "/".
%
% PROPOSAL: Eliminate. Use strip(... , 'right', '/') instead.
%   PROPOSAL: Use "filesep".

%assert(isempty(regexp('', '^/+$')))

path = regexprep(path, '/*$', '');
end
