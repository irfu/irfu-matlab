%
% regexpf = regexp (MATLAB builtin) + f (=full match)
%
% Roughly like MATLAB's "regexp", except
% (1) it ONLY MATCHES ENTIRE strings (not substrings). In practice, it surrounds
%     the submitted regular expressions with ^ and $.
% (2) it can match empty strings (option in regexp),
% (3) it returns a more usable logic array.
% (4) it can iterate over both (a) strings, and (b) regular expressions, albeit
%     not simultaneously.
%
%
% ARGUMENTS
% =========
% str
%       String or cell array of strings. Empty string must be of class char,
%       i.e. e.g. not [] (due to function "regexp").
%       NOTE: Will permit empty strings to match a regular expression.
% regexPattern
%       String or cell array of strings. Each string is a regexp which may or
%       may not be surrounded by ^ and $.
% --
% NOTE: str and regexPattern must not simultaneously be non-scalar cell arrays.
%
%
% RETURN VALUE
% ============
% isMatch
%       Logical array. True iff str is matched by corresponding regexp. Same
%       size and indices as any non-scalar argument cell array.
%       NOTE: This is different from "regexp" which returns a cell array of
%       arrays.
%
%
% Initially created 2018-07-12 by Erik P G Johansson, IRF, Uppsala, Sweden.
%
function isMatch = regexpf(str, regexpPattern)
%
% PROPOSAL: Allow both str and regexpPattern to simultaneously be cell arrays of strings.
%   CON: Using both simultaneously means having to return a 2D matrix and the
%        caller having to keep track of indices with different meanings.
%       PRO: Bad for backwards compatibility.
%       CON: Breaks existing functionality.
%           CON-PROPOSAL: Allow both str and regexpPattern to be arbitrary
%                         dimension cell arrays, as long as dimensions do not
%                         conflict (must not have any overlapping non-size-one
%                         dimensions/indices).
%               CON/PROBLEM: Difficult to implement.




%====================================================
% ASSERTIONS:
% NORMALIZE input: Turn strings into 1x1 cell arrays
%====================================================
% IMPLEMENTATION NOTE: Empty string must be of class char, i.e. not [] (due
% to function "regexp"). This is also more consistent.
if ischar(str)
  str = {str};
end
if ischar(regexpPattern)
  regexpPattern = {regexpPattern};
end
assert(iscell(str))
assert(iscell(regexpPattern))



if  numel(regexpPattern) == 1
  % NOTE: Handles case of both str and regexpPattern being scalar.

  isMatch = cellfun(...
    @(s) (~isempty(regexp(...
    s, ['^', regexpPattern{1}, '$'], 'once', 'emptymatch'))), ...
    str);

elseif numel(str) == 1

  % IMPLEMENTATION NOTE: regexp option "emptymatch" is required to match
  % empty strings.
  isMatch = cellfun(...
    @(re) (~isempty(regexp(...
    str{1}, ['^', re, '$'], 'once', 'emptymatch'))), ...
    regexpPattern);
else
  % CASE: Neither "str" or "regexpPattern" is scalar.
  % ASSERTION
  error(...
    ['Both arguments are non-scalar (more than one element)', ...
    ' cell arrays. Function can not handle that case.'])
end

% IMPLEMENTATION NOTE: Empirically, cellfun() returns empty DOUBLE array for
% empty cell arrays, but the size is correct. ==> Must correct type.
if isempty(isMatch)
  isMatch = logical(isMatch);
end

assert(islogical(isMatch))
end
