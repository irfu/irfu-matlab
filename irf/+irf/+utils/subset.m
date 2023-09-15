%
% NOTE: Can also be used for checking supersets.
%
%
% ARGUMENTS
% =========
% set1, set2 : (1) Cell arrays of strings.
%              (2) Numeric arrays.
%              NOTE: May be empty arrays.
%
%
% RETURN VALUE
% ============
% isSubset : Whether set1 is a subset of set2 (including equal sets).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-25.
%
function isSubset = subset(set1, set2)
% NOTE: all({}) == all([]) == true
isSubset = all(ismember(set1, set2));
end
