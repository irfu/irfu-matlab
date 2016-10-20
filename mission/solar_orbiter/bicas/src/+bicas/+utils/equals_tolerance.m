function eq = equals_tolerance(A,B, epsilon)
% Check if two numerical arrays are identical.
% Intended as a utility function for automatic test code.
%
% Equality requires:
% - Same MATLAB class.
% - Same array size, where
%   * Empty array == Empty array
% - Same values, where
%   * NaN == NaN
%   * Inf == Inf
%   * -Inf == -Inf
%   * Finite values are equal if within a max difference epsilon.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-14
%

% NOTE: Check for NaN, Inf, class, size/Ndimensions, empty matrix.

if ~strcmp(class(A), class(B))
    eq = 0;
    return
elseif ~isequal(size(A), size(B))
    eq = 0;
    return
elseif ~isequal(isnan(A), isnan(B))
    eq = 0;
    return
end
A(isnan(A)) = [];   % Seems like this produces a 1D array, but only if something is removed.
B(isnan(B)) = [];
if isempty(A)  % NOTE: Must check for empty matrix after removing NaN since that changes the size..
    eq = 1;
    return
end
max_diff = max(max(abs(double(A-B))));

eq = (max_diff <= epsilon);
end
