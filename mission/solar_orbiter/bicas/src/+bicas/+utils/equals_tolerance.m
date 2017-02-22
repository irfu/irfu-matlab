function eq = equals_tolerance(A, B, epsilon)
% function eq = equals_tolerance(A, B, epsilon)   Check if two numerical arrays are identical.
% 
% Primarily intended as a utility function for automatic test code.
%
% Equality requires:
% - Same MATLAB class.
% - Same array size, where
%   * Empty array == Empty array, if arrays also have the same array sizes.
%     NOTE: MATLAB permits multiple different empty (zero-component) array sizes, e.g. 0x2 and 3x0 which in this
%           function count as different arrays.
% - Same values, where
%   * NaN == NaN
%   * Inf == Inf
%   * -Inf == -Inf
%   * Finite values are equal if within a max difference epsilon.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-10-14
%
%
% ARGUMENTS AND RETURN VALUE
% ==========================
% A, B    : The two numerical arrays to be compared.
% epsilon : The maximum difference (absolute value) permitted between corresponding elements in A and B.
% eq      : 1 iff A and B are deemed equal.

% PROPOSAL: Return difference value.
% PROPOSAL: Create test code.
% NOTE: Checks for NaN, Inf, class, size/Ndimensions, empty matrix.

% ASSERTIONS
if ~isnumeric(A) || ~isnumeric(B)
    error('BICAS:equals_tolerance:Assertion:IllegalArgument', 'Illegal argument class')
end


% Check equality of sizes and occurrances of NaN.
if ~strcmp(class(A), class(B))
    eq = 0;
    return
elseif ~isequal(size(A), size(B))
    eq = 0;
    return
elseif ~isequal(isnan(A), isnan(B))
    % A and B do NOT have NaN in the same places.
    eq = 0;
    return
end

% Convert to 1D ROW arrays. 1D arrays make it easier to construct code which works with all dimensionalities.
A = A(:)';
B = B(:)';

% ASSUMES: 1D row arrays.
% NOTE: For a 2D array, this syntax produces a 1D ROW array, but only if something is removed. Therefore best to work only 1D ROW arrays.
A(isnan(A)) = [];
B(isnan(B)) = [];

% Check for empty matrix (removing NaN changes the size).
if isempty(A)    
    eq = 1;
    return
end

% ASSUMES: Non-empty 1D arrays ("max" requires 1D vectors).
% NOTE: This should alos work for -inf/+inf values.
maxDiff = max(abs(double(A-B)));    % Greatest difference between any pair of components.

eq = (maxDiff <= epsilon);
end
