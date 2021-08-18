%
% Find out whether numeric matrices are approximately equal.
% Dissimilar dimensions count as being unequal.
%
% NOTE: There might be a function for this in newer versions of MATLAB. Not sure
% which.
% NOTE: Inf is always dissimilar from Inf, reagrdless of tolerance.
%
%
% ARGUMENTS
% =========
% a,b
%       Numeric arbitrary dimension matrices.
% 
%
% Initially created 2018-07-18 by Erik P G Johansson.
%
function result = approx_equals(a, b, epsilon, nanPolicy)
% PROPOSAL: Assertions for equal dimensions instead of result=false?!!
%   PRO: Assertion useful automatic tests.
%   CON: equals_recursive uses approx_equals. result=false is appropriate there.
%   CON: equals_recursive is probably the function that one should use for those
%        case. It is ~made for automatic tests.
%   --
%   PROPOSAL: Separate function for testing same-size, same class.
%   PROPOSAL: Setting for same-size assertion.
% 
% PROPOSAL: More sophisticated comparison?! log(a/b) < epsilon, for |a|,|b| > omega
%   CON: So far not needed.
%   CON: More arguments.
%
% PROPOSAL: New name implying numeric comparison.
    
    % ASSERTIONS
    assert(isnumeric(a))
    assert(isnumeric(b))
    assert(epsilon >= 0, 'Not non-negative epsilon argument.')
    
    % Check that dimensions.
    % NOTE: Should not be assertion.
    if ndims(a) ~= ndims(b)
        result = false;
        return
    end
    if ~all(size(a)  == size(b))
        result = false;
        return
    end
    
    % Convert to 1D ROW arrays. 1D arrays make it easier to construct code which
    % works with all dimensionalities.
    a = a(:);
    b = b(:);
    
    switch nanPolicy
        case 'NaN equal to itself'
            
            iA = isnan(a);
            iB = isnan(b);
            if all(iA == iB)
                % CASE: All NaN are in the exact same locations.
                a(iA) = [];
                b(iB) = [];
            else
                result = false;
                return
            end                
                        
        case 'NaN unequal to itself'
            % Do nothing
            
        otherwise
            error('Illegal nanPolicy')
    end
    
    % Comparing with Nan with any value (incuding NaN) ==> false.
    result = all(abs(a-b) <= epsilon);
end
