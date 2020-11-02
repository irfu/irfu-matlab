%
% Calculate "modified standard deviation", but using an arbitrary
% reference rather than the average.
%
%
% ARGUMENTS
% =========
% v    : Arbitrary number of dimensions.
% iDim : Dimension along which to calculate.
%        
%
%
%
% RETURN VALUES
% =============
% s : NaN, for zero values
%     0,   for one value
%     Modified standard deviation otherwise.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-10-22.
%
function s = modif_std_deviation(v, ref, iDim)
    % PROPOSAL: Automatic test code.
    % PROPOSAL: Return NaN for zero samples (zero-length dimension).
    
    assert(isscalar(ref))
    assert(isscalar(iDim))
    
    N = size(v, iDim);
    if N == 0 || N == 1
        vSize = size(v);
        vSize(iDim) = 1;
        if N == 0
            s = NaN(vSize);
        else
            % CASE: N == 1
            s = zeros(vSize);
        end
    else
        s = sqrt( sum((v-ref).^2, iDim) / (N-1) );
    end
    
end
