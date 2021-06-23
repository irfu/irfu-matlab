%
% Calculate "modified standard deviation", i.e. standard deviation but using an
% arbitrary reference (e.g. median) instead of the conventional mean.
%
%
% ARGUMENTS
% =========
% v    : Aray with arbitrary size and number of dimensions.
% ref  :
% iDim : Dimension along which to calculate MSTD.
%        
%
% RETURN VALUES
% =============
% mstd : Modified STandard Deviation.
%        Same size as "v", except for that size(mstd, iDim) == 1.
%        NaN, if size(v, iDim) == 0 or 1,
%             if "v" contains any NaN (for the particular.
%        Modified standard deviation otherwise.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-10-22.
%
function mstd = modif_std_deviation(v, ref, iDim)
    % PROPOSAL: Return NaN for zero samples (zero-length dimension).
    % TODO-DEC: Handle NaN?
    
    assert(isscalar(ref))
    assert(isscalar(iDim))
    
    N = size(v, iDim);
    if N == 0 || N == 1
        mstdSize       = size(v);
        mstdSize(iDim) = 1;
        if N == 0
            % CASE: N==0
            mstd = NaN(mstdSize);
        else
            % CASE: N == 1
            % Standard deviation is undefined / ill-defined.
            %   ref<>v ==> mstd ~ 1/0   (       1 / (N-1) )
            %   ref==v ==> mstd ~ 0/0   ( (v-ref) / (N-1) )
            mstd = NaN(mstdSize);
        end
    else
        % CASE: N>=2
        mstd = sqrt( sum((v-ref).^2, iDim) / (N-1) );
    end
    
    % ASSERTION
    % NOTE: Requires mstd to always be set (outside if-then).
    %EJ_library.assert.sizes(mstd, mstdSize)   % Disable?
end
