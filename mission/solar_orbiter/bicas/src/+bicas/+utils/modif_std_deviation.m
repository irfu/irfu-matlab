%
% Calculate "modified standard deviation" (MSTD). Calculate it over one
% arbitrary dimension.
%
%
% ARGUMENTS
% =========
% v
%       Array with arbitrary size and number of dimensions. May be NaN.
% ref
%       Scalar value that replaces the mean in standard deviation. May be NaN.
% iDim
%       Scalar value. Dimension along which to calculate MSTD.
%        
%
% RETURN VALUES
% =============
% mstd
%       MSTD.
%       Same size as "v", except for that size(mstd, iDim) == 1.
%       --
%       NaN, if size(v, iDim) == 0 or 1, or
%            if "v" contains at least one NaN for the particular 1D-sub-vector
%            for which one scalar MSTD values is calculated.
%       MSTD otherwise.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-10-22.
%
function mstd = modif_std_deviation(v, ref, iDim)
    % PROPOSAL: Rename using abbreviation MSTD.
    % PROPOSAL: Ignore NaN in the calculation of MSTD.
    %           ==> Be able to return non-NaN if there are sufficiently many non-NaN values.
    % PROPOSAL: Argument for smallest number of non-NaN values separately per
    %           1D-subvector.
    % PROPOSAL: Non-scalar "ref".
    %   CON: Complicates assertion on size.
    
    assert(isscalar(ref))    % Permit NaN.
    assert(isscalar(iDim))
    
    N = size(v, iDim);
    if N == 0 || N == 1
        % NOTE:
        % N==0 ==> MSTD undefined.
        % N==1 ==> Standard deviation is undefined / ill-defined.
        %   ref<>v ==> mstd ~ 1/0   (       1 / (N-1) )
        %   ref==v ==> mstd ~ 0/0   ( (v-ref) / (N-1) )
        
        mstdSize       = size(v);
        mstdSize(iDim) = 1;
        mstd           = NaN(mstdSize);
    else
        % CASE: N>=2
        mstd = sqrt( sum((v-ref).^2, iDim) / (N-1) );
        
        % EXPERIMENTAL
        % Try ignore NaN.
        % ==> Handle absense of non-NaN.
        %mstd = sqrt( sum((v-ref).^2, iDim, 'ignorenan') / (N-1) );
        %nNaN = sum(isnan(v), iDim);
        %mstd(nNaN == N) = NaN;
        
    end
    
    % ASSERTION
    % NOTE: Requires mstd to always be set (outside if-then).
    %irf.assert.sizes(mstd, mstdSize)   % Disable?
end
