%
% Calculate "modified standard deviation" (MSTD). Calculate it over one
% arbitrary dimension.
%
%
% IMPLEMENTATION NOTE
% ===================
% The implementation and interface only works with double, using NaN as a "fill
% value" which is allowed to propagate through the calculations and make the
% result NaN. This is use instead of using e.g. arbitrary fill value/positions
% (which are not counted/ignored) or bicas.utils.FPArray since:
% (1) bicas.utils.FPArray would probably be too slow in realistic applications,
% (2) the operation only makes sense for floats (square root, division), and
% (3) NaN propagates naturally and meaningfully through the mathematical
%     operations.
% Due to (3), it also makes sense to use all values in the calculation rather
% than ignoring (fill) values. It is therefore the caller's responsibility to
% omit values, e.g. NaN values/FPs if it wishes them to be omitted from the
% calculation.
%
%
% ARGUMENTS
% =========
% v
%       Double. Array with arbitrary size and number of dimensions. May be NaN.
% ref
%       Double. Scalar value that replaces the mean in standard deviation. May
%       be NaN.
% iDim
%       Scalar value. Dimension along which to calculate MSTD.
%
%
% RETURN VALUES
% =============
% mstd
%       MSTD. Double.
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
function mstd = mstd(v, ref, iDim)

    % Permit ref=NaN.
    assert(isscalar(ref))
    assert(isscalar(iDim))
    assert(isa(v,   'double'))
    assert(isa(ref, 'double'))

    N = size(v, iDim);
    if N == 0 || N == 1
        % NOTE:
        % N==0 ==> MSTD undefined.
        % N==1 ==> MSTD is undefined / ill-defined.
        %   ref<>v ==> mstd ~ 1/0        (       1 / (N-1) )
        %   ref==v ==> mstd ~ 1/0 or 0/0 ( (v-ref) / (N-1) )

        mstdSize       = size(v);
        mstdSize(iDim) = 1;
        mstd           = NaN(mstdSize);
    else
        % CASE: N>=2
        mstd = sqrt( sum((v-ref).^2, iDim) / (N-1) );
    end
end
