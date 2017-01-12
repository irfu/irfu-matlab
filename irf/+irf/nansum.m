function y = nansum(x,dim)
%IRF.NANMEAN Sum ignoring NaNs.
%   Y = IRF.NANSUM(X) returns the sum of X, treating NaNs as missing
%   values.  For vector input, Y is the sum of the non-NaN elements
%   in X.  For matrix input, Y is a row vector containing the sum value of
%   non-NaN elements in each column.  For N-D arrays, IRF.NANSUM operates
%   along the first non-singleton dimension.
%
%   IRF.NANSUM(X,DIM) takes the sum along dimension DIM of X.
%
%   See also SUM, IRF.NANMEAN
%
%	NANSUM is also part of the Signal Processing Toolbox since Matlab R2013a.


% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    y = sum(x);
else
    % Count up non-NaNs.
    y = sum(x,dim);
end
