%
% Given a function consisting of piecewise constant intervals and described by the data points (x_i, y_i) at the
% beginning of each such segment, calculate
% (1) the values of the function at arbitrary positions, and
% (2) the integral from the first data point (x,y).
% Function extends to positive infinity but is undefined for x<x_1 (x_1=first data point).
%
%
% ARGUMENTS AND RETURN VALUES
% ===========================
% xp, yp : Same-sized arrays describing a function y=y(x) that is piecewise constant.
%           y=yp(i),   if xp(i) =< x < xp(i+1), i<numel(xp)
%           y=yp(end), if xp(end) < x.
%           xp must be monotonically increasing.
% X      : x values in which the function should be evaluated. Does not have to be sorted.
% Y      : y values corresponding to the x values in X.               NaN if X(i) < xp(1).
% YI     : YI(i) = The integral over the function from xp(1) to X(i). NaN if X(i) < xp(1).
%
%
% NOTE: The location of this function is preliminary since it is a generic function.
%
%
% Created 2018-01-08 by Erik Johansson, IRF Uppsala.
%
function [Y, YI] = step_func(xp, yp, X)
% PROPOSAL: Move to directory of generic util functions for package.
% NOTE: Can implement the evaluation of the function (not integral) using: Y = interp1(x, y, X, 'previous', 'extrap');
% PROPOSAL: Describe how to retrieve Y from YI, with ~diff somehow.

% ASSERTIONS
if ~all(size(xp) == size(yp))
    error('step_func:Assertion', 'xp and yp have dissimilar sizes.')
end
if ~(diff(xp) > 0)   % Check for monotonically increasing xp. Function "issorted" is not enough.
    error('step_func:Assertion', 'xp is not sorted.')
end

yip = 0;    % Integral value at x=xp(i), initial value for i=1. Initial value is integration constant.
Y  = zeros(size(X)) * NaN;
YI = zeros(size(X)) * NaN;
for i = 1:numel(xp)    % ASSUMES: xp monotonically increasing.
    if i >= 2
        yip = yip + (xp(i) - xp(i-1)) * yp(i-1);
    end
    k = (xp(i) <= X);    % ASSUMES: xp monotonically increasing. NOTE: Does not require xp(i+1) to exist.
    Y(k)  = yp(i);    % NOTE: Works for scalar k (1x1 array) since k is class "logical".
    YI(k) = yip + (X(k) - xp(i)) * yp(i);
end



% ALTERNATIVE OLDER IMPLEMENTATION. OBSOLETED. FASTER? DELETE?
%======================================================================================================
% if numel(xp) == 1
% 
%     Y  = yp * ones(size(X));   % NOTE: yp is scalar, but has to use X to produce size(Y)==size(X).
%     YI = (X - xp) * yp;
% 
%     Y( X < xp) = NaN;
%     YI(X < xp) = NaN;
% 
% elseif numel(xp) >= 2
% 
%     Y = sw_func(xp, yp, X);
% 
%     % yiAll(i) = Integral from beginning of function, to end of interval i
%     % (omitting last interval which goes to infinity).
%     yiInterval = diff(xp) .* yp(1:end-1);
%     yiAll = [0, cumsum(yiInterval)];
% 
%     YI = sw_func(xp, yiAll, X) + sw_func(xp, yp, X) .* (X - sw_func(xp, xp, X));
% 
% else
%     error('xp has illegal size (empty).')
% end
%======================================================================================================
end


% % Stepwise (SW) function. 
% function Y = sw_func(x, y, X)
% % Requires numel(x) >= 2, numel(y) >= 2 although one could expect it to work for numel == 1.
% Y = interp1(x, y, X, 'previous', 'extrap');   
% end
