%
% Nearest-point interpolation from (x1,y1) to (x2,y2) which only works within a
% certain distance xMargin of min(x1) & max(x1). Outside of that interval,
% y2==NaN.
%
% Basically an extension of interp1 with EXTRAPVAL=NaN.
%
%
% ARGUMENTS
% =========
% xMargin : Scalar numeric. Non-negative.
% x1      : 1D array. Finite. May be empty, unsorted.
% y1      : 1D array. Same size as x1.
% x2      : Numeric. Any size.
%
%
% RETURN VALUE
% ============
% y2      : Same size as x2.
%           y2(i)=NaN if    x2(i) < min(x2) - xmargin,
%                        or x2(i) > min(x2) + xmargin
%           NOTE: Always double (if numel(x1) >= 2).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-28
%
function y2 = interpolate_nearest(xMargin, x1, y1, x2)
    % NOTE: In principle a weakness that y2 is always a double. Should ideally
    % be same MATLAB class as y1.
    
    assert(xMargin >= 0)
    
    
    nX1 = numel(x1);
    if     nX1 == 0
        y2 = NaN(size(x2));
        return
        
    elseif nX1 == 1
        % NOTE: Special treatment of scalar x1 since "interp1" does not support
        % it.
        assert(isfinite(x1))
        
        y2 = NaN(size(x2));
        y2(x2 == x1) = y1;
        
    else
        % Required by interp1.
        x1 = double(x1);
        y1 = double(y1);
        x2 = double(x2);
    
        % NOTE: For interp1
        %   x1 : Does NOT have to be sorted.
        %        Must
        %         - be a vector
        %         - have unique values
        %         - have at least 2 values
        %         - be finite
        %         - be single or double
        %   y1 : May be NaN, Inf.
        %        Must
        %         - be single or double
        %   x2 : May be NaN, Inf.
        %        Must
        %         - be single or double
        % Therefore no such assertions.
        y2 = interp1(x1, y1, x2, 'nearest', NaN);
    end
    
    [x1a, iA] = min(x1);
    [x1b, iB] = max(x1);
    
    bLowerMargin = (x1a - xMargin <= x2) & (x2 <= x1a);
    bUpperMargin = (x1b           <= x2) & (x2 <= x1b + xMargin);
    y2(bLowerMargin) = y1(iA);
    y2(bUpperMargin) = y1(iB);
end
