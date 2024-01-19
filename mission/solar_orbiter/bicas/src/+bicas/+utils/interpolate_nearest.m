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
% xMargin
%       Scalar numeric. Non-negative.
% xArray1
%       1D array. Finite. May be empty, unsorted.
% yArray1
%       1D array. Same size as x1.
% xArray2
%       Numeric. Any size.
%
%
% RETURN VALUE
% ============
% yArray2
%       Same size as xArray2.
%       y2(i)=NaN if    xArray2(i) < min(xArray2) - xMargin,
%                    or xArray2(i) > min(xArray2) + xMargin
%       NOTE: Always double (if numel(xArray1) >= 2).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-04-28
%
function yArray2 = interpolate_nearest(xMargin, xArray1, yArray1, xArray2)
    % NOTE: In principle a weakness that y2 is always a double. Should ideally
    %       be same MATLAB class as y1.
    %
    % PROPOSAL: Make y values work for any MATLAB class: input and output (same type).
    %   PRO: Useful for integer+logical ZVs while caring about fill values.
    %   CON: interp1() requires floats. Would have to convert integers and
    %        logicals to floats and back.
    %       PRO: Can not represent fill values.

    assert(isscalar(xMargin) && (xMargin >= 0))



    nX1 = numel(xArray1);
    if     nX1 == 0
        yArray2 = NaN(size(xArray2));
        return

    elseif nX1 == 1
        % NOTE: Special treatment of scalar x1 since interp1() does not support
        % it.
        assert(isfinite(xArray1))

        yArray2 = NaN(size(xArray2));
        yArray2(xArray2 == xArray1) = yArray1;

    else
        % IMPLEMENTATION NOTE: interp1() requires floats.
        % Required by interp1.
        xArray1 = double(xArray1);
        yArray1 = double(yArray1);
        xArray2 = double(xArray2);

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
        yArray2 = interp1(xArray1, yArray1, xArray2, 'nearest', NaN);
    end

    [x1a, iA] = min(xArray1);
    [x1b, iB] = max(xArray1);

    bLowerMargin = (x1a - xMargin <= xArray2) & (xArray2 <= x1a);
    bUpperMargin = (x1b           <= xArray2) & (xArray2 <= x1b + xMargin);
    yArray2(bLowerMargin) = yArray1(iA);
    yArray2(bUpperMargin) = yArray1(iB);
end
