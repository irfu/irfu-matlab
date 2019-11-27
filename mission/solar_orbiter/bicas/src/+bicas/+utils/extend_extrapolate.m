%
% Extend a tabulated function y=y(x) by extrapolation.
% NOTE: This is different from just extrapolating a functions since it also extrapolates the x values.
%
% NOTE: It is unclear what is the best ("natural") way of specifying where the tabulated function should be extended (at
% positive/negative x; how much). Might change.
%
%
% ARGUMENTS
% =========
% xMinDelta : Distance on the x axis by which the table should be extended at a minimum.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-11-27
%
function [x1,y1] = extend_extrapolate(x0, y0, xMinDelta, xExtrapolationDir, xExtrapolationMethod, yExtrapolationMethod)
    % PROPOSAL: Instead specify where extra margin should be via new (after extrapolation) x limit, which may be below/above/within the range of x0.
    %   PRO: xExtrapolationDir is unnecessary.
    %   CON: Must search for min-max x.
    % PROPOSAL: Use the sign in xMinDelta to specify extrapolation direction.
    % PROPOSAL: Recursive call to handle negative extrapolation direction.
    %   CON: Asserting x order hinders.
    
    assert(issorted(x0))
    EJ_library.utils.assert.vector(x0);
    EJ_library.utils.assert.vector(y0);
    %assert(xMinDelta > 0)
    
    % NOTE: Column vectors are needed for (1) producing xn and yn column vectors and being able to merge them with x0,
    % y1.
    x0 = x0(:);
    y0 = y0(:);

    switch(xExtrapolationDir)
        case 'positive'
            [x1,y1] = extend_extrapolate_end_index(x0, y0, xMinDelta, xExtrapolationMethod, yExtrapolationMethod);

        case 'negative'
            % NOTE: Assumes that all x & y extrapolation methods are invariant under sign change for x.
            x0 = x0(end:-1:1);
            y0 = y0(end:-1:1);
            [x1,y1] = extend_extrapolate_end_index(x0, y0, xMinDelta, xExtrapolationMethod, yExtrapolationMethod);
            x1 = x1(end:-1:1);
            y1 = y1(end:-1:1);

        otherwise
            error('BICAS:extend_extrapolate:Assertion:IllegalArgument', 'Illegal argument xExtrapolationDir="%s"', xExtrapolationDir)
    end

end



function [x1,y1] = extend_extrapolate_end_index(x0, y0, xMinDelta, xExtrapolationMethod, yExtrapolationMethod)
    
    switch(xExtrapolationMethod)
        case 'linear'
            dx = x0(end) - x0(end-1);
            %xn = (x0(end)+dx) : dx : (x0(end)+xMinDelta+dx);
            
            xp = x0(end);
            xn = [];
            while true
                xp = xp + dx;
                xn(end+1, 1) = xp;
                
                if abs(xp-x0(end)) >= xMinDelta
                    break
                end
            end

        case 'exponential'
            kx = x0(end) / x0(end-1);
            assert(kx>=0)

            xp = x0(end);
            xn = [];
            assert( (kx>1) || ((kx<1) && (abs(xp) > xMinDelta)), 'Can not fill requested x interval with new x values using exponential extrapolation of x values.')
            while true
                xp = xp * kx;
                xn(end+1, 1) = xp;
                
                if abs(xp-x0(end)) >= xMinDelta
                    break
                end
            end

        otherwise
            error('BICAS:extend_extrapolate:Assertion:IllegalArgument', 'Illegal argument xExtrapolationMethod="%s"', xExtrapolationMethod)

    end

    switch(yExtrapolationMethod)
        case 'linear'
            dx = x0(end) - x0(end-1);
            dy = y0(end) - y0(end-1);
            yn = y0(end) + dy/dx * (xn-x0(end));

        case 'exponential'
            assert(y0(end) / y0(end-1) > 0)

            k  = log(y0(end) / y0(end-1)) / (x0(end) - x0(end-1));
            yn = y0(end) * exp(k*(xn-x0(end)));

        otherwise
            error('BICAS:extend_extrapolate:Assertion:IllegalArgument', 'Illegal argument yExtrapolationMethod="%s"', yExtrapolationMethod)

    end
    
    x1 = [x0; xn];
    y1 = [y0; yn];    
end
