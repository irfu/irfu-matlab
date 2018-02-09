classdef engine_utils
    
    methods(Static, Access=public)
        
        function [xp, y1, yArray] = linear_extrapolate_limit(x0, x1, y0, dydx, xArray, yArray, yMin, yMax)
            % Assume function y=y(x). Linearly extrapolate y values from x0,y0 to x1,y1 but keep track of when exceeding bounds.
            % Intended to be used for evolving variables which can not exceed certain thresholds over time.
            %
            %
            %
            % ARGUMENTS
            % =========
            % x0, y0          : Starting point for extrapolation.
            % dydx            : Constant derivative used for the extrapolation.
            % x1              : End point for extrapolation.
            % yMin, yMax      : Interval of y values which the extrapolation from x0,y0 to x1,y1 must stay within.
            % xArray, yArray  : Arrays with x and (potentially) y values.
            % --
            % NOTE: x0 < x1
            % NOTE: x0==yMax => dydx =< 0
            %       x0==yMin => dydx => 0
            %
            %
            % RETURN VALUES
            % =============
            % xp     : Point where extrapolation reached min or max limit.
            % y1     : y value for x=x1. Is strictly within the closed interval yMin--yMax, if xp == x1.
            % yArray : Argument yArray, but with the y values for the range x0-x1 (in xArray) set to yArray(i)==y(xArray(i)).
            %
            %
            % IMPLEMENTATION NOTE: Many subtle and hard to resolve bugs in "engine" have been associated with this
            % function. It therefore has many assertions and its own test code.


            % ASSERTIONS
            if x0 >= x1    % NOTE: Modifying yArray requires this check. Preventing zero length extrapolations requires check.
                error('x0 >= x1')
            end
            %if x0 > x1    % NOTE: Modifying yArray requires this check.
            %    error('x0 > x1')
            %end
            if (yMin > yMax)
                error('yMin > yMax')
            end
            [~, belowMin, aboveMax] = TM_power_budget.engine_utils.outside_boundaries(y0, yMin, yMax);
            if (belowMin || aboveMax)
                error('y0=%g out of bounds yMin=%g, yMax=%g.', y0, yMin, yMax)
            end
            if (y0 == yMin) && (dydx < 0)
                error('x0=yMin=%g + illegal dydx=%e ==> Can not extrapolate positive x distance.', x0, dydx)
            elseif (y0 == yMax) && (dydx > 0)
                error('x0=yMax=%g + illegal dydx=%e ==> Can not extrapolate positive x distance.', x0, dydx)
            end
            
            %============
            % Set xp, y1
            %============
            % Check for y0 < y <= y1 being out of bounds via y1 (knowing that we are working with a linear function y=y(x)).
            y1 = y0 + (x1 - x0) * dydx;
            [yp, belowMin, aboveMax] = TM_power_budget.engine_utils.outside_boundaries(y1, yMin, yMax);
            if (belowMin || aboveMax)
                % CASE: y1 out of bounds.
                
                % NOTE: Does not work for dydx==0, but that should never happen here, since
                % y0 within bounds, dydx==0   ==>   y0==y1   ==>  y1 within bounds   ==>   Can not be here.
                xp = x0 + (yp - y0) / dydx;
                
                % IMPLEMENTATION NOTE: On rare occasions, the function will return yp~=y0 and xp==x0 du to numerical rounding of the latter.
                % One does NOT want to return xp==x0 in order to change y. We want y=y(x) to describe a one-valued
                % function. Therefore, return xp==x0+epsilon.
                if xp == x0
                    xp = xp + eps(xp);
                    
                    % ASSERTION
                    if xp == x0
                        error('linear_extrapolate:Assertion', 'xp==x0')
                    end
                end
                
                % IMPLEMENTATION NOTE: If x1 == xp, then y1 can still be marginally out of bounds due to numerical
                % error. Therefore, use yp (bounded) instead.
                y1 = yp;
            else
                    
                % CASE: y1 within bounds.
                xp = x1;
                y1 = yp;
                % NOTE: Using yp == y1 always absolutely within closed interval yMin--yMax. No risk of slightly
                % exceeding it due to numerical inaccuracies.
                
            end
            
            % IMPLEMENTATION NOTE: Should in principle only set this if y1 within bounds, but that does not follow the
            % if-then-else above due to numerical errors and seems harder to implement than at first sight.
            i = (x0 <= xArray) & (xArray <= x1);
            yArray(i) = y0 + (xArray(i) - x0) * dydx;
            
            
            
            % ASSERTION
            if xp > x1
                error('linear_extrapolate_limit:Assertion', 'xp > x1')
            elseif xp < x0
                error('linear_extrapolate_limit:Assertion', 'xp < x0')
            end
        end
        
        
        
        % yp : If inside yMin-yMax, then equal to y, otherwise that value yMin, yMax which is closest/exceeded.
        function [yp, belowMin, aboveMax] = outside_boundaries(y, yMin, yMax)
            belowMin = false;
            aboveMax = false;
            
            if y < yMin
                belowMin = true;
                yp       = yMin;
            elseif y > yMax
                aboveMax = true;
                yp       = yMax;
            else
                yp = y;
            end
        end
        
    end
end
