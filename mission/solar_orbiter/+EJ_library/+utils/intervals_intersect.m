%
% Determine whether pairs of intervals, u1 <= x <= u2 and v1 <= x <= v2, intersect (the way sets do).
%
% NOTE: Assumes that boundaries are included.
% NOTE: Upper limit < lower limit is equivalent to empty interval (set).
%
%
% ARGUMENTS
% =========
% u1, u2 : Same-sized numeric arrays.
% v1, v2 : Same-sized numeric arrays.
% If all of u1,u2,v1,v2 are non-scalar, then they must all have the same size.
% NOTE: 
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-05-20.
%
function intervalsIntersect = intervals_intersect(u1, u2, v1, v2)
% PROPOSAL: Return intersection (interval) instead.
%   PRO: Does not care about whether boundaries are included, as long as same interpretation for both input intervals.
%   Return result follows same convention as arguments.
% PROPOSAL: Implement using function that returns intersection instead.
% PROPOSAL: Additionally return intersection AND whether intersection is empty.
%   ~CON: Empty intersection depends on whether boundaries are included or not.
% NOTE: Inconsistent boundary inclusion in arguments sometimes leads to inconsistent (varying) boundary inclusion in
% return value.
%
%     w1 = max(u1, v1);
%     w2 = min(u2, v2);
%
% PROPOSAL: Flag for whether boundaries are included.
%   NOTE: Could have different policy for upper & lower boundary.
%       upperBoundaryIncluded
%       lowerBoundaryIncluded
%       CON: No. If one or two boundaries are included, then intersection w1==w2 must be interpreted as a non-intersection.

    intervalsIntersect = (u1 <= v2) & (v1 <= u2);
    
end