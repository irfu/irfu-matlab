%
% Determine whether a pairs of intervals on the real number line ([a,b] or
% (a,b)), intersect (intersect the way sets do).
%
% NOTE: Assumes that boundaries are included.
% NOTE: Upper limit < lower limit is equivalent to empty interval (set).
%
%
% ARGUMENTS
% =========
% u1, u2     : Same-sized numeric arrays.
% v1, v2     : Same-sized numeric arrays.
% edgePolicy : String constant for whether to use open or closed intervals.
% --
% If all of u1,u2,v1,v2 are all non-scalar, then they must all have the same
% size.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-05-20.
%
function intervalsIntersect = intervals_intersect(u1, u2, v1, v2, edgePolicy)
% PROPOSAL: Implement using function that returns intersection instead.
%   NOTE: This also depends on whether boundaries are included.
%
% PROPOSAL: Return intersection (interval) instead.
%   PRO: Does not care about whether boundaries are included, as long as one uses the same interpretation for both input intervals.
%        Return result follows same convention as arguments.
%       CON: NOT TRUE for when intervals touch and boundaries are excluded.
%           Ex: [0,1] and [1,2] intersect, but
%               (0,1) and (1,2) do NOT intersect, but use same boundary-included convention.
%
% PROPOSAL: Additionally return intersection AND whether intersection is empty.
%   ~CON: Empty intersection depends on whether boundaries are included or not.
% NOTE: Inconsistent boundary inclusion in arguments sometimes leads to inconsistent (varying) boundary inclusion in
% return value.
%
%     w1 = max(u1, v1);
%     w2 = min(u2, v2);
%
% PROPOSAL: Different policy for upper & lower boundary.
%   upperBoundaryIncluded
%   lowerBoundaryIncluded
%   CON: No. If one or two boundaries are included, then intersection w1==w2 must be interpreted as a non-intersection.
%
% PROPOSAL: Forbid negativ length interval. (zero-length OK).
%   NOTE: Mathematical expressions work also for negative length intervals.

switch(edgePolicy)
  case 'closed intervals'
    intervalsIntersect = (u1 <= v2) & (v1 <= u2);
  case 'open intervals'
    intervalsIntersect = (u1 < v2) & (v1 < u2);
  otherwise
    error('Illegal argument edgePolicy="%s".', edgePolicy)
end

end
