%
% If one regards a 1D vector of numeric values as defining a range (min to max), return whether v1 and v2 intersect
% (have a non-empty intersection).
%
%
% NOTE: Ranges that only have exactly one number in common count as overlapping, i.e. the boundaries are included.
% 
%
% Author: Erik P G Johansson
% Initially created <2020-04-09.
%
function rangesIntersect = ranges_intersect(v1, v2)

    EJ_library.assert.vector(v1)
    EJ_library.assert.vector(v2)
    
    %rangesIntersect = (min(v2) <= max(v1)) && (min(v1) <= max(v2));   % NOTE: Equality.
    rangesIntersect = EJ_library.utils.intervals_intersect(min(v1), max(v1), min(v2), max(v2));   % NOTE: Includes boundaries.
end
