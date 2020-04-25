%
% If one regards a 1D vector of numeric values as defining a range (min to max), return whether v1 and v2 intersect
% (have a non-empty intersection).
%
%
% NOTE: Ranges that only have exactly one number in common count as overlapping.
% 
%
% Author: Erik P G Johansson
% Initially created <2020-04-09.
%
function rangesOverlap = ranges_intersect(v1, v2)

    EJ_library.assert.vector(v1)
    EJ_library.assert.vector(v2)
    
    rangesOverlap = (min(v2) <= max(v1)) && (min(v1) <= max(v2));   % NOTE: 
end
