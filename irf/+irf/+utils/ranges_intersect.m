%
% If one regards a 1D vector of numeric values as defining a range (min to max),
% return whether v1 and v2 intersect (have a non-empty intersection).
%
%
% NOTE: Ranges that only have exactly one number in common count as overlapping,
% i.e. the boundaries are included.
%
% ABOLISH AND DELETE?
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% Initially created <2020-04-09.
%
function rangesIntersect = ranges_intersect(v1, v2)

irf.assert.vector(v1)
irf.assert.vector(v2)

rangesIntersect = irf.utils.intervals_intersect(...
  min(v1), max(v1), min(v2), max(v2), 'closed intervals');   % NOTE: Includes boundaries.
end
