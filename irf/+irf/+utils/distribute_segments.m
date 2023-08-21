%
% Split a fixed-length 1D "resource" into pieces according to rules.
%
% Should be useful for e.g. automatically placing and re-sizing graphical
% objects.
%
%
% ARGUMENTS
% =========
% totalSize      : Scalar number. Total size to be divided.
% fixedSizeArray : 1D array. Fixed (minimum/starting) size for each segment.
% weightArray    : 1D array. Weights for how to distribute the remaining "size",
%                  not in fixedSize, between segments. only relative magnitudes
%                  matter.
%
%
% RETURN VALUES
% =============
% sizeArray : 1D array. Sums up to totalSize (barring rounding errors).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-07-08.
%
function [sizeArray] = distribute_segments(totalSize, fixedSizeArray, weightArray)
% PROPOSAL: Automatic test code.
%
% PROPOSAL: Better name. Needs to be easier to find.
%   PROPOSAL: ~split
%       CON: Only determines exact locations, sizes; not number of segments.
%   PROPOSAL: Imply "scarse resource", fixed "resource".
%   PROPOSAL: ~autoscale
%   PROPOSAL: autoscale_segments
%       NOTE: There was an old function autoscale" (deleted 2020-07-16) that had subset of functionality.
%
% TODO-DEC: Which package? utils? graph?
%
% TODO-DEC: How deal with failure, when totalSize<sum(fixedArraySize)?
%   Ex: When using it for placing graphical objects it is natural to resort to some backup scaling.
%   PROPOSAL: Assertion policy.
%   PROPOSAL: Return success flag
%   PROPOSAL: Caller should check for this case by inspecting arguments.
%
% PROBLEM: Bad terms: fixed size, weights
%   PRO: Two and one words respectively.
%   PROPOSAL: A=fixed size, B=weights

% ASSERTIONS: Whether to assert not using negative values. Algorithm should
% work with negative values, but unsure what the application would be.
DISALLOW_NEGATIVE_SIZES_WEIGHTS = 1;



% ASSERTIONS
assert(isscalar(totalSize))
irf.assert.vector(fixedSizeArray)
irf.assert.vector(weightArray)
assert(numel(fixedSizeArray) == numel(weightArray))

% Normalize argument variable sizes.
fixedSizeArray   = fixedSizeArray(:);
weightArray      = weightArray(:);

% Algorithm help variables
totalFixedSize   = sum(fixedSizeArray);
totalUnfixedSize = totalSize - totalFixedSize;
totalWeights     = sum(weightArray);

% ASSERTIONS
% NOTE: Algorithm is undetermined without totalWeights<>0 (multiple
% solutions; practically divide-by-zero.)
assert(totalWeights ~= 0)
if DISALLOW_NEGATIVE_SIZES_WEIGHTS
  assert(totalSize >= 0)
  assert(all(fixedSizeArray >= 0))
  assert(all(weightArray >= 0))
  assert(totalFixedSize < totalSize, ...
    'Must have sum(fixedSizeArray) < totalSize, but that is not true.')
end



% "Algorithm"
sizeArray = fixedSizeArray + weightArray/totalWeights * totalUnfixedSize;
end
