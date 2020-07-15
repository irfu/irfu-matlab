%
% Split a fixed-length 1D "resource" into pieces according to rules.
%
% Should be useful for e.g. automatically placing and re-sizing graphical objects.
%
%
% ARGUMENTS
% =========
% totalSize      : Scalar number. Total size to be divided.
% fixedSizeArray : 1D array. Fixed (minimum/starting) size for each segment.
% weightArray    : 1D array. Weights for how to distribute the remaining "size", not in fixedSize, between segments.
%                  only relative magnitudes matter.
%
%
% RETURN VALUES
% =============
% sizeArray : 1D array. Sums up to totalSize (barring rounding errors).
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-07-08.
%
function [sizeArray] = distribute_segments(totalSize, fixedSizeArray, weightArray)
    % PROPOSAL: Automatic test code.
    % PROPOSAL: Better name
    %   PROPOSAL: ~split
    %   PROPOSAL: Imply "scarse resource", fixed "resource".
    %
    % TODO-DEC: Which package? utils? graph?
    
    DISALLOW_NEGATIVE_SIZES_WEIGHTS = 1;
    
    
    % ASSERTIONS:
    assert(isscalar(totalSize))
    EJ_library.assert.vector(fixedSizeArray)
    EJ_library.assert.vector(weightArray)
    
    % Normalize
    fixedSizeArray = fixedSizeArray(:);
    weightArray    = weightArray(:);
    
    totalFixedSize   = sum(fixedSizeArray);
    totalUnfixedSize = totalSize - totalFixedSize;
    totalWeights     = sum(weightArray);
    
    % ASSERTIONS
    assert(numel(fixedSizeArray) == numel(weightArray))
    assert(totalWeights ~= 0)
    
    if DISALLOW_NEGATIVE_SIZES_WEIGHTS
        % ASSERTIONS: Assert not using negative values. Algorithm should work anyway, but unsure what the application is.
        assert(totalSize >= 0)
        assert(all(fixedSizeArray >= 0))
        assert(all(weightArray >= 0))
        assert(totalFixedSize <= totalSize)
    end

    
    sizeArray      = fixedSizeArray + weightArray/totalWeights * totalUnfixedSize;
end
