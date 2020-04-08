% ARGUMENTS
% =========
% x    : Array of arbitrary numeric values.
% xNew : Scalar, or array of the same size as iSet. New value(s) that should be set at indices iSet.
% iSet : Array of indices.
%
%
% RETURN VALUE
% ============
% x : Modified x so that x(iSet) = xNew, and that other x components are uniformly scaled to maintain the same sum.
%
function x = reweight(x, xNew, iSet)
    % Variable naming convention
    % S = Set     = x components which should be set.
    % N = Not set = x components which are not explicitly set, but indirectly.
    %
    % PROPOSAL: Better name?
    % PROPOSAL: Move to ~utils?
    
    assert(numel(unique(iSet)) == numel(iSet))
    
    xSum1 = sum(x);
    iS = iSet;   clear iSet    % Change variable name.
    iN = setdiff(1:numel(x), iS);
    
    x(iS) = xNew;
    x(iN) = x(iN) * (xSum1 - sum(x(iS))) / sum(x(iN));
    
    %assert(sum(x) == xSum1)    % Does not always work due to rounding errors. Good enough for testing.
end