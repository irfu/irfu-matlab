%
% UNFINISHED / EXPERIMENTAL
%
% For a given function x-->y, replace every NaN with the nearest value (in x distance), unless it exceeds some
% threeshold.
%
%
% ARGUMENTS
% =========
% x : 1D numeric array. Must be finite.
% y : 1D numeric array. Same length as x. Must be finite or NaN.
% x,y represents a function x-->y.
%
%
% RETURN VALUES
% =============
%
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-08-14.
%
function [y] = fill_NaN(x, y, maxDist)
    % PROPOSAL: Automatic test code.
    %
    % PROPOSAL: Better name.
    %   ~data gaps, ~replace, ~nearest, ~fill
    %   fill_data_gaps_with_nearest
    
    % ASSERTIONS
    EJ_library.assert.vector(x)
    EJ_library.assert.vector(y)
    assert(isnumeric(x) && isnumeric(y))
    assert(numel(x) == numel(y), 'x and y have different lengths.')
    assert(isscalar(maxDist))
    assert(all(~isnan(x)))
    
    
    
    bf = ~isnan(y);
    xf  = x(bf);
    y1f = y(bf);
    
    for i = 1:numel(x)
        if isnan(y(i))
            [xDist, jNearest] = min(abs(x(i) - xf));
            
            if xDist <= maxDist
                y(i) = y1f(jNearest);
            end
        end
    end
end
