%
% Collection of trivial ~utility functions shared by summary plots.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-03-15
%
classdef utils   % < handle
    % PROPOSAL: Automatic test code.
    
    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)
        
        
        
        % Merge two ~zVars, with identical timestamps, but with non-overlapping
        % regions of non-NaN indices.
        function zv = merge_zvs(zv1, zv2)
            EJ_library.assert.sizes(...
                zv1, [-1], ...
                zv2, [-1])
            
            b1 = ~isnan(zv1);
            b2 = ~isnan(zv2);
            
            assert(~all(b1 & b2))
            
            zv     = zv1;
            zv(b2) = zv2(b2);
            
        end
    end
    
end