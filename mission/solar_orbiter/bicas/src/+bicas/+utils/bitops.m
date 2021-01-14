%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-10-22
%
classdef bitops
    % PROPOSAL: Automatic test code.
    % PROPOSAL: Convert to simple function in proc_utils.
    

    
    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)

        % Repeated bit-wise OR on all values in vector.
        %
        % ARGUMENTS
        % =========
        % a : Non-empty vector of non-negative integers.
        function b = or(a)
            assert(~isempty(a))
            
            b = a(1);
            for i = 2:numel(a) 
                b = bitor(b, a(i));
            end
        end



    end
    
end