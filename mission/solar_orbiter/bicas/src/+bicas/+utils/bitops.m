%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-10-22
%
classdef bitops
    % PROPOSAL: Automatic test code.

    methods(Static)

        % Repeated bit-wise OR on all values in vector.
        %
        % ARGUMENTS
        % =========
        % a : Non-empty vector of non-negative integers.
        %
        function b = or(a)
            % PROPOSAL: Convert to simple function in bicas.proc.utils.
            % PROPOSAL: Extend to empty vectors. Return all bits=0.
            %   NOTE: An empty vector contains the integer type, so the function
            %         knows which integer type (and how many bits) to return.
            % PROPOSAL: Can be replaced by any().
            %   CON: No it can't. any() operates on logical, not the multiple
            %        bits in integers.

            assert(~isempty(a))

            b = a(1);
            for i = 2:numel(a)
                b = bitor(b, a(i));
            end
        end



    end

end
