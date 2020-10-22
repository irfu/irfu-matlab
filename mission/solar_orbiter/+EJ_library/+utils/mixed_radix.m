%
% Functions for working with mixed radix numbers.
%
% Ex: Can interpret day+hour+minute+second as digits in mixed radix number (if
% not using leap seconds).
%
%
% NOTE: Compare EJ_library.utils.extract_subinteger.
%
%
% CONVENTIONS
% ===========
% MRD = Mixed Radix Digits, i.e. the digits in mixed radix number.
% --
% bases : (iDigit)
% mrd   : (iNbr, iDigit). iDigit=1 <==> least significant digit.
%         NOTE: This means that digits are in the opposite order to what is
%         convention when printing out digits as a literal.
% n     : (iNbr). Regular integer. May be negative.
% --
% It is implicit that all arguments will be rounded to integers.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-10-21
%
classdef mixed_radix   % < handle
    % PROPOSAL: Shorter class name.
    % PROPOSAL: Outside of utils package.
    % PROPOSAL: Automatic test code.
    

    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)
        % PROPOSAL: Permit infinitely high higest base.
        %   CON: Caller can simply set the highest base to a high value.
        %   CON: Ill-defined result for negative integers.
        % PROPOSAL: Permit negative integers (not MRD).
        %   CON: Caller can do mod before calling instead.
        %       CON: Extra work for caller.
        
        
        
        % NOTE: Permits mrd components to be outside of range as defined by
        % bases, including negative numbers.
        % RATIONALE: This is useful for in particular the highest component
        % where it is often useful to think of the highest base as
        % "infinite", or for negative values which can not be emulated with a
        % very high highest base.
        %   Ex: Day numbers+time of day.
        function n = MRD_to_integer(mrd, bases)
        
            % ASSERTIONS
            [nNbrs, nBases] = EJ_library.assert.sizes(mrd, [-1, -2], bases, [-2]);
            assert(all(bases >= 1))
            %assert(all(mrd >= 0, 'all'))
            for i = 1:nBases
                assert(all(mrd(:, i) < bases(i)), 'Mixed radix digits must be lower than digit base.')
            end
            
            mrd   = int64(mrd);
            bases = int64(bases);

            n = int64(zeros(nNbrs, 1));
            B = int64(1);
            for i = 1:nBases
                n = n + mrd(:, i) * B;
                B = B * bases(i);
            end
        end
        
        
        
        % NOTE: Permits n to be outside of range defined by bases. Result can be
        % interpreted as being used for n := mod(n, prod(bases)).
        % RATIONALE: Useful for in particular negative values where one can not
        % emulate the behaviour by setting a very large highest base.
        %   Ex: Negative day numbers+time of day.
        function mrd = integer_to_MRD(n, bases)
            % ASSERTIONS
            [nNbrs, nBases] = EJ_library.assert.sizes(n, [-1, 1], bases, [-2]);
            assert(all(bases >= 1))
            
            n     = int64(n);
            bases = int64(bases);
            
            mrd = int64(zeros(nNbrs, nBases));
            B   = int64(1);
            for i = 1:nBases
                mrd(:, i) = mod(idivide(int64(n), B, 'floor'), bases(i));
                B = B*bases(i);
            end
            
            % ASSERTION
            %assert(all(n<B), 'n is too large for the specified bases.')
        end
        
        
        
    end
    
end
