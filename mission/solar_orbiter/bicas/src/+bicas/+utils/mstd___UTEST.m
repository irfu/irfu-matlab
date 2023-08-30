%
% matlab.unittest automatic test code for bicas.utils.mstd().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-10 from older test code.
%
classdef mstd___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)

        
        
        function test0(testCase)
            
            function test(v, ref, iDim, mstd)
                actOutput = bicas.utils.mstd(v, ref, iDim);
                testCase.verifyEqual(actOutput, mstd)
            end
            %===================================================================
            
            % Function declaration: mstd = mstd(v, ref, iDim)

            % Empty data.
            test(zeros(0,3), 5, 1, NaN(1,3));
            test(zeros(0,3), 5, 2, NaN(0,1));
            test(zeros(0,3), 5, 3, NaN(0,3,1));

            % 1D vector data.
            test(ones(1,3), 5, 1, NaN( 1,3,1));
            test(ones(1,3), 5, 2, ones(1,1,1) * sqrt(3*4^2 / 2));
            test(ones(1,3), 5, 3, NaN( 1,3,1));

            % One data point.
            test([0], 0, 1, NaN);
            test([5], 5, 1, NaN);
            test([5], 4, 1, NaN);
            test([5], 4, 2, NaN);

            test([5,6,7; 2,3,4], 4, 1, sqrt([(1^2+2^2),     2^2+1^2,    3^2+0^2] / 1));
            test([5,6,7; 2,3,4], 4, 2, sqrt([(1^2+2^2+3^2); 2^2+1^2+0^2]         / 2));
            test([5,6,7; 2,3,4], 4, 3, NaN(2,3,1));

            % Test NaN sample.
            test([NaN,6,7; 2,3,4], 4, 1, sqrt([NaN, 2^2+1^2,     3^2+0^2] / 1));
            test([NaN,6,7; 2,3,4], 4, 2, sqrt([NaN; 2^2+1^2+0^2]          / 2));
            test([NaN,6,7; 2,3,4], 4, 3, NaN(2,3,1));

            % Test ref=NaN
            test([5,6,7; 2,3,4], NaN, 1, NaN(1,3,1));
            test([5,6,7; 2,3,4], NaN, 2, NaN(2,1,1));
            test([5,6,7; 2,3,4], NaN, 3, NaN(2,3,1));

        end
        
        
        
    end    % methods(Test)

    
    
end
