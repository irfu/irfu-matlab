%
% matlab.unittest automatic test code for solo.qli.ensure_data_tick_margins().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ensure_data_tick_margins___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)
            function test(inputsCa, expOutput)
                actOutput = solo.qli.ensure_data_tick_margins(inputsCa{:});
                testCase.verifyEqual(actOutput, expOutput)
            end

            %===================================================================

            %=============
            % scale = lin
            %=============
            % Zero data interval (uses tick interval for deriving margin).
            test({[2, 3], [2.5,   2.5  ], 'linear'}, [ 1.9;  3.1])
            test({[2, 3], [2,     3    ], 'linear'}, [ 1.9;  3.1])
            test({[2, 3], [1.5,   3.5  ], 'linear'}, [ 1.5;  3.5])

            %=============
            % scale = log
            %=============
            % IMPLEMENTATION NOTE: Using constants makes it easier to update the
            % tests if the corresponding constants in the tested function are
            % modified.
            C = 0.9;
            D = 1.1;

            % Only positive values
            test({[ 2,  3], [ 2.5,  2.5], 'log'}, [ 2*C;  3*D])
            test({[ 2,  3], [ 2*C,  3*D], 'log'}, [ 2*C;  3*D])
            test({[ 2,  3], [ 2,    3  ], 'log'}, [ 2*C;  3*D])
            test({[ 2,  3], [ 1,    4  ], 'log'}, [ 1;    4  ])

            % Only negative values
            test({[-3, -2], [-2.5, -2.5], 'log'}, [-3*D; -2*C])
            test({[-3, -2], [-2*D, -2*C], 'log'}, [-3*D; -2*C])
            test({[-3, -2], [-3,   -2  ], 'log'}, [-3*D; -2*C])
            test({[-3, -2], [-4,   -1  ], 'log'}, [-4;   -1  ])

            % Interval that covers both negative and positive numbers, with tick
            % limits<>0.
            test({[-2,  3], [ -1,   2  ], 'log'}, [-2*D;  3*D])
            test({[-2,  3], [ -3,   4  ], 'log'}, [-3;    4  ])

            % Intervals that have tick limit=zero.
            test({[ 0,  4], [  1,   2  ], 'log'}, [-0.4;  4*D])
            test({[ 0,  4], [  0,   4  ], 'log'}, [-0.4;  4*D])
            test({[ 0,  4], [ -1,   2  ], 'log'}, [-1;    4*D])
            test({[-4,  0], [ -2,  -1  ], 'log'}, [-4*D;  0.4])
            test({[-4,  0], [ -4,   0  ], 'log'}, [-4*D;  0.4])
            test({[-4,  0], [ -2,   1  ], 'log'}, [-4*D;  1  ])
        end



    end    % methods(Test)



end
