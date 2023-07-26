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
            % Zero data interval (uses data interval for deriving margin).
            test({[2, 5], [3      3    ], 'linear'}, [   3;  3  ])
            test({[],     [3      3    ], 'linear'}, [   3;  3  ])

            test({[],     [2      4    ], 'linear'}, [   2;  4  ])

            test({[1, 2, 3, 4, 5, 6], [3    4  ], 'linear'}, [2.95; 4.05])
            test({[1, 2, 3, 4, 5, 6], [2.1  4.9], 'linear'}, [2.1;  4.9])
            test({[1, 2, 3, 4, 5, 6], [7    8  ], 'linear'}, [7  ;  8  ])
            test({[1, 2, 3, 4, 5, 6], [0    8  ], 'linear'}, [0  ;  8  ])

            %=============
            % scale = log
            %=============
            % IMPLEMENTATION NOTE: Using constants makes it easier to update the
            % tests if the corresponding constants in the tested function are
            % modified.
            C = 0.9;
            D = 1.1;

            % Only positive values
            test({[ 2,  5], [ 3,    4  ], 'log'}, [ 3  ;  4  ])
            test({[ 2,  5], [ 2,    5  ], 'log'}, [ 2*C;  5*D])
            test({[ 2,  5], [ 2*C,  5*D], 'log'}, [ 2*C;  5*D])
            test({[ 2,  5], [ 1,    6  ], 'log'}, [ 1;    6  ])

            % Only negative values
            test({[-5, -2], [-4,   -3  ], 'log'}, [-4  ; -3  ])
            test({[-5, -2], [-5,   -2  ], 'log'}, [-5*D; -2*C])
            test({[-5, -2], [-5*D, -2*C], 'log'}, [-5*D; -2*C])
            test({[-5, -2], [-6,   -1  ], 'log'}, [-6;   -1  ])

            % Interval that covers both negative and positive numbers, with tick
            % limits<>0.
            test({[-2,  3], [ -1,   2  ], 'log'}, [-1;    2  ])
            test({[-2,  3], [ -3,   4  ], 'log'}, [-3;    4  ])

            % Intervals that have tick limit=zero.
            % Worked in old implementation.
            %test({[ 0,  4], [  1,   2  ], 'log'}, [-0.4;  4*D])
            %test({[ 0,  4], [  0,   4  ], 'log'}, [-0.4;  4*D])
            %test({[ 0,  4], [ -1,   2  ], 'log'}, [-1;    4*D])
            %test({[-4,  0], [ -2,  -1  ], 'log'}, [-4*D;  0.4])
            %test({[-4,  0], [ -4,   0  ], 'log'}, [-4*D;  0.4])
            %test({[-4,  0], [ -2,   1  ], 'log'}, [-4*D;  1  ])
        end



    end    % methods(Test)



end
