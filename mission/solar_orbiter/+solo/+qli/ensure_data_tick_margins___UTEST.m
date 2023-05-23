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

           % Only positive values
            test({[2, 3], [2.5,   2.5  ]}, [ 2*0.9;  3*1.1])
            test({[2, 3], [2*0.9, 3*1.1]}, [ 2*0.9;  3*1.1])
            test({[2, 3], [2,     3    ]}, [ 2*0.9;  3*1.1])
            test({[2, 3], [1,     4    ]}, [ 1;      4    ])

            % Only negative values
            test({[-3, -2], [-2.5,   -2.5  ]}, [ -3*1.1;  -2*0.9])
            test({[-3, -2], [-2*1.1, -2*0.9]}, [ -3*1.1;  -2*0.9])
            test({[-3, -2], [-3,     -2    ]}, [ -3*1.1;  -2*0.9])
            test({[-3, -2], [-4,     -1    ]}, [ -4;      -1    ])

            % Interval that covers both negative and positive numbers, with tick
            % limits<>0.
            test({[-2,  3], [-1,    2  ]}, [-2*1.1;  3*1.1])
            test({[-2,  3], [-3,    4  ]}, [-3;      4    ])

            % Intervals that have tick limit=zero.
            test({[ 0, 4], [  1,  2]}, [ 1;       4*1.1])      % BAD?!!
            test({[ 0, 4], [  0,  4]}, [ 0;       4*1.1])
            test({[ 0, 4], [ -1,  2]}, [-1;       4*1.1])
            test({[-4, 0], [ -2, -1]}, [-4*1.1;  -1    ])      % BAD?!!
            test({[-4, 0], [ -4,  0]}, [-4*1.1;   0    ])
            test({[-4, 0], [ -2,  1]}, [-4*1.1;   1    ])
        end



    end    % methods(Test)



end
