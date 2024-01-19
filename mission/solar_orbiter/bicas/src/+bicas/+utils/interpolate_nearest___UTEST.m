%
% matlab.unittest automatic test code for bicas.utils.interpolate_nearest().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-09-09 from older test code.
%
classdef interpolate_nearest___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)

            function test(inputsCa, expOutputsCa)
                % Pre-allocate correct size for later assignment via function
                actOutputs = cell(size(expOutputsCa));

                [actOutputs{:}] = bicas.utils.interpolate_nearest(inputsCa{:});
                testCase.verifyEqual(actOutputs, expOutputsCa)
            end
            %===================================================================

            test({1,   [], [], []        },  {[]});
            test({1,   [], [], ones(13,0)},  {ones(13,0)});

            test({0,   [], [], [1 2 3 4 5]}, {NaN(1,5)});

            test({0,   [3],   [13],    [1 2 3 4 5]},   {[NaN NaN 13  NaN NaN]});
            test({1,   [3],   [13],    [1 2 3 4 5]},   {[NaN 13  13  13  NaN]});
            test({2,   [3],   [13],    [1 2 3 4 5]},   {[13  13  13  13  13 ]});

            % Test xMargin
            test({0,   [3 4], [13 14], [1 2 3 4 5 6]}, {[NaN NaN 13 14 NaN NaN]});
            test({0.9, [3 4], [13 14], [1 2 3 4 5 6]}, {[NaN NaN 13 14 NaN NaN]});
            test({1,   [3 4], [13 14], [1 2 3 4 5 6]}, {[NaN 13  13 14 14  NaN]});
            test({Inf, [3 4], [13 14], [1 2 3 4 5 6]}, {[13  13  13 14 14  14 ]});

            test({1,   [3 4], [13 14], []}, {[]});

            % Non-incrementing xArray1,xArray2.
            test({1,   [3 5 4], [13 15 14], [3 4 5]}, {[13 14 15]});
            test({1,   [3 5 4], [13 15 14], [5 4 3]}, {[15 14 13]});

            % Row/column vectors xArray2, yArray2.
            test({1,   [3 4],   [13 14],    [3 4] }, {[13 14] });
            test({1,   [3 4],   [13 14],    [3 4]'}, {[13 14]'});

        end



    end    % methods(Test)



end
