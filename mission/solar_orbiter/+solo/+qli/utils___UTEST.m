%
% UNFINISHED
%
% matlab.unittest automatic test code for solo.qli.utils.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_generate_data_source_info(testCase)

            % Merely test that function (1) does not crash, and (2) returns string.
            function test(inputsCa)
                actOutput = solo.qli.utils.generate_data_source_info(inputsCa{:});
                testCase.verifyInstanceOf(actOutput, 'char')
                testCase.verifyTrue(isrow(actOutput))
            end

            %===================================================================

            test({})
        end



        function test_get_plot_filename(testCase)

            function test(inputsCa, expOutput)
                actOutput = solo.qli.utils.get_plot_filename(inputsCa{:});
                testCase.verifyEqual(actOutput, expOutput)
            end

            %===================================================================

            Tint = EpochTT( ...
                [ ...
                    '2024-01-10T02:09:04.900000009Z'; ...
                    '2024-01-11T04:04:09.900000004Z'; ...
                ] ...
            );
            test({Tint}, '20240110T02_20240111T04.png')
        end


    end    % methods(Test)



end
