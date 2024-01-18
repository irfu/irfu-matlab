%
% matlab.unittest automatic test code for
% bicas.settings_value_to_display_str___UTEST().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-12
%
classdef settings_value_to_display_str___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)
            %=========
            % STRINGS
            %=========
            bicas.settings_value_to_display_str___UTEST.test(testCase, {...
                '',   '""';
                'text', '"text"';
            })

            %======================
            % 1D VECTOR OF NUMBERS
            %======================
            % 0x0 array ==> ERROR
            testCase.verifyError(...
                @() bicas.settings_value_to_display_str(zeros(0, 0)), ...
                ?MException)

            % 0x1, 1x0 array ==> OK
            bicas.settings_value_to_display_str___UTEST.test(testCase, {
                zeros(0, 1), '[]';
                zeros(1, 0), '[]';
                1.0, '1';
                int64(1), '1'
                3.14, '3.14';
            })

            % Row and column vector.
            bicas.settings_value_to_display_str___UTEST.test(testCase, {
                [3.14, 29],  '[3.14, 29]';
                [3.14, 29]', '[3.14, 29]';
            })

            %=================
            % LOGICAL/BOOLEAN
            %=================
            bicas.settings_value_to_display_str___UTEST.test(testCase, {...
                false, 'false';
                true, 'true';
            })

            %============
            % CELL ARRAY
            %============
            % 0x0 array ==> ERROR
            testCase.verifyError(...
                @() bicas.settings_value_to_display_str(zeros(0, 0)), ...
                ?MException)
            % 0x1, 1x0 ==> OK
            bicas.settings_value_to_display_str___UTEST.test(testCase, {
                cell(1,0), '{}';
                cell(0,1), '{}';
            })

            bicas.settings_value_to_display_str___UTEST.test(testCase, {
                {1}, '{1}'
            })

            % Row and column vector.
            bicas.settings_value_to_display_str___UTEST.test(testCase, {
                {3.14, 29, 'text'},  '{3.14, 29, "text"}';
                {3.14, 29}', '{3.14, 29}';
            })
        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        function test(testCase, table)
            assert(size(table, 2) == 2)

            for i = 1:size(table, 1)
                value     = table{i, 1};
                expOutput = table{i, 2};

                testCase.verifyEqual(...
                    bicas.settings_value_to_display_str(value), ...
                    expOutput)
            end
        end



    end    % methods(Static, Access=private)



end
