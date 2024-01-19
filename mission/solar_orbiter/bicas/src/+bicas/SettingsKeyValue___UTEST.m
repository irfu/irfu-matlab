%
% matlab.unittest automatic test code for bicas.SettingsKeyValue class.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef SettingsKeyValue___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test_constructor(testCase)
            Skv = bicas.SettingsKeyValue(99, 'default');

            testCase.assertEqual(Skv.valuesCa,       {99})
            testCase.assertEqual(Skv.valueSourcesCa, {'default'})
        end



        function test_override(testCase)
            Skv = bicas.SettingsKeyValue(99,  'default');
            Skv = Skv.override(          123, 'override');

            % NOTE: Column arrays.
            testCase.assertEqual(Skv.valuesCa,       {99;        123})
            testCase.assertEqual(Skv.valueSourcesCa, {'default'; 'override'})

            % Test assertion agains reusing valueSource.
            testCase.verifyError(...
                @() Skv.override(111, 'default'), ...
                ?MException)
            testCase.verifyError(...
                @() Skv.override(222, 'override'), ...
                ?MException)
        end



        % NOTE: Operators == and ~= are not defined.
        function test_equality(testCase)
            Skv1 = bicas.SettingsKeyValue(NaN, 'default');
            Skv2 = bicas.SettingsKeyValue(NaN, 'default');
            Skv3 = bicas.SettingsKeyValue(1,   'default');

            testCase.assertTrue( isequaln(Skv1, Skv2))
            testCase.assertFalse(isequaln(Skv1, Skv3))
        end



        function test0(testCase)

%             % Arbitrary number output variables.
%             function test(inputsCa, expOutputsCa)
%                 % Pre-allocate correct size for later assignment via function.
%                 actOutputs = cell(size(expOutputsCa));
%
%                 [actOutputs{:}] = FUNCTION_TO_TEST(inputsCa{:});
%                 testCase.verifyEqual(actOutputs, expOutputsCa)
%             end

%             % One output variable.
%             function test(inputsCa, expOutput)
%                 actOutput = FUNCTION_TO_TEST(inputsCa{:});
%                 testCase.verifyEqual(actOutput, expOutput)
%             end

%             function test_exc(varargin)
%                 testCase.verifyError(...
%                     @() FUNCTION_TO_TEST(varargin{:}), ...
%                     ?MException)
%             end
            %===================================================================


        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
    end    % methods(Static, Access=private)



end
