%
% matlab.unittest automatic test code for bicas.main().
%
% Since bicas.main() is a very top-level function, this code does not try to
% test everything, but only that which is easy.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef main___UTEST < matlab.unittest.TestCase



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)

            % Test successful call (error code).
            function test(varargin)
                errorCode = bicas.main(varargin{:});

                testCase.verifyEqual(errorCode, 0)
            end

            % Test unsuccessful call (non-zero error code)
            % NOTE: bicas.main() raises exception but catches it itself.
            function test_error(varargin)
                errorCode = bicas.main(varargin{:});

                testCase.verifyEqual(errorCode, 1)
            end
            %===================================================================

           test('--help')
           test('--version')
           test('--identification')
           test('--swdescriptor')

           test('--help', '--set', 'SWM.L1-L2_ENABLED', '0')
           test('--help', '--set', 'SWM.L1-L2_ENABLED', '1')
           test('--version', '--log',        '/ignored_path_to_log_file')

           % ----------------------------------
           % --log-matlab : Log file is created
           % ----------------------------------
           fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
           logFilePath = fullfile(fixture.Folder, 'bicas.log');
           testCase.verifyFalse(isfile(logFilePath))
           test('--version', '--log-matlab', logFilePath)
           testCase.verifyTrue(isfile(logFilePath))

           test_error()    % No CLI argument.
           test_error('illegal_argument')
           test_error('--illegal_argument')
        end



    end    % methods(Test)



end
