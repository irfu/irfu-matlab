%
% matlab.unittest automatic test code for bicas.utils.Timekeeper.
%
% NOTE: Testing of output strings "has to" be manual (not really, but writing
% tests for it are harder).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Timekeeper___UTEST < matlab.unittest.TestCase



    %#####################
    %#####################
    % CONSTANT PROPERTIES
    %#####################
    %#####################
    properties(Constant)
        DELAY_SEC = 0.0
    end



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        % Zero quantities
        function test_stop_log_0(testCase)
            L = bicas.Logger('human-readable', false);

            Tmk = bicas.utils.Timekeeper('CODE_NAME', L);
            pause(bicas.utils.Timekeeper___UTEST.DELAY_SEC)
            Tmk.stop_log()
        end



        % One quantity
        function test_stop_log_1(testCase)
            L = bicas.Logger('human-readable', false);

            Tmk = bicas.utils.Timekeeper('CODE_NAME', L);
            pause(bicas.utils.Timekeeper___UTEST.DELAY_SEC)
            Tmk.stop_log(10, 'gadget')

            testCase.verifyError(...
                @() Tmk.stop_log(), ...
                ?MException)
            testCase.verifyError(...
                @() Tmk.stop_log(10, 'gadget'), ...
                ?MException)
        end



        % Multiple quantities
        function test_stop_log_N(testCase)
            L = bicas.Logger('human-readable', false);

            Tmk = bicas.utils.Timekeeper('CODE_NAME', L);
            pause(bicas.utils.Timekeeper___UTEST.DELAY_SEC)
            Tmk.stop_log(10, 'gadget', 20, 'bin', 1000, 'byte')
        end



    end    % methods(Test)



end
