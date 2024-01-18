%
% matlab.unittest automatic test code for bicas.tf.drt___UTEST().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-12
%
classdef drt___UTEST < matlab.unittest.TestCase



    properties(TestParameter)
        % "All" legal combinations of detDegreeOf, retEnabled which returns the
        % original signal if unity TF.
        %
        % NOTE: Can get "Warning: Polynomial is badly conditioned" if using too
        % high detDegreeOf.
        UNITY_DET_RET_SETTINGS = {{-2, false}, {0, true}, {1, true}, {2, true}, {4, true}}
        %UNITY_DET_RET_SETTINGS = {{5, true}}
    end



    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        % De- & retrending disabled.
        function test_disabled(testCase)
            t = [0:0.01:1]';

            % Complex, arbitrary function.
            y1a = exp(t) .* cos(t) + t.^3;

            [y12b, y2a] = bicas.tf.drt___UTEST.run_drt(y1a, {-2, false}, 29);

            testCase.verifyEqual(y1a, y12b)
            testCase.verifyEqual(y12b, y2a)
        end



        % Restore arbitrary signal completely.
        function test_unity(testCase, UNITY_DET_RET_SETTINGS)

            t = [0:0.01:1]';

            % Complex, arbitrary function.
            y1a = exp(t) .* cos(5*t) + t.^3;

            % NOTE: Must use scale=1 to obtain original function.
            [y12b, y2a] = bicas.tf.drt___UTEST.run_drt(y1a, UNITY_DET_RET_SETTINGS, 1);

            % ASSERTIONS
            testCase.verifyEqual(y1a, y2a, 'AbsTol', 1e-15)
            if UNITY_DET_RET_SETTINGS{1} >= 0    % DRT not disabled entirely.
                testCase.assertTrue(max(abs(y1a-y12b)) > 0.1)
                testCase.assertTrue(max(abs(y2a-y12b)) > 0.1)
            end
        end



        function test0(testCase)
            x = [-1:0.1:1]';    % "Normalized".

            y1a = 2 + 3*x;

            [y12b, y2a] = bicas.tf.drt___UTEST.run_drt(y1a, {0, true}, 1);

            testCase.verifyEqual(y12b, 3*x,   'AbsTol', 1e-15)
            testCase.verifyEqual(y2a,  2+3*x, 'AbsTol', 1e-15)



            % NOTE: Test retrending scaling.
            [y12b, y2a] = bicas.tf.drt___UTEST.run_drt(y1a, {1, true}, 2);

            % 1e-15 works on irony but fails in GitHub CI.
            testCase.verifyEqual(y12b,   0*x, 'AbsTol', 1e-14)
            testCase.verifyEqual(y2a,  4+6*x, 'AbsTol', 1e-14)



            [y12b, y2a] = bicas.tf.drt___UTEST.run_drt(y1a, {0, false}, 1);

            testCase.verifyEqual(y12b,   3*x, 'AbsTol', 1e-15)
            testCase.verifyEqual(y2a,    3*x, 'AbsTol', 1e-15)

        end



    end    % methods(Test)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % NOTE: Does not modify signal between detrending and retrending
        % (y1b-->y2b).
        function [y12b, y2a] = run_drt(y1a, drtInitArgsCa, retScaleFactor)

            Drt = bicas.tf.drt(drtInitArgsCa{:});
            y12b = Drt.detrend(y1a);

            y2a = Drt.retrend(y12b, retScaleFactor);
        end



    end    % methods(Static, Access=private)



end
