%
% matlab.unittest automatic test code for bicas.gamods.VersionEntry.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef VersionEntry___UTEST < matlab.unittest.TestCase

    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)


        function test_constructor(testCase)
            test_constructor_exc('2020-20-01', '1.2.3', {})
            test_constructor_exc('2020-02-01', 'V1.2.3', {})
            test_constructor_exc('2020-02-01', '1.2.3', {123, 'Comment1'})

            function test_constructor_exc(varargin)
                testCase.verifyError(...
                    @() bicas.gamods.VersionEntry(varargin{:}), ...
                    ?MException)
            end
        end



        function test_add_comments(testCase)
            ve = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {'Comment1.'});
            ve2 = ve.add_comments({'Comment2.'});

            testCase.verifyEqual(ve.bicasVersionStr, ve2.bicasVersionStr)
            testCase.verifyEqual(ve.dateStr,         ve2.dateStr)

            testCase.verifyEqual(ve.commentsCa,  {'Comment1.'})
            testCase.verifyEqual(ve2.commentsCa, {'Comment1.'; 'Comment2.'})
        end



        function test_get_str(testCase)
            ve = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {});
            testCase.verifyError(@() ve.getstr(), ?MException)

            ve = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
                'A first comment.'});
            actStr = ve.get_str();
            expStr = '2020-01-01 -- V1.2.3 -- A first comment.';
            testCase.verifyEqual(actStr, expStr)

            ve = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
                'A first comment.', 'A second comment.'});
            actStr = ve.get_str();
            expStr = '2020-01-01 -- V1.2.3 -- A first comment. | A second comment.';
            testCase.verifyEqual(actStr, expStr)
        end



    end    % methods(Test)

end
