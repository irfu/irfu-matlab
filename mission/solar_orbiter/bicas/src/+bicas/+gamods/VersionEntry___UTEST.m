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
            Gmve1 = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {'Comment1.'});
            Gmve2 = Gmve1.add_comments({'Comment2.'});

            testCase.verifyEqual(Gmve1.bicasVersionStr, Gmve2.bicasVersionStr)
            testCase.verifyEqual(Gmve1.dateStr,         Gmve2.dateStr)

            testCase.verifyEqual(Gmve1.commentsCa,  {'Comment1.'})
            testCase.verifyEqual(Gmve2.commentsCa, {'Comment1.'; 'Comment2.'})

            testCase.verifyError(@() bicas.gamods.VersionEntry('2020-01-01', '1.2.3', ...
                {}), ...
                ?MException)
            testCase.verifyError(@() bicas.gamods.VersionEntry('2020-01-01', '1.2.3', ...
                {'Comment without trailing period'}), ...
                ?MException)
        end



        function test_get_str(testCase)
            Gmve = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
                'A first comment.'});
            actStr = Gmve.get_str();
            expStr = '2020-01-01 -- V1.2.3 -- A first comment.';
            testCase.verifyEqual(actStr, expStr)

            Gmve = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
                'A first comment.', 'A second comment.'});
            actStr = Gmve.get_str();
            expStr = '2020-01-01 -- V1.2.3 -- A first comment. | A second comment.';
            testCase.verifyEqual(actStr, expStr)
        end



    end    % methods(Test)



end
