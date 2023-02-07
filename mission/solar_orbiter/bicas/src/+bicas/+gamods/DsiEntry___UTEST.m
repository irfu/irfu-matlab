%
% matlab.unittest automatic test code for bicas.gamods.DsiEntry.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef DsiEntry___UTEST < matlab.unittest.TestCase

    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)

            de = bicas.gamods.DsiEntry();

            ve = bicas.gamods.VersionEntry('2020-01-01', '1.0.0', {'Comment1.'});
            de.add_version_entry(ve)

            % Add entry with reused date.
            ve = bicas.gamods.VersionEntry('2020-01-01', '2.0.0', {'Comment2.'});
            testCase.verifyError(...
                @() de.add_version_entry(ve), ...
                ?MException)

            % Add entry with reused BICAS version.
            ve = bicas.gamods.VersionEntry('2021-01-01', '1.0.0', {'Comment3.'});
            % testCase.verifyError(...
            %     @() de.add_version_entry(ve), ...
            %     ?MException)
            de.add_version_entry(ve)

            % Add entry with no comments.
            ve = bicas.gamods.VersionEntry('2020-01-01', '2.0.0', {});
            testCase.verifyError(...
                @() de.add_version_entry(ve), ...
                ?MException)

            % Add second entry (successfully).
            ve = bicas.gamods.VersionEntry('2022-01-01', '2.0.0', {'Comment4.'});
            de.add_version_entry(ve)

            actStrCa = de.get_MODS_strings_CA();
            expStrCa = {...
                '2020-01-01 -- V1.0.0 -- Comment1.'; ...
                '2021-01-01 -- V1.0.0 -- Comment3.'; ...
                '2022-01-01 -- V2.0.0 -- Comment4.'; ...
            };
            testCase.verifyEqual(actStrCa, expStrCa)
        end



    end    % methods(Test)

end
