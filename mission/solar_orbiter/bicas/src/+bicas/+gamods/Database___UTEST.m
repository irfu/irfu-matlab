%
% matlab.unittest automatic test code for bicas.gamods.Database.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Database___UTEST < matlab.unittest.TestCase

    %##############
    %##############
    % TEST METHODS
    %##############
    %##############
    methods(Test)



        function test0(testCase)
            % ===========
            % Constructor
            % ===========
            Gmdb = bicas.gamods.Database({});
            Gmdb = bicas.gamods.Database({'DSI_1', 'DSI_2'});



            % ===================
            % add_version_entry()
            % ===================
            
            % Add VE to zero DSIs.
            ve = bicas.gamods.VersionEntry('2020-01-01', '1.0.0', {'Comment for zero DSIs.'});
            Gmdb.add_version_entry({}, ve)

            % Add to one DSI.
            ve1 = bicas.gamods.VersionEntry('2020-01-01', '3.0.0', {...
                'Comment1.'});
            Gmdb.add_version_entry({'DSI_1'}, ve1)

            % Add to two DSIs.
            ve2 = bicas.gamods.VersionEntry('2021-01-01', '4.0.0', {...
                'Comment2.'});
            Gmdb.add_version_entry({'DSI_1', 'DSI_2'}, ve2)



            % =====================
            % get_MODS_strings_CA()
            % =====================
            
            actGaModsStrCa = Gmdb.get_MODS_strings_CA('DSI_1');
            expGaModsStrCa = {ve1.get_str(); ve2.get_str()};
            testCase.verifyEqual(actGaModsStrCa, expGaModsStrCa)

            actGaModsStrCa = Gmdb.get_MODS_strings_CA('DSI_2');
            expGaModsStrCa = {ve2.get_str()};
            testCase.verifyEqual(actGaModsStrCa, expGaModsStrCa)
            
            testCase.verifyError(...
                @() Gmdb.get_MODS_strings_CA('DSI_UNKNOWN'), ...
                ?MException)

        end



    end    % methods(Test)

end
