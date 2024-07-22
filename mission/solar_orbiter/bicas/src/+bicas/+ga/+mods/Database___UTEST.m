%
% matlab.unittest automatic test code for bicas.ga.mods.Database.
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

      Gmdb = bicas.ga.mods.Database({});
      Gmdb = bicas.ga.mods.Database({'DSI_1', 'DSI_2'});



      % ==========
      % add_GMVE()
      % ==========

      % Add GMVE to zero DSIs.
      Gmve0 = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment for zero DSIs.'});
      Gmdb.add_GMVE({}, Gmve0)

      % Add to one DSI.
      Gmve1 = bicas.ga.mods.VersionEntry('2020-01-01', '3.0.0', {...
        'Comment1.'});
      Gmdb.add_GMVE({'DSI_1'}, Gmve1)

      % Add to two DSIs.
      Gmve2 = bicas.ga.mods.VersionEntry('2021-01-01', '4.0.0', {...
        'Comment2.'});
      Gmdb.add_GMVE({'DSI_1', 'DSI_2'}, Gmve2)



      % =====================
      % get_MODS_strings_CA()
      % =====================

      actGaModsStrCa = Gmdb.get_MODS_strings_CA('DSI_1');
      expGaModsStrCa = {Gmve1.get_str(); Gmve2.get_str()};
      testCase.assertEqual(actGaModsStrCa, expGaModsStrCa)

      actGaModsStrCa = Gmdb.get_MODS_strings_CA('DSI_2');
      expGaModsStrCa = {Gmve2.get_str()};
      testCase.assertEqual(actGaModsStrCa, expGaModsStrCa)

      testCase.assertError(...
        @() Gmdb.get_MODS_strings_CA('DSI_UNKNOWN'), ...
        ?MException)
    end



  end    % methods(Test)



end
