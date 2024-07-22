%
% matlab.unittest automatic test code for bicas.ga.mods.DsiEntry.
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



    function test_add_GMVE___one(testCase)
      Gmde = bicas.ga.mods.DsiEntry();

      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment1.'});
      Gmde.add_GMVE(Gmve)

      actStrCa = Gmde.get_MODS_strings_CA();
      expStrCa = {...
        '2020-01-01 -- V1.0.0 -- Comment1.'; ...
        };
      testCase.assertEqual(actStrCa, expStrCa)
    end



    function test_add_GMVE___reuse_date(testCase)
      Gmde = bicas.ga.mods.DsiEntry();

      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment 1.'});
      Gmde.add_GMVE(Gmve)

      % Add entry with reused date (different BICAS version).
      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '2.0.0', {'Comment 2.'});
      Gmde.add_GMVE(Gmve)

      actStrCa = Gmde.get_MODS_strings_CA();
      expStrCa = {...
        '2020-01-01 -- V1.0.0 -- Comment 1.'; ...
        '2020-01-01 -- V2.0.0 -- Comment 2.'; ...
        };
      testCase.assertEqual(actStrCa, expStrCa)
    end



    function test_add_GMVE___reuse_BICAS_version(testCase)
      Gmde = bicas.ga.mods.DsiEntry();

      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment 1.'});
      Gmde.add_GMVE(Gmve)

      % Add entry with reused BICAS version (different date).
      Gmve = bicas.ga.mods.VersionEntry('2021-01-01', '1.0.0', {'Comment 2.'});
      Gmde.add_GMVE(Gmve)

      actStrCa = Gmde.get_MODS_strings_CA();
      expStrCa = {...
        '2020-01-01 -- V1.0.0 -- Comment 1.'; ...
        '2021-01-01 -- V1.0.0 -- Comment 2.'; ...
        };
      testCase.assertEqual(actStrCa, expStrCa)
    end



    % First and second GMDE have identical dates and BICAS version.
    function test_add_GMVE___reuse_date_BICAS_version_1(testCase)
      Gmde = bicas.ga.mods.DsiEntry();

      % Add entry.
      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment 1.'});
      Gmde.add_GMVE(Gmve)

      % Add entry with reused date and BICAS version.
      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment 2.'});
      Gmde.add_GMVE(Gmve)

      actStrCa = Gmde.get_MODS_strings_CA();
      expStrCa = {...
        '2020-01-01 -- V1.0.0 -- Comment 1. | Comment 2.'; ...
        };
      testCase.assertEqual(actStrCa, expStrCa)
    end


    % second and third GMDE have identical dates and BICAS version, i.e.
    % there is another GMDE not involved in the merger of GMDEs.
    function test_add_GMVE___reuse_date_BICAS_version_2(testCase)
      Gmde = bicas.ga.mods.DsiEntry();

      % Add entry.
      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment 1.'});
      Gmde.add_GMVE(Gmve)
      Gmve = bicas.ga.mods.VersionEntry('2021-01-01', '2.0.0', {'Comment 2.'});
      Gmde.add_GMVE(Gmve)

      % Add entry with reused date and BICAS version.
      Gmve = bicas.ga.mods.VersionEntry('2021-01-01', '2.0.0', {'Comment 3.'});
      Gmde.add_GMVE(Gmve)

      actStrCa = Gmde.get_MODS_strings_CA();
      expStrCa = {...
        '2020-01-01 -- V1.0.0 -- Comment 1.'; ...
        '2021-01-01 -- V2.0.0 -- Comment 2. | Comment 3.'; ...
        };
      testCase.assertEqual(actStrCa, expStrCa)
    end



    function test_add_GMVE___reuse_date_illegal(testCase)

      Gmde = bicas.ga.mods.DsiEntry();

      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '1.0.0', {'Comment 1.'});
      Gmde.add_GMVE(Gmve)

      % Add entry with reused date (different BICAS version).
      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '2.0.0', {'Comment 2.'});
      Gmde.add_GMVE(Gmve)

      % Add entry with reused BICAS version (different date).
      Gmve = bicas.ga.mods.VersionEntry('2021-01-01', '2.0.0', {'Comment 3.'});
      Gmde.add_GMVE(Gmve)

      % Add entry with reused date (non-last) and BICAS version -- ILLEGAL
      Gmve = bicas.ga.mods.VersionEntry('2020-01-01', '2.0.0', {'Comment 4.'});
      testCase.assertError(...
        @() Gmde.add_GMVE(Gmve), ...
        ?MException)



      actStrCa = Gmde.get_MODS_strings_CA();
      expStrCa = {...
        '2020-01-01 -- V1.0.0 -- Comment 1.'; ...
        '2020-01-01 -- V2.0.0 -- Comment 2.'; ...
        '2021-01-01 -- V2.0.0 -- Comment 3.'; ...
        };
      testCase.assertEqual(actStrCa, expStrCa)
    end



  end    % methods(Test)



end
