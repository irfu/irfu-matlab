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



    function test_plus_legal(testCase)
      % NOTE: Class has assertions against duplicated comment strings and
      % submitting zero comment strings. Can therefore not that in tests.

      GMVE_1 = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
        'Comment 1.'});
      GMVE_2 = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
        'Comment 2.'});
      GMVE_3 = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
        'Comment 3a.', 'Comment 3b.'});

      GMVE_4 = GMVE_1 + GMVE_2;
      testCase.assertEqual(...
        GMVE_4.get_str(), ...
        '2020-01-01 -- V1.2.3 -- Comment 1. | Comment 2.')

      GMVE5 = GMVE_1 + GMVE_3;
      testCase.assertEqual(...
        GMVE5.get_str(), ...
        '2020-01-01 -- V1.2.3 -- Comment 1. | Comment 3a. | Comment 3b.')

      GMVE_6 = GMVE_3 + GMVE_1;
      testCase.assertEqual(...
        GMVE_6.get_str(), ...
        '2020-01-01 -- V1.2.3 -- Comment 3a. | Comment 3b. | Comment 1.')
    end



    % Test incompatible date strings and BICAS version numbers.
    function test_plus_compatibility(testCase)
      function test_exc(GmveA, GmveB)
        testCase.verifyError(...
          @() (GmveA + GmveB), ...
          ?MException)
      end

      % Incompatible date strings.
      GMVE_1 = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
        'Comment 1.'});
      GMVE_2 = bicas.gamods.VersionEntry('2020-01-02', '1.2.3', {...
        'Comment 2.'});
      test_exc(GMVE_1, GMVE_2)

      % Incompatible BICAS version numbers.
      GMVE_3 = bicas.gamods.VersionEntry('2020-01-01', '1.2.3', {...
        'Comment 1.'});
      GMVE_4 = bicas.gamods.VersionEntry('2020-01-01', '1.2.9', {...
        'Comment 2.'});
      test_exc(GMVE_3, GMVE_4)

      % Incompatible date strings and BICAS version numbers.
      test_exc(GMVE_1, GMVE_1)
      test_exc(GMVE_3, GMVE_3)
    end



  end    % methods(Test)



end
