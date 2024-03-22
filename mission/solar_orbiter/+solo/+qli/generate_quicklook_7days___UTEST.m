%
% matlab.unittest automatic test code for solo.qli.generate_quicklook_7days().
%
%
% NOTES
% =====
% * solo.qli.const contains some constants which can be used to disable certain
%   aspects of plotting for speeding it up during manual testing. These need to
%   be considered when running this test code manually.
% * Tests use data from solo.qli.testdata.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef generate_quicklook_7days___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    outputPath
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods (TestMethodSetup)



    function create_output_directories(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.outputPath = F.Folder;
    end



  end



  %##########
  %##########
  % TEARDOWN
  %##########
  %##########
  methods (TestMethodTeardown)



    function close_figures(testCase)
      close all
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_no_data(testCase)
      QuicklooksTint = irf.tint(...
        '2023-12-27T00:00:00.00Z', ...
        '2024-01-03T00:00:00.00Z');
      SpacePosTint = irf.tint(...
        '2023-12-01T00:00:00.00Z', ...
        '2024-02-02T00:00:00.00Z');

      Data = solo.qli.testdata.generate_empty_test_data(SpacePosTint);

      solo.qli.generate_quicklook_7days(Data, testCase.outputPath, QuicklooksTint, [])
    end



    function test_all_data(testCase)
      QuicklooksTint = irf.tint(...
        '2023-12-27T00:00:00.00Z', ...
        '2024-01-03T00:00:00.00Z');
      SpacePosTint = irf.tint(...
        '2023-12-01T00:00:00.00Z', ...
        '2024-02-02T00:00:00.00Z');
      % Short time interval with B data to speed up test.
      BTint = irf.tint(...
        '2024-01-01T01:00:00', ...
        '2024-01-01T01:02:00');

      Data = solo.qli.testdata.generate_test_data(...
        QuicklooksTint, SpacePosTint, BTint);

      solo.qli.generate_quicklook_7days(Data, testCase.outputPath, QuicklooksTint, [])
    end



  end    % methods(Test)



end
