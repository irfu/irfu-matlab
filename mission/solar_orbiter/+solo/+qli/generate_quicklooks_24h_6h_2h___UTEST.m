%
% matlab.unittest automatic test code for
% solo.qli.generate_quicklooks_24h_6h_2h().
%
%
% NOTES
% =====
% * solo.qli.const contains some constants which can be used to disable
%   certain aspects of plotting for speeding it up during manual testing. These
%   need to be considered when running this test code manually.
% * Tests use data from solo.qli.testdata.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef generate_quicklooks_24h_6h_2h___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    OutputPaths
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods (TestMethodSetup)



    function create_output_directories(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.OutputPaths = solo.qli.utils.create_output_directories(F.Folder);
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



    % Plotting without data is important for testing many special cases.
    function test_no_data(testCase)
      QuicklooksTint = irf.tint(...
        '2024-01-01T00:00:00.00Z', ...
        '2024-01-02T00:00:00.00Z');
      SpacePosTint = irf.tint(...
        '2023-12-01T00:00:00.00Z', ...
        '2024-02-02T00:00:00.00Z');

      Data = solo.qli.testdata.generate_empty_test_data(SpacePosTint);

      solo.qli.generate_quicklooks_24h_6h_2h(Data, testCase.OutputPaths, QuicklooksTint, [])
    end



    function test_all_data(testCase)
      QuicklooksTint = irf.tint(...
        '2024-01-01T00:00:00.00Z', ...
        '2024-01-02T00:00:00.00Z');
      SpacePosTint = irf.tint(...
        '2023-12-01T00:00:00.00Z', ...
        '2024-02-02T00:00:00.00Z');
      % Short time interval with B data to speed up test.
      BTint = irf.tint(...
        '2024-01-01T01:00:00', ...
        '2024-01-01T01:02:00');

      Data = solo.qli.testdata.generate_test_data(...
        QuicklooksTint, SpacePosTint, BTint);
      irfLogoPath = solo.qli.testdata.get_test_logo_path();

      solo.qli.generate_quicklooks_24h_6h_2h(Data, testCase.OutputPaths, QuicklooksTint, irfLogoPath)
    end



  end    % methods(Test)



end
