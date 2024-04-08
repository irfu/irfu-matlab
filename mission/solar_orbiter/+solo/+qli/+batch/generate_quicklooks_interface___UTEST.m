%
% matlab.unittest automatic test code for
% solo.qli.batch.generate_quicklooks_interface() and indirectly
% solo.qli.batch.interface.get_days_from_time_interval().
%
%
% NOTES
% =====
% Does not test the details of separate commands, since that is the
% purpose of other test code which is closer to the relevant parts of the code:
%   * solo.qli.batch.generate_quicklooks___UTEST
%   * solo.qli.batch.extract_dataset_dates_from_logs___UTEST
%   * solo.qli.batch.fmd___UTEST
%   * solo.qli.batch.interface_*_UTEST
% --
% Only testing the TIME_INTERVAL algorithm since it is
% the easiest one to use here. Other algorithms are tested through the
% corresponding subfunctions.
% --
% Indirectly testing the TIME_INTERVAL subfunction
% (solo.qli.batch.interface.get_days_from_time_interval()) which thus does
% not really need any other tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef generate_quicklooks_interface___UTEST < matlab.unittest.TestCase



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  properties(TestParameter)
    % Technically, additional properties of testCase objects with cell array
    % default values. Test methods with arguments with the same name will be
    % called once for every element in the cell arrays.
    OPERATION_ID = {'LIST', 'GENERATE'}
  end



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    Settings
    inputLogDir
    outputDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods (TestMethodSetup)



    function create_output_directories(testCase)
      VhtFixture      = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      OutputFixture   = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      %InputLogFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      Settings             = [];
      Settings.irfLogoPath = solo.qli.testdata.get_test_logo_path();
      Settings.vhtDir      = VhtFixture.Folder;
      Settings.Gql         = solo.qli.batch.GenerateQuicklooksTest();
      % Settings which are required but should not be relevant.
      Settings.datasetDirsCa         = {};
      Settings.LogFileDirPatternDict = dictionary();

      testCase.Settings    = Settings;
      %testCase.inputLogDir = InputLogFixture.Folder;
      testCase.outputDir   = OutputFixture.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_time_interval_zero_days(testCase, OPERATION_ID)
      Settings = testCase.Settings; %#ok<*PROP>

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        OPERATION_ID, 'TIME_INTERVAL', '2024-01-01', '2024-01-01' ...
        )

      testCase.assertEqual(Settings.Gql.Dt24h6h2hArray, solo.qli.const.EMPTY_DT_ARRAY)
      testCase.assertEqual(Settings.Gql.Dt7daysArray,   solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test_time_interval_two_days(testCase, OPERATION_ID)
      Settings = testCase.Settings;

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        OPERATION_ID, 'TIME_INTERVAL', '2024-01-01', '2024-01-04' ...
        )

      switch(OPERATION_ID)
        case 'LIST'
          ExpDt24h6h2hArray = solo.qli.const.EMPTY_DT_ARRAY;
          ExpDt7daysArray   = solo.qli.const.EMPTY_DT_ARRAY;

        case 'GENERATE'
          ExpDt24h6h2hArray = solo.qli.utils.umdt('2024-01-01') + caldays([0:2]');
          ExpDt7daysArray   = solo.qli.utils.umdt('2023-12-27') + calweeks([0:1]');

        otherwise
          error('')
      end

      testCase.assertEqual(Settings.Gql.Dt24h6h2hArray, ExpDt24h6h2hArray)
      testCase.assertEqual(Settings.Gql.Dt7daysArray,   ExpDt7daysArray)
    end



  end    % methods(Test)



end
