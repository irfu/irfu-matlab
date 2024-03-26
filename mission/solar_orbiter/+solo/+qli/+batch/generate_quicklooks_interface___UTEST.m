%
% matlab.unittest automatic test code for
% solo.qli.batch.generate_quicklooks_interface().
%
% NOTE: Does not test the details of separate commands, since that is the
% purpose of other test code which is closer to the relevant parts of the code:
%   * solo.qli.batch.generate_quicklooks___UTEST
%   * solo.qli.batch.extract_dataset_dates_from_logs___UTEST
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef generate_quicklooks_interface___UTEST < matlab.unittest.TestCase



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
      InputLogFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      Settings             = [];
      Settings.irfLogoPath = solo.qli.testdata.get_test_logo_path();
      Settings.vhtDir      = VhtFixture.Folder;
      Settings.Gql         = solo.qli.batch.GenerateQuicklooksTest();

      testCase.Settings    = Settings;
      testCase.inputLogDir = InputLogFixture.Folder;
      testCase.outputDir   = OutputFixture.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)
    % ======================
    % TIME_INTERVAL commands
    % ======================



    function test_time_interval_zero_days(testCase)
      Settings = testCase.Settings; %#ok<*PROP>

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'TIME_INTERVAL', '2024-01-01', '2024-01-01' ...
      )

      testCase.assertEqual(Settings.Gql.Dt24h6h2hArray, solo.qli.const.DT_EMPTY_ARRAY)
      testCase.assertEqual(Settings.Gql.Dt7daysArray,   solo.qli.const.DT_EMPTY_ARRAY)
    end



    function test_time_interval_two_days(testCase)
      Settings = testCase.Settings;

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'TIME_INTERVAL', '2024-01-01', '2024-01-03' ...
      )

      testCase.assertEqual(...
        Settings.Gql.Dt24h6h2hArray, ...
        datetime('2024-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') + caldays([0; 1]))
      testCase.assertEqual(...
        Settings.Gql.Dt7daysArray, ...
        datetime('2023-12-27T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'))
    end



  end    % methods(Test)
  methods(Test)
    % ==========================
    % GENERATE_FROM_LOG commands
    % ==========================



    function test_log_zero_days(testCase)
      Settings = testCase.Settings;

      % Empty log files / without matches.
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'LESIA_2024-01-01.log'), {})
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'SOAR_2024-01-01.log'), {})

      Settings.lesiaLogFileDirPattern = fullfile(testCase.inputLogDir, 'LESIA_*.log');
      Settings.soarLogFileDirPattern  = fullfile(testCase.inputLogDir, 'SOAR_*.log');

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'GENERATE_FROM_LOGS', 'LESIA', 'SOAR' ...
      )

      testCase.assertEqual(Settings.Gql.Dt24h6h2hArray, solo.qli.const.DT_EMPTY_ARRAY)
      testCase.assertEqual(Settings.Gql.Dt7daysArray,   solo.qli.const.DT_EMPTY_ARRAY)
    end



    function test_log_two_days(testCase)
      Settings = testCase.Settings;

      % Empty log files / without matches.
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'LESIA_2024-01-01.log'), ...
        {'solo_L2_rpw-tnr-surv_20230101_V03.cdf'})
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'SOAR_2024-01-01.log'), ...
        {'solo_L2_mag-rtn-normal_20230202_V99.cdf'})

      Settings.lesiaLogFileDirPattern = fullfile(testCase.inputLogDir, 'LESIA_*.log');
      Settings.soarLogFileDirPattern  = fullfile(testCase.inputLogDir, 'SOAR_*.log');

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'GENERATE_FROM_LOGS', 'LESIA', 'SOAR' ...
      )

      ExpDt24h6h2h = datetime({'2023-01-01T00:00:00Z'; '2023-02-02T00:00:00Z'}, 'TimeZone', 'UTCLeapSeconds');
      ExpDt7days   = datetime({'2022-12-28T00:00:00Z'; '2023-02-01T00:00:00Z'}, 'TimeZone', 'UTCLeapSeconds');
      testCase.assertEqual(Settings.Gql.Dt24h6h2hArray, ExpDt24h6h2h)
      testCase.assertEqual(Settings.Gql.Dt7daysArray,   ExpDt7days)
    end



  end    % methods(Test)



end
