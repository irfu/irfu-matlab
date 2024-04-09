%
% matlab.unittest automatic test code for
% solo.qli.batch.generate_quicklooks_interface().
%
%
% NOTE: Does not test the details of separate commands, since that is the
% purpose of other test code which is closer to the relevant parts of the code:
%   * solo.qli.batch.generate_quicklooks___UTEST
%   * solo.qli.batch.extract_dataset_dates_from_logs___UTEST
%   * solo.qli.batch.fmd___UTEST
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef generate_quicklooks_interface___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Reorg. to have separate test classes for the resp. functions.
  %   PRO: Can check return value (DaysDtArray). Current implementation can only
  %        check for non-crashes.
  %     PRO: Important for LOGS, FMDS.



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
      InputLogFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      QliFixture      = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      Settings             = [];
      Settings.irfLogoPath = solo.qli.testdata.get_test_logo_path();
      Settings.vhtDir      = VhtFixture.Folder;
      Settings.Gql         = solo.qli.batch.GenerateQuicklooksTest();
      Settings.qliDir      = QliFixture.Folder;

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
        'GENERATE', 'TIME_INTERVAL', '2024-01-01', '2024-01-03' ...
        )

      testCase.assertEqual(...
        Settings.Gql.Dt24h6h2hArray, ...
        solo.qli.utils.umdt('2024-01-01') + caldays([0; 1]))
      testCase.assertEqual(...
        Settings.Gql.Dt7daysArray, ...
        solo.qli.utils.umdt('2023-12-27'))
    end



  end    % methods(Test)
  methods(Test)
    % =============
    % LOGS commands
    % =============



    function test_log_zero_days(testCase)
      Settings = testCase.Settings;

      % Empty log files / without matches.
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'LESIA_2024-01-01.log'), {})
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'SOAR_2024-01-01.log'), {})


      LOG_FILE_DIR_PATTERN_DICT          = dictionary();
      LOG_FILE_DIR_PATTERN_DICT('LESIA') = fullfile(testCase.inputLogDir, 'LESIA_*.log');
      LOG_FILE_DIR_PATTERN_DICT('SOAR')  = fullfile(testCase.inputLogDir, 'SOAR_*.log');
      Settings.LogFileDirPatternDict     = LOG_FILE_DIR_PATTERN_DICT;

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'GENERATE', 'LOGS', 'LESIA', 'SOAR' ...
        )

      testCase.assertEqual(Settings.Gql.Dt24h6h2hArray, solo.qli.const.EMPTY_DT_ARRAY)
      testCase.assertEqual(Settings.Gql.Dt7daysArray,   solo.qli.const.EMPTY_DT_ARRAY)
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

      LOG_FILE_DIR_PATTERN_DICT          = dictionary();
      LOG_FILE_DIR_PATTERN_DICT('LESIA') = fullfile(testCase.inputLogDir, 'LESIA_*.log');
      LOG_FILE_DIR_PATTERN_DICT('SOAR')  = fullfile(testCase.inputLogDir, 'SOAR_*.log');
      Settings.LogFileDirPatternDict     = LOG_FILE_DIR_PATTERN_DICT;

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'GENERATE', 'LOGS', 'LESIA', 'SOAR' ...
        )

      ExpDt24h6h2h = solo.qli.utils.umdt({'2023-01-01'; '2023-02-02'});
      ExpDt7days   = solo.qli.utils.umdt({'2022-12-28'; '2023-02-01'});
      testCase.assertEqual(Settings.Gql.Dt24h6h2hArray, ExpDt24h6h2h)
      testCase.assertEqual(Settings.Gql.Dt7daysArray,   ExpDt7days)
    end



  end    % methods(Test)
  methods(Test)
    % =============
    % FMDS commands
    % =============



    function test_empty(testCase)
      Settings = testCase.Settings;
      Settings.datasetDirsCa = cell(0, 1);

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'GENERATE', 'FMDS' ...
        )
    end



    function test_nonempty(testCase)
      % NOTE: Would be more meaningful if could assert the returned dates.
      Settings = testCase.Settings;

      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir1 = F.Folder;
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir2 = F.Folder;

      qliPath      = irf.fs.create_empty_file({Settings.qliDir, '20240101T00_20240102T00.png'});
      datasetPath1 = irf.fs.create_empty_file({dir1, 'solo_L2_swa-pas-eflux_20240101_V02.cdf'});
      datasetPath2 = irf.fs.create_empty_file({dir2, 'solo_L2_mag-rtn-normal_20240101_V02.cdf'});

      Settings.datasetDirsCa = {dir1; dir2};

      solo.qli.batch.generate_quicklooks_interface(...
        Settings, testCase.outputDir, '1', '1', ...
        'GENERATE', 'FMDS' ...
        )
    end



  end    % methods(Test)



end
