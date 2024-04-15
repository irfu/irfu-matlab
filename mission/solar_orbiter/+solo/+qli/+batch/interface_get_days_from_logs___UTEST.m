%
% matlab.unittest automatic test code for
% solo.qli.batch.interface.get_days_from_logs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef interface_get_days_from_logs___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and
    % teardown methods which store/read their own data from the testCase
    % object.
    inputLogDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      InputLogFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      testCase.inputLogDir = InputLogFixture.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_log_zero_days(testCase)
      % Empty log files / without matches.
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'LESIA_2024-01-01.log'), {})
      solo.qli.batch.utils.write_file(...
        fullfile(testCase.inputLogDir, 'SOAR_2024-01-01.log'), {})

      LOG_FILE_DIR_PATTERN_DICT          = dictionary();
      LOG_FILE_DIR_PATTERN_DICT('LESIA') = fullfile(testCase.inputLogDir, 'LESIA_*.log');
      LOG_FILE_DIR_PATTERN_DICT('SOAR')  = fullfile(testCase.inputLogDir, 'SOAR_*.log');
      LogFileDirPatternDict              = LOG_FILE_DIR_PATTERN_DICT;

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_logs(...
        LogFileDirPatternDict, {'LESIA', 'SOAR'} ...
        );

      testCase.assertEqual(ActDaysDtArray, solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test_log_two_days(testCase)
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
      LogFileDirPatternDict              = LOG_FILE_DIR_PATTERN_DICT;

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_logs(...
        LogFileDirPatternDict, {'LESIA', 'SOAR'} ...
        );

      testCase.assertEqual(ActDaysDtArray, solo.qli.utils.umddt({'2023-01-01'; '2023-02-02'}))
    end



  end    % methods(Test)



end
