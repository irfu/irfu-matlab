
% matlab.unittest automatic test code for
% solo.qli.batch.extract_dataset_dates_from_logs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef extract_dataset_dates_from_logs___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects.
    testDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods (TestMethodSetup)



    function close_figures(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.testDir = F.Folder;
    end



  end



  %################
  %################
  % HELPER METHODS
  %################
  %################
  methods



    % Create path to file (existing/non-existing) in test directory.
    function path = fullfile(testCase, filename)
      path = fullfile(testCase.testDir, filename);
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_logs(testCase)
      function call_raise_exception(varargin)
        [~, ~] = solo.qli.batch.extract_dataset_dates_from_logs(...
          '/nonexisting_dir/nonexisting_log.log', {'SOLO_L3_RPW-BIA-DENSITY'});
      end

      testCase.verifyError(...
        @() call_raise_exception(), ...
        ?MException)
    end



    function test___one_matching_log___empty_log_file(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, {})

      [ActDtArray, actLogFilePath] = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

      testCase.assertEqual(ActDtArray, solo.qli.const.EMPTY_DT_ARRAY)
      testCase.assertEqual(actLogFilePath, filePath)
    end



    function test___one_log___no_dataset(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'ABC', ...
        'IRRELEVANT TEXT', ...
        })

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

      testCase.assertEqual(sort(ActDtArray), solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test___one_log___one_dataset_match_simple(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...
        })

      ExpDtArray = irf.dt.um('2024-01-01');

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test___one_log___zero_DSIs_simple(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...
        })

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {});

      testCase.assertEqual(sort(ActDtArray), solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test___one_log___two_identical_datasets(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...    % Same filename.
        })

      ExpDtArray = irf.dt.um('2024-01-01');

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    % Test that can detect both CDAG and non-CDAG filenames.
    function test___one_log___CDAG_and_nonCDAG_datasets(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...         % Same filename, non-CDAG.
        'solo_L3_rpw-bia-density-cdag_20240201_V01.cdf', ...    % Same filename, CDAG.
        })

      ExpDtArray = irf.dt.um({'2024-01-01'; '2024-02-01'});

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test___one_log___one_DSI___multiple_datasets___complex(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        '/data/solo_L3_rpw-bia-density_20240101_V01.cdfQWE', ...              % Matching filename
        '      solo_L3_rpw-bia-density-10-seconds_20240102_V01.cdf   ', ...   % Not matching filename.
        '      solo_L3_rpw-bia-density_20240201_V01.cdf   ', ...              % Matching filename.
        'IRRELEVANT TEXT', ...
        '      solo_L3_rpw-bia-density_20240101_V01.cdf   ', ...              % Same matching filename again.
        })

      ExpDtArray = irf.dt.um({'2024-01-01'; '2024-02-01'});

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test___one_log___multiple_DSIs___complex(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      % IMPLEMENTATION NOTE: Two density, three MAG to test concatenation of
      % different-length cell arrays inside called function.
      solo.qli.batch.utils.write_file(filePath, ...
        {
        'ABCsolo_L3_rpw-bia-density_20240101_V01.cdfQWE', ...   % Matching filename
        '   solo_L3_rpw-bia-density_20240201_V01.cdf   ', ...   % Matching filename.
        'IRRELEVANT TEXT', ...
        'solo_L2_mag-rtn-normal-1-minute_20240101_V03.cdf', ... % Matching filename. Identical date as other match.
        'solo_L2_mag-rtn-normal-1-minute_20240301_V03.cdf', ... % Matching filename
        'solo_L2_mag-rtn-normal-1-minute_20240401_V03.cdf', ... % Matching filename
        })

      ExpDtArray = irf.dt.um('2024-01-01') + calmonths([0;1;2;3]);
      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY', 'SOLO_L2_MAG-RTN-NORMAL-1-MINUTE'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test___one_log___multiple_matches_on_same_row(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath          = testCase.fullfile('processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf  solo_L3_rpw-bia-density_20240201_V01.cdf', ...   % Two matching filenames.
        })

      ExpDtArray = irf.dt.um({'2024-01-01'; '2024-02-01'});

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test___multiple_logs___multiple_DSIs_complex(testCase)
      logFileDirPattern = testCase.fullfile('processing*.log');
      filePath1         = testCase.fullfile('processing2024-01-01T12.00.00.log');
      filePath2         = testCase.fullfile('processing2024-01-01T18.00.00.log');

      % NOTE: Log file which will be ignored.
      solo.qli.batch.utils.write_file(filePath1, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...         % Matching filename
        'solo_L2_mag-rtn-normal-1-minute_20240201_V03.cdf', ... % Matching filename
        })

      solo.qli.batch.utils.write_file(filePath2, ...
        {
        'solo_L3_rpw-bia-density_20250101_V01.cdf', ...         % Matching filename.
        'solo_L2_mag-rtn-normal-1-minute_20250201_V03.cdf', ... % Matching filename
        })

      ExpDtArray = irf.dt.um({'2025-01-01'; '2025-02-01'});

      [ActDtArray, actLogFilePath] = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY', 'SOLO_L2_MAG-RTN-NORMAL-1-MINUTE'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
      testCase.assertEqual(actLogFilePath, filePath2)
    end



    % Test that only recognizes specified log files.
    %
    % NOTE: Test creates two log files, and checks if code can recognize either
    % file separately.
    function test___matching_nonmatching_logs___one_dataset(testCase)
      logFileDirPattern1 = testCase.fullfile('processing*.log');
      logFileDirPattern2 = testCase.fullfile('PROCESSING*.LOG');
      logFilePath1       = testCase.fullfile('processing2024-01-01T12.00.00.log');
      logFilePath2       = testCase.fullfile('PROCESSING2024-01-01T12.00.00.LOG');

      % NOTE: Log file which will be ignored.
      solo.qli.batch.utils.write_file(logFilePath1, ...
        {'solo_L3_rpw-bia-density_20230101_V01.cdf'})         % Matching filename.

      solo.qli.batch.utils.write_file(logFilePath2, ...
        {'solo_L3_rpw-bia-density_20240101_V01.cdf'})         % Matching filename.

      EXP_DT_ARRAY_1 = irf.dt.um({'2023-01-01'});
      EXP_DT_ARRAY_2 = irf.dt.um({'2024-01-01'});

      function test(logFileDirPattern, ExpDtArray, expLogFilePath)
        [ActDtArray, actLogFilePath] = solo.qli.batch.extract_dataset_dates_from_logs(...
          logFileDirPattern, {'SOLO_L3_RPW-BIA-DENSITY'});

        testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
        testCase.assertEqual(actLogFilePath, expLogFilePath)
      end

      test(logFileDirPattern1, EXP_DT_ARRAY_1, logFilePath1)
      test(logFileDirPattern2, EXP_DT_ARRAY_2, logFilePath2)
    end



  end    % methods(Test)



end
