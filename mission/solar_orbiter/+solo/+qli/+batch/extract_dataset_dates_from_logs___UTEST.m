
% matlab.unittest automatic test code for
% solo.qli.batch.extract_dataset_dates_from_logs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef extract_dataset_dates_from_logs___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Use unit test setup for temporary directory.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_logs(testCase)
      function call_raise_exception(varargin)
        [~, ~] = solo.qli.batch.extract_dataset_dates_from_logs(...
          '/nonexisting_dir/nonexisting_log.log', {'solo_L3_rpw-bia-density'});
      end

      testCase.verifyError(...
        @() call_raise_exception(), ...
        ?MException)
    end



    function test_one_log_empty(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, {})

      ExpDtArray = solo.qli.const.DT_EMPTY_ARRAY;

      [ActDtArray, actLogFilePath] = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density'});

      testCase.assertEqual(ActDtArray, ExpDtArray)
      testCase.assertEqual(actLogFilePath, filePath)
    end



    function test_one_log_no_match(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'ABC', ...
        'IRRELEVANT TEXT', ...
        })

      ExpDtArray = solo.qli.const.DT_EMPTY_ARRAY;

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test_one_log_one_match_simple(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...
        })

      ExpDtArray = [
        datetime(2024, 1, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        ];

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test_one_log_zero_DSIs_simple(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...
        })

      ExpDtArray = solo.qli.const.DT_EMPTY_ARRAY;

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {});

      testCase.assertEqual(sort(ActDtArray), ExpDtArray)
    end



    function test_one_log_two_identical_matches(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...    % Same filename.
        })

      ExpDtArray = [
        datetime(2024, 1, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        ];

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    % Test that can detect both CDAG and non-CDAG filenames.
    function test_one_log_CDAG_nonCDAG(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf', ...         % Same filename, non-CDAG.
        'solo_L3_rpw-bia-density-cdag_20240201_V01.cdf', ...    % Same filename, CDAG.
        })

      ExpDtArray = [
        datetime(2024, 1, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        datetime(2024, 2, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        ];

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test_one_log_one_DSI_complex(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        '/data/solo_L3_rpw-bia-density_20240101_V01.cdfQWE', ...   % Matching filename
        '      solo_L3_rpw-bia-density-10-seconds_20240102_V01.cdf   ', ...   % Potentially mistakenly matching filename.
        '      solo_L3_rpw-bia-density_20240201_V01.cdf   ', ...   % Matching filename.
        'IRRELEVANT TEXT', ...
        '      solo_L3_rpw-bia-density_20240101_V01.cdf   ', ...   % Same matching filename again.
        })

      ExpDtArray = [
        datetime(2024, 1, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        datetime(2024, 2, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        ];

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test_one_log_multiple_DSIs_complex(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

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

      ExpDtArray = [
        datetime(2024, 1, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        datetime(2024, 2, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        datetime(2024, 3, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        datetime(2024, 4, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        ];

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density', 'solo_L2_mag-rtn-normal-1-minute'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test_one_log_multiple_matches_on_same_row(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath          = fullfile(F.Folder, 'processing2024-01-01.log');

      solo.qli.batch.utils.write_file(filePath, ...
        {
        'solo_L3_rpw-bia-density_20240101_V01.cdf  solo_L3_rpw-bia-density_20240201_V01.cdf', ...   % Two matching filenames.
        })

      ExpDtArray = [
        datetime(2024, 1, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        datetime(2024, 2, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        ];

      ActDtArray = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
    end



    function test_multiple_logs_multiple_DSIs_complex(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      logFileDirPattern = fullfile(F.Folder, 'processing*.log');
      filePath1         = fullfile(F.Folder, 'processing2024-01-01T12.00.00.log');
      filePath2         = fullfile(F.Folder, 'processing2024-01-01T18.00.00.log');

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

      ExpDtArray = [
        datetime(2025, 1, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        datetime(2025, 2, 1, 'TimeZone', 'UTCLeapSeconds'); ...
        ];

      [ActDtArray, actLogFilePath] = solo.qli.batch.extract_dataset_dates_from_logs(...
        logFileDirPattern, {'solo_L3_rpw-bia-density', 'solo_L2_mag-rtn-normal-1-minute'});

      testCase.assertEqual(sort(ActDtArray), sort(ExpDtArray))
      testCase.assertEqual(actLogFilePath, filePath2)
    end



  end    % methods(Test)



end
