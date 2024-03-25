%
% UNFINISHED
%
% matlab.unittest automatic test code for solo.qli.batch.generate_quicklooks().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef generate_quicklooks___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    vhtDir
    outputDir
  end



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  properties(TestParameter)
    EXC_24H6H2H = {...
      solo.qli.const.DT_EMPTY_ARRAY, ...
      datetime('2024-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') ...
      }
    EXC_7DAYS = {...
      solo.qli.const.DT_EMPTY_ARRAY, ...
      datetime('2023-12-27T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') ...
      }
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods (TestMethodSetup)



    function create_output_directories(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.vhtDir = F.Folder;

      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.outputDir = F.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_dates(testCase)
      Gql      = solo.qli.batch.GenerateQuicklooksTest();
      logoPath = solo.qli.testdata.get_test_logo_path();

      % Everything disabled.
      solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        false, false, solo.qli.const.DT_EMPTY_ARRAY, Gql)

      % Logo, quicklooks enabled.
      solo.qli.batch.generate_quicklooks(...
        logoPath, testCase.vhtDir, testCase.outputDir, ...
        true, true, solo.qli.const.DT_EMPTY_ARRAY, Gql)
    end



    function test_one_day_24h6h2h(testCase)
      DaysDtArray = datetime('2024-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');

      Gql = solo.qli.batch.GenerateQuicklooksTest();

      solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        true, false, DaysDtArray, Gql)

      testCase.assertEqual(Gql.Dt24h6h2hArray, DaysDtArray)
      testCase.assertEqual(Gql.Dt7daysArray,   solo.qli.const.DT_EMPTY_ARRAY)
    end



    function test_one_day_week_7days(testCase)
      DaysDtArray     = datetime('2024-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      ExpDt7daysArray = datetime('2023-12-27T00:00:00Z', 'TimeZone', 'UTCLeapSeconds');
      Gql = solo.qli.batch.GenerateQuicklooksTest();

      solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        false, true, DaysDtArray, Gql)

      testCase.assertEqual(Gql.Dt24h6h2hArray, solo.qli.const.DT_EMPTY_ARRAY)
      testCase.assertEqual(Gql.Dt7daysArray,   ExpDt7daysArray)
    end



    % Generate two days, but raise exception for the first one.
    function test_two_days_24h6h2h_exception(testCase)
      DaysDtArray = datetime('2024-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') + caldays([0; 1]);

      Gql = solo.qli.batch.GenerateQuicklooksTest(...
        DaysDtArray(1), ...
        solo.qli.const.DT_EMPTY_ARRAY);

      testCase.assertError(...
        @() solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        true, false, DaysDtArray, Gql), ...
        ?MException)

      testCase.assertEqual(Gql.Dt24h6h2hArray, DaysDtArray)
      testCase.assertEqual(Gql.Dt7daysArray,   solo.qli.const.DT_EMPTY_ARRAY)
    end



    % Generate two weeks, but raise exception for the first one.
    function test_two_weeks_7days_exception(testCase)
      DaysDtArray     = datetime('2024-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') + caldays([0; 2]);  % Mon, Wed
      ExpDt7daysArray = datetime('2023-12-27T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') + calweeks([0; 1]);
      Gql = solo.qli.batch.GenerateQuicklooksTest(...
        solo.qli.const.DT_EMPTY_ARRAY, ...
        ExpDt7daysArray(1));

      testCase.assertError(...
        @() solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        false, true, DaysDtArray, Gql), ...
        ?MException)

      testCase.assertEqual(Gql.Dt24h6h2hArray, solo.qli.const.DT_EMPTY_ARRAY)
      testCase.assertEqual(Gql.Dt7daysArray,   ExpDt7daysArray)
    end



    function test_complex(testCase)
      DaysDtArray     = datetime('2024-01-01T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') + caldays([0; 1; 10; 11]); % Mon-Tue, Wed-Thu
      ExpDt7daysArray = datetime('2023-12-27T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') + calweeks([0; 2]);        % Skip middle week.

      Gql = solo.qli.batch.GenerateQuicklooksTest(...
        DaysDtArray(2), ...
        ExpDt7daysArray(1));

      testCase.assertError(...
        @() solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        true, true, DaysDtArray, Gql), ...
        ?MException)

      testCase.assertEqual(Gql.Dt24h6h2hArray, DaysDtArray)
      testCase.assertEqual(Gql.Dt7daysArray,   ExpDt7daysArray)
    end



  end    % methods(Test)



end
