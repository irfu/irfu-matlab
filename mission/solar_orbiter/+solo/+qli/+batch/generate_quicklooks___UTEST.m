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
      solo.qli.const.EMPTY_DT_ARRAY, ...
      solo.qli.utils.umddt('2024-01-01') ...
      }
    EXC_7DAYS = {...
      solo.qli.const.EMPTY_DT_ARRAY, ...
      solo.qli.utils.umddt('2023-12-27') ...
      }
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods (TestMethodSetup)



    function create_output_directories(testCase)
      VhtFixture    = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      OutputFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      testCase.vhtDir    = VhtFixture.Folder;
      testCase.outputDir = OutputFixture.Folder;
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
        false, false, solo.qli.const.EMPTY_DT_ARRAY, Gql)

      % Logo & quicklooks enabled.
      solo.qli.batch.generate_quicklooks(...
        logoPath, testCase.vhtDir, testCase.outputDir, ...
        true, true, solo.qli.const.EMPTY_DT_ARRAY, Gql)
    end



    function test_one_day_24h6h2h(testCase)
      UmdDtArray = solo.qli.utils.umddt('2024-01-01');

      Gql = solo.qli.batch.GenerateQuicklooksTest();

      solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        true, false, UmdDtArray, Gql)

      testCase.assertEqual(Gql.UmdDt24h6h2hArray, UmdDtArray)
      testCase.assertEqual(Gql.UmdDt7daysArray,   solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test_one_day_week_7days(testCase)
      UmdDtArray         = solo.qli.utils.umddt('2024-01-01');
      ExpUmdDt7daysArray = solo.qli.utils.umddt('2023-12-27');
      Gql = solo.qli.batch.GenerateQuicklooksTest();

      solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        false, true, UmdDtArray, Gql)

      testCase.assertEqual(Gql.UmdDt24h6h2hArray, solo.qli.const.EMPTY_DT_ARRAY)
      testCase.assertEqual(Gql.UmdDt7daysArray,   ExpUmdDt7daysArray)
    end



    % Generate two days, but raise exception for the first one.
    function test_two_days_24h6h2h_exception(testCase)
      UmdDtArray = solo.qli.utils.umddt('2024-01-01') + caldays([0; 1]);

      Gql = solo.qli.batch.GenerateQuicklooksTest(...
        UmdDtArray(1), ...
        solo.qli.const.EMPTY_DT_ARRAY);

      testCase.assertError(...
        @() solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        true, false, UmdDtArray, Gql), ...
        ?MException)

      testCase.assertEqual(Gql.UmdDt24h6h2hArray, UmdDtArray)
      testCase.assertEqual(Gql.UmdDt7daysArray,   solo.qli.const.EMPTY_DT_ARRAY)
    end



    % Generate two weeks, but raise exception for the first one.
    function test_two_weeks_7days_exception(testCase)
      UmdDtArray         = solo.qli.utils.umddt('2024-01-01') + caldays( [0; 2]);  % Mon, Wed
      ExpUmdDt7daysArray = solo.qli.utils.umddt('2023-12-27') + calweeks([0; 1]);
      Gql = solo.qli.batch.GenerateQuicklooksTest(...
        solo.qli.const.EMPTY_DT_ARRAY, ...
        ExpUmdDt7daysArray(1));

      testCase.assertError(...
        @() solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        false, true, UmdDtArray, Gql), ...
        ?MException)

      testCase.assertEqual(Gql.UmdDt24h6h2hArray, solo.qli.const.EMPTY_DT_ARRAY)
      testCase.assertEqual(Gql.UmdDt7daysArray,   ExpUmdDt7daysArray)
    end



    function test_complex(testCase)
      UmdDtArray         = solo.qli.utils.umddt('2024-01-01') + caldays([0; 1; 10; 11]); % Mon-Tue, Wed-Thu
      ExpUmdDt7daysArray = solo.qli.utils.umddt('2023-12-27') + calweeks([0; 2]);        % Skip middle week.

      Gql = solo.qli.batch.GenerateQuicklooksTest(...
        UmdDtArray(2), ...
        ExpUmdDt7daysArray(1));

      testCase.assertError(...
        @() solo.qli.batch.generate_quicklooks(...
        [], testCase.vhtDir, testCase.outputDir, ...
        true, true, UmdDtArray, Gql), ...
        ?MException)

      testCase.assertEqual(Gql.UmdDt24h6h2hArray, UmdDtArray)
      testCase.assertEqual(Gql.UmdDt7daysArray,   ExpUmdDt7daysArray)
    end



  end    % methods(Test)



end
