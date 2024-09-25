%
% matlab.unittest automatic test code for inverse functions
% solo.adm.dsfn.time_interval_str.create() and
% solo.adm.dsfn.time_interval_str.parse().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef time_interval_str___UTEST < matlab.unittest.TestCase



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  properties(TestParameter)
    % Technically, additional properties of testCase objects with cell array
    % default values. Test methods with arguments with the same name will be
    % called once for every element in the cell arrays.
    DATE_SEPARATOR = {'', '-'};
    TIME_SEPARATOR = {'', '.', ':'};
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_DAY(testCase)
      testCase.test_OK(        [2021, 2, 3, 0, 0, 0], [2021, 2, 4, 0, 0, 0], 'DAY', '20210203')

      testCase.test_create_exc([2021, 2, 3, 0, 0, 1], [2021, 2, 4, 0, 0, 0], 'DAY')
      testCase.test_create_exc([2021, 2, 3, 0, 0, 0], [2021, 2, 4, 0, 0, 1], 'DAY')
      testCase.test_create_exc([2021, 2, 3, 0, 0, 1], [2021, 2, 4, 0, 0, 1], 'DAY')

      testCase.test_parse_fail('210203')
      testCase.test_parse_fail('2021203')
    end



    function test_DAY_exc_separators(testCase, DATE_SEPARATOR)
      if isempty(DATE_SEPARATOR)
        return
      end

      s = sprintf('2021%s02%s03', DATE_SEPARATOR, DATE_SEPARATOR);
      testCase.test_parse_fail(s)
    end



    function test_DAY_TO_DAY(testCase)
      testCase.test_OK(        [2021, 2, 3, 0, 0, 0], [ 2021,  3,  4, 0, 0, 0], 'DAY_TO_DAY', '20210203-20210303')
      testCase.test_OK(        [2020, 1, 1, 0, 0, 0], [10000, 01, 01, 0, 0, 0], 'DAY_TO_DAY', '20200101-99991231')

      testCase.test_create_exc([2021, 2, 3, 0, 0, 1], [2021, 3, 4, 0, 0, 0], 'DAY_TO_DAY')
      testCase.test_create_exc([2021, 2, 3, 0, 0, 0], [2021, 3, 4, 0, 0, 1], 'DAY_TO_DAY')
      testCase.test_create_exc([2021, 2, 3, 0, 0, 1], [2021, 3, 4, 0, 0, 1], 'DAY_TO_DAY')
    end



    function test_SECOND_TO_SECOND(testCase)
      testCase.test_OK([2020,12,31, 23,58,59], [2021, 1, 2, 3, 4, 5], 'SECOND_TO_SECOND', '20201231T235859-20210102T030405')
      testCase.test_OK([2020, 1, 1,  0, 0, 0], [9999,12,31,23,59,59], 'SECOND_TO_SECOND', '20200101T000000-99991231T235959')
    end



    % Test timestamps using illegal timestamp field separators.
    function test_SECOND_TO_SECOND_separator_parse_fail(testCase, DATE_SEPARATOR, TIME_SEPARATOR)
      if isempty(DATE_SEPARATOR) && isempty(TIME_SEPARATOR)
        % Skip, since there are no field separators (which is legal/required).
        return
      end

      s = sprintf('2021%s02%s03T23%s59%s59', ...
        DATE_SEPARATOR, DATE_SEPARATOR, ...
        TIME_SEPARATOR, TIME_SEPARATOR);
      s2 = [s, '-', s];
      testCase.test_parse_fail(s2)
    end



    function test_SECOND_TO_SECOND_length_parse_fail(testCase)
      testCase.test_parse_fail('2020121T235859-20210102T030405')
      testCase.test_parse_fail('20201231T23585-20210102T030405')
      testCase.test_parse_fail('20201231T235859-2021102T030405')
      testCase.test_parse_fail('20201231T235859-20210102T03040')
    end



    % Same time interval on different formats.
    % Over newyear.
    function test_same_timestamps_multiple_formats_newyear(testCase)
      testCase.test_OK([2020,12,31, 0,0,0], [2021,1,1,0,0,0], 'DAY',              '20201231')
      testCase.test_OK([2020,12,31, 0,0,0], [2021,1,1,0,0,0], 'DAY_TO_DAY',       '20201231-20201231')
      testCase.test_OK([2020,12,31, 0,0,0], [2021,1,1,0,0,0], 'SECOND_TO_SECOND', '20201231T000000-20210101T000000')
    end



    % Leap second as timestamp.
    function test_leap_second(testCase)
      testCase.test_OK([2016,12,31, 23,59,60], [2017,01,01,00,00,00], 'SECOND_TO_SECOND', '20161231T235960-20170101T000000')
      testCase.test_OK([2016,12,31, 23,59,59], [2016,12,31,23,59,60], 'SECOND_TO_SECOND', '20161231T235959-20161231T235960')
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Test that conversions in both directions are consistent.
    function test_OK(testCase, dateVec1, dateVec2, timeIntervalFormat, timeIntervalStr)
      Dt1 = datetime(dateVec1, 'TimeZone', 'UTCLeapSeconds');
      Dt2 = datetime(dateVec2, 'TimeZone', 'UTCLeapSeconds');

      [actDt1, actDt2, actTimeIntervalFormat, actBSuccess] = ...
        solo.adm.dsfn.time_interval_str.parse(timeIntervalStr);
      actTimeIntervalStr                      = ...
        solo.adm.dsfn.time_interval_str.create(Dt1, Dt2, timeIntervalFormat);

      testCase.assertEqual(actDt1, Dt1)
      testCase.assertEqual(actDt2, Dt2)
      testCase.assertEqual(actTimeIntervalFormat, timeIntervalFormat)
      testCase.assertTrue( actBSuccess)
      testCase.assertEqual(actTimeIntervalStr,    timeIntervalStr)
    end



    function test_create_exc(testCase, dateVec1, dateVec2, timeIntervalFormat)
      Dt1 = datetime(dateVec1, 'TimeZone', 'UTCLeapSeconds');
      Dt2 = datetime(dateVec2, 'TimeZone', 'UTCLeapSeconds');

      testCase.assertError(...
        @() solo.adm.dsfn.time_interval_str.create(Dt1, Dt2, timeIntervalFormat), ...
        ?MException)
    end



    function test_parse_fail(testCase, timeIntervalStr)
      [actDt1, actDt2, actTimeIntervalFormat, actBSuccess] = ...
          solo.adm.dsfn.time_interval_str.parse(timeIntervalStr);

      testCase.assertEqual(actDt1, irf.dt.UTC('NaT'))
      testCase.assertEqual(actDt2, irf.dt.UTC('NaT'))
      testCase.assertEqual(actTimeIntervalFormat, [])
      testCase.assertFalse(actBSuccess)
    end



  end    % methods(Access=private)



end
