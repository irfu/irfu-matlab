%
% matlab.unittest automatic test code for solo.qli.utils.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_assert_UTC_midnight_datetime(testCase)

      function test_OK(Dt)
        % Check that does not crash/raise exception.
        solo.qli.utils.assert_UTC_midnight_datetime(Dt);
      end

      function test_exc(Dt)
        testCase.verifyError(...
          @() solo.qli.utils.assert_UTC_midnight_datetime(Dt), ...
          ?MException)
      end

      EDT          = datetime.empty(0, 0);
      EDT.TimeZone = 'UTCLeapSeconds';

      test_OK( EDT)
      test_OK( datetime('2024-03-06T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'))
      test_OK( [
        datetime('2024-03-06T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-07T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'); ...
        datetime('2024-03-08T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-09T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'); ...
        ])

      test_exc(datetime('2024-03-06T00:00:00.001Z', 'TimeZone', 'UTCLeapSeconds'))
      test_exc(datetime('2024-03-06T23:59:59.999Z', 'TimeZone', 'UTCLeapSeconds'))
      test_exc(datetime('2024-03-06 00:00:00'))
      test_exc( [
        datetime('2024-03-06T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-07T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'); ...
        datetime('2024-03-08T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-09T00:00:01Z', 'TimeZone', 'UTCLeapSeconds'); ...
        ])
    end



    function test_scalar_datetime_to_EpochTT(testCase)

      function test(utcStr)
        Dt          = datetime(utcStr, 'TimeZone', 'UTCLeapSeconds');
        ExpEpochTT  = EpochTT(utcStr);

        ActEpochTT = solo.qli.utils.scalar_datetime_to_EpochTT(Dt);

        testCase.assertEqual(ActEpochTT, ExpEpochTT)
      end

      test('2024-03-05T01:02:03.004Z')
      test('2024-03-05T01:02:03.123456789Z')

      % On and around leap second.
      test('2016-12-31T23:59:58.000Z')
      test('2016-12-31T23:59:59.000Z')
      test('2016-12-31T23:59:60.000Z')
      test('2016-12-31T23:59:60.500Z')
      test('2017-01-01T00:00:00.000Z')
      test('2017-01-01T00:00:01.000Z')
    end



    function test_derive_weeks(testCase)

      function test(DayBeginDtArray, firstDayOfWeek, ExpWeekBeginDtArray)
        DayBeginDtArray.TimeZone     = 'UTCLeapSeconds';
        ExpWeekBeginDtArray.TimeZone = 'UTCLeapSeconds';

        ActWeekBeginDtArray = solo.qli.utils.derive_weeks(DayBeginDtArray, firstDayOfWeek);

        testCase.assertEqual(ActWeekBeginDtArray, ExpWeekBeginDtArray)
      end

      SUNDAY   = 1;
      SATURDAY = 7;
      EDT          = datetime.empty(0, 1);
      EDT.TimeZone = 'UTCLeapSeconds';

      for dayOfWeek = 1:7
        test(...
          EDT, dayOfWeek, ...
          EDT)
      end

      % 2024-03-01: Fri
      %         02: Sat
      %         03: Sun
      %         04: Mon
      %         05: Tue
      %         06: Wed
      %         07: Thu
      %         08: Fri
      %         09: Sat
      test(...
        datetime('2024-03-03'), SUNDAY, ...
        datetime('2024-03-03'))
      test(...
        datetime('2024-03-09'), SUNDAY, ...
        datetime('2024-03-03'))
      test(...
        datetime('2024-03-10'), SUNDAY, ...
        datetime('2024-03-10'))

      test(...
        datetime('2024-03-02'), SATURDAY, ...
        datetime('2024-03-02'))
      test(...
        datetime('2024-03-08'), SATURDAY, ...
        datetime('2024-03-02'))
      test(...
        datetime('2024-03-09'), SATURDAY, ...
        datetime('2024-03-09'))

      % Multiple input days for same week
      % Non-incrementing order.
      % 1 week.
      test(...
        [
        datetime('2024-03-10');
        datetime('2024-03-09'); ...
        datetime('2024-03-15'); ...
        ], ...
        SATURDAY, ...
        datetime('2024-03-09'))

      % Duplicate timestamps.
      % Non-incrementing order.
      % 2 weeks.
      test(...
        [
        datetime('2024-03-08'); ...
        datetime('2024-03-10'); ...
        datetime('2024-03-10'); ...
        datetime('2024-03-09'); ...
        datetime('2024-03-15'); ...
        ], ...
        SATURDAY, ...
        [
        datetime('2024-03-02'); ...
        datetime('2024-03-09'); ...
        ] ...
        )
    end



    % Merely test that function (1) does not crash, and (2) returns string.
    %
    % NOTE: Function returns time-dependent string.
    function test_generate_data_source_info_string(testCase)

      actOutput = solo.qli.utils.generate_data_source_info_string();

      testCase.verifyInstanceOf(actOutput, 'char')
      testCase.verifyTrue(isrow(actOutput))
    end



    function test_context_info_strings(testCase)

      % Arbitrary number output variables.
      function test(...
          soloPosTSeries, earthPosTSeries, Tint, ...
          expSoloStr, expEarthStr)

        [actSoloStr, actEarthStr] = solo.qli.utils.context_info_strings(soloPosTSeries, earthPosTSeries, Tint);
        testCase.verifyEqual(actSoloStr,  expSoloStr)
        testCase.verifyEqual(actEarthStr, expEarthStr)
      end

      %===================================================================

      Units = irf_units;
      AU_KM = Units.AU / Units.km;   % Astronomical unit [km]

      ETT = EpochTT( ...
        [ ...
        '2024-01-10T00:00:00.000000000Z'; ...
        '2024-01-11T00:00:00.000000000Z'; ...
        '2024-01-12T00:00:00.000000000Z' ...
        ] ...
        );
      POSITION_TS = TSeries( ...
        ETT, ...
        [ ...
        1*AU_KM, 3, 4; ...
        2*AU_KM, 5, 6; ...
        3*AU_KM, 7, 8; ...
        ] ...
        );

      % In-range time interval.
      TI_1 = EpochTT(['2024-01-09T00:00:00.000000000Z'; '2024-03-13T00:00:00.000000000Z']);
      test(...
        POSITION_TS, POSITION_TS, TI_1, ...
        'SolO:  1.00 AU,  EcLat 229\circ,  EcLon 172\circ', ...
        'Earth:  EcLon 172\circ')

      % Out-of-range time interval.
      TI_2 = EpochTT(['2024-01-01T00:00:00.000000000Z'; '2024-01-02T00:00:00.000000000Z']);
      test(...
        POSITION_TS, POSITION_TS, TI_2, ...
        '', ...
        '')
    end



    function test_get_plot_filename(testCase)

      function test(inputsCa, expOutput)
        actOutput = solo.qli.utils.get_plot_filename(inputsCa{:});
        testCase.verifyEqual(actOutput, expOutput)
      end

      %===================================================================

      Tint = EpochTT( ...
        [ ...
        '2024-01-10T02:09:04.900000009Z'; ...
        '2024-01-11T04:04:09.900000004Z'; ...
        ] ...
        );
      test({Tint}, '20240110T02_20240111T04.png')
    end


  end    % methods(Test)



end
