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



    function test_umdt(testCase)
      function test(strCa, expStrCa)
        assert(iscell(expStrCa))
        ActDt = solo.qli.utils.umdt(strCa);

        ExpDt = datetime(expStrCa, 'TimeZone', 'UTCLeapSeconds');

        testCase.assertEqual(ActDt, ExpDt)
        testCase.assertEqual(size(ActDt), size(expStrCa))
      end

      test(...
        '2024-02-03', ...
        {'2024-02-03T00:00:00Z'})
      test(...
        cell(0, 1), ...
        cell(0, 1))
      test(...
        {'2024-02-03'}, ...
        {'2024-02-03T00:00:00Z'})
      test(...
        {'2024-02-03';           '2024-03-04'}, ...
        {'2024-02-03T00:00:00Z'; '2024-03-04T00:00:00Z'})

      testCase.assertError(...
        @() solo.qli.utils.umdt('2024-02-03T00:00:00Z'), ...
        ?MException)
      testCase.assertError(...
        @() solo.qli.utils.umdt(datetime('2024-02-03')), ...
        ?MException)
    end



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

      test_OK( EDT )
      test_OK( datetime('2024-03-06T00:00:00Z', 'TimeZone', 'UTCLeapSeconds') )
      test_OK( [
        datetime('2024-03-06T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-07T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'); ...
        datetime('2024-03-08T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-09T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'); ...
        ])

      test_exc(datetime('2024-03-06T00:00:00.001Z', 'TimeZone', 'UTCLeapSeconds'))
      test_exc(datetime('2024-03-06T23:59:59.999Z', 'TimeZone', 'UTCLeapSeconds'))
      test_exc(datetime('2024-03-06 00:00:00'))   % Not UTC.
      test_exc( [
        datetime('2024-03-06T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-07T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'); ...
        datetime('2024-03-08T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'), ...
        datetime('2024-03-09T00:00:01Z', 'TimeZone', 'UTCLeapSeconds'); ...
        ])
    end



    function test_derive_weeks(testCase)

      function test(DayBeginStrCa, firstDayOfWeek, ExpWeekBeginStrCa)
        DayBeginDtArray     = solo.qli.utils.umdt(DayBeginStrCa);
        ExpWeekBeginDtArray = solo.qli.utils.umdt(ExpWeekBeginStrCa);

        ActWeekBeginDtArray = solo.qli.utils.derive_weeks(DayBeginDtArray, firstDayOfWeek);

        testCase.assertEqual(ActWeekBeginDtArray, ExpWeekBeginDtArray)
      end

      SUNDAY   = 1;
      SATURDAY = 7;

      for dayOfWeek = 1:7
        test(...
          cell(0, 1), dayOfWeek, ...
          cell(0, 1))
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
        '2024-03-03', SUNDAY, ...
        '2024-03-03')
      test(...
        '2024-03-09', SUNDAY, ...
        '2024-03-03')
      test(...
        '2024-03-10', SUNDAY, ...
        '2024-03-10')

      test(...
        '2024-03-02', SATURDAY, ...
        '2024-03-02')
      test(...
        '2024-03-08', SATURDAY, ...
        '2024-03-02')
      test(...
        '2024-03-09', SATURDAY, ...
        '2024-03-09')

      % Multiple input days for same week
      % Non-incrementing order.
      % 1 week.
      test(...
        {
        '2024-03-10';
        '2024-03-09'; ...
        '2024-03-15'; ...
        }, ...
        SATURDAY, ...
        '2024-03-09')

      % Duplicate timestamps.
      % Non-incrementing order.
      % 2 weeks.
      test(...
        {
        '2024-03-08'; ...
        '2024-03-10'; ...
        '2024-03-10'; ...
        '2024-03-09'; ...
        '2024-03-15'; ...
        }, ...
        SATURDAY, ...
        {
        '2024-03-02'; ...
        '2024-03-09'; ...
        } ...
        )
    end



    % Merely test that function (1) does not crash, and (2) returns string.
    %
    % NOTE: Function returns time-dependent string.
    function test_get_data_source_info_string(testCase)

      actOutput = solo.qli.utils.get_data_source_info_string();

      testCase.verifyInstanceOf(actOutput, 'char')
      testCase.verifyTrue(isrow(actOutput))
    end



    function test_get_context_info_strings(testCase)

      % Arbitrary number of output variables.
      function test(...
          soloPosTSeries, earthPosTSeries, Tint, ...
          expSoloStr, expEarthStr)

        [actSoloStr, actEarthStr] = solo.qli.utils.get_context_info_strings(soloPosTSeries, earthPosTSeries, Tint);
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



    function test_create_quicklook_filename(testCase)
      function test(Tint, expFilename)
        actFilename = solo.qli.utils.create_quicklook_filename(Tint);
        testCase.assertEqual(actFilename, expFilename)
      end

      %===================================================================

      Tint = EpochTT( ...
        [ ...
        '2024-01-10T02:09:04.900000009Z'; ...
        '2024-01-11T04:04:09.900000004Z'; ...
        ] ...
        );
      test(Tint, '20240110T02_20240111T04.png')
    end



  end    % methods(Test)



end
