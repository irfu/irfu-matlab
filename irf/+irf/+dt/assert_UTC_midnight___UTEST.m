%
% matlab.unittest automatic test code for irf.dt.assert_UTC_midnight().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef assert_UTC_midnight___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test_OK(Dt)
        % Check that does not crash/raise exception.
        irf.dt.assert_UTC_midnight(Dt);
      end

      function test_exc(Dt)
        testCase.verifyError(...
          @() irf.dt.assert_UTC_midnight(Dt), ...
          ?MException)
      end



      EDT          = datetime.empty(0, 0);
      EDT.TimeZone = 'UTCLeapSeconds';

      test_OK( EDT )
      test_OK(datetime('2024-03-06T00:00:00Z', 'TimeZone', 'UTCLeapSeconds'))
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



  end    % methods(Test)



end
