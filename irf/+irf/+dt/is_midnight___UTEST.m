%
% matlab.unittest automatic test code for irf.dt.is_midnight().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef is_midnight___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test(expIsMidnight, varargin)
        Dt = datetime(varargin{:});

        actIsMidnight = irf.dt.is_midnight(Dt);

        testCase.assertEqual(actIsMidnight, expIsMidnight)
      end

      % Scalars, 1x1
      test(false, [2016,12,31,23,59,59])
      test(false, [2016,12,31,23,59,60], 'TimeZone', 'UTCLeapSeconds')   % Leap second
      test(true,  [2017,01,01,00,00,00])
      test(false, [2023,01,01,00,00,01])

      % 0x0
      test(false(0,0), cell(0,0))
      % 2x1
      test(...
        [false; true], ...
        {'2023-03-03T01:01:01'; '2024-04-04T00:00:00'})
      % 1x2
      test(...
        [false, true], ...
        {'2023-03-03T01:01:01', '2024-04-04T00:00:00'})
    end



  end    % methods(Test)



end
