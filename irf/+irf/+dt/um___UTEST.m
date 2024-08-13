%
% matlab.unittest automatic test code for irf.dt.um().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef um___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test(strCa, expStrCa)
        assert(iscell(expStrCa))

        ActDt = irf.dt.um(strCa);

        ExpDt = datetime(expStrCa, 'TimeZone', 'UTCLeapSeconds');

        testCase.assertEqual(ActDt, ExpDt)
        testCase.assertEqual(size(ActDt), size(expStrCa))
      end



      test(...
        '2024-02-03', ...
        {'2024-02-03T00:00:00Z'})

      % Empty arrays.
      test(...
        cell(0, 0), ...
        cell(0, 0))
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
        @() irf.dt.um('2024-02-03T00:00:00Z'), ...
        ?MException)
      testCase.assertError(...
        @() irf.dt.um(datetime('2024-02-03')), ...
        ?MException)

    end



  end    % methods(Test)
end
