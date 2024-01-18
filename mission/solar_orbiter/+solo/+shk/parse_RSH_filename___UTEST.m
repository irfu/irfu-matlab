%
% matlab.unittest automatic test code for
% solo.shk.parse_RSH_filename().
%
%
% Author: Erik P G Johansson, IRF Uppsala, Sweden
% First created 2021-09-08
%
classdef parse_RSH_filename___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)

      function test(filename, expOutputsCa)
        % Pre-allocate correct size for later assignment via function
        actOutputs = cell(size(expOutputsCa));

        [actOutputs{:}] = solo.shk.parse_RSH_filename(filename);
        testCase.verifyEqual(actOutputs, expOutputsCa)
      end
      %===================================================================

      test('solo_HK_platform_20210901_V2.xml',   {false, NaT, NaN})
      test('solo_HK_platform_20210901_V02.xml',  {true, datetime('2021-09-01'), 2})
      test('solo_HK_platform_20200219_V023.xml', {true, datetime('2020-02-19'), 23})
    end



  end    % methods(Test)



end
