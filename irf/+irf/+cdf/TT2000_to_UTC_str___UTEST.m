%
% matlab.unittest automatic test code for irf.cdf.TT2000_to_UTC_str().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef TT2000_to_UTC_str___UTEST < matlab.unittest.TestCase



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  properties(TestParameter)
    % Technically, additional properties of testCase objects with cell array
    % default values. Test methods with arguments with the same name will be
    % called once for every element in the cell arrays.

    TEST_INPUT_OUTPUT = {
      {     0, 9, '2000-01-01T11:58:55.816000000Z'}, ...
      {456789, 9, '2000-01-01T11:58:55.816456789Z'} ...
      {456789, 1, '2000-01-01T11:58:55.8Z'} ...
      {456789, 0, '2000-01-01T11:58:56Z'} ...
      }
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase, TEST_INPUT_OUTPUT)
      % NOTE: Converting argument to int64. This should maybe required in the
      % future.
      tt2000    = int64(TEST_INPUT_OUTPUT{1});
      nDecimals = TEST_INPUT_OUTPUT{2};
      expUtcStr = TEST_INPUT_OUTPUT{3};

      actUtcStr = irf.cdf.TT2000_to_UTC_str(tt2000, nDecimals);

      testCase.assertEqual(actUtcStr, expUtcStr)
    end



  end    % methods(Test)



end
