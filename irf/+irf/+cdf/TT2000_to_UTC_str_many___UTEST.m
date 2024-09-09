%
% matlab.unittest automatic test code for irf.cdf.TT2000_to_UTC_str_many().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef TT2000_to_UTC_str_many___UTEST < matlab.unittest.TestCase



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
      {zeros(0, 0), 9, cell(0, 0)}, ...
      {zeros(1, 0), 9, cell(1, 0)}, ...
      {zeros(0, 1), 9, cell(0, 1)}, ...
      ...
      ...    % Test decimals, rounding thereof, special case zero decimals.
      {1, 0,  {'2000-01-01T11:58:56Z'}}, ...
      {1, 1,  {'2000-01-01T11:58:55.8Z'}}, ...
      {1, 2,  {'2000-01-01T11:58:55.82Z'}}, ...
      {1, 9,  {'2000-01-01T11:58:55.816000001Z'}}, ...
      {1, 11, {'2000-01-01T11:58:55.81600000100Z'}}, ...
      ...
      ...    % Test array size.
      {
      [0, 1, 2; 3, 4, 5], 9, ...
      {
      '2000-01-01T11:58:55.816000000Z', ...
      '2000-01-01T11:58:55.816000001Z', ...
      '2000-01-01T11:58:55.816000002Z'; ...
      '2000-01-01T11:58:55.816000003Z', ...
      '2000-01-01T11:58:55.816000004Z', ...
      '2000-01-01T11:58:55.816000005Z' ...
      } ...
      } ...
      ...
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
      tt2000Array = int64(TEST_INPUT_OUTPUT{1});
      nDecimals   =       TEST_INPUT_OUTPUT{2};
      expUtcStrCa =       TEST_INPUT_OUTPUT{3};

      actUtcStrCa = irf.cdf.TT2000_to_UTC_str_many(tt2000Array, nDecimals);

      testCase.assertEqual(actUtcStrCa, expUtcStrCa)
    end



  end    % methods(Test)



end
