%
% matlab.unittest automatic test code for bicas.tools.rct.create_RCT_filename().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef create_RCT_filename___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test(testCase)
      [actDestFilename, actGaCALIBRATION_VERSION] = bicas.tools.rct.create_RCT_filename();
    end



  end    % methods(Test)



end
