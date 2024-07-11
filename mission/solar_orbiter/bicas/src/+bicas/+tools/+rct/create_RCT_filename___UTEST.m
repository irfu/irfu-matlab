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
      beginDt = datetime('2020-01-01T00:00:00');
      endDt   = datetime('2099-12-31T00:00:00');

      [actDestFilename, actGaCALIBRATION_VERSION] = bicas.tools.rct.create_RCT_filename(beginDt, endDt, 3);
    end



  end    % methods(Test)



end
