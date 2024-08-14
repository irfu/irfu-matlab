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
      Dt1 = irf.dt.UTC('2020-01-01T00:00:00Z');
      Dt2   = irf.dt.UTC('2099-12-31T00:00:00Z');

      actDestFilename = bicas.tools.rct.create_RCT_filename(Dt1, Dt2, 3);

      irf.assert.castring(actDestFilename)
      testCase.assertTrue(numel(actDestFilename) > 4)
    end



  end    % methods(Test)



end
