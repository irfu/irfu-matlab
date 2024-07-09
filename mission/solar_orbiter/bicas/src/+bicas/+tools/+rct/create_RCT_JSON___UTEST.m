%
% matlab.unittest automatic test code for bicas.tools.rct.create_RCT_JSON().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef create_RCT_JSON___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and
    % teardown methods which store/read their own data from the testCase
    % object.
    testDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      Fixture          = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.testDir = Fixture.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test(testCase)
      actRctJsonPath = bicas.tools.rct.create_RCT_JSON(testCase.testDir, 'biasRctFilename.cdf');
      irf.assert.file_exists(actRctJsonPath)
    end



  end    % methods(Test)



end
