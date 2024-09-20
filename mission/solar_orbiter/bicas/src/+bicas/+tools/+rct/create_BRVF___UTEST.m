%
% matlab.unittest automatic test code for bicas.tools.rct.create_BRVF().
%
% NOTE: bicas.proc.L1L2.cal.rct.findread___UTEST also indirectly tests
% bicas.tools.rct.create_BRVF()
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef create_BRVF___UTEST < matlab.unittest.TestCase



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
      BEGIN_DT = datetime('2020-01-02T03:04:05Z', 'TimeZone', 'UTCLeapSeconds');
      END_DT   = datetime('2099-12-31T23:59:59Z', 'TimeZone', 'UTCLeapSeconds');

      actRctJsonPath = bicas.tools.rct.create_BRVF(testCase.testDir, 'biasRctFilename.cdf', BEGIN_DT, END_DT);

      irf.assert.file_exists(actRctJsonPath)
    end



  end    % methods(Test)



end
