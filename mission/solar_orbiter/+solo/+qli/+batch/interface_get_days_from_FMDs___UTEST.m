%
% matlab.unittest automatic test code for
% solo.qli.batch.interface.get_days_from_FMDs().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef interface_get_days_from_FMDs___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and
    % teardown methods which store/read their own data from the testCase
    % object.
    qliDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      InputLogFixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      QliFixture      = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      testCase.qliDir      = QliFixture.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_empty(testCase)
      datasetDirsCa = cell(0, 1);

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_FMDs(...
        datasetDirsCa, testCase.qliDir);

      testCase.assertEqual(ActDaysDtArray, solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test_nonempty(testCase)
      F    = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir1 = F.Folder;
      F    = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir2 = F.Folder;

      % NOTE: Quicklook (created first) is older than datasets.
      qliPath      = irf.fs.create_empty_file({testCase.qliDir, '20240101T00_20240102T00.png'});
      pause(0.1)    % Does not appear to be needed. Added for safety.
      datasetPath1 = irf.fs.create_empty_file({dir1, 'solo_L2_swa-pas-eflux_20240101_V02.cdf'});
      datasetPath2 = irf.fs.create_empty_file({dir2, 'solo_L2_mag-rtn-normal_20240101_V02.cdf'});

      datasetDirsCa = {dir1; dir2};

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_FMDs(...
        datasetDirsCa, testCase.qliDir);

      testCase.assertEqual(ActDaysDtArray, ...
        solo.qli.utils.umdt('2024-01-01'))
    end



  end    % methods(Test)



end
