%
% matlab.unittest automatic test code for
% solo.qli.batch.interface.get_days_from_IDMRQ().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef interface_get_days_from_IDMRQ___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and
    % teardown methods which store/read their own data from the testCase
    % object.
    fmdQliDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      QliFixture         = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.fmdQliDir = QliFixture.Folder;
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

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_IDMRQ(...
        datasetDirsCa, testCase.fmdQliDir, {'999', '2000-01-01', '2099-01-01'});

      testCase.assertEqual(ActDaysDtArray, solo.qli.const.EMPTY_DT_ARRAY)
    end



    function test_nonempty(testCase)
      F    = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir1 = F.Folder;
      F    = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      dir2 = F.Folder;



      %================================
      % Create quicklooks and datasets
      %================================
      % NOTE: Quicklook (created first) is older than datasets.
      qliPath1     = irf.fs.create_empty_file({testCase.fmdQliDir, '20240101T00_20240102T00.png'});

      % Delay does not appear to be needed. Added for safety.
      pause(0.1)

      datasetPath1 = irf.fs.create_empty_file({dir1, 'solo_L2_swa-pas-eflux_20240101_V02.cdf'});
      datasetPath2 = irf.fs.create_empty_file({dir2, 'solo_L2_mag-rtn-normal_20240102_V02.cdf'});

      % Should be needed since FMDs appear to have a time-resolution of 1 s (?)
      % (Linux).
      % NOTE: This delay is a potential source of mistaken test failures if the
      % local file modification time or clock does not support a time resolution
      % smaller than this delay.
      pause(1.1)

      % NOTE: Quicklook (created last) is more recent than datasets.
      qliPath2     = irf.fs.create_empty_file({testCase.fmdQliDir, '20240102T00_20240103T00.png'});



      %==============
      % Execute test
      %==============
      datasetDirsCa = {dir1; dir2};

      ActDaysDtArray = solo.qli.batch.interface.get_days_from_IDMRQ(...
        datasetDirsCa, testCase.fmdQliDir, {'9999', '2000-01-01', '2099-01-01'});

      testCase.assertEqual(ActDaysDtArray, ...
        solo.qli.utils.umdt({'2024-01-01'}))
    end



  end    % methods(Test)



end
