%
% matlab.unittest automatic test code for solo.db_list_files().
%
%
% NOTES
% =====
% * In reality, the tested code solo.db_list_files(), is a wrapper around method
%   solo_db.list_files(). Code should maybe be rewritten for this.
% * The tested code underpins solo.db_get_ts(). Therefore, this test code
%   implicitly tests that function.
% * Tests have to begin with "clear global" since the tested code uses
%   global variables.
% * Tests do not cover all behaviour. They were written long after the tested
%   code was written by other authors in anticipation of an update that was
%   never implemented. In particular, there are no tests for having multiple
%   datasets of the same type (2024-03-12).
% * Tests do not require access to SolO datasets. The tests create their own
%   CDFs (sic!).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef db_list_files___UTEST < matlab.unittest.TestCase



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



    function init(testCase)
      % Tests have to begin with "clear global" since the tested code uses
      % global variables.
      clear global

      F = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.testDir = F.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Zero CDFs.
    function test_empty(testCase)
        Tint = irf.tint('2023-12-31T00:00:00/2024-01-02T00:00:00');
        solo.db_init('local_file_db', testCase.testDir);

        FileList = solo.db_list_files('solo_L2_nonexisting-dataset', Tint);
        testCase.assertEqual(FileList, [])
    end



    % ================================
    % DATASET WITH MONTHLY DIRECTORIES
    % ================================

    function test_one_file_M_1(testCase)
      % Only 24h, UTC midnights. Everything aligns.
      CDF_TIME_CA = {'2024-01-01T00:00:00', '2024-01-02T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_month(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA)
    end

    function test_one_file_M_2(testCase)
      % File contains MUCH more (superset) than 24h in filename.
      CDF_TIME_CA = {'2023-01-01T00:00:00', '2025-01-02T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_month(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA)
    end

    % -----------------------------------------------------------------------
    % Check for mismatch between filename and file data
    % -------------------------------------------------
    % NOTE: Maybe unnecessary to check the difference in behaviour between
    % files with data not in the day and month. The difference is not too
    % important for the user and is probably an artefact from implementation.
    % Can however not check behaviour without explicitly making that
    % difference.
    % -----------------------------------------------------------------------

    % File does not contain data for the DAY in the filename (but has data for
    % the same month).
    function test_one_file_M_3(testCase)
      CDF_TIME_CA = {'2024-01-31T00:00:00', '2024-02-01T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_month(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA)
    end

    function test_one_file_M_4(testCase)
      % File does not contain data for the MONTH in the filename.
      CDF_TIME_CA = {'2024-02-03T00:00:00', '2024-02-04T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_month(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        [])
    end

    function test_one_file_M_5(testCase)
      CDF_TIME_CA = {'2023-12-03T00:00:00', '2023-12-31T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_month(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        [])
    end

    % File does not contain the data in the filename (not same month) and
    % not the data that is searched for.
    function test_one_file_M_6(testCase)
      solo.db_list_files___UTEST.test_20240101_V03_month(...
        testCase, ...
        {'2024-02-03T00:00:00', '2024-02-04T00:00:00'}, ...
        {'2024-01-01T00:00:00', '2024-01-02T00:00:00'}, ...
        [])
    end

    % File does not contain the data in the filename (same month) and not
    % the data that is searched for.
    function test_one_file_M_7(testCase)
      solo.db_list_files___UTEST.test_20240101_V03_month(...
        testCase, ...
        {'2024-01-03T00:00:00', '2024-01-04T00:00:00'}, ...
        {'2024-01-01T00:00:00', '2024-01-02T00:00:00'}, ...
        [])
    end



    % ==============================
    % DATASET WITH DAILY DIRECTORIES
    % ==============================

    % Only 24h, UTC midnights. Everything aligns.
    function test_one_file_D_1(testCase)
      CDF_TIME_CA = {'2024-01-01T00:00:00', '2024-01-02T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_day(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA)
    end

    % File contains MUCH more (superset) than 24h in filename.
    function test_one_file_D_2(testCase)
      CDF_TIME_CA = {'2023-01-01T00:00:00', '2025-01-02T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_day(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA)
    end



    % -----------------------------------------------------------------------
    % File data overlaps with day in the filename (part inside, part outside)
    % -----------------------------------------------------------------------

    function test_one_file_D_3(testCase)
      CDF_TIME_CA = {'2024-01-01T12:00:00', '2024-01-02T12:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_day(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA)
    end

    function test_one_file_D_4(testCase)
      CDF_TIME_CA = {'2023-12-31T12:00:00', '2024-01-01T12:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_day(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA)
    end



    % ------------------------------------------------------
    % File does not contain data for the DAY in the filename
    % ------------------------------------------------------

    function test_one_file_D_5(testCase)
      CDF_TIME_CA = {'2024-01-31T00:00:00', '2024-02-01T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_day(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        [])
    end

    function test_one_file_D_6(testCase)
      CDF_TIME_CA = {'2023-12-30T00:00:00', '2023-12-31T00:00:00'};
      solo.db_list_files___UTEST.test_20240101_V03_day(...
        testCase, ...
        CDF_TIME_CA, ...
        CDF_TIME_CA, ...
        [])
    end



    % =======================
    % Test CDAG (Mis)matching
    % =======================

    function test_CDAG_prefix_1(testCase)
      solo.db_list_files___UTEST.test_differing_prefixes(...
        testCase, ...
        'solo_L2_rpw-lfr-surv-cwf-e-cdag', ...
        'solo_L2_rpw-lfr-surv-cwf-e-cdag', true)
    end

    function test_CDAG_prefix_2(testCase)
      solo.db_list_files___UTEST.test_differing_prefixes(...
        testCase, ...
        'solo_L2_rpw-lfr-surv-cwf-e', ...
        'solo_L2_rpw-lfr-surv-cwf-e', true)
    end

    function test_CDAG_prefix_3(testCase)
      solo.db_list_files___UTEST.test_differing_prefixes(...
        testCase, ...
        'solo_L2_rpw-lfr-surv-cwf-e', ...
        'solo_L2_rpw-lfr-surv-cwf-e-cdag', false)
    end

    function test_CDAG_prefix_4(testCase)
      solo.db_list_files___UTEST.test_differing_prefixes(...
        testCase, ...
        'solo_L2_rpw-lfr-surv-cwf-e-cdag', ...
        'solo_L2_rpw-lfr-surv-cwf-e', false)
    end



    function test_simultaneous_CDAG_nonCDAG_files(testCase)
      % Use same time interval for everything.
      TIME_CA = {'2024-01-01T00:00:00', '2024-01-02T00:00:00'};

      FILENAME_PREFIX_CA = {...
        'solo_L2_rpw-lfr-surv-cwf-e'; ...
        'solo_L2_rpw-lfr-surv-cwf-e-cdag'};

      Tint = irf.tint(spdfparsett2000(TIME_CA));

      % Create CDFs on disk.
      testFileCa = {};
      for i = 1:numel(FILENAME_PREFIX_CA)
        testFileCa{i} = fullfile(testCase.testDir, ...
          'remote', 'data', 'L2', 'lfr_wf_e', ...
          '2024', '01', sprintf('%s_20240101_V03.cdf', FILENAME_PREFIX_CA{i}));

        solo.db_list_files___UTEST.write_CDF(...
          testFileCa{i}, TIME_CA{:})
      end

      % List CDFs.
      solo.db_init('local_file_db', testCase.testDir);
      for i = 1:numel(FILENAME_PREFIX_CA)
        expFilename = irf.fs.get_name(testFileCa{i});

        FileList = solo.db_list_files(FILENAME_PREFIX_CA{i}, Tint);

        % Assert only get the file explicitly searched for
        % ------------------------------------------------
        % NOTE: Simplified test of return value.
        % NOTE: Both filename prefixes are treated as separate types of
        % datasets.
        testCase.assertTrue(isscalar(FileList))
        testCase.assertEqual(FileList.name, expFilename)
      end

    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Test with one CDF file.
    %
    % ARGUMENTS
    % =========
    % argFilenamePrefix
    %     Actual filename prefix that is search for (used as argument).
    % fileFilenamePrefix
    %     Prefix of actual file.
    function test_differing_prefixes(testCase, argFilenamePrefix, fileFilenamePrefix, actFound)
      assert(islogical(actFound))

      % Use same time interval for everything.
      CDF_TIME_CA    = {'2024-01-01T00:00:00', '2024-01-02T00:00:00'};
      TINT_TIME_CA   = CDF_TIME_CA;
      SOLODB_TIME_CA = CDF_TIME_CA;

      Tint = irf.tint(spdfparsett2000(TINT_TIME_CA));

      testFile = fullfile(testCase.testDir, ...
        'remote', 'data', 'L2', 'lfr_wf_e', ...
        '2024', '01', sprintf('%s_20240101_V03.cdf', fileFilenamePrefix));

      testFileName      = irf.fs.get_name(testFile);
      testFileParentDir = fileparts(testFile);

      solo.db_list_files___UTEST.write_CDF(testFile, CDF_TIME_CA{:})

      solo.db_init('local_file_db', testCase.testDir);
      FileList = solo.db_list_files(argFilenamePrefix, Tint);

      if actFound
        testCase.assertTrue(isscalar(FileList))
        testCase.assertEqual(FileList.ver,        3)
        testCase.assertEqual(FileList.name,       testFileName)
        testCase.assertEqual(FileList.path,       testFileParentDir)
        testCase.assertEqual(FileList.start.ttns, spdfparsett2000(SOLODB_TIME_CA{1}))
        testCase.assertEqual(FileList.stop.ttns,  spdfparsett2000(SOLODB_TIME_CA{2}))
      else
        testCase.assertEqual(FileList, [])
      end
    end



    % Test with one CDF file.
    %
    % ARGUMENTS
    % =========
    % 1x2 cell arrays with timestamps as UTC strings.
    % --
    % cdfUtcCa
    %     Range of timestamps to use in CDF file.
    % tintTimeCa
    %     Min-max timestamps for search interval (solo.db_list_files() tint
    %     argument).
    % expSoloDbTimeCa
    %     (1) Expected time range specified for file in return value from
    %         solo.db_list_files().
    %     (2) Empty <=> Expect no recognized file.
    function test_20240101_V03(testCase, fileRpath, cdfUtcCa, tintTimeCa, expSoloDbTimeCa)
      Tint = irf.tint(spdfparsett2000(tintTimeCa));

      testFile          = fullfile(testCase.testDir, 'remote', 'data', fileRpath);
      testFileName      = irf.fs.get_name(testFile);
      testFileParentDir = fileparts(testFile);
      R = solo.adm.dsfn.parse_dataset_filename(testFileName);
      filePrefix        = R.fnDatasetIdCdag;

      solo.db_list_files___UTEST.write_CDF(testFile, cdfUtcCa{:})

      solo.db_init('local_file_db', testCase.testDir);
      FileList = solo.db_list_files(filePrefix, Tint);

      if ~isempty(expSoloDbTimeCa)
        testCase.assertTrue(isscalar(FileList))
        testCase.assertEqual(FileList.ver,        3)
        testCase.assertEqual(FileList.name,       testFileName)
        testCase.assertEqual(FileList.path,       testFileParentDir)
        testCase.assertEqual(FileList.start.ttns, spdfparsett2000(expSoloDbTimeCa{1}))
        testCase.assertEqual(FileList.stop.ttns,  spdfparsett2000(expSoloDbTimeCa{2}))
      else
        testCase.assertEqual(FileList, [])
      end
    end



    % Type of dataset where the standard directory structure covers one day
    % per directory.
    function test_20240101_V03_day(testCase, cdfUtcCa, tintTimeCa, expSoloDbTimeCa)
      % NOTE: File can not be chosen freely. Has to be consistent with test().
      FILE_RPATH = fullfile(...
        'HK', '2024', '01', '01', 'solo_HK_rpw-bia_20240101_V03.cdf');

      solo.db_list_files___UTEST.test_20240101_V03(...
        testCase, FILE_RPATH, cdfUtcCa, tintTimeCa, expSoloDbTimeCa)
    end



    % Type of dataset where the standard directory structure covers one month
    % per directory.
    function test_20240101_V03_month(testCase, cdfUtcCa, tintTimeCa, expSoloDbTimeCa)
      % NOTE: File can not be chosen freely. Has to be consistent with test().
      FILE_R_PATH = fullfile(...
        'L2', 'lfr_wf_e', ...
        '2024', '01', 'solo_L2_rpw-lfr-surv-cwf-e-cdag_20240101_V03.cdf');

      solo.db_list_files___UTEST.test_20240101_V03(...
        testCase, FILE_R_PATH, cdfUtcCa, tintTimeCa, expSoloDbTimeCa)
    end



    % Create minimal CDF with information which SolO DB will read (can not
    % just write empty files). Also creates parent directory (recursively) if
    % needed.
    function write_CDF(filePath, epochStartTimeUtcStr, epochTimeUtcStr)
      epochStartTt2000 = spdfparsett2000(epochStartTimeUtcStr);
      epochStopTt2000  = spdfparsett2000(epochTimeUtcStr);

      [parentDir, ~, ~] = fileparts(filePath);
      if ~exist(parentDir, 'dir')
          mkdir(parentDir)
      end

      spdfcdfwrite(filePath, ...
        {'Epoch', int64([epochStartTt2000; epochStopTt2000])}, ...
        'Vardatatypes', {'Epoch', 'CDF_TIME_TT2000'}, ...
        'RecordBound',  {'Epoch'})
    end



  end    % methods(Static, Access=private)



end
