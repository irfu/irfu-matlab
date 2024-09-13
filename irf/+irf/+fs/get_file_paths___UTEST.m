%
% matlab.unittest automatic test code for irf.fs.get_file_paths().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_file_paths___UTEST < matlab.unittest.TestCase



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



    function some_method_name1(testCase)
      F = testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
      testCase.testDir = F.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % NOTE: Mostly for manual inspection of error messages (by inspecting
    % variables; error messages are are not printed on screen when using
    % runtest()). Otherwise not a great test.
    function test_nonexisting(testCase)
      existingFilePath = irf.fs.write_empty_file({testCase.testDir, 'existing_file'});

      function test(fileDirPathsCa)
        assert(iscolumn(fileDirPathsCa))

        testCase.assertError(...
          @() irf.fs.get_file_paths(fileDirPathsCa), ...
          ?MException)
      end

      test({'/nonexisting_path'});
      test({ ...
        existingFilePath; ...
        '/nonexisting_path1'; ...
        '/nonexisting_path2'; ...
        testCase.testDir});
    end



    % Calls which return zero file paths.
    function test_empty(testCase)
      testDir = testCase.testDir;

      ExpFsoiArray1 = dir('~');
      ExpFsoiArray1 = ExpFsoiArray1([], 1);    % Column array.

      emptyDir1 = irf.fs.get_file_paths___UTEST.create_directory(testDir, {'empty_dir1'});
      emptyDir2 = irf.fs.get_file_paths___UTEST.create_directory(testDir, {'empty_dir2'});

      function test(fileDirPathsCa)
        [actFilePathsCa1, ActFsoiArray1] = irf.fs.get_file_paths___UTEST.test_call(...
          testCase, fileDirPathsCa);

        testCase.assertEqual(actFilePathsCa1, cell(0, 1))
        testCase.assertEqual(ActFsoiArray1,   ExpFsoiArray1)
      end

      % Zero directories.
      test(cell(0, 1))
      % One directory with empty subdirectories.
      test({testDir})
      % Two empty directories.
      test({emptyDir1; emptyDir2})
    end



    % Calls which return one file path.
    function test_one_file(testCase)
      testDir = testCase.testDir;

      dir1     = fullfile(testDir, 'dir1');
      emptyDir = irf.fs.get_file_paths___UTEST.create_directory( testDir, {'empty_dir'});
      file1    = irf.fs.write_empty_file({testDir, 'dir1', 'dir1', 'file1'});

      function test(fileDirPathsCa)
        [actFilePathsCa1, ActFsoiArray1] = irf.fs.get_file_paths___UTEST.test_call(...
          testCase, fileDirPathsCa);

        testCase.assertEqual(actFilePathsCa1,      {file1})
      end

      test({file1})
      test({dir1});
      test({testDir});
      test({dir1; emptyDir});
    end



    % Empty directory, non-empty directory, subdirectories, explicit file paths
    function test_mixed_complex(testCase)
      testDir = testCase.testDir;

      dirWithEmptySubdir = irf.fs.get_file_paths___UTEST.create_directory( testDir, {'dir_with_empty_subdir'});
      [~]                = irf.fs.get_file_paths___UTEST.create_directory( testDir, {'dir_with_empty_subdir', 'dir1'});

      dir1     = fullfile(                 testDir, 'dir1');
      file11   = irf.fs.write_empty_file({testDir, 'dir1', 'file1'});
      file111  = irf.fs.write_empty_file({testDir, 'dir1', 'dir1', 'file11'});
      file112  = irf.fs.write_empty_file({testDir, 'dir1', 'dir1', 'file12'});

      file21   = irf.fs.write_empty_file({testDir, 'dir2', 'file1'});
      file211  = irf.fs.write_empty_file({testDir, 'dir2', 'dir1', 'file11'});

      [actFilePathsCa1, ~] = irf.fs.get_file_paths___UTEST.test_call(...
        testCase, {dirWithEmptySubdir; dir1; file21; file211} ...
        );

      testCase.assertEqual(...
        sort(actFilePathsCa1), ...
        sort({file11; file111; file112; file21; file211})...
        )
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function [actFilePathsCa1, ActFsoiArray1] = test_call(testCase, fileDirPathsCa)
      expFieldNamesCa = fieldnames(dir('~'));

      [actFilePathsCa1, ActFsoiArray1] = irf.fs.get_file_paths(fileDirPathsCa);

      % Check existence of struct fields.
      testCase.assertEqual(...
        sort(fieldnames(ActFsoiArray1)), ...
        sort(expFieldNamesCa))

      actFsoiPathsCa = arrayfun(...
        @(Fsoi) (fullfile(Fsoi.folder, Fsoi.name)), ActFsoiArray1, ...
        'UniformOutput', false);
      testCase.assertEqual(actFilePathsCa1, actFsoiPathsCa)
    end



    function dirPath = create_directory(parentDir, dirRpathPartsCa)
      % PROPOSAL: Make into generic function.
      %   Cf. irf.fs.write_empty_file().

      dirPath = fullfile(parentDir, dirRpathPartsCa{:});
      mkdir(dirPath)
    end



  end    % methods(Static, Access=private)



end
