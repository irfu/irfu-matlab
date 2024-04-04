%
% matlab.unittest automatic test code for bicas.tools.batch.get_file_paths().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_file_paths___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Calls which return zero file paths.
    function test_empty(testCase)
      ExpFsoiArray1 = dir('~');
      ExpFsoiArray1 = ExpFsoiArray1([], 1);    % Column array.

      f = testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);

      testDir = f.Folder;
      emptyDir1 = bicas.tools.batch.get_file_paths___UTEST.create_directory(testDir, 'empty_dir1');
      emptyDir2 = bicas.tools.batch.get_file_paths___UTEST.create_directory(testDir, 'empty_dir2');

      function test(fileDirPathsCa)
        [actFilePathsCa1, ActFsoiArray1] = bicas.tools.batch.get_file_paths___UTEST.test_call(...
          testCase, fileDirPathsCa);

        testCase.assertEqual(actFilePathsCa1, cell(0, 1))
        testCase.assertEqual(ActFsoiArray1,   ExpFsoiArray1)
      end

      test(cell(0, 1))
      test({testDir})
      test({emptyDir1; emptyDir2})
    end



    % Calls which return one file path.
    function test_one_file(testCase)
      f = testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);

      testDir  = f.Folder;
      dir1     = fullfile(testDir, 'dir1');
      emptyDir = bicas.tools.batch.get_file_paths___UTEST.create_directory( testDir, 'empty_dir');
      file1    = irf.fs.create_empty_file({testDir, 'dir1/dir1/file1'});

      function test(fileDirPathsCa)
        [actFilePathsCa1, ActFsoiArray1] = bicas.tools.batch.get_file_paths___UTEST.test_call(...
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
      f = testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);

      testDir  = f.Folder;

      dirE    = bicas.tools.batch.get_file_paths___UTEST.create_directory( testDir, 'empty_dir');
      [~]     = bicas.tools.batch.get_file_paths___UTEST.create_directory( testDir, 'empty_dir/dir1');
      dir1    =                                                   fullfile(testDir, 'dir1');
      file11  = irf.fs.create_empty_file({testDir, 'dir1/file1'});
      file111 = irf.fs.create_empty_file({testDir, 'dir1/dir1/file11'});
      file112 = irf.fs.create_empty_file({testDir, 'dir1/dir1/file12'});

      file21  = irf.fs.create_empty_file({testDir, 'dir2/file1'});
      file211 = irf.fs.create_empty_file({testDir, 'dir2/dir1/file11'});

      [actFilePathsCa1, ~] = bicas.tools.batch.get_file_paths___UTEST.test_call(...
        testCase, {dirE; dir1; file21; file211} ...
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

      [actFilePathsCa1, ActFsoiArray1] = bicas.tools.batch.get_file_paths(fileDirPathsCa);

       % Only check presence of fields.
      testCase.assertEqual(...
        sort(fieldnames(ActFsoiArray1)), ...
        sort(expFieldNamesCa))

      actFsoiPathsCa = arrayfun(@(Fsoi) (fullfile(Fsoi.folder, Fsoi.name)), ActFsoiArray1, 'UniformOutput', false);
      testCase.assertEqual(actFilePathsCa1, actFsoiPathsCa)
    end



    function dirPath = create_directory(parentDir, dirRpath)
      % PROPOSAL: Make into generic function.
      %   Cf. irf.fs.create_empty_file().

      dirPath = fullfile(parentDir, dirRpath);
      mkdir(dirPath)
    end



  end    % methods(Static, Access=private)



end
