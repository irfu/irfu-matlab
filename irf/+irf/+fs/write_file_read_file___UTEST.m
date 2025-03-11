%
% matlab.unittest automatic test code for irf.fs.write_file() and
% irf.fs.read_file().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef write_file_read_file___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and teardown
    % methods which store/read their own data from the testCase object.
    filePath
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      Fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testDir = Fixture.Folder;

      testCase.filePath = fullfile(testDir, 'file.txt');
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_overwrite(testCase)
      irf.fs.write_file(testCase.filePath, uint8((0:15)'))

      irf.assert.file_exists(testCase.filePath)
      fileSize = dir(testCase.filePath).bytes;
      testCase.assertEqual(fileSize, 16)

      irf.fs.write_file(testCase.filePath, uint8((0:255)'))

      fileSize = dir(testCase.filePath).bytes;
      testCase.assertEqual(fileSize, 256)

      irf.fs.write_file(testCase.filePath, uint8((0:31)'))

      fileSize = dir(testCase.filePath).bytes;
      testCase.assertEqual(fileSize, 32)
    end



    function test_write_read(testCase)

      function test(doubleArrayWrite)
        assert(iscolumn(doubleArrayWrite))

        byteArrayWrite = uint8(doubleArrayWrite);

        irf.fs.write_file(testCase.filePath, byteArrayWrite);
        actByteArray = irf.fs.read_file(testCase.filePath);

        testCase.assertTrue(iscolumn(actByteArray))
        testCase.assertEqual(actByteArray, byteArrayWrite)
      end

      test(zeros(0, 1))
      test(1)
      test((0:255)')
    end



  end    % methods(Test)



end
