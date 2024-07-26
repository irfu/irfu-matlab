%
% UNFINISHED
%
% matlab.unittest automatic test code for ?().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef generate_official_SWD_file___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and teardown
    % methods which store/read their own data from the testCase object.
    testDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(testCase)
      Fixture = testCase.applyFixture(...
        matlab.unittest.fixtures.TemporaryFolderFixture);
      % NOTE: The same fixture should always return the same directory.
      testCase.testDir = Fixture.Folder;
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_write_compare(testCase)

      %==========================
      % Test generating SWF file
      %==========================
      officialSwdFile = bicas.utils.get_SWD_file();
      testSwdFile     = fullfile(testCase.testDir, 'descriptor.json');

      bicas.tools.generate_official_SWD_file(testSwdFile)

      % ASSERT: File exists and is non-empty.
      irf.assert.file_exists(testSwdFile)
      testCase.assertTrue(dir(testSwdFile).bytes > 0)



      if 1
        %================================================================
        % Compare generated SWD file with SWD file stored with BICAS (as
        % specified by RCS ICD).
        %================================================================
        officialBytesArray = irf.fs.read_file(officialSwdFile);
        testBytesArray     = irf.fs.read_file(testSwdFile);

        if ~isequaln(officialBytesArray, testBytesArray)
          error( ...
            'The official SWD ("%s") stored in BICAS root directory is not identical to the generated SWD ("%s").', ...
            officialBytesArray, testBytesArray)
        end
      end
    end



  end    % methods(Test)



end
