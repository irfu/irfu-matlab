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
      officialSwdFile  = bicas.utils.get_SWD_file();
      generatedSwdFile = fullfile(testCase.testDir, bicas.const.SWD_FILENAME);

      bicas.tools.generate_official_SWD_file(generatedSwdFile)

      % ASSERT: File exists and is non-empty.
      irf.assert.file_exists(generatedSwdFile)
      testCase.assertTrue(dir(generatedSwdFile).bytes > 0)



      if 1
        %=======================================================================
        % Compare (1) generated SWD file with (2) SWD file stored with BICAS (as
        % specified by RCS ICD).
        %=======================================================================
        officialSwdBytesArray  = irf.fs.read_file(officialSwdFile);
        generatedSwdBytesArray = irf.fs.read_file(generatedSwdFile);

        if ~isequaln(officialSwdBytesArray, generatedSwdBytesArray)
          error( ...
            ['The official SWD ("%s") stored in BICAS root directory is', ...
            ' not identical to the generated SWD ("%s").'], ...
            officialSwdFile, generatedSwdFile)
        end
      end
    end



  end    % methods(Test)



end
