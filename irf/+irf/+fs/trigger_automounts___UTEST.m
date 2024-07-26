%
% matlab.unittest automatic test code for
% irf.fs.trigger_automounts().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef trigger_automounts___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    dirPath
    filePath
    nonexistentFilePath
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function some_method_name1(testCase)
      Fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      testCase.dirPath  = Fixture.Folder;
      testCase.filePath = irf.fs.write_empty_file(...
        {testCase.dirPath, 'empty_test_file.dat'} ...
      );
      testCase.nonexistentFilePath = fullfile(testCase.dirPath, 'nonexistent_file.dat');
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      irf.fs.trigger_automounts(cell(0, 1))
    end



    function test1(testCase)
      irf.fs.trigger_automounts({testCase.dirPath})
    end



    function test2(testCase)
      irf.fs.trigger_automounts({testCase.filePath})
    end



    function test3(testCase)
      irf.fs.trigger_automounts({testCase.dirPath; testCase.filePath})
    end



    function test_nonexistent_file(testCase)
      function call_function()
        irf.fs.trigger_automounts(...
          {testCase.dirPath; testCase.filePath, testCase.nonexistentFilePath}...
          )
      end

      testCase.assertError(@() call_function(), ?MException)
    end



  end    % methods(Test)



end
