%
% matlab.unittest automatic test code for bicas.main().
%
% Since bicas.main() is a very top-level function, this code does not try to
% test everything, but only that which is easy.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef main___UTEST < matlab.unittest.TestCase
  % NOTE: The irfu-matlab git repo does not contain any BICAS config file. That
  %       "has to" be created by the tests (and ideally should be).
  % PROBLEM: Not clear how to test the use of config file: Explicit path .default path, content
  % PROBLEM: How be able to both (1) run the tests locally, and (2) test loading
  %          from the default config file path, without overwriting a
  %          potentially pre-existing config file in the default location?
  %   PROPOSAL: Refuse to run tests if there is a config file in the default
  %             location. -- IMPLEMENTED
  %       PRO: Safe.
  %       CON: Local user has to manually temporarily remove the
  %            pre-existing config file on default path.
  %   PROPOSAL: Tests move a pre-existing default path config file out of the way,
  %             and move it back after the tests.
  %       CON: Dangerous w.r.t. bugs and failed tests(!).
  %           PROPOSAL: Backup default config file.
  %               CON: Pollutes local file system.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      % IMPLEMENTATION NOTE: Code makes sure to always remove the created
      % default path config file so that it does not pollute the local
      % file system, if tests are run locally.

      % Test successful call (error code).
      function test(varargin)
        errorCode = bicas.main(varargin{:});
        testCase.verifyEqual(errorCode, 0)
      end

      % Test unsuccessful call (non-zero error code)
      % NOTE: bicas.main() raises exception but catches it itself.
      function test_error(varargin)
        errorCode = bicas.main(varargin{:});

        testCase.verifyEqual(errorCode, 1)
      end
      %===================================================================

      % Use default path config file.
      configFileAPath = bicas.main___UTEST.setup_default_config_file();
      test('--help')
      test('--version')
      test('--identification')
      test('--swdescriptor')
      %
      test('--help', '--set', 'SWM.L1-L2_ENABLED', '0')
      test('--help', '--set', 'SWM.L1-L2_ENABLED', '1')
      test('--version', '--log', '/ignored_path_to_log_file')
      delete(configFileAPath)

      % Test specifying config file path.
      configFileAPath = bicas.main___UTEST.setup_specified_config_file(testCase, 'test.conf');
      test('--help', '--config', configFileAPath, '--set', 'SWM.L1-L2_ENABLED', '1')

      % ----------------------------------
      % --log-matlab : Log file is created
      % ----------------------------------
      configFileAPath = bicas.main___UTEST.setup_default_config_file();
      logFilePath = fullfile(bicas.main___UTEST.get_temp_dir(testCase), 'bicas.log');
      testCase.verifyFalse(isfile(logFilePath))
      test('--version', '--log-matlab', logFilePath)
      testCase.verifyTrue(isfile(logFilePath))
      delete(configFileAPath)

      configFileAPath = bicas.main___UTEST.setup_default_config_file();
      test_error()    % No CLI argument.
      test_error('illegal_argument')
      test_error('--illegal_argument')
      delete(configFileAPath)
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function tempDir = get_temp_dir(testCase)
      fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      tempDir = fixture.Folder;
    end



    function configFileAPath = setup_default_config_file()
      % Create empty BICAS config file in the default location.
      % (Asserts that no such file pre-exists.
      bicasRootPath = bicas.utils.get_BICAS_root_path();
      configFileAPath = fullfile(bicasRootPath, 'config', 'bicas.conf');

      % ASSERT: Default path config file does not pre-exists.
      % IMPLEMENTATION NOTE: Refuse to overwrite pre-existing file to
      % make sure that the test does not overwrite a legitimate custom
      % config file when running the tests locally.
      if isfile(configFileAPath)
        error('"%s" already pre-exists. Aborting to avoid overwriting presumed local user-defined config file.', configFileAPath)
      end

      fclose(fopen(configFileAPath, 'w'));
    end



    function configFileAPath = setup_specified_config_file(testCase, configFileRPath)
      configFileAPath = fullfile(bicas.main___UTEST.get_temp_dir(testCase), configFileRPath);
      fclose(fopen(configFileAPath, 'w'));
    end



  end    % methods(Static, Access=private)



end
