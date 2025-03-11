%
% matlab.unittest automatic test code for bicas.main().
%
% Since bicas.main() is a very top-level function, this code does not try to
% test everything, but only that which is easy.
%
%
% PROBLEMS WITH TESTS WHICH RELY ON THE DEFAULT CONFIG FILE
% =========================================================
% BICAS uses a config file which it looks for in a default location (relative to
% the source code), unless the config file location is explicitly specified with
% a special (unofficial) CLI option. This creates problems for testing since:
% (1) It is difficult to set up a truly and completely independent "test
%     environment" with its own default config file location: It would need at
%     least (a) its own copy of source code (since the default location is
%     relative to the source code) and (b) set its own MATLAB paths.
% (2) A developer who runs the tests likely has his/her own default config file
%     in the default location already.
% (3) It is dangerous for a test to move away and then move back a valid config
%     file in the default location to enable tests to themselves specify the
%     config file in the default location (e.g. no default cofig file).
%
% NOTE: The tests (2024-07-26) do temporarily replace any (optinally)
%       pre-existing BICAS default config file (or directory, symlink), and then
%       moves it back when the tests have completed (setup+teardown).
%
% NOTE: The irfu-matlab git repo does not contain any BICAS config file. That
%       "has to" be created by the tests (and ideally should be).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef main___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    testDir

    % Path to config file which was found in the default location before the
    % tests launched. Empty if there was no such file.
    oldDefaultConfigFile
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestClassSetup)



    function setup_class(testCase)
      defaultConfigFile = bicas.utils.get_BICAS_default_config_file();

      %==============================================================
      % Temporarily move/rename default config file, if there is one
      %==============================================================
      if ismember(exist(defaultConfigFile), [2, 7])
        Dt        = datetime();
        Dt.Format = 'yyyy-MM-dd''T''hhmmss';
        timestampStr = char(Dt);

        filename = sprintf('%s.%s.%s', ...
          bicas.const.DEFAULT_CONFIG_FILENAME, ...
          mfilename(), timestampStr);

        testCase.oldDefaultConfigFile = fullfile(...
          bicas.utils.get_BICAS_config_dir(), filename);

        [success, errorMessage, ~] = movefile(defaultConfigFile, testCase.oldDefaultConfigFile);
        assert(success, errorMessage)
      else
        testCase.oldDefaultConfigFile = [];
      end
    end



  end
  methods(TestMethodSetup)



    function setup_method(testCase)
      Fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      testCase.testDir = Fixture.Folder;
    end



  end



  %##########
  %##########
  % TEARDOWN
  %##########
  %##########
  methods(TestClassTeardown)



    % NOTE: Empirically, this method is executed also if
    % (1) tests are manually interrupted (Ctrl-C), or
    % (2) tests raise exception.
    function teardown_class(testCase)
      defaultConfigFile = bicas.utils.get_BICAS_default_config_file();

      if isfile(defaultConfigFile)
        % NOTE: delete() does not return any error info!
        delete(defaultConfigFile)
      end

      %===================================================
      % Restore old default config file, if there was one
      %===================================================
      if ~isempty(testCase.oldDefaultConfigFile)
        [success, errorMessage, ~] = movefile(testCase.oldDefaultConfigFile, defaultConfigFile);
        assert(success, ...
          'Failed to restore old default config file from saved copy "%s".', ...
          testCase.oldDefaultConfigFile, errorMessage)
      end
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Use default path config file (location).
    function test_no_error_default_config_file(testCase)

      function test(varargin)
        errorCode = bicas.main(varargin{:});
        testCase.assertEqual(errorCode, 0)
      end

      configFileAPath = bicas.main___UTEST.write_default_config_file();

      test('--help')
      test('--version')
      test('--identification')
      test('--swdescriptor')

      test('--swdescriptor', ...
        '--set', 'SWM.L1-L2_ENABLED',         'true', ...
        '--set', 'SWM.L2-L2_CWF-DSR_ENABLED', 'true', ...
        '--set', 'SWM.L2-L3_ENABLED',         'true')
      test('--help', '--set', 'SWM.L1-L2_ENABLED', 'true')

      test('--version', '--log', '/ignored_path_to_log_file')

      % ----------------------------------
      % --log-matlab : Log file is created
      % ----------------------------------
      logFilePath     = fullfile(testCase.testDir, 'bicas.log');
      testCase.assertFalse(isfile(logFilePath))

      test('--version', '--log-matlab', logFilePath)
      testCase.assertTrue(isfile(logFilePath))

      delete(configFileAPath)
    end



    % Test specifying config file path.
    function test_no_error_explicit_config_file(testCase)

      function test(varargin)
        errorCode = bicas.main(varargin{:});
        testCase.assertEqual(errorCode, 0)
      end

      configFileAPath = bicas.main___UTEST.write_specified_config_file(testCase, 'test.conf');

      test('--help', '--config', configFileAPath, '--set', 'SWM.L1-L2_ENABLED', 'true')
    end



    function test_error(testCase)

      % Test unsuccessful call (non-zero error code)
      % NOTE: bicas.main() raises exception but catches it itself.
      function test_error(varargin)
        errorCode = bicas.main(varargin{:});

        testCase.assertEqual(errorCode, 1)
      end

      configFileAPath = bicas.main___UTEST.write_default_config_file();
      test_error()    % Zero CLI arguments.
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



    function defaultConfigFileAPath = write_default_config_file()
      % Create empty BICAS config file in the default location.
      % Asserts that no such file pre-exists.
      defaultConfigFileAPath = irf.fs.write_empty_file({...
        bicas.utils.get_BICAS_default_config_file()});
    end



    function configFileAPath = write_specified_config_file(testCase, configFileRPath)
      configFileAPath = irf.fs.write_empty_file({testCase.testDir, configFileRPath});
    end



  end    % methods(Static, Access=private)



end
