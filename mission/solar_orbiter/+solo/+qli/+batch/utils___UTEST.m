%
% matlab.unittest automatic test code for solo.qli.batch.utils.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils___UTEST < matlab.unittest.TestCase



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    configFilePath
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function init_test(testCase)
      Fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);

      testCase.configFilePath = fullfile(Fixture.Folder, 'config.json');
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function read_config_file(testCase)
      IRF_LOGO_PATH = '/home/erjo/so_irfu-matlab_qli_IRFpl/mission/solar_orbiter/+solo/+qli/+batch/irf_logo.png';

      solo.qli.batch.utils___UTEST.write_config_file(testCase.configFilePath, IRF_LOGO_PATH)

      ActConfig = solo.qli.batch.utils.read_config_file(testCase.configFilePath);

      % NOTE: Only checking some values.
      testCase.assertEqual(ActConfig.irfLogoPath, IRF_LOGO_PATH)
      testCase.assertEqual(ActConfig.vhtDir, '/data/solo/data_yuri/')
      testCase.assertEqual( ...
        sort(ActConfig.datasetDirsCa), ...
        sort({'/data/solo/remote/data'; '/data/solo/soar'}) ...
        )

      ExpLogFileDirPatternDict = dictionary();
      ExpLogFileDirPatternDict('SOAR')  = "/home/erjo/logs/so_soar_irfu_mirror_sync.*.log";
      ExpLogFileDirPatternDict('LESIA') = "/home/erjo/logs/pull.so.data.cron.brain.*.log";
      testCase.assertEqual(ActConfig.LogFileDirPatternDict, ExpLogFileDirPatternDict)
    end



    function read_config_file_no_logo(testCase)
      solo.qli.batch.utils___UTEST.write_config_file(testCase.configFilePath, [])

      ActConfig = solo.qli.batch.utils.read_config_file(testCase.configFilePath);

      % NOTE: Only checking some values.
      testCase.assertEqual(ActConfig.irfLogoPath, [])
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function write_config_file(configFilePath, irfLogoPath)
      % '  "irfLogoPath": "/home/erjo/so_irfu-matlab_qli_IRFpl/mission/solar_orbiter/+solo/+qli/+batch/irf_logo.png",'

      if isempty(irfLogoPath)
        irfLogoPathJsonStr = 'null';
      else
        irfLogoPathJsonStr = sprintf('"%s"', irfLogoPath);
      end

      ROWS_CA = {
        '{'
        '  "irfLogoPath": "/home/erjo/so_irfu-matlab_qli_IRFpl/mission/solar_orbiter/+solo/+qli/+batch/irf_logo.png",'
        '  "vhtDir":       "/data/solo/data_yuri/",'
        '  "solo_db_init_local_file_db": "/data/solo/",'
        '  "datasetDirs": ['
        '    "/data/solo/remote/data",'
        '    "/data/solo/soar"'
        '  ],'
        '  "fmdQliDir":    "/data/juice/EJ_temp/quicklooks_SOLAR_ORBITER/www/24h",'
        '  "logFileDirPatterns": {'
        '    "LESIA": "/home/erjo/logs/pull.so.data.cron.brain.*.log",'
        '    "SOAR":  "/home/erjo/logs/so_soar_irfu_mirror_sync.*.log"'
        '  }'
        '}'
        };

      % Ugly hack?
      ROWS_CA{2} = sprintf('  "irfLogoPath": %s,', irfLogoPathJsonStr);

      solo.qli.batch.utils.write_file(...
        configFilePath, ...
        ROWS_CA)
    end



  end    % methods(Static, Access=private)



end
