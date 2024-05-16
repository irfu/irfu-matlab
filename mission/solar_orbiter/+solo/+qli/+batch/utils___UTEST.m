%
% matlab.unittest automatic test code for solo.qli.batch.utils.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Replace hardcoded JSON config file with actual "example config
  %           file" on disk.
  %   PROPOSAL: Use as OFFGEN config file.



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function read_config_file(testCase)
      Fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
      configFilePath = fullfile(Fixture.Folder, 'config.json');

      solo.qli.batch.utils.write_file(...
        configFilePath, ...
        {
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
        } ...
      )

      ActConfig = solo.qli.batch.utils.read_config_file(configFilePath);      

      testCase.verifyEqual(ActConfig.vhtDir, '/data/solo/data_yuri/')
    end



  end    % methods(Test)



end
