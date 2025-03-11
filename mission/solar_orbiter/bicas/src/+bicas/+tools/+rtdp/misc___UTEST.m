%
% matlab.unittest automatic test code for bicas.tools.rtdp.misc().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef misc___UTEST < matlab.unittest.TestCase
% PROPOSAL: Verify existence of all .txt files.



  %############
  %############
  % PROPERTIES
  %############
  %############
  properties
    % Additional properties of testCase objects. Needed for setup and teardown
    % methods which store/read their own data from the testCase object.
    dir
    configFile
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
      testCase.dir = Fixture.Folder;

      testCase.configFile = fullfile(testCase.dir, 'config.json');
      testCase.create_config_file(testCase.configFile)

      testCase.create_empty_config_file_datasets(testCase.configFile)
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test0(testCase)
      % TEMP_CONFIG_FILE = '/home/erjo/temp/temp/roctestpkg.json';
      TEMP_CONFIG_FILE = testCase.configFile;

      [actRtpdDir, actZippedRtdp] = bicas.tools.rtdp.misc.create_RCS_test_data_package(...
        testCase.dir, 'A', TEMP_CONFIG_FILE, true, true);

      irf.assert.dir_exists(actRtpdDir)
      irf.assert.file_exists(actZippedRtdp)
    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    function create_config_file(testCase, configFile)
      % PROPOSAL: Put files in subdirectories.
      %   PRO: Avoids filenaming collisions.
      %   PRO: Better test.

      % NOTE: Input filenames (versions) have been modified to not be duplicated
      %       (avoid naming collisions).
      ROWS_CA = {
      '{'
      '  "expectedBicasRootDir": "<BICAS_ROOT_DIR>",'
      '  "inputDatasets": {'
      '    "LFR-SBM1-CWF-E": {'
      '      "in_hk":  "<PARENT_DIR>/solo_HK_rpw-bia_20210715_V05.cdf",'
      '      "in_sci": "<PARENT_DIR>/solo_L1R_rpw-lfr-sbm1-cwf-e-cdag_20210715T234148-20210715T235548_V02.cdf",'
      '      "in_cur": "<PARENT_DIR>/solo_L1_rpw-bia-current-cdag_20210701-20210731_V34.cdf"'
      '    },'
      '    "LFR-SBM2-CWF-E": {'
      '      "in_hk":  "<PARENT_DIR>/solo_HK_rpw-bia_20220922_V06.cdf",'
      '      "in_sci": "<PARENT_DIR>/solo_L1R_rpw-lfr-sbm2-cwf-e-cdag_20220922T134335-20220922T154536_V01.cdf",'
      '      "in_cur": "<PARENT_DIR>/solo_L1_rpw-bia-current-cdag_20220901-20220930_V37.cdf"'
      '    },'
      '    "LFR-SURV-CWF-E": {'
      '      "in_hk":  "<PARENT_DIR>/solo_HK_rpw-bia_20200213_V06.cdf",'
      '      "in_sci": "<PARENT_DIR>/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20200213_V10.cdf",'
      '      "in_cur": "<PARENT_DIR>/solo_L1_rpw-bia-current-cdag_20200211-20200229_V02.cdf"'
      '    },'
      '    "LFR-SURV-SWF-E": {'
      '      "in_hk":  "<PARENT_DIR>/solo_HK_rpw-bia_20200213_V07.cdf",'
      '      "in_sci": "<PARENT_DIR>/solo_L1R_rpw-lfr-surv-swf-e-cdag_20200213_V10.cdf",'
      '      "in_cur": "<PARENT_DIR>/solo_L1_rpw-bia-current-cdag_20200211-20200229_V03.cdf"'
      '    },'
      '    "TDS-LFM-CWF-E": {'
      '      "in_hk":  "<PARENT_DIR>/solo_HK_rpw-bia_20200225_V06.cdf",'
      '      "in_sci": "<PARENT_DIR>/solo_L1R_rpw-tds-lfm-cwf-e-cdag_20200225_V07.cdf",'
      '      "in_cur": "<PARENT_DIR>/solo_L1_rpw-bia-current-cdag_20200211-20200229_V04.cdf"'
      '    },'
      '    "TDS-LFM-RSWF-E": {'
      '      "in_hk":  "<PARENT_DIR>/solo_HK_rpw-bia_20200409_V08.cdf",'
      '      "in_sci": "<PARENT_DIR>/solo_L1R_rpw-tds-lfm-rswf-e-cdag_20200409_V08.cdf",'
      '      "in_cur": "<PARENT_DIR>/solo_L1_rpw-bia-current-cdag_20200401-20200430_V03.cdf"'
      '    }'
      '  }'
      '}'
      };

      s = strjoin(ROWS_CA, '\n');
      s = replace(s, '<PARENT_DIR>',     testCase.dir);
      s = replace(s, '<BICAS_ROOT_DIR>', bicas.utils.get_BICAS_root_dir());

      irf.fs.write_file(configFile, uint8(s)')
    end



    % Read config file and create empty files for the specified datasets.
    function create_empty_config_file_datasets(testCase, configFile)
      uint8Array       = irf.fs.read_file(configFile);
      FileJson         = jsondecode(char(uint8Array)');
      InputDatasetJson = FileJson.inputDatasets;

      fnCa1 = fieldnames(InputDatasetJson);
      for i1 = 1:numel(fnCa1)
        SwmJson = InputDatasetJson.(fnCa1{i1});
        fnCa2 = fieldnames(SwmJson);

        for i2 = 1:numel(fnCa2)
          filePath = SwmJson.(fnCa2{i2});

          irf.assert.path_is_available(filePath)
          irf.fs.write_empty_file({filePath});
        end
      end
    end



  end    % methods(Access=private)



end
