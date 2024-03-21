%
% matlab.unittest automatic test code for
% bicas.tools.batch.get_BPCI_output_path2().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef get_BPCI_output_path2___UTEST < matlab.unittest.TestCase



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  properties(TestParameter)
    FILENAME_VERSION_ALGORITHM = {'HIGHEST_USED', 'ABOVE_HIGHEST_USED'}

    % YMD = Year-Month-Day
    BPCI_INPUT_YMD_FILES = {
      {
      'solo_L1R_rpw-lfr-surv-cwf-e_20220101_V03.cdf', ...  % No CDAG
      }, ...
      { ...
      'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20220101_V03.cdf' ...   % CDAG
      }, ...
      {
      % Different time intervals.
      'solo_L1_rpw-bia-current-cdag_20220101-20220131_V15.cdf', ...
      'solo_L1R_rpw-lfr-surv-cwf-e_20220101_V03.cdf', ...  % No CDAG
      }, ...
      }

    % YMDHMS = Year-Month-Day-Hour-Minute-Second
    % No CDAG
    BPCI_INPUT_YMDHMS_FILES = {
      {
      'solo_L1R_rpw-lfr-sbm1-cwf-e_20220101T012345-20220101T123456_V02.cdf', ...
      }, ...
      {
      % Different time intervals.
      'solo_L1_rpw-bia-current_20220101-20220131_V15.cdf', ...
      'solo_L1R_rpw-lfr-sbm1-cwf-e_20220101T012345-20220101T123456_V02.cdf', ...
      }, ...
      }
  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_YMD_CDAG_output_dir_no_preexisting(testCase, FILENAME_VERSION_ALGORITHM, BPCI_INPUT_YMD_FILES)
      % CDAG output filename
      % Non-empty-string output directory.
      % No pre-existing datasets.
      bicas.tools.batch.get_BPCI_output_path2___UTEST.test(...
        testCase, ...
        BPCI_INPUT_YMD_FILES, ...
        {...
        }, ...
        'SOLO_L2_RPW-LFR-SURV-CWF-E', FILENAME_VERSION_ALGORITHM, 'out', true, ...
        fullfile('out', 'solo_L2_rpw-lfr-surv-cwf-e-cdag_20220101_V01.cdf') ...
        )
    end



    function test_YMD_preexisting_output(testCase, BPCI_INPUT_YMD_FILES, FILENAME_VERSION_ALGORITHM)

      switch(FILENAME_VERSION_ALGORITHM)
        case 'HIGHEST_USED'
          expOutputFilename = 'solo_L2_rpw-lfr-surv-cwf-e_20220101_V04.cdf';
        case 'ABOVE_HIGHEST_USED'
          expOutputFilename = 'solo_L2_rpw-lfr-surv-cwf-e_20220101_V05.cdf';
      end

      bicas.tools.batch.get_BPCI_output_path2___UTEST.test(...
        testCase, ...
        BPCI_INPUT_YMD_FILES, ...
        {...
        fullfile('ref', 'solo_L2_rpw-lfr-surv-cwf-e-cdag_20211231_V01.cdf'), ...  % Irrelevant adjacent date.
        fullfile('ref', 'solo_L2_rpw-lfr-surv-cwf-e-cdag_20220101_V04.cdf'), ...
        fullfile('ref', 'solo_L2_rpw-lfr-surv-cwf-e-cdag_20220102_V02.cdf'), ...  % Irrelevant adjacent date.
        }, ...
        'SOLO_L2_RPW-LFR-SURV-CWF-E', FILENAME_VERSION_ALGORITHM, '', false, ...
        expOutputFilename ...
        )
    end



    function test_YMDHMS_output_dir_no_preexisting(testCase, FILENAME_VERSION_ALGORITHM, BPCI_INPUT_YMDHMS_FILES)
      % CDAG output filename
      % Non-empty-string output directory.
      % No pre-existing datasets.
      bicas.tools.batch.get_BPCI_output_path2___UTEST.test(...
        testCase, ...
        BPCI_INPUT_YMDHMS_FILES, ...
        {...
        }, ...
        'SOLO_L2_RPW-LFR-SBM1-CWF-E', FILENAME_VERSION_ALGORITHM, ...
        'out', true, ...
        fullfile('out', 'solo_L2_rpw-lfr-sbm1-cwf-e-cdag_20220101T012345-20220101T123456_V01.cdf') ...
        )
    end



    function test_YMDHMS_preexisting_output(testCase, BPCI_INPUT_YMDHMS_FILES, FILENAME_VERSION_ALGORITHM)

      switch(FILENAME_VERSION_ALGORITHM)
        case 'HIGHEST_USED'
          expOutputFilename = 'solo_L2_rpw-lfr-sbm1-cwf-e_20220101T012345-20220101T123456_V04.cdf';
        case 'ABOVE_HIGHEST_USED'
          expOutputFilename = 'solo_L2_rpw-lfr-sbm1-cwf-e_20220101T012345-20220101T123456_V05.cdf';
      end

      bicas.tools.batch.get_BPCI_output_path2___UTEST.test(...
        testCase, ...
        BPCI_INPUT_YMDHMS_FILES, ...
        {...
        'solo_L2_rpw-lfr-sbm1-cwf-e_20220101_V01.cdf', ...  % Irrelevant time interval. Unexpected time format.
        'solo_L2_rpw-lfr-sbm1-cwf-e_20220101T012345-20220101T123450_V01.cdf', ...  % Different but similar time interval.
        'solo_L2_rpw-lfr-sbm1-cwf-e_20220101T012345-20220101T123456_V04.cdf', ...
        'solo_L2_rpw-lfr-sbm1-cwf-e_20220101T012340-20220101T123456_V01.cdf', ...  % Different but similar time interval.
        }, ...
        'SOLO_L2_RPW-LFR-SBM1-CWF-E', FILENAME_VERSION_ALGORITHM, ...
        '', false, ...
        expOutputFilename ...
        )
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function test(...
        testCase, ...
        bpciInputPathCa, ...        % Converted to DSMD.
        preexistingOutputLvCa, ...  % Converted to DSMD.
        outputDsi, fnVerAlgorithm, outputDir, outputIsCdag, expFilePath)

      BpciInputDsmdArray           = solo.adm.paths_to_DSMD_array(bpciInputPathCa(:));
      PreexistingOutputLvDsmdArray = solo.adm.paths_to_DSMD_array(preexistingOutputLvCa(:));

      % CALL TESTED CODE
      actFilePath = bicas.tools.batch.get_BPCI_output_path2(...
        BpciInputDsmdArray, PreexistingOutputLvDsmdArray, ...
        outputDsi, fnVerAlgorithm, outputDir, outputIsCdag ...
        );

      testCase.assertEqual(actFilePath, expFilePath)
    end



  end    % methods(Static, Access=private)



end
