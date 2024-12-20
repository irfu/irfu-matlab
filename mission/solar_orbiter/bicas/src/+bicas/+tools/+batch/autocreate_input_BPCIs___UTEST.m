%
% matlab.unittest automatic test code for
% bicas.tools.batch.autocreate_input_BPCIs().
%
% NOTE: Not (yet) testing
% indirect use of solo.adm.convert_DSMD_DATASET_ID_to_SOLO()
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef autocreate_input_BPCIs___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_one_BPCIs(testCase)

      function test(inputDatasetsPathsCa, SwmArray, ExpBpciArray)
        bicas.tools.batch.autocreate_input_BPCIs___UTEST.test(...
          testCase, inputDatasetsPathsCa, SwmArray, 0, ExpBpciArray)
      end



      %==================
      % Define constants
      %==================

      FILE_1           = 'L1R/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240101_V01.cdf';
      FILE_NONMATCHING = 'HK_/solo_HK_rpw-bia_20240101_V05.cdf';

      % SWM that matches one file.
      SWMP = bicas.tools.batch.TestSwmProcessing();
      SWM_1 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI-TEST-1', 'Human readable purpose', ...
        bicas.swm.InputDataset(...
        'cli_input', 'SOLO_L1R_RPW-LFR-SURV-CWF-E', 'IN_cdf' ...
        ), ...
        bicas.swm.OutputDataset(...
        'cli_output', 'SOLO_L2_RPW-LFR-SURV-CWF-E', 'OUT_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        );
      % SWM that does not match any file.
      SWM_NONMATCHING = bicas.swm.SoftwareMode(...
        SWMP, 'CLI-TEST-NONMATCHING', 'Human readable purpose', ...
        bicas.swm.InputDataset(...
        'cli_input', 'SOLO_L1_RPW-BIA-CURRENT', 'IN_cdf' ...
        ), ...
        bicas.swm.OutputDataset(...
        'cli_output', 'SOLO_L2_RPW-LFR-SURV-CWF-E', 'OUT_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        );

      EXP_BPCI = bicas.tools.batch.BicasProcessingCallInfo(...
        SWM_1.cliOption, ...
        bicas.tools.batch.BpciInput(...
        SWM_1.inputsList.cliOptionHeaderBody, ...
        SWM_1.inputsList.dsi, FILE_1), ...
        bicas.tools.batch.BpciOutput(...
        SWM_1.outputsList.cliOptionHeaderBody, ...
        SWM_1.outputsList.dsi, ...
        'SOLO_L2_RPW-LFR-SURV-CWF-E___output_dataset.cdf'));



      ZERO_BPCI = bicas.tools.batch.BicasProcessingCallInfo.empty(0, 1);

      %============
      % Zero BPCIs
      %============
      if 1
        test(...
          {}, ...
          bicas.swm.SoftwareMode.empty(0, 1), ...
          ZERO_BPCI)
        test(...
          {}, ...
          [SWM_1, SWM_NONMATCHING], ...
          ZERO_BPCI)

        test(...
          {FILE_NONMATCHING}, ...
          [SWM_1], ...
          ZERO_BPCI)
        test(...
          {FILE_NONMATCHING}, ...
          [SWM_1, SWM_NONMATCHING], ...
          ZERO_BPCI)
        test(...
          {FILE_1}, ...
          [SWM_NONMATCHING], ...
          ZERO_BPCI)
      end

      %================
      % Non-zero BPCIs
      %================
      test(...
        {FILE_1}, ...
        [SWM_1], ...
        EXP_BPCI)
      test(...
        {FILE_1, FILE_NONMATCHING}, ...
        [SWM_1,  SWM_NONMATCHING], ...
        EXP_BPCI)
    end



    % Test
    % (1) can find BPCIs only for datasets that overlap in time, and
    % (2) that currentDatasetExtensionDays works.
    function test_time_overlap_currentDatasetExtensionDays(testCase)

      % NOTE: Files DO NOT overlap in time.
      FILE_1       = 'L1R/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240101_V01.cdf';
      FILE_CURRENT = 'BIA/solo_L1_rpw-bia-current-cdag_20231201-20231231_V05.cdf';
      ZERO_BPCI    = bicas.tools.batch.BicasProcessingCallInfo.empty(0, 1);

      SWMP = bicas.tools.batch.TestSwmProcessing();
      SWM_SCI_CURRENT = bicas.swm.SoftwareMode(...
        SWMP, 'CLI-TEST-SCI-CURRENT', 'Human readable purpose', ...
        [...
        bicas.swm.InputDataset('cli_sci',     'SOLO_L1R_RPW-LFR-SURV-CWF-E', 'SCI_cdf'); ...
        bicas.swm.InputDataset('cli_current', 'SOLO_L1_RPW-BIA-CURRENT',     'CUR_cdf') ...
        ], ...
        bicas.swm.OutputDataset(...
        'cli_output', 'SOLO_L2_RPW-LFR-SURV-CWF-E', 'OUT_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        );

      EXP_BPCI = bicas.tools.batch.BicasProcessingCallInfo(...
        SWM_SCI_CURRENT.cliOption, ...
        [...
        bicas.tools.batch.BpciInput(...
        SWM_SCI_CURRENT.inputsList(1).cliOptionHeaderBody, ...
        SWM_SCI_CURRENT.inputsList(1).dsi, ...
        FILE_1); ...
        bicas.tools.batch.BpciInput(...
        SWM_SCI_CURRENT.inputsList(2).cliOptionHeaderBody, ...
        SWM_SCI_CURRENT.inputsList(2).dsi, ...
        FILE_CURRENT), ...
        ], ...
        bicas.tools.batch.BpciOutput(...
        SWM_SCI_CURRENT.outputsList(1).cliOptionHeaderBody, ...
        SWM_SCI_CURRENT.outputsList(1).dsi, ...
        'SOLO_L2_RPW-LFR-SURV-CWF-E___output_dataset.cdf'));

      % No BPCI
      bicas.tools.batch.autocreate_input_BPCIs___UTEST.test(...
        testCase, ...
        {FILE_1; FILE_CURRENT}, ...
        [SWM_SCI_CURRENT], 0, ...
        ZERO_BPCI)

      % One BPCI
      bicas.tools.batch.autocreate_input_BPCIs___UTEST.test(...
        testCase, ...
        {FILE_1; FILE_CURRENT}, ...
        [SWM_SCI_CURRENT], 1, ...
        EXP_BPCI)
    end



    % Check that function only uses latest version of input datasets.
    function test_latest_version(testCase)

      function test(inputDatasetsPathsCa, SwmArray, ExpBpciArray)
        bicas.tools.batch.autocreate_input_BPCIs___UTEST.test(...
          testCase, inputDatasetsPathsCa, SwmArray, 0, ExpBpciArray)
      end

      FILE_1_V1 = 'L1R/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240101_V01.cdf';
      FILE_1_V2 = 'L1R/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240101_V02.cdf';
      FILE_2_V1 = 'L1R/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240102_V01.cdf';
      FILE_3_V2 = 'L1R/solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240103_V02.cdf';

      SWMP = bicas.tools.batch.TestSwmProcessing();
      SWM_1 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI-TEST-1', 'Human readable purpose', ...
        bicas.swm.InputDataset(...
        'cli_input', 'SOLO_L1R_RPW-LFR-SURV-CWF-E', 'IN_cdf' ...
        ), ...
        bicas.swm.OutputDataset(...
        'cli_output', 'SOLO_L2_RPW-LFR-SURV-CWF-E', 'OUT_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        );

      function ExpBpci = exp_BPCI(outputPath)
        ExpBpci = bicas.tools.batch.BicasProcessingCallInfo(...
          SWM_1.cliOption, ...
          bicas.tools.batch.BpciInput(...
          SWM_1.inputsList.cliOptionHeaderBody, ...
          SWM_1.inputsList.dsi, ...
          outputPath), ...
          bicas.tools.batch.BpciOutput(...
          SWM_1.outputsList.cliOptionHeaderBody, ...
          SWM_1.outputsList.dsi, ...
          'SOLO_L2_RPW-LFR-SURV-CWF-E___output_dataset.cdf'));
      end

      % NOTE: ORDER IS NOT GUARANTEED BUT TEST ASSUMES ORDER FOR
      % SIMPLICITY.
      test(...
        {FILE_1_V1, FILE_1_V2, FILE_2_V1, FILE_3_V2}, ...
        SWM_1, ...
        [exp_BPCI(FILE_1_V2); ...
        exp_BPCI(FILE_2_V1); ...
        exp_BPCI(FILE_3_V2)]...
        )

    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function test(testCase, inputDatasetsPathsCa, SwmArray, currentDatasetExtensionDays, ExpBpciArray)

      function filename = get_BPCI_output_path(...
          outputDsi, BpciInputDsmdArray)

        filename = sprintf('%s___output_dataset.cdf', ...
          outputDsi);
      end



      InputDsmdArray = solo.adm.paths_to_DSMD_array(inputDatasetsPathsCa(:));

      % Assert that all paths could be converted to DSMDs (e.g. that
      % filenames conformed to filenaming conventions).
      assert(numel(InputDsmdArray) == numel(inputDatasetsPathsCa))

      % CALL TESTED FUNCTION
      ActBpciArray = bicas.tools.batch.autocreate_input_BPCIs(...
        InputDsmdArray, @get_BPCI_output_path, ...
        SwmArray, currentDatasetExtensionDays);

      testCase.assertEqual(ActBpciArray, ExpBpciArray)
    end



  end    % methods(Static, Access=private)



end
