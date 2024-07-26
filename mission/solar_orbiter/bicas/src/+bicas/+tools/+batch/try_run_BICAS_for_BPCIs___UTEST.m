%
% matlab.unittest automatic test code for
% bicas.tools.batch.try_run_BICAS_for_BPCIs().
%
% The code works, but is too complicated to be worth it in hindsight?!
%
% NOTE: Does not test BICAS failing (returning on-zero error code).
% NOTE: Does not really test argument automountTriggerPathsCa.
%
% IMPLEMENTATION NOTE: Code has to use (almost) real DSIs due to BICAS data
% structure assertions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef try_run_BICAS_for_BPCIs___UTEST < matlab.unittest.TestCase



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    function test_zero_BPCIs(testCase)

      BPA = bicas.tools.batch.BicasProcessingAccessTest(...
        bicas.swm.SoftwareMode.empty(0, 1), []);

      bicas.tools.batch.try_run_BICAS_for_BPCIs___UTEST.test(...
        testCase, BPA, ...
        bicas.tools.batch.BicasProcessingCallInfo.empty(0, 1), ...
        bicas.tools.batch.BicasProcessingCallSummary.empty(0, 1))

    end



    function test_one_BPCI(testCase)

      % NOTE: Changes current working directory.
      testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture)

      % NOTE: Dependence on BICAS.
      SWMP = bicas.tools.batch.TestSwmProcessing();

      DSI_1 = 'SOLO_L1R_RPW-LFR-SURV-CWF-E';
      DSI_2 = 'SOLO_L2_RPW-LFR-SURV-CWF-E';

      SWM = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_1', 'SWD purpose', ...
        bicas.swm.InputDataset(...
        'cli_input', DSI_1, 'IN_cdf' ...
        ), ...
        bicas.swm.OutputDataset(...
        'cli_output', DSI_2, 'OUT_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        );
      BPA = bicas.tools.batch.BicasProcessingAccessTest(SWM, []);

      OUTPUT_DIR = fullfile(pwd, 'dir2');
      mkdir(OUTPUT_DIR)

      INPUT_FILE  = fullfile(pwd, 'dir1', 'input_dataset.cdf');
      OUTPUT_FILE = fullfile(OUTPUT_DIR, 'output_dataset.cdf');

      BPCI = bicas.tools.batch.BicasProcessingCallInfo(...
        'CLI_SWM_1', ...
        bicas.tools.batch.BpciInput( 'cli_input',  DSI_1, INPUT_FILE), ...
        bicas.tools.batch.BpciOutput('cli_output', DSI_2, OUTPUT_FILE));

      % NOTE: BPCS should only contain filenames.
      EXP_BPCS = bicas.tools.batch.BicasProcessingCallSummary(...
        BPCI, 0);



      bicas.tools.batch.try_run_BICAS_for_BPCIs___UTEST.test(...
        testCase, BPA, ...
        BPCI, ...
        EXP_BPCS)
    end



    function test_complex(testCase)

      % NOTE: Dependence on BICAS.
      SWMP = bicas.tools.batch.TestSwmProcessing();

      DSI_1 = 'SOLO_L1R_RPW-LFR-SURV-CWF-E';
      DSI_2 = 'SOLO_L1R_RPW-LFR-SURV-SWF-E';
      DSI_3 = 'SOLO_L1R_RPW-LFR-SBM1-CWF-E';

      SWM_1 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_1', 'SWD purpose', ...
        [
        bicas.swm.InputDataset('cli_input_1', DSI_1, 'IN_1_cdf'), ...
        bicas.swm.InputDataset('cli_input_2', DSI_2, 'IN_2_cdf') ...
        ], ...
        bicas.swm.OutputDataset(...
        'cli_output', DSI_3, 'OUT_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        );
      SWM_2 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_2', 'Human readable purpose', ...
        bicas.swm.InputDataset('cli_input', DSI_1, 'IN_cdf'), ...
        [...
        bicas.swm.OutputDataset(...
        'cli_output_1', DSI_2, 'OUT_1_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        bicas.swm.OutputDataset(...
        'cli_output_2', DSI_3, 'OUT_2_cdf', ...
        'SWD name', 'SWD description', '02' ...
        ) ...
        ] ...
        );

      BPA = bicas.tools.batch.BicasProcessingAccessTest([SWM_1; SWM_2], []);



      function [Bpci, Bpcs] = get_BPCI_1(inputFile1, inputFile2, outputFile)
        Bpci = bicas.tools.batch.BicasProcessingCallInfo(...
          'CLI_SWM_1', ...
          [
          bicas.tools.batch.BpciInput('cli_input_1',  DSI_1, inputFile1); ...
          bicas.tools.batch.BpciInput('cli_input_2',  DSI_2, inputFile2) ...
          ], ...
          bicas.tools.batch.BpciOutput('cli_output', DSI_3, outputFile));

        Bpcs = bicas.tools.batch.BicasProcessingCallSummary(...
          Bpci, 0);
      end

      function [Bpci, Bpcs] = get_BPCI_2(inputFile, outputFile1, outputFile2)
        Bpci = bicas.tools.batch.BicasProcessingCallInfo(...
          'CLI_SWM_2', ...
          bicas.tools.batch.BpciInput('cli_input',  DSI_1, inputFile), ...
          [...
          bicas.tools.batch.BpciOutput('cli_output_1', DSI_2, outputFile1); ...
          bicas.tools.batch.BpciOutput('cli_output_2', DSI_3, outputFile2) ...
          ]...
          );
        Bpcs = bicas.tools.batch.BicasProcessingCallSummary(...
          Bpci, 0);
      end

      function C = prepare_test()
        % Create temporary directory.
        % NOTE: Changes current working directory.
        testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture)

        % NOTE: Only relative paths.
        % Using subdirectories just to test paths (not just filenames).
        [C.BPCI_1a, C.EXP_BPCS_1a] = get_BPCI_1('11/i1', '11/i2', '11/o1'      );
        [C.BPCI_1b, C.EXP_BPCS_1b] = get_BPCI_1('12/i1', '12/i2', '12/o1'      );
        [C.BPCI_2a, C.EXP_BPCS_2a] = get_BPCI_2('21/i1',          '21/o1', '21/o2');
        [C.BPCI_2b, C.EXP_BPCS_2b] = get_BPCI_2('22/i1',          '22/o1', '22/o2');
      end

      % NOTE: BPCS should only contain filenames.
      % NOTE: Comparisons rely on that the order in arrays stays the same
      %       even though the code does not really guarantee it. Could
      %       sort CAs of strings in relevant data structures to make
      %       comparisons easier.
      % IMPLEMENTATION NOTE: Makes sure to have a new temporary
      % directory for every call to
      % bicas.tools.batch.try_run_BICAS_for_BPCIs___UTEST.test()
      % so that no files from previous call can affect the new call.

      C = prepare_test();
      bicas.tools.batch.try_run_BICAS_for_BPCIs___UTEST.test(...
        testCase, BPA, ...
        C.BPCI_1a, ...
        C.EXP_BPCS_1a)

      C = prepare_test();
      bicas.tools.batch.try_run_BICAS_for_BPCIs___UTEST.test(...
        testCase, BPA, ...
        C.BPCI_2a, ...
        C.EXP_BPCS_2a)

      C = prepare_test();
      bicas.tools.batch.try_run_BICAS_for_BPCIs___UTEST.test(...
        testCase, BPA, ...
        [    C.BPCI_1a;     C.BPCI_2a;     C.BPCI_1b;     C.BPCI_2b], ...
        [C.EXP_BPCS_1a; C.EXP_BPCS_2a; C.EXP_BPCS_1b; C.EXP_BPCS_2b])
    end



  end    % methods(Test)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Make test call to function. Create input files before, check existence
    % of output files after.
    %
    % NOTE: Can not create test directory in this function, since the BPCIs
    % contain the full paths to datasets.
    function test(testCase, Bpa, BpciArray, ExpBpcsArray)

      % ===========================================
      % Create empty input files specified in BPCIs
      % ===========================================
      function create_input_files(BpciArray)
        for iBpci = 1:numel(BpciArray)
          Bpci = BpciArray(iBpci);
          for iFile = 1:numel(Bpci.inputsArray)

            path = Bpci.inputsArray(iFile).path;

            irf.fs.write_empty_file({path});
          end
        end
      end

      function main()
        assert(iscolumn(ExpBpcsArray))

        CONFIG_FILE            = 'NO_CONFIG_FILE';
        BICAS_SETTINGS_ARGS_CA = {};

        create_input_files(BpciArray)



        %======================
        % CALL TESTED FUNCTION
        %======================
        ActBpcsArray = bicas.tools.batch.try_run_BICAS_for_BPCIs(...
          Bpa, BpciArray, ...
          CONFIG_FILE, cell(0, 1), BICAS_SETTINGS_ARGS_CA);



        testCase.assertEqual(...
          ActBpcsArray, ...
          ExpBpcsArray)

        %=====================================
        % Assert that BPCI output files exist
        %=====================================
        for i = 1:numel(BpciArray)
          Bpci         = BpciArray(i);
          outputPathCa = {Bpci.outputsArray(:).path};

          for iFile = 1:numel(outputPathCa)
            irf.assert.file_exists(outputPathCa{iFile})
          end
        end
      end

      main()
    end



  end    % methods(Static, Access=private)



end
