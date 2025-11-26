%
% matlab.unittest automatic test code for
% bicas.tools.batch.run_BICAS_all_passes().
%
% NOTE: Some tests use both possible values for OUTPUT_VER_ALGORITHM_ID since
% the result should be the same meaning it is easy to implement both tests. This
% slows down the testing while only adding little value. Could change this to
% speed up tests.
%
%
% NLV = Not Latest (Dataset) Version
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef run_BICAS_all_passes___UTEST < matlab.unittest.TestCase
  % PROPOSAL: Replace by tests on bicas.tools.batch.main()
  %   NOTE: Do not want tests on both since much would be duplicated (if being thorough).
  %   PRO: More complete test.
  %   CON: Must change interface of bicas.tools.batch.main() to have BPA object.
  %       PRO: Makes it less convenient as manual command-line command.
  %
  % PROPOSAL: Qualitatively different cases:
  %   outputVerAlgorithmId = 'ABOVE_HIGHEST_USED';
  %   outputVerAlgorithmId = 'HIGHEST_USED';
  %   Match/non-match in reference directory.
  %       In the presence/absence of V01.
  %   BICAS reads latest version of input datasets.
  %       In the presence/absence of V01.
  %   Specifying SWMs to include/exclude
  %       (functionality not yet implemented 2024-01-12)
  %   Using output datasets to create new datasets.
  %   BICAS returns non-zero error code.
  %   TODO: BICAS raises exception.
  %
  % PROBLEM: Can not test for input files being added/deleted during/between passes.
  %   PROPOSAL: Convert iteration inside
  %       bicas.tools.batch.run_BICAS_all_passes() into separate
  %       function. Test code for it can call it multiple times to emulate
  %       calling bicas.tools.batch.run_BICAS_all_passes(), but
  %       add/delete files between calls.
  %       CON: There is communication & "processing" between iterations.
  %           Ex: BpcsAllArray is accumulated.
  %               Next pass does not use this information. Is just a sum over
  %               iteration results.
  %           Ex: tpdFilenamesCa is accumulated.
  %
  % PROBLEM: Can not test for using the correct input files (latest version).
  %   PROPOSAL: Have function return BPCSs or BPCIs.



  %#################
  %#################
  % TEST PARAMETERS
  %#################
  %#################
  properties(TestParameter)
    OUTPUT_VER_ALGORITHM_ID = {'HIGHEST_USED', 'ABOVE_HIGHEST_USED'}
  end



  %############
  %############
  % PROPERTIES
  %############
  %############
  % Additional properties of testCase objects. Needed for setup and teardown
  % methods which store/read their own data from the testCase object.
  properties
    inDir
    outDir
    refDir
  end



  %#######
  %#######
  % SETUP
  %#######
  %#######
  methods(TestMethodSetup)



    function setup(T)
      Fixture = T.applyFixture(...
        matlab.unittest.fixtures.TemporaryFolderFixture);
      testDir = Fixture.Folder;

      T.inDir  = fullfile(testDir, 'in');
      T.outDir = fullfile(testDir, 'out');
      T.refDir = fullfile(testDir, 'ref');

      mkdir(T.inDir)
      mkdir(T.outDir)
      mkdir(T.refDir)

      cd('~')    % Move to any OTHER unrelated directory.
    end



  end



  %##############
  %##############
  % TEST METHODS
  %##############
  %##############
  methods(Test)



    % Zero relevant input files
    %
    function test1_zero_relevant_input(T)
      INPUT_1 = irf.fs.write_empty_file({T.inDir, 'NOT_DATASET.cdf'});
      INPUT_2 = irf.fs.write_empty_file({T.inDir, 'solo_L1_rpw-bia-current_20240101-20240131_V02.cdf'});

      % ===============================
      % Test specifying input directory
      % ===============================
      T.test1(...
        {T.inDir}, '', T.outDir, 'HIGHEST_USED');

      % ====================================================
      % Test specifying explicit (irrelevant) input datasets
      % ====================================================
      T.test1(...
        {INPUT_1, INPUT_2}, '', T.outDir, 'HIGHEST_USED');
    end


    % Test if it can figure out which is the latest version input file, and only
    % use it and not the other input file.
    %
    % 1 LV  input file
    % 1 NLV input file
    % No ref. dir.
    % (There is no V01 input file.)
    % --> 1 output file
    function test1_1LV_1NLV_to_1(T, OUTPUT_VER_ALGORITHM_ID)
      irf.fs.write_empty_file(          {T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e_20240101_V02.cdf'});
      INPUT_2 = irf.fs.write_empty_file({T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e_20240101_V03.cdf'});

      ActBpcsArray = T.test1(...
        {T.inDir}, '', T.outDir, OUTPUT_VER_ALGORITHM_ID);

      assert(numel(ActBpcsArray) == 1)
      % Assert used correct version of input file.
      assert(strcmp(ActBpcsArray(1).Bpci.inputsArray.path, INPUT_2))
      irf.assert.file_exists(fullfile(T.outDir, 'solo_L2_rpw-lfr-surv-cwf-e_20240101_V01.cdf'))
    end



    % 1 input file
    % Ref. dir. file collision depending on output filename versioning.
    % HIGHEST_USED       ==> No output
    % ABOVE_HIGHEST_USED ==> Increment version number
    %
    % NOTE: Ref. dir. file. is not V01! Still blocks output.
    function test1_1_to_01_ref_collision(T, OUTPUT_VER_ALGORITHM_ID)
      irf.fs.write_empty_file({T.inDir,  'solo_L1R_rpw-lfr-surv-cwf-e_20240101_V02.cdf'});
      irf.fs.write_empty_file({T.refDir, 'solo_L2_rpw-lfr-surv-cwf-e_20240101_V05.cdf' });    % Not V01.

      ActBpcsArray = T.test1(...
        {T.inDir}, T.refDir, T.outDir, OUTPUT_VER_ALGORITHM_ID);

      % T.disp_dir_tree(testDir)    % DEBUG

      switch(OUTPUT_VER_ALGORITHM_ID)
        case 'HIGHEST_USED'
          % NOTE: No output file.
          assert(numel(ActBpcsArray) == 0)
        case 'ABOVE_HIGHEST_USED'
          assert(numel(ActBpcsArray) == 1)
          irf.assert.file_exists(fullfile(T.outDir, 'solo_L2_rpw-lfr-surv-cwf-e_20240101_V06.cdf'))
        otherwise
          error('')
      end
    end



    % 1x L1
    % --> 1x L2
    % --> 1x L3
    %
    % NOTE: Output directory is also input directory.
    function test1_1_to_1_to_1(T, OUTPUT_VER_ALGORITHM_ID)
      function assert_actual_result()
        assert(numel(ActBpcsArray) == 2)
        irf.assert.file_exists(fullfile(T.outDir, 'solo_L2_rpw-lfr-surv-cwf-e_20240101_V01.cdf'))
        irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density_20240101_V01.cdf'))
      end

      INPUT_1 = irf.fs.write_empty_file({T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240101_V02.cdf'});

      % ===============================
      % Test specifying input directory
      % ===============================
      ActBpcsArray = T.test1(...
        {T.inDir, T.outDir}, '', T.outDir, OUTPUT_VER_ALGORITHM_ID);
      assert_actual_result()

      % =================================================
      % Test specifying explicit (relevant) input dataset
      % =================================================
      % NOTE: Overwrites old output datasets?
      ActBpcsArray = T.test1(...
        {INPUT_1, T.outDir}, '', T.outDir, OUTPUT_VER_ALGORITHM_ID);
      assert_actual_result()
    end



    % Test the usage of filenames specifying timeintervals on another format
    % than DAY (in this case, SECOND_TO_SECOND).
    function test1_filename_time_interval(T, OUTPUT_VER_ALGORITHM_ID)
      function assert_actual_result()
        T.assertEqual(numel(ActBpcsArray), 2)
        T.assertEqual(ActBpcsArray(1).Bpci.inputsArray(1).path, INPUT_2)

        irf.assert.file_exists(fullfile(T.outDir, 'solo_L2_rpw-lfr-surv-cwf-e_20240102T030405-20240908T070605_V01.cdf'))
        irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density_20240102T030405-20240908T070605_V01.cdf'))
      end

      INPUT_1 = irf.fs.write_empty_file({T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240102T030405-20240908T070605_V02.cdf'});
      INPUT_2 = irf.fs.write_empty_file({T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240102T030405-20240908T070605_V03.cdf'});

      % ===============================
      % Test specifying input directory
      % ===============================
      ActBpcsArray = T.test1(...
        {T.inDir, T.outDir}, '', T.outDir, OUTPUT_VER_ALGORITHM_ID);
      assert_actual_result()
    end



    % 2x L1
    % --> 1x L2 + 1 non-zero error
    % --> 1x L3 (2x L3 without error)
    %
    % IMPLEMENTATION NOTE: Want to test that code continues after BICAS
    % returning non-zero error code. Therefore good to do this when the code
    % is guaranteed to run one more pass after the pass where BICAS failed.
    % This a failsafe against the code executing BPCIs in any order (within
    % a given pass).
    function test1_2_to_1_and_crash_to_1(T, OUTPUT_VER_ALGORITHM_ID)
      INPUT_FILE_1 = irf.fs.write_empty_file(...
        {T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e_20240101_V02.cdf'});  % Crashes
      irf.fs.write_empty_file(...
        {T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e_20240102_V02.cdf'});

      % RUN TESTED CODE
      ActBpcsArray = T.test1(...
        {T.inDir, T.outDir}, '', T.outDir, OUTPUT_VER_ALGORITHM_ID, ...
        {{INPUT_FILE_1}});

      assert(numel(ActBpcsArray) == 3)
      irf.assert.file_exists(fullfile(T.outDir, 'solo_L2_rpw-lfr-surv-cwf-e_20240102_V01.cdf'))
      irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density_20240102_V01.cdf'))

      % Assert one error, one non-error (without assuming order).
      errorCodeArray = [ActBpcsArray.errorCode];
      T.assertEqual(sort(errorCodeArray), [0, 0, 1])

      % Assert correct BPCI failed.
      iError = find(errorCodeArray);
      T.assertTrue(strcmp(ActBpcsArray(iError).Bpci.inputsArray(1).path, INPUT_FILE_1))
    end



    % Process
    % 2x L1 --> 1x L2 --> 2x L3
    % Empty ref. dir..
    function test2_2_to_1_to_2(T, OUTPUT_VER_ALGORITHM_ID)
      irf.fs.write_empty_file({T.inDir, 'solo_L1R_rpw-lfr-surv-cwf-e-cdag_20240101_V02.cdf'});
      irf.fs.write_empty_file({T.inDir, 'solo_L1_rpw-bia-current_20240101-20240131_V02.cdf'});

      ActBpcsArray = T.test2(...
        {T.inDir, T.outDir}, '', T.outDir, OUTPUT_VER_ALGORITHM_ID);

      assert(numel(ActBpcsArray) == 2)
      irf.assert.file_exists(fullfile(T.outDir, 'solo_L2_rpw-lfr-surv-cwf-e_20240101_V01.cdf'))
      irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density_20240101_V01.cdf'))
      irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density-10-seconds_20240101_V01.cdf'))
    end



    % Process
    % 1x L2 --> 2x L3
    % 1x L3 ref. dir. collision.
    % ==> Ref. dir. does not block since there is still one output dataset
    % which is not in the ref. dir..
    function test2_1_to_2_ref_collision(T, OUTPUT_VER_ALGORITHM_ID)
      irf.fs.write_empty_file({T.inDir,  'solo_L2_rpw-lfr-surv-cwf-e_20240101_V01.cdf'});
      irf.fs.write_empty_file({T.refDir, 'solo_L3_rpw-bia-density_20240101_V01.cdf'});

      ActBpcsArray = T.test2(...
        {T.inDir, T.outDir}, T.refDir, T.outDir, OUTPUT_VER_ALGORITHM_ID);

      assert(numel(ActBpcsArray) == 1)
      irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density-10-seconds_20240101_V01.cdf'))
      switch(OUTPUT_VER_ALGORITHM_ID)
        case 'HIGHEST_USED'
          irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density_20240101_V01.cdf'))

        case 'ABOVE_HIGHEST_USED'
          irf.assert.file_exists(fullfile(T.outDir, 'solo_L3_rpw-bia-density_20240101_V02.cdf'))

        otherwise
          error('')
      end

    end



  end    % methods(Test)



  %##########################
  %##########################
  % PRIVATE INSTANCE METHODS
  %##########################
  %##########################
  methods(Access=private)



    % Call BICAS for predefined SWMs, but with datasets specified by caller.
    %
    % Hardcoded SWMs:
    %   1x L1R in --> 1x L2 out
    %   1x L2  in --> 1x L3 out
    %
    % NOTE: Hardcoded DSIs.
    % NOTE: The function does not verify the result. The caller has to create
    %       input files and verify output files and BPCSs.
    %
    % NOTE: Test functions which use this function should be prefixed "test1".
    %
    function ActBpcsArray = test1(T, ...
        inputPathsCa, referenceDir, outputDir, outputVerAlgorithmId, varargin)

      switch(numel(varargin))
        case 0
          errorInputFilesCaCa = cell(0, 1);
        case 1
          errorInputFilesCaCa = varargin{1};
        otherwise
          error('Wrong number of arguments.')
      end

      BICAS_SETTINGS_ARGS_CA = {};
      BICAS_CONFIG_FILE      = 'NO_CONFIG_FILE.conf';
      SETTINGS = [];
      SETTINGS.currentDatasetExtensionDays = 0;

      DSI_1 = 'SOLO_L1R_RPW-LFR-SURV-CWF-E';
      DSI_2 = 'SOLO_L2_RPW-LFR-SURV-CWF-E';
      DSI_3 = 'SOLO_L3_RPW-BIA-DENSITY';

      SWMP = bicas.tools.batch.TestSwmProcessing();

      SWM_1 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_1', 'SWD purpose', ...
        bicas.swm.InputDataset( 'cli_in',  DSI_1, 'IN_cdf'), ...
        bicas.swm.OutputDataset('cli_out', DSI_2, 'OUT_cdf', 'SWD name', 'SWD descr.', '02') ...
        );
      SWM_2 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_2', 'SWD purpose', ...
        bicas.swm.InputDataset( 'cli_in',  DSI_2, 'IN_cdf'), ...
        bicas.swm.OutputDataset('cli_out', DSI_3, 'OUT_cdf', 'SWD name', 'SWD descr.', '03') ...
        );

      BPA = bicas.tools.batch.BicasProcessingAccessTest(...
        [SWM_1; SWM_2], errorInputFilesCaCa);

      % CALL TESTED CODE
      ActBpcsArray = bicas.tools.batch.run_BICAS_all_passes(...
        BPA, BICAS_SETTINGS_ARGS_CA, ...
        BICAS_CONFIG_FILE, outputDir, referenceDir, inputPathsCa(:), ...
        outputVerAlgorithmId, false, [SWM_1; SWM_2], SETTINGS);
    end



    % Call BICAS for predefined SWMs, but with datasets specified by caller.
    %
    % SWMs:
    %   2x L1 in --> 1x L2 out
    %   1x L2 in --> 2x L3 out
    %
    % NOTE: Test functions which use this function should be prefixed "test2".
    %
    function ActBpcsArray = test2(T, inputPathsCa, referenceDir, outputDir, outputVerAlgorithmId)
      BICAS_SETTINGS_ARGS_CA = {};
      BICAS_CONFIG_FILE      = 'NO_CONFIG_FILE.conf';
      SETTINGS = [];
      SETTINGS.currentDatasetExtensionDays = 0;

      DSI_1a = 'SOLO_L1R_RPW-LFR-SURV-CWF-E';
      DSI_1b = 'SOLO_L1_RPW-BIA-CURRENT';
      DSI_2  = 'SOLO_L2_RPW-LFR-SURV-CWF-E';
      DSI_3a = 'SOLO_L3_RPW-BIA-DENSITY';
      DSI_3b = 'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS';

      SWMP = bicas.tools.batch.TestSwmProcessing();

      SWM_1 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_1', 'SWD purpose', ...
        [...
        bicas.swm.InputDataset( 'cli_in1',  DSI_1a, 'IN_1_cdf'), ...
        bicas.swm.InputDataset( 'cli_in2',  DSI_1b, 'IN_2_cdf') ...
        ], ...
        bicas.swm.OutputDataset('cli_out', DSI_2, 'OUT_cdf', 'SWD ', 'SWD ', '02') ...
        );
      SWM_2 = bicas.swm.SoftwareMode(...
        SWMP, 'CLI_SWM_2', 'SWD purpose', ...
        bicas.swm.InputDataset( 'cli_in',  DSI_2, 'IN_cdf'), ...
        [
        bicas.swm.OutputDataset('cli_out1', DSI_3a, 'OUT_1_cdf', 'SWD ', 'SWD ', '03'), ...
        bicas.swm.OutputDataset('cli_out2', DSI_3b, 'OUT_2_cdf', 'SWD ', 'SWD ', '03') ...
        ] ...
        );

      BPA = bicas.tools.batch.BicasProcessingAccessTest([SWM_1; SWM_2], cell(0, 1));

      % CALL TESTED CODE
      ActBpcsArray = bicas.tools.batch.run_BICAS_all_passes(...
        BPA, BICAS_SETTINGS_ARGS_CA, ...
        BICAS_CONFIG_FILE, outputDir, referenceDir, inputPathsCa(:), ...
        outputVerAlgorithmId, false, [SWM_1; SWM_2], SETTINGS);
    end



    % For debugging.
    %
    % Print files in directory tree. For debugging failed tests.
    %
    function disp_dir_tree(T, dirPath)
      ROW_LINE = [repmat('=', 1, 120), '\n'];
      fprintf(ROW_LINE)

      FsoiArray = dir(fullfile(dirPath, '**'));
      for i = 1:numel(FsoiArray)
        Fsoi = FsoiArray(i);
        if ~ismember(Fsoi.name, {'.', '..'})
          path = fullfile(Fsoi.folder, Fsoi.name);
          fprintf('%s\n', path);
        end
      end

      fprintf(ROW_LINE)
    end



  end    % methods(Access=private)



end
