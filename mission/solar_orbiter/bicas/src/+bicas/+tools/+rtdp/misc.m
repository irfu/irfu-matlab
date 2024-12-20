%
% Collection of various functions. Class name is likely temporary. Exact set of
% functions is likely temporary until a better organization of code is found.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef misc
  % PROPOSAL: Support wildcards in config file input dataset paths.
  %   PRO: Easier to link to ROC-mirrored datasets.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Constant)
    % PROPOSAL: Copy from bicas.const.SWD_METADATA('SWD.release.author').
    CONTACT_PERSON = bicas.const.SWD_METADATA('SWD.release.author');
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % De facto top-level function for the bicas.tools.rtdp package, with the
    % added argument automatedTestRun for tests. The nominal user is supposed
    % to call bicas.tools.rtdp.create_RCS_test_data_package() (a trivial
    % wrapper) but tests should call this function.
    %
    %
    % ARGUMENTS
    % =========
    % See bicas.tools.rtdp.create_RCS_test_data_package().
    % automatedTestRun
    %       Whether the function is called by an automated test or not.
    %       NOTE: This is substitute for submitting a class for calling BICAS
    %       (or a function handle), which would then be mocked for tests, which
    %       would be overkill for this application.
    %
    function [rtdpDir, rdtpZipFile] = create_RCS_test_data_package( ...
        outputParentDir, letterVersion, configFile, automatedTestRun)
      %
      % PROPOSAL: Separate function file.
      %   CON: Name will likely be confusing compare to the actual top-level
      %        function file for nominal users.
      %
      % PROPOSAL: Zip package. -- IMPLEMENTED
      %   CON: Can not manually update readme.txt, release_notes.txt
      %   PROPOSAL: Separate command.
      %     CON: ~Superfluous?
      %   PROPOSAL: Flag
      %   PROPOSAL: Only zip for letterVersion="A".
      %
      % PROPOSAL: Use bicas.tools.batch functionality.
      %   Ex: bicas.tools.batch.autocreate_input_BPCIs()
      %
      % PROPOSAL: Check that using the correct directory with source code (bicas_ROC
      %           git repo). Specify in config file.
      %   NOTE: Must be able to run the TEST code in arbitrary irfu-matlab directory.
      %         Tests can specify the current irfu-matlab directory when generating
      %         the config file.
      %
      % PROBLEM: Not checking BICAS bash file.
      %
      % PROPOSAL: Should verify that the config file input datasets do exist.

      assert(islogical(automatedTestRun))


      Bso = bicas.create_default_BSO();
      Bso.make_read_only();

      Swml = bicas.swm.get_SWML(Bso);

      Config = bicas.tools.rtdp.Config(configFile, Swml);

      %====================================
      % ASSERT: Expected BICAS source code
      %====================================
      % NOTE: Only checked for non-automated tests, since automated tests may
      % be run on other MATLAB versions (though in theory they should not).
      % This degrades the value of the test (very) slightly.
      actBicasRootDir = irf.fs.remove_trailing_slash(bicas.utils.get_BICAS_root_dir());
      expBicasRootDir = irf.fs.remove_trailing_slash(Config.get_BICAS_root_dir());
      if ~strcmp(actBicasRootDir, expBicasRootDir)
        error( ...
          ['The actual BICAS root directory and the expected BICAS root' ...
          ' directory are not the same. You might be using the wrong git repo.\n' ...
          'Actual:   "%s"\nExpected: "%s"'], ...
          actBicasRootDir, expBicasRootDir)
      end

      %========================
      % ASSERT: MATLAB version
      %========================
      actMatlabVersion = version('-release');
      expMatlabVersion = bicas.const.OFFICIAL_MATLAB_VERSION;
      if ~strcmp(actMatlabVersion, expMatlabVersion) && ~automatedTestRun
        error( ...
          ['The actual MATLAB version (%s) and the expected MATLAB version' ...
          ' (%s) are not the same.'], ...
          actMatlabVersion, expMatlabVersion)
      end

      %===========================
      % ASSERT: SWD is up-to-date
      %===========================
      % IMPLEMENTATION NOTE: Not really related to generating an RCS test
      % package, but since test packages are generated for deliveries, this is
      % still useful to run.
      TestResult = runtests('bicas.tools.generate_official_SWD_file___UTEST');
      if TestResult.Failed
        error(['The SWD file is not up-to-date (does not correspond to' ...
          ' the BICAS implementation).'])
      end

      % Create root directory.
      rtdpDirName = bicas.tools.rtdp.misc.create_RTDP_directory_name(letterVersion);
      rtdpDir     = bicas.tools.rtdp.misc.mkdir(outputParentDir, rtdpDirName);

      bicas.tools.rtdp.misc.create_readme_file(rtdpDir)
      bicas.tools.rtdp.misc.create_release_notes_file(rtdpDir, letterVersion)

      for iSwm = 1:numel(Swml.List)
        Swm = Swml.List(iSwm);

        bicas.tools.rtdp.misc.create_SWM_directory(rtdpDir, Swm, Config, automatedTestRun)
      end

      rdtpZipFile = bicas.tools.rtdp.misc.zip_RTDP_directory(rtdpDir);
    end



    function zippedRtdp = zip_RTDP_directory(rtdpDir)
      zippedRtdp = [irf.fs.remove_trailing_slash(rtdpDir), '.zip'];

      irf.assert.path_is_available(zippedRtdp)

      zip(zippedRtdp, rtdpDir)
    end



    function create_SWM_directory(parentDir, Swm, Config, automatedTestRun)
      swmDir = bicas.tools.rtdp.misc.mkdir(parentDir, Swm.cliOption);

      inputsDir  = bicas.tools.rtdp.misc.mkdir(swmDir, 'inputs');
      outputsDir = bicas.tools.rtdp.misc.mkdir(swmDir, 'expected_outputs');

      bicasArgsCa        = {Swm.cliOption};

      BpciInputDsmdArray = solo.adm.DSMD.empty(0, 1);
      inputFilenamesCa   = cell(0, 1);
      inputCohbCa        = cell(0, 1);
      for i = 1:numel(Swm.inputsList)
        InputDataset = Swm.inputsList(i);
        cohb         = InputDataset.cliOptionHeaderBody;

        inputSrcFile = Config.get_input_dataset(Swm.cliOption, cohb);
        irf.assert.file_exists(inputSrcFile)

        bicasArgsCa = [bicasArgsCa; {['--', cohb]}; {inputSrcFile}];

        inputFilenamesCa{end+1, 1} = irf.fs.get_name(inputSrcFile);
        inputCohbCa     {end+1, 1} = cohb;

        % Store information for later being able to automatically construct the
        % output filename.
        InputDsmd = solo.adm.paths_to_DSMD_array({inputSrcFile});
        % Check that filename can be parsed.
        assert(isscalar(InputDsmd), ...
          'Can not parse the filename of file "%s".', inputSrcFile)
        BpciInputDsmdArray(end+1, 1) = InputDsmd;

        copyfile(inputSrcFile, inputsDir)
      end
      bicas.tools.rtdp.misc.create_manifest_file(...
        inputsDir, inputFilenamesCa, inputCohbCa)



      outputFilenamesCa = cell(0, 1);
      outputCohbCa      = cell(0, 1);
      for i = 1:numel(Swm.outputsList)
        OutputDataset = Swm.outputsList(i);
        cohb          = OutputDataset.cliOptionHeaderBody;

        % Automatically construct output dataset filename.
        % NOTE: Uses bicas.tools.batch function.
        outputFile    = bicas.tools.batch.get_BPCI_output_path2(...
          BpciInputDsmdArray, solo.adm.DSMD.empty(0, 1), ...
          OutputDataset.dsi, 'HIGHEST_USED', outputsDir, false);

        irf.assert.path_is_available(outputFile)

        bicasArgsCa = [bicasArgsCa; {['--', cohb]}; {outputFile}];

        outputFilenamesCa{end+1, 1} = irf.fs.get_name(inputSrcFile);
        outputCohbCa     {end+1, 1} = cohb;
      end
      bicas.tools.rtdp.misc.create_manifest_file(...
        outputsDir, outputFilenamesCa, outputCohbCa)



      if ~automatedTestRun
        %============
        % CALL BICAS
        %============
        errorCode = bicas.main(bicasArgsCa{:});
        assert(errorCode == 0)
      end
    end



    function rtdpDirName = create_RTDP_directory_name(letterVersion)
      irf.assert.castring_regexp(letterVersion, '[A-Z]')

      % Ex: TESTDATA_RODP_BICAS_V8.2.1A
      bicasVerStr = bicas.const.SWD_METADATA('SWD.release.version');

      rtdpDirName = sprintf('TESTDATA_RODP_BICAS_V%s%s', bicasVerStr, letterVersion);
    end



    function newDirPath = mkdir(parentDir, dirName)
      newDirPath = fullfile(parentDir, dirName);

      % IMPLEMENTATION NOTE: Checking for pre-existing directory for preventing
      % the user from overwriting(?) a pre-existing RTDP directory.
      irf.assert.path_is_available(newDirPath)

      [success, message, errorMessageId] = mkdir(parentDir, dirName);
      assert(success, errorMessageId, message)
    end



    function create_manifest_file(datasetDir, inputFilenamesCa, cohbCa)
      assert(iscell(inputFilenamesCa))
      assert(iscell(cohbCa))
      nDatasets = irf.assert.sizes(...
        inputFilenamesCa, [-1, 1], ...
        cohbCa,           [-1, 1]);

      manifestFile = fullfile(datasetDir, "manifest.txt");

      fid = fopen(manifestFile, "w");
      for i = 1:nDatasets
        inputFilename = inputFilenamesCa{i};
        cohb          = cohbCa{i};
        fprintf(fid, "%s\t%s\n", inputFilename, cohb);
      end
      fclose(fid);
    end



    function create_readme_file(parentDir)
      s = sprintf('Contact person; %s\n', bicas.tools.rtdp.misc.CONTACT_PERSON);

      readmePath = fullfile(parentDir, "readme.txt");
      irf.fs.write_file(readmePath, uint8(s)')
    end



    function create_release_notes_file(parentDir, letterVersion)
      dateStr = string(datetime('now','TimeZone','local','Format','yyyy-mm-d'));
      s = sprintf('Version %s; %s; %s:\n', letterVersion, dateStr, bicas.tools.rtdp.misc.CONTACT_PERSON);

      readmePath = fullfile(parentDir, "release_notes.txt");
      irf.fs.write_file(readmePath, uint8(s)')
    end



  end    % methods(Static)



end
